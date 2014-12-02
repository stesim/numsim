#include "communication.h"
#include "metamath/metamath.h"
#include <cassert>

using namespace mm;

const MPI_Datatype Communication::REAL_MPI = ( sizeof(real) == sizeof(double) ? MPI_DOUBLE : MPI_FLOAT );
const MPI_Comm     Communication::comm = MPI_COMM_WORLD;

MultiIndex         Communication::s_Size;
MultiIndex         Communication::s_Rank;
int                Communication::s_iSize;
int                Communication::s_iRank;

void Communication::init( int& argc, char**& argv )
{
	MPI_Init( &argc, &argv );

	MPI_Comm_size( comm, &s_iSize );
	MPI_Comm_rank( comm, &s_iRank );

	s_Size = sizeTo2D( s_iSize );

	s_Rank.y = s_iRank / s_Size.x;
	s_Rank.x = s_iRank % s_Size.x;
}

void Communication::finalize()
{
	MPI_Finalize();
}

const MultiIndex& Communication::size()
{
	return s_Size;
}

const MultiIndex& Communication::rank()
{
	return s_Rank;
}

int Communication::linearRank()
{
	return s_iRank;
}

void Communication::exchangePressure( GridFunction& p )
{
	exchangeFunctionBoundary( p, MultiIndex::ONE, 1 );
}

void Communication::exchangeVelocities( GridFunction& u, GridFunction& v )
{
	exchangeFunctionBoundary( u, MultiIndex( 2, 1 ), 2 );
	exchangeFunctionBoundary( v, MultiIndex( 1, 2 ), 3 );
}

real Communication::min( const real& local )
{
	real res;
	MPI_Allreduce( &local, &res, 1, REAL_MPI, MPI_MIN, comm );
	return res;
}

real Communication::sum( const real& local )
{
	real res;
	MPI_Allreduce( &local, &res, 1, REAL_MPI, MPI_SUM, comm );
	return res;
}

MultiIndex Communication::sizeTo2D( int size )
{
	MultiIndex res;
	MPI_Dims_create( size, 2, res.coord );
	return res;
}

int Communication::rankTo1D( const MultiIndex& rank, bool cyclic )
{
	MultiIndex rank2D = rank;
	if( rank2D.x < 0 )
	{
		rank2D.x = cyclic ? rank2D.x + s_Size.x : -1;
	}
	else if( rank2D.x >= s_Size.x )
	{
		rank2D.x = cyclic ? rank2D.x - s_Size.x : -1;
	}

	if( rank2D.y < 0 )
	{
		rank2D.y = cyclic ? rank2D.y + s_Size.y : -1;
	}
	else if( rank2D.y >= s_Size.y )
	{
		rank2D.y = cyclic ? rank2D.y - s_Size.y : -1;
	}

	if( rank2D.x < 0 || rank2D.y < 0 )
	{
		return -1;
	}
	else
	{
		return ( rank2D.y * s_Size.x + rank2D.x );
	}
}

// alternate between send and receive in successive processes
// e.g. Proc1 sends to Proc2 while Proc2 receives from Proc1
#define ALTERNATE_SEND_RECV(arg,sendCode,recvCode) if( arg ) {\
	sendCode\
	recvCode\
} else {\
	recvCode\
	sendCode\
}

void Communication::exchangeFunctionBoundary( GridFunction& f,
		const MultiIndex& readOffset, int tag )
{
	static const int rankLeft = rankTo1D( s_Rank - MultiIndex( 1, 0 ) );
	static const int rankRight = rankTo1D( s_Rank + MultiIndex( 1, 0 ) );
	static const int rankTop = rankTo1D( s_Rank + MultiIndex( 0, 1 ) );
	static const int rankBottom = rankTo1D( s_Rank - MultiIndex( 0, 1 ) );
	static const bool sendRecvOrder = ( ( ( s_Rank.x + s_Rank.y ) % 2 ) == 0 );

	MPI_Status stat;

	// itermediate (linear) buffer for vertical boundaries
	GridFunction transposedBoundary( MultiIndex( f.size().y - 2, 1 ) );

	ALTERNATE_SEND_RECV( sendRecvOrder,
		// send left
		if( rankLeft >= 0 )
		{
			FunctionView<GridFunction> fViewRead = view( f,
				MultiIndex( readOffset.x, 1 ),
				MultiIndex( readOffset.x + 1, f.size().y - 1 ) );
			transposedBoundary = transpose( fViewRead );
			MPI_Send( &transposedBoundary( 0, 0 ),
				transposedBoundary.size().x,
				REAL_MPI,
				rankLeft,
				tag,
				comm );
		}
	,
		// receive right
		if( rankRight >= 0 )
		{
			FunctionView<GridFunction> fViewWrite = view( f,
				MultiIndex( f.size().x - 1, 1 ),
				f.size() - MultiIndex( 0, 1 ) );
			MPI_Recv( &transposedBoundary( 0, 0 ),
				transposedBoundary.size().x,
				REAL_MPI,
				rankRight,
				tag,
				comm,
				&stat );
			fViewWrite = transpose( transposedBoundary );
		}
	)

	ALTERNATE_SEND_RECV( sendRecvOrder,
		// send right
		if( rankRight >= 0 )
		{
			FunctionView<GridFunction> fViewRead = view( f,
				MultiIndex( f.size().x - readOffset.x - 1, 1 ),
				f.size() - MultiIndex( readOffset.x, 1 ) );
			transposedBoundary = transpose( fViewRead );
			MPI_Send( &transposedBoundary( 0, 0 ),
				transposedBoundary.size().x,
				REAL_MPI,
				rankRight,
				tag,
				comm );
		}
	,
		// receive left
		if( rankLeft >= 0 )
		{
			FunctionView<GridFunction> fViewWrite = view( f,
				MultiIndex( 0, 1 ),
				MultiIndex( 1, f.size().y - 1 ) );
			MPI_Recv( &transposedBoundary( 0, 0 ),
				transposedBoundary.size().x,
				REAL_MPI,
				rankLeft,
				tag,
				comm,
				&stat );
			fViewWrite = transpose( transposedBoundary );
		}
	)

	ALTERNATE_SEND_RECV( sendRecvOrder,
		// send bottom
		if( rankBottom >= 0 )
		{
			MPI_Send( &f( 1, readOffset.y ),
				f.size().x - 2,
				REAL_MPI,
				rankBottom,
				tag,
				comm );
		}
	,
		// receive top
		if( rankTop >= 0 )
		{
			MPI_Recv( &f( 1, f.size().y - 1 ),
				f.size().x - 2,
				REAL_MPI,
				rankTop,
				tag,
				comm,
				&stat );
		}
	)

	ALTERNATE_SEND_RECV( sendRecvOrder,
		// send top
		if( rankTop >= 0 )
		{
			MPI_Send( &f( 1, f.size().y - readOffset.y - 1 ),
				f.size().x - 2,
				REAL_MPI,
				rankTop,
				tag,
				comm );
		}
	,
		// receive bottom
		if( rankBottom >= 0 )
		{
			MPI_Recv( &f( 1, 0 ),
				f.size().x - 2,
				REAL_MPI,
				rankBottom,
				tag,
				comm,
				&stat );
		}
	)
}

#undef ALTERNATE_SEND_RECV

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

int                Communication::s_iRankLeft;
int                Communication::s_iRankRight;
int                Communication::s_iRankTop;
int                Communication::s_iRankBottom;
bool               Communication::s_SendRecvOrder;

void Communication::init( int& argc, char**& argv )
{
	MPI_Init( &argc, &argv );

	MPI_Comm_size( comm, &s_iSize );
	MPI_Comm_rank( comm, &s_iRank );

	s_Size = sizeTo2D( s_iSize );

	s_Rank.y = s_iRank / s_Size.x;
	s_Rank.x = s_iRank % s_Size.x;
	
	s_iRankLeft = rankTo1D( s_Rank - MultiIndex( 1, 0 ) );
	s_iRankRight = rankTo1D( s_Rank + MultiIndex( 1, 0 ) );
	s_iRankTop = rankTo1D( s_Rank + MultiIndex( 0, 1 ) );
	s_iRankBottom = rankTo1D( s_Rank - MultiIndex( 0, 1 ) );
	s_SendRecvOrder = ( ( ( s_Rank.x + s_Rank.y ) % 2 ) == 0 );
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
	sendCode;\
	recvCode;\
} else {\
	recvCode;\
	sendCode;\
}

void Communication::exchangeFunctionBoundary( GridFunction& f,
		const MultiIndex& readOffset )
{
	MPI_Status stat;

	// itermediate linear buffer for vertical boundaries
	GridFunction transposedBoundary( MultiIndex( f.size().y - 2, 1 ) );

	ALTERNATE_SEND_RECV( s_SendRecvOrder,
			sendLeft( f, readOffset.x, transposedBoundary, true ),
			recvRight( f, stat, transposedBoundary, true ) );

	ALTERNATE_SEND_RECV( s_SendRecvOrder,
			sendRight( f, readOffset.x, transposedBoundary, true ),
			recvLeft( f, stat, transposedBoundary, true ) );

	ALTERNATE_SEND_RECV( s_SendRecvOrder,
			sendTop( f, readOffset.y, true ),
			recvBottom( f, stat, true ) );

	ALTERNATE_SEND_RECV( s_SendRecvOrder,
			sendBottom( f, readOffset.y, true ),
			recvTop( f, stat, true ) );
}

#undef ALTERNATE_SEND_RECV

void Communication::sendBoundary( GridFunction& f, Boundary boundary,
		int readOffset, bool excludeCorners )
{
	switch( boundary )
	{
		case BOUNDARY_LEFT:
		{
			GridFunction transposedBoundary(
					MultiIndex( f.size().y - 2 * excludeCorners, 1 ) );
			sendLeft( f, readOffset, transposedBoundary, excludeCorners );
			break;
		}
		case BOUNDARY_TOP:
		{
			sendTop( f, readOffset, excludeCorners );
			break;
		}
		case BOUNDARY_RIGHT:
		{
			GridFunction transposedBoundary(
					MultiIndex( f.size().y - 2 * excludeCorners, 1 ) );
			sendRight( f, readOffset, transposedBoundary, excludeCorners );
			break;
		}
		case BOUNDARY_BOTTOM:
		{
			sendBottom( f, readOffset, excludeCorners );
			break;
		}
		default:
			break;
	}
}

void Communication::receiveBoundary(
		GridFunction& f, Boundary boundary, bool excludeCorners )
{
	MPI_Status stat;
	switch( boundary )
	{
		case BOUNDARY_LEFT:
		{
			GridFunction transposedBoundary(
					MultiIndex( f.size().y - 2 * excludeCorners, 1 ) );
			recvLeft( f, stat, transposedBoundary, excludeCorners );
			break;
		}
		case BOUNDARY_TOP:
		{
			recvTop( f, stat, excludeCorners );
			break;
		}
		case BOUNDARY_RIGHT:
		{
			GridFunction transposedBoundary(
					MultiIndex( f.size().y - 2 * excludeCorners, 1 ) );
			recvRight( f, stat, transposedBoundary, excludeCorners );
			break;
		}
		case BOUNDARY_BOTTOM:
		{
			recvBottom( f, stat, excludeCorners );
			break;
		}
		default:
			break;
	}
}

void Communication::sendLeft( const GridFunction& f, int offset,
		GridFunction& transposedBoundary, bool excludeCorners )
{
	if( s_iRankLeft >= 0 )
	{
		FunctionView<const GridFunction> fViewRead = view( f,
			MultiIndex( offset, excludeCorners ),
			MultiIndex( offset + 1, f.size().y - excludeCorners ) );
		transposedBoundary = transpose( fViewRead );
		MPI_Send( &transposedBoundary( 0, 0 ),
			f.size().y - 2 * excludeCorners,
			REAL_MPI,
			s_iRankLeft,
			1,
			comm );
	}
}

void Communication::sendTop(
		const GridFunction& f, int offset, bool excludeCorners )
{
	if( s_iRankTop >= 0 )
	{
		MPI_Send( &f( excludeCorners, f.size().y - offset - 1 ),
			f.size().x - 2 * excludeCorners,
			REAL_MPI,
			s_iRankTop,
			1,
			comm );
	}
}

void Communication::sendRight( const GridFunction& f, int offset,
		GridFunction& transposedBoundary, bool excludeCorners )
{
	if( s_iRankRight >= 0 )
	{
		FunctionView<const GridFunction> fViewRead = view( f,
			MultiIndex( f.size().x - offset - 1, excludeCorners ),
			f.size() - MultiIndex( offset, excludeCorners ) );
		transposedBoundary = transpose( fViewRead );
		MPI_Send( &transposedBoundary( 0, 0 ),
			f.size().y - 2 * excludeCorners,
			REAL_MPI,
			s_iRankRight,
			1,
			comm );
	}
}

void Communication::sendBottom(
		const GridFunction& f, int offset, bool excludeCorners )
{
	if( s_iRankBottom >= 0 )
	{
		MPI_Send( &f( excludeCorners, offset ),
			f.size().x - 2 * excludeCorners,
			REAL_MPI,
			s_iRankBottom,
			1,
			comm );
	}
}

void Communication::recvLeft( GridFunction& f, MPI_Status& stat,
		GridFunction& transposedBoundary, bool excludeCorners )
{
	if( s_iRankLeft >= 0 )
	{
		FunctionView<GridFunction> fViewWrite = view( f,
			MultiIndex( 0, excludeCorners ),
			MultiIndex( 1, f.size().y - excludeCorners ) );
		MPI_Recv( &transposedBoundary( 0, 0 ),
			f.size().y - 2 * excludeCorners,
			REAL_MPI,
			s_iRankLeft,
			1,
			comm,
			&stat );
		fViewWrite = transpose( transposedBoundary );
	}
}

void Communication::recvTop(
		GridFunction& f, MPI_Status& stat, bool excludeCorners )
{
	if( s_iRankTop >= 0 )
	{
		MPI_Recv( &f( excludeCorners, f.size().y - 1 ),
			f.size().x - 2 * excludeCorners,
			REAL_MPI,
			s_iRankTop,
			1,
			comm,
			&stat );
	}
}

void Communication::recvRight( GridFunction& f, MPI_Status& stat,
		GridFunction& transposedBoundary, bool excludeCorners )
{
	if( s_iRankRight >= 0 )
	{
		FunctionView<GridFunction> fViewWrite = view( f,
			MultiIndex( f.size().x - 1, excludeCorners ),
			f.size() - MultiIndex( 0, excludeCorners ) );
		MPI_Recv( &transposedBoundary( 0, 0 ),
			f.size().y - 2 * excludeCorners,
			REAL_MPI,
			s_iRankRight,
			1,
			comm,
			&stat );
		fViewWrite = transpose( transposedBoundary );
	}
}

void Communication::recvBottom(
		GridFunction& f, MPI_Status& stat, bool excludeCorners )
{
	if( s_iRankBottom >= 0 )
	{
		MPI_Recv( &f( excludeCorners, 0 ),
			f.size().x - 2 * excludeCorners,
			REAL_MPI,
			s_iRankBottom,
			1,
			comm,
			&stat );
	}
}

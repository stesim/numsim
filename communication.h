#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include "typedef.h"
#include "mpi.h"

class Communication
{
public:
	//! Initialize MPI
	static void init( int& argc, char**& argv );

	//! Finalize MPI
	static void finalize();

	//! Return 2D MPI size
	static const MultiIndex& size();

	//! Return 2D MPI rank of the current process
	static const MultiIndex& rank();

	//! Return linear MPI rank
	static int linearRank();

	//! Exchange pressure boundary values
	static void exchangePressure( GridFunction& p );

	//! Exchange (intermediate) velocity boundary values
	static void exchangeVelocities( GridFunction& u, GridFunction& v );

	//! Determine and return the minimum of the local values provided by all processes
	static real min( const real& local );

	//! Determine and return the sum of the local values provided by all processes
	static real sum( const real& local );

private:
	//! Determine 2D grid dimensions from number of processes
	static MultiIndex sizeTo2D( int size );

	//! Calculate linear rank from 2D coordinates
	static inline int rankTo1D( const MultiIndex& rank, bool cyclic = false );

	//! Exchange boundary values of \pm{f}
	//! \param f Function to operate on
	//! \param readOffset Distance of the data that is read from the boundary in the respective direction
	//! \param tag MPI tag
	static inline void exchangeFunctionBoundary( GridFunction& f,
			const MultiIndex& readOffset, int tag );

private:
	static const MPI_Datatype REAL_MPI;
	static const MPI_Comm comm;
	static MultiIndex s_Size;
	static MultiIndex s_Rank;
	static int s_iSize;
	static int s_iRank;
};

#endif

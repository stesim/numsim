#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include "typedef.h"
#include "mpi.h"

class Communication
{
public:
	enum Boundary
	{
		BOUNDARY_NONE   = 0,
		BOUNDARY_LEFT   = 1,
		BOUNDARY_TOP    = 2,
		BOUNDARY_RIGHT  = 4,
		BOUNDARY_BOTTOM = 8,
		BOUNDARY_ALL    = BOUNDARY_LEFT | BOUNDARY_TOP | BOUNDARY_RIGHT | BOUNDARY_RIGHT
	};

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

	//! Exchange boundary values of \pm{f}
	//! \param f Function to operate on
	//! \param readOffset Distance of the data that is read from the boundary in the respective direction
	static void exchangeFunctionBoundary( GridFunction& f,
			const MultiIndex& readOffset );

	//! Send boundary to neighbor \pm{f}
	//! \param f Function to operate on
	//! \param boundary Which boundary to propagate
	//! \param readOffset Distance of the data that is read from the boundary in the respective direction
	static void sendBoundary( GridFunction& f, Boundary boundary,
			int readOffset, bool excludeCorners );

	//! Receive boundary from neighbor \pm{f}
	//! \param f Function to operate on
	//! \param boundary Which boundary to propagate
	//! \param readOffset Distance of the data that is read from the boundary in the respective direction
	static void receiveBoundary(
			GridFunction& f, Boundary boundary, bool excludeCorners );

	//! Determine and return the minimum of the local values provided by all processes
	static real min( const real& local );

	//! Determine and return the sum of the local values provided by all processes
	static real sum( const real& local );

private:
	//! Determine 2D grid dimensions from number of processes
	static MultiIndex sizeTo2D( int size );

	//! Calculate linear rank from 2D coordinates
	static inline int rankTo1D( const MultiIndex& rank, bool cyclic = false );

	static inline void sendLeft( const GridFunction& f, int offset,
			GridFunction& transposedBoundary, bool excludeCorners );
	static inline void sendTop(
			const GridFunction& f, int offset, bool excludeCorners );
	static inline void sendRight( const GridFunction& f, int offset,
			GridFunction& transposedBoundary, bool excludeCorners );
	static inline void sendBottom(
			const GridFunction& f, int offset, bool excludeCorners );

	static inline void recvLeft( GridFunction& f, MPI_Status& stat,
			GridFunction& transposedBoundary, bool excludeCorners );
	static inline void recvTop(
			GridFunction& f, MPI_Status& stat, bool excludeCorners );
	static inline void recvRight( GridFunction& f, MPI_Status& stat,
			GridFunction& transposedBoundary, bool excludeCorners );
	static inline void recvBottom(
			GridFunction& f, MPI_Status& stat, bool excludeCorners );

private:
	static const MPI_Datatype REAL_MPI;
	static const MPI_Comm comm;

	static MultiIndex s_Size;
	static MultiIndex s_Rank;
	static int s_iSize;
	static int s_iRank;

	static int s_iRankLeft;
	static int s_iRankTop;
	static int s_iRankRight;
	static int s_iRankBottom;
	static bool s_SendRecvOrder;
};

#endif

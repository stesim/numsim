#include "IO.h"
#include "computation.h"
#include "solver.h"
#include "communication.h"
#include <iostream>
#include <vector>
#include "metamath/mmfunction.h"
#include "metamath/metamath.h"
#include <cmath>
#include <cassert>

#ifdef _CPP11
#include <thread>
#include <chrono>
inline std::chrono::high_resolution_clock::time_point getTime()
{
	return std::chrono::high_resolution_clock::now();
}
#define tick() std::chrono::high_resolution_clock::time_point _timeStart = getTime()
#define tock() ( std::chrono::duration_cast<std::chrono::duration<real>>( getTime() - _timeStart ) ).count()
#else
#define tick()
#define tock() 0
#endif

#define DO_ON_RANK(rnk,todo) if( Communication::linearRank() == rnk ) { todo }
//#define DO_ON_RANK(rnk,todo) todo

#define PRINT_ON_RANK 1

using namespace mm;

void simulate( const char* filename, int instance )
{
	MultiIndex commSize = Communication::size();
	MultiIndex commRank = Communication::rank();

	bool leftBoundary = ( commRank.x == 0 );
	bool rightBoundary = ( commRank.x == commSize.x - 1 );
	bool bottomBoundary = ( commRank.y == 0 );
	bool topBoundary = ( commRank.y == commSize.y - 1 );

	const char* paramFile = "inputvals";
	if( filename != NULL )
	{
		paramFile = filename;
	}

	DO_ON_RANK( PRINT_ON_RANK,
		std::cout << "[Initialization] Loading parameters from file '" << paramFile << "'." << std::endl;
	);

	Params params = IO::readInputFile( paramFile );

	DO_ON_RANK( PRINT_ON_RANK,
		std::cout << "[Initialization] Parameter overview:" << std::endl;
		params.print( std::cout );
	);

	MultiIndex localGridSize = params.gridSize;
	localGridSize.x /= commSize.x;
	localGridSize.y /= commSize.y;

	params.gridSize.x = localGridSize.x * commSize.x;
	params.gridSize.y = localGridSize.y * commSize.y;

	DO_ON_RANK( PRINT_ON_RANK,
		std::cout << "[Parallelization]" << std::endl
			<< "  Number of processes:   " << commSize.x << " x " << commSize.y << std::endl
			<< "  Current process:       " << "[" << commRank.x << "," << commRank.y << "]" << std::endl
			<< "  Local grid size:       " << localGridSize.x << " x " << localGridSize.y << std::endl;
	);

	Point h( params.domainSize.x / params.gridSize.x,
			params.domainSize.y / params.gridSize.y );

	Point localGridOffset( localGridSize.x * commRank.x * h.x,
			localGridSize.y * commRank.y * h.y );

	GridFunction u( localGridSize + MultiIndex( 3, 2 ) );
	GridFunction v( localGridSize + MultiIndex( 2, 3 ) );
	GridFunction f( localGridSize + MultiIndex( 3, 2 ) );
	GridFunction g( localGridSize + MultiIndex( 2, 3 ) );
	GridFunction p( localGridSize + MultiIndex( 2, 2 ) );
	GridFunction rhs( localGridSize );

	u = constant( params.initialVelocity.x );
	v = constant( params.initialVelocity.y );
	p = constant( params.initialPressure );

	f = constant( params.initialVelocity.x );
	g = constant( params.initialVelocity.y );

	tick();

	real t = 0.0;
	index_t step = 0;
	real tOut = 0.0;
	while( t < params.T )
	{
		real timeBeginStep = tock();

		real dt = Computation::computeTimeStep(
				u, v, h, params.Re, params.tau );

		Computation::setVelocityBoundary( u, v,
				leftBoundary, topBoundary, rightBoundary, bottomBoundary );
		
		// compute intermediate velocities
		Computation::computeMomentumEquations(
				f, g, u, v, h, dt, params.Re, params.alpha );
		Computation::setVelocityBoundary( f, g,
				leftBoundary, topBoundary, rightBoundary, bottomBoundary );
		Communication::exchangeVelocities( f, g );

		Computation::computeRightHandSide( rhs, f, g, h, dt );

		real timeBeginSOR = tock();

		index_t iter = 0;
		real res = params.eps + 1.0;
		while( iter < params.maxIter && res > params.eps )
		{
			Computation::setPressureBoundary( p,
					leftBoundary, topBoundary, rightBoundary, bottomBoundary );

			Solver::SORSubcycle( p, rhs, h, params.omega, true );
			Communication::exchangePressure( p );
			Solver::SORSubcycle( p, rhs, h, params.omega, false );
			Communication::exchangePressure( p );

			res = Solver::computeResidual( p, rhs, params.gridSize, h );

			++iter;
		}

		real timeEndSOR = tock();

		Computation::computeNewVelocities( u, v, f, g, p, h, dt );
		Communication::exchangeVelocities( u, v );

		real timeEndStep = tock();

		t += dt;
		++step;

		tOut += dt;
		if( step <= 1 || tOut > params.deltaVec || t >= params.T )
		{
			tOut -= params.deltaVec;
			IO::writeRawOutput( u, v, p, h, localGridOffset, step, commRank );
			/*
			IO::writeVTKFile( params.gridSize, u, v,
					p, h, params, instance, step );
			*/

			DO_ON_RANK( PRINT_ON_RANK,
				std::cout << "[Progress] Instance " << instance << ", step " << step << ": " << (int)( t / params.T * 100.0 ) << "%" << std::endl
					<< "  t:            " << t                                  << std::endl
					<< "  dt:           " << dt                                 << std::endl
					<< "  iter:         " << iter                               << std::endl
					<< "  res:          " << res                                << std::endl
					<< "  SOR time:     " << timeEndSOR - timeBeginSOR << 's'   << std::endl
					<< "  Step time:    " << timeEndStep - timeBeginStep << 's' << std::endl
					<< "  Elapsed time: " << tock() << 's'                      << std::endl;
			);
		}
	}

	DO_ON_RANK( PRINT_ON_RANK,
		std::cout << "[Profiling]: Total execution time: " << tock() << 's' << std::endl;
	);
}

int main( int argc, char* argv[] )
{
	assert( argc <= 2 );
	
	Communication::init( argc, argv );

	if( argc > 2 )
	{
#ifdef _CPP11
		std::vector<std::thread> threads( argc - 1 );
		for( int i = 1; i < argc; ++i )
		{
			threads[ i - 1 ] = std::thread( simulate, argv[ i ], i - 1 );
		}
		for( int i = 1; i < argc; ++i )
		{
			threads[ i - 1 ].join();
		}
#else
		for( int i = 1; i < argc; ++i )
		{
			simulate( argv[ i ], i - 1 );
		}
#endif
	}
	else if( argc == 2 )
	{
		simulate( argv[ 1 ], 0 );
	}
	else
	{
		simulate( NULL, 0 );
	}

	Communication::finalize();

	return 0;
}

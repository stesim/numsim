#include "IO.h"
#include "computation.h"
#include "solver.h"
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

using namespace mm;

void simulate( const char* filename, int instance )
{
	const char* paramFile = "inputvals";
	if( filename != NULL )
	{
		paramFile = filename;
	}

	std::cout << "[Initialization] Loading parameters from file '" << paramFile << "'." << std::endl;

	Params params = IO::readInputFile( paramFile );

	std::cout << "[Initialization] Parameter overview:" << std::endl;
	params.print( std::cout );

	Point h( params.domainSize.x / params.gridSize.x,
			params.domainSize.y / params.gridSize.y );

	GridFunction u( params.gridSize + MultiIndex( 1, 2 ) );
	GridFunction v( params.gridSize + MultiIndex( 2, 1 ) );
	GridFunction f( params.gridSize + MultiIndex( 1, 2 ) );
	GridFunction g( params.gridSize + MultiIndex( 2, 1 ) );
	GridFunction p( params.gridSize + MultiIndex( 2, 2 ) );
	GridFunction rhs( params.gridSize );

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
		//real dt = 0.001;

		Computation::setVelocityBoundary( u, v );
		// compute intermediate velocities
		Computation::computeMomentumEquations(
				f, g, u, v, h, dt, params.Re, params.alpha );
		Computation::setVelocityBoundary( f, g );
		Computation::computeRightHandSide( rhs, f, g, h, dt );

		real timeBeginSOR = tock();

		index_t iter = 0;
		real res = params.eps + 1.0;
		while( iter < params.maxIter && res > params.eps )
		{
			Computation::setPressureBoundary( p );
			Solver::SORCycle( p, rhs, h, params.omega );

			++iter;
			res = Solver::computeResidual( p, rhs, h );
		}

		real timeEndSOR = tock();

		Computation::computeNewVelocities( u, v, f, g, p, h, dt );

		real timeEndStep = tock();

		t += dt;
		++step;

		tOut += dt;
		if( tOut > params.deltaVec || t >= params.T )
		{
			tOut -= params.deltaVec;
			IO::writeRawOutput( params.gridSize, u, v,
					p, h, instance, step );
			IO::writeVTKFile( params.gridSize, u, v,
					p, h, params, instance, step );

			std::cout << "[Progress] Instance " << instance << ", step " << step << ": " << (int)( t / params.T * 100.0 ) << "%" << std::endl
				<< "  t:            " << t                                  << std::endl
				<< "  dt:           " << dt                                 << std::endl
				<< "  iter:         " << iter                               << std::endl
				<< "  res:          " << res                                << std::endl
				<< "  SOR time:     " << timeEndSOR - timeBeginSOR << 's'   << std::endl
				<< "  Step time:    " << timeEndStep - timeBeginStep << 's' << std::endl
				<< "  Elapsed time: " << tock() << 's'                      << std::endl;
		}
	}

	std::cout << "[Profiling]: Total execution time: " << tock() << 's' << std::endl;
}

int main( int argc, char* argv[] )
{
	if( argc > 1 )
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
	else
	{
		simulate( NULL, 0 );
	}
	return 0;
}

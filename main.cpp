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
#include <thread>
#include <chrono>

inline std::chrono::high_resolution_clock::time_point getTime()
{
	return std::chrono::high_resolution_clock::now();
}

inline double getElapsedTime( std::chrono::high_resolution_clock::time_point t )
{
	return ( std::chrono::duration_cast<std::chrono::duration<real>>(
				getTime() - t ) ).count();
}

#define DO_ON_RANK(rnk,todo) if( Communication::linearRank() == rnk ) { todo }
//#define DO_ON_RANK(rnk,todo) todo

#define PRINT_ON_RANK 0

using namespace mm;

void simulate( Params params, int instance )
{
	MultiIndex commSize = Communication::size();
	MultiIndex commRank = Communication::rank();

	bool leftBoundary = ( commRank.x == 0 );
	bool rightBoundary = ( commRank.x == commSize.x - 1 );
	bool bottomBoundary = ( commRank.y == 0 );
	bool topBoundary = ( commRank.y == commSize.y - 1 );

	DO_ON_RANK( PRINT_ON_RANK,
		std::cout << "[Initialization] Parameter overview:" << std::endl;
		std::cout << params;
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

	GridFunction psi( localGridSize + MultiIndex::ONE );
	GridFunction zeta( localGridSize - MultiIndex::ONE );

	// set boundary values of psi once
	if( leftBoundary )
	{
		set( psi, MultiIndex::ZERO, MultiIndex( 1, psi.size().y ),
				constant( 0.0 ) );
	}

	u = constant( params.initialVelocity.x );
	v = constant( params.initialVelocity.y );
	p = constant( params.initialPressure );

	f = constant( params.initialVelocity.x );
	g = constant( params.initialVelocity.y );

	auto startTime = getTime();

	Computation::setVelocityBoundary( u, v,
			leftBoundary, topBoundary, rightBoundary, bottomBoundary );

	real t = 0.0;
	index_t step = 0;
	real tOut = 0.0;
	while( t < params.T )
	{
		real timeBeginStep = getElapsedTime( startTime );

		real dt = Computation::computeTimeStep(
				u, v, h, params.Re, params.tau );
		
		// compute intermediate velocities
		Computation::computeMomentumEquations(
				f, g, u, v, h, dt, params.Re, params.alpha );
		Computation::setVelocityBoundary( f, g,
				leftBoundary, topBoundary, rightBoundary, bottomBoundary );
		Communication::exchangeFunctionBoundary( f, MultiIndex( 2, 1 ) );
		Communication::exchangeFunctionBoundary( g, MultiIndex( 1, 2 ) );

		Computation::computeRightHandSide( rhs, f, g, h, dt );

		real timeBeginSOR = getElapsedTime( startTime );

		index_t iter = 0;
		real res = params.eps + 1.0;
		while( iter < params.maxIter && res > params.eps )
		{
			Computation::setPressureBoundary( p,
					leftBoundary, topBoundary, rightBoundary, bottomBoundary );

			Solver::SORSubcycle( p, rhs, h, params.omega, true );
			Communication::exchangeFunctionBoundary( p, MultiIndex::ONE );
			Solver::SORSubcycle( p, rhs, h, params.omega, false );
			Communication::exchangeFunctionBoundary( p, MultiIndex::ONE );

			res = Solver::computeResidual( p, rhs, params.gridSize, h );

			++iter;
		}

		real timeEndSOR = getElapsedTime( startTime );

		// compute new velocities
		Computation::computeNewVelocities( u, v, f, g, p, h, dt );
		Communication::exchangeFunctionBoundary( u, MultiIndex( 2, 1 ) );
		Communication::exchangeFunctionBoundary( v, MultiIndex( 1, 2 ) );

		Computation::setVelocityBoundary( u, v,
				leftBoundary, topBoundary, rightBoundary, bottomBoundary );

		real timeEndStep = getElapsedTime( startTime );

		t += dt;
		++step;

		tOut += dt;
		if( step <= 1 || tOut > params.deltaVec || t >= params.T )
		{
			tOut -= params.deltaVec;

			// compute velocity potential and vorticity
			Communication::receiveBoundary(
					psi, Communication::BOUNDARY_LEFT, false );
			Computation::computeVelocityPotential( psi, u, v, h );
			Communication::sendBoundary(
					psi, Communication::BOUNDARY_RIGHT, 0, false );

			Computation::computeVorticity( zeta, u, v, h );

			IO::writeRawOutput( u, v, p, psi, zeta,
					h, localGridOffset, step, commRank );
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
					<< "  SOR time:     " << timeEndSOR - timeBeginSOR   << 's' << std::endl
					<< "  Step time:    " << timeEndStep - timeBeginStep << 's' << std::endl
					<< "  Elapsed time: " << getElapsedTime( startTime ) << 's' << std::endl;
			);
		}
	}

	DO_ON_RANK( PRINT_ON_RANK,
		std::cout << "[Profiling]: Total execution time: " << getElapsedTime( startTime ) << 's' << std::endl;
	);
}

Params defaultParams()
{
	Params p;
	p.domainSize = Point( 1.0, 1.0 );
	p.gridSize = MultiIndex( 128, 128 );
	p.T = 16.5;
	p.dt = 0;
	p.tau = 0.5;
	p.deltaVec = 1.4;
	p.maxIter = 100;
	p.eps = 0.001;
	p.omega = 1.7;
	p.alpha = 0.9;
	p.Re = 1000;
	p.initialVelocity = Point( 0.0, 0.0 );
	p.initialPressure = 0;

	return p;
}

int main( int argc, char* argv[] )
{
	Communication::init( argc, argv );

	std::vector<Params> params;

	--argc;
	++argv;

	do
	{
		Params p = defaultParams();
		int endIdx = p.parseCmdArgs( argc, argv );

		params.push_back( p );

		argc -= endIdx;
		argv += endIdx;
	} while( argc > 0 );

#ifndef COMM_MPI
	if( params.size() > 1 )
	{
		std::vector<std::thread> threads( params.size() );
		for( std::size_t i = 0; i < params.size(); ++i )
		{
			threads[ i ] = std::thread( simulate, params[ i ], i );
		}
		for( std::size_t i = 0; i < params.size(); ++i )
		{
			threads[ i ].join();
		}
	}
	else
	{
		simulate( params[ 0 ], 0 );
	}
#else
	simulate( params[ 0 ], 0 );
#endif

	Communication::finalize();

	return 0;
}

#include "IO.h"
#include "computation.h"
#include "solver.h"
#include <iostream>

int main( int argc, char* argv[] )
{
	const char* paramFile = "inputvals";
	if( argc > 1 )
	{
		paramFile = argv[ 1 ];
	}

	std::cout << "[Initialization] Loading parameters from file '" << paramFile << "'." << std::endl;

	Params params = IO::readInputFile( paramFile );

	std::cout << "[Initialization] Parameter overview:" << std::endl
		<< "  Domain size:           " << params.domainSize.x << " x " << params.domainSize.y << std::endl
		<< "  Grid size:             " << params.gridSize.x << " x " << params.gridSize.y << std::endl
		<< "  T:                     " << params.T << std::endl
		<< "  dt:                    " << params.dt << std::endl
		<< "  tau:                   " << params.tau << std::endl
		<< "  Output time step:      " << params.deltaVec << std::endl
		<< "  Max SOR iterations:    " << params.maxIter << std::endl
		<< "  SOR error tolerance:   " << params.eps << std::endl
		<< "  Relaxation factor:     " << params.omega << std::endl
		<< "  Upwinding factor:      " << params.alpha << std::endl
		<< "  Re:                    " << params.Re << std::endl
		<< "  Initial velocity:      (" << params.initialVelocity.x << "," << params.initialVelocity.y << ")" << std::endl
		<< "  Initial pressure:      " << params.initialPressure << std::endl;

	Point h( params.domainSize.x / params.gridSize.x,
			params.domainSize.y / params.gridSize.y );

	Computation computation( params );
	Solver solver( params.gridSize, h );

	real t = 0.0;
	index_t step = 0;
	real tOut = 0.0;
	while( t < params.T )
	{
		real dt = computation.computeTimeStep();

		computation.setVelocityBoundary( computation.u, computation.v );
		computation.computeMomentumEquations( dt );
		computation.setVelocityBoundary( computation.f, computation.g );
		computation.computeRightHandSide( dt );

		index_t iter = 0;
		real res = params.eps + 1.0;
		while( iter < params.maxIter && res > params.eps )
		{
			computation.setPressureBoundary();
			solver.SORCycle( computation.p, computation.rhs, params.omega );

			++iter;
			res = solver.computeResidual( computation.p, computation.rhs );
		}

		computation.computeNewVelocities( dt );

		t += dt;
		++step;

		tOut += dt;
		if( tOut > params.deltaVec )
		{
			tOut -= params.deltaVec;
			IO::writeRawOutput( params.gridSize, computation.u, computation.v,
					computation.p, h, step );

			std::cout << "[Progress] Step " << step << ": " << (int)( t / params.T * 100.0 ) << "%" << std::endl
				<< "  t:    " << t << std::endl
				<< "  iter: " << iter << std::endl
				<< "  res:  " << res << std::endl;
		}
	}

	return 0;
}

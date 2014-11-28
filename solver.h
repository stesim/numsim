#ifndef SOLVER_H
#define SOLVER_H

#include "typedef.h"

class Solver
{
public:
	//! Compute SOR residual
	//! \param p Pressure function
	//! \param rhs Right hand side
	static real computeResidual( const GridFunction& p,
			const GridFunction& rhs,
			const Point& h );

	//! Compute one SOR cycle
	//! \param p Pressure function
	//! \param rhs Right hand side
	//! \param omega Relaxation factor
	static void SORCycle( GridFunction& p,
			const GridFunction& rhs,
			const Point& h,
			real omega );
};

#endif

#ifndef SOLVER_H
#define SOLVER_H

#include "gridfunction.h"
#include "stencil.h"

class Solver
{
public:
	Solver( const MultiIndex& dim,
			const Point& h );

	real computeResidual( const GridFunction& p,
			const GridFunction& rhs );

	void SORCycle( GridFunction& p,
			const GridFunction& rhs,
			real omega );

private:
	MultiIndex   m_miSize;
	Point        m_ptH;
	GridFunction m_funcTemp;
	Stencil      m_stSOR;
};

#endif

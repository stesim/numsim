#ifndef COMPUTATION_H
#define COMPUTATION_H

#include "typedef.h"
#include "gridfunction.h"
#include "stencil.h"
#include "params.h"

class Computation
{
public:
	Computation( const Params& params );

	real computeTimeStep();

	void computeNewVelocities( real dt );

	void computeMomentumEquations( real dt );

	void computeRightHandSide( real dt );

	void setVelocityBoundary( GridFunction& u, GridFunction& v );

	void setPressureBoundary();

private:
	Params m_Params;
	Point  h;

public:
	GridFunction u;
	GridFunction v;
	GridFunction f;
	GridFunction g;
	GridFunction p;
	GridFunction rhs;

private:
	GridFunction m_funcTemp1;
	GridFunction m_funcTemp2;
	GridFunction m_funcTemp3;

	Stencil m_stDxBw;
	Stencil m_stDyBw;
	Stencil m_stDxx;
	Stencil m_stDyy;
	Stencil m_stInterpx;
	Stencil m_stInterpy;
};

#endif

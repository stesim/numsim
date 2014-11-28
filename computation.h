#ifndef COMPUTATION_H
#define COMPUTATION_H

#include "typedef.h"
#include "params.h"

class Computation
{
public:
	//! Compute time step
	static real computeTimeStep( const GridFunction& u,
			const GridFunction& v,
			const Point& h,
			real Re,
			real tau );

	//! Compute velocities from intermediate velocities
	//! \param dt Time step
	static void computeNewVelocities( GridFunction& u,
			GridFunction& v,
			const GridFunction& f,
			const GridFunction& g,
			const GridFunction& p,
			const Point& h,
			real dt );

	//! Compute intermediate velocities using the momentum equations
	//! \param dt Time step
	static void computeMomentumEquations( GridFunction& f,
			GridFunction& g,
			const GridFunction& u,
			const GridFunction& v,
			const Point& h,
			real dt,
			real Re,
			real alpha );

	//! Compute the right hand side of the linear system for the pressure
	//! \param dt Time step
	static void computeRightHandSide( GridFunction& rhs,
			const GridFunction& f,
			const GridFunction& g,
			const Point& h,
			real dt );

	//! Set (intermediate) velocity boundary values
	//! \param u First velocity component
	//! \param u Second velocity component
	static void setVelocityBoundary( GridFunction& u, GridFunction& v );

	//! Set pressure boundary values
	static void setPressureBoundary( GridFunction& p );
};

#endif

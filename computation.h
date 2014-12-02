#ifndef COMPUTATION_H
#define COMPUTATION_H

#include "typedef.h"

class Computation
{
public:
	//! Compute time step
	//! \param u First velocity component
	//! \param v Second velocity component
	//! \param h Spacial step size
	//! \param Re Reynolds number
	static real computeTimeStep( const GridFunction& u,
			const GridFunction& v,
			const Point& h,
			real Re,
			real tau );

	//! Compute velocities from intermediate velocities
	//! \param u First velocity component
	//! \param v Second velocity component
	//! \param f First intermediate velocity component
	//! \param g Second intermediate velocity component
	//! \param p Pressure function
	//! \param h Spacial step size
	//! \param dt Time step
	static void computeNewVelocities( GridFunction& u,
			GridFunction& v,
			const GridFunction& f,
			const GridFunction& g,
			const GridFunction& p,
			const Point& h,
			real dt );

	//! Compute intermediate velocities using the momentum equations
	//! \param f First intermediate velocity component
	//! \param g Second intermediate velocity component
	//! \param u First velocity component
	//! \param v Second velocity component
	//! \param h Spacial step size
	//! \param dt Time step
	//! \param Re Reynolds number
	static void computeMomentumEquations( GridFunction& f,
			GridFunction& g,
			const GridFunction& u,
			const GridFunction& v,
			const Point& h,
			real dt,
			real Re,
			real alpha );

	//! Compute the right hand side of the linear system for the pressure
	//! \param f First intermediate velocity component
	//! \param g Second intermediate velocity component
	//! \param h Spacial step size
	//! \param dt Time step
	static void computeRightHandSide( GridFunction& rhs,
			const GridFunction& f,
			const GridFunction& g,
			const Point& h,
			real dt );

	//! Compute the velocity potential function
	//! \param psi Velocity potential
	//! \param u First velocity component
	//! \param v Second velocity component
	//! \param h Spacial step size
	static void computeVelocityPotential( GridFunction& psi,
			const GridFunction& u,
			const GridFunction& v,
			const Point& h );

	//! Compute the vorticity function
	//! \param zeta Vorticity
	//! \param u First velocity component
	//! \param v Second velocity component
	//! \param h Spacial step size
	static void computeVorticity( GridFunction& zeta,
			const GridFunction& u,
			const GridFunction& v,
			const Point& h );

	//! Set (intermediate) velocity boundary values
	//! \param u First velocity component
	//! \param v Second velocity component
	//! \param left Determines if left boundary is set
	//! \param top Determines if top boundary is set
	//! \param right Determines if right boundary is set
	//! \param bottom Determines if bottom boundary is set
	static void setVelocityBoundary( GridFunction& u, GridFunction& v,
			bool left, bool top, bool right, bool bottom );

	//! Set pressure boundary values
	//! \param p Pressure function
	//! \param left Determines if left boundary is set
	//! \param top Determines if top boundary is set
	//! \param right Determines if right boundary is set
	//! \param bottom Determines if bottom boundary is set
	static void setPressureBoundary( GridFunction& p,
			bool left, bool top, bool right, bool bottom );
};

#endif

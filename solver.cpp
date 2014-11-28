#include "solver.h"
#include <cmath>
#include <cassert>
#include "metamath/metamath.h"
#include "metamath/mmutils.h"

using namespace mm;

real Solver::computeResidual( const GridFunction& p,
		const GridFunction& rhs,
		const Point& h )
{
	assert( p.size() - MultiIndex( 2, 2 ) == m_miSize &&
			rhs.size() == m_miSize );

	real sum = mm::sum( sqr( utils::diffXX_YY( eval<+1,+1>( p ), h ) - rhs ),
			MultiIndex::ZERO, rhs.size() );
	
	return std::sqrt( sum / ( rhs.size().y * rhs.size().x ) );
}

void Solver::SORCycle( GridFunction& p,
		const GridFunction& rhs,
		const Point& h,
		real omega )
{
	assert( p.size() - MultiIndex( 2, 2 ) == m_miSize &&
			rhs.size() == m_miSize );

	real dx2 = h.x * h.x;
	real dy2 = h.y * h.y;

	// "HACK": use p as source and image function
	// -> both old and new values are used in the calculations (Gauss-Seidel)
	set( p, MultiIndex::ONE, p.size() - MultiIndex::ONE,
			( 1.0 - omega ) % p + 0.5 * omega * dx2 * dy2 / ( dx2 + dy2 ) % (
				1 / dx2 % ( eval<-1,0>( p ) + eval<+1,0>( p ) ) +  1 / dy2 % ( eval<0,-1>( p ) + eval<0,+1>( p ) )
				- eval<-1,-1>( rhs ) ) );
}

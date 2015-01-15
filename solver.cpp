#include "solver.h"
#include <cmath>
#include <cassert>
#include "metamath/metamath.h"
#include "metamath/mmutils.h"
#include "communication.h"

#include <cstdio>

using namespace mm;

real Solver::computeResidual( const GridFunction& p,
		const GridFunction& rhs,
		const MultiIndex& gridSize,
		const Point& h )
{
	assert( p.size() - MultiIndex( 2, 2 ) == m_miSize &&
			rhs.size() == m_miSize );

	real sum_local =
		mm::sum( sqr( utils::diffXX_YY( eval<+1,+1>( p ), h ) - rhs ),
				MultiIndex::ZERO, rhs.size() );

	real sum = Communication::sum( sum_local );
	return std::sqrt( sum / ( gridSize.y * gridSize.x ) );
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

	// "HACK": use p as source _and_ image function
	// -> both old and new values are used in the calculations (Gauss-Seidel)
	set( p, MultiIndex::ONE, p.size() - MultiIndex::ONE,
			( 1.0 - omega ) * p + 0.5 * omega * dx2 * dy2 / ( dx2 + dy2 ) * (
				1 / dx2 * ( eval<-1,0>( p ) + eval<+1,0>( p ) ) +  1 / dy2 * ( eval<0,-1>( p ) + eval<0,+1>( p ) )
				- eval<-1,-1>( rhs ) ) );
}

void Solver::SORSubcycle( GridFunction& p,
		const GridFunction& rhs,
		const Point& h,
		real omega,
		bool color )
{
	assert( p.size() - MultiIndex( 2, 2 ) == m_miSize &&
			rhs.size() == m_miSize );

	real dx2 = h.x * h.x;
	real dy2 = h.y * h.y;

	utils::setCheckered( p, MultiIndex::ONE, p.size() - MultiIndex::ONE, color,
			( 1.0 - omega ) * p + 0.5 * omega * dx2 * dy2 / ( dx2 + dy2 ) * (
				1 / dx2 * ( eval<-1,0>( p ) + eval<+1,0>( p ) ) +  1 / dy2 * ( eval<0,-1>( p ) + eval<0,+1>( p ) )
				- eval<-1,-1>( rhs ) ) );
}

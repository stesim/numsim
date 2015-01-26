#include "solver.h"
#include <cmath>
#include <cassert>
#include "metamath/metamath.h"
#include "metamath/mmutils.h"
#include "communication.h"

#include <cstdio>

using namespace mm;

template<typename Top, typename Tbegin, typename Tend, typename Tmask>
op_dtype<Top> sumMasked( const Top& op, const Tbegin& begin, const Tend& end,
		const Tmask& mask )
{
	op_dtype<Top> sum = 0;
	for( int j = begin.y; j < end.y; ++j )
	{
		for( int i = begin.x; i < end.x; ++i )
		{
			if( mask( i, j ) )
			{
				sum += op( i, j );
			}
		}
	}
	return sum;
}

template<typename Tfunc, typename Tbegin, typename Tend,
	typename Tmask, typename Top>
void setCheckeredMasked( Tfunc& func, const Tbegin& begin,
		const Tend& end, const Tmask& mask, int color, const Top& op )
{
	for( int j = begin.y; j < end.y; ++j )
	{
		for( int i = begin.x + ( j + color ) % 2; i < end.x; i += 2 )
		{
			if( mask( i, j ) )
			{
				func( i, j ) = op( i, j );
			}
		}
	}
}

real Solver::computeResidual( const GridFunction& p,
		const GridFunction& rhs,
		const MaskFunction& pMask,
		const MultiIndex& gridSize,
		const Point& h )
{
	assert( p.size() - MultiIndex( 2, 2 ) == m_miSize &&
			rhs.size() == m_miSize );

	real sum_local =
		sumMasked( sqr( utils::diffXX_YY( eval<+1,+1>( p ), h ) - rhs ),
				MultiIndex::ZERO, rhs.size(), eval<+1,+1>( pMask ) );

	real sum = Communication::sum( sum_local );
	return std::sqrt( sum / ( gridSize.y * gridSize.x ) );
}

void Solver::SORCycle( GridFunction& p,
		const GridFunction& rhs,
		const MaskFunction& pMask,
		const Point& h,
		real omega )
{
	assert( p.size() - MultiIndex( 2, 2 ) == m_miSize &&
			rhs.size() == m_miSize );

	real dx2 = h.x * h.x;
	real dy2 = h.y * h.y;

	// "HACK": use p as source _and_ image function
	// -> both old and new values are used in the calculations (Gauss-Seidel)
	utils::setMasked( p, MultiIndex::ONE, p.size() - MultiIndex::ONE, pMask,
			( 1.0 - omega ) * p + 0.5 * omega * dx2 * dy2 / ( dx2 + dy2 ) * (
				1 / dx2 * ( eval<-1,0>( p ) + eval<+1,0>( p ) ) +  1 / dy2 * ( eval<0,-1>( p ) + eval<0,+1>( p ) )
				- eval<-1,-1>( rhs ) ) );
}

void Solver::SORSubcycle( GridFunction& p,
		const GridFunction& rhs,
		const MaskFunction& pMask,
		const Point& h,
		real omega,
		bool color )
{
	assert( p.size() - MultiIndex( 2, 2 ) == m_miSize &&
			rhs.size() == m_miSize );

	real dx2 = h.x * h.x;
	real dy2 = h.y * h.y;

	setCheckeredMasked( p, MultiIndex::ONE, p.size() - MultiIndex::ONE, pMask, color,
			( 1.0 - omega ) * p + 0.5 * omega * dx2 * dy2 / ( dx2 + dy2 ) * (
				1 / dx2 * ( eval<-1,0>( p ) + eval<+1,0>( p ) ) +  1 / dy2 * ( eval<0,-1>( p ) + eval<0,+1>( p ) )
				- eval<-1,-1>( rhs ) ) );
}

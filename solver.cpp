#include "solver.h"
#include <cmath>

Solver::Solver( const MultiIndex& dim,
		const Point& h )
	: m_miSize( dim ),
	m_ptH( h ),
	m_funcTemp( dim ),
	m_stSOR( Stencil::sor( h ) )
{
}

real Solver::computeResidual( const GridFunction& p,
		const GridFunction& rhs )
{
	assert( p.getSize() - MultiIndex( 2, 2 ) == m_miSize &&
			rhs.getSize() == m_miSize );

	real dx2 = m_ptH.x * m_ptH.x;
	real dy2 = m_ptH.y * m_ptH.y;
	real sum = 0.0;
	for( index_t j = 1; j < p.getSize().y - 1; ++j )
	{
		for( index_t i = 1; i < p.getSize().x - 1; ++i )
		{
			real res_loc = ( p( i + 1, j ) - 2.0 * p( i, j ) + p( i - 1, j ) ) / dx2 +
				( p( i, j + 1 ) - 2.0 * p( i, j ) + p( i, j - 1 ) ) / dy2 - rhs( i - 1, j - 1 );
			sum += res_loc * res_loc;
		}
	}
	sum = std::sqrt( sum / ( ( p.getSize().y - 1 ) * ( p.getSize().x - 1 ) ) );

	return sum;
}

void Solver::SORCycle( GridFunction& p,
		const GridFunction& rhs,
		real omega )
{
	assert( p.getSize() - MultiIndex( 2, 2 ) == m_miSize &&
			rhs.getSize() == m_miSize );

	real dx2 = m_ptH.x * m_ptH.x;
	real dy2 = m_ptH.y * m_ptH.y;

	m_stSOR.apply( p, MultiIndex::ONE, p.getSize() - MultiIndex::ONE,
			m_funcTemp, MultiIndex::ZERO, m_funcTemp.getSize() );
	m_funcTemp.addScaledCopy( MultiIndex::ZERO, m_funcTemp.getSize(),
			rhs, -1.0 );
	m_funcTemp.scale( MultiIndex::ZERO, m_funcTemp.getSize(),
			omega * dx2 * dy2 / ( 2.0 * ( dx2 + dy2 ) ) );

	p.scale( MultiIndex::ONE, p.getSize() - MultiIndex::ONE, 1.0 - omega );
	p.addOffsetCopy( MultiIndex::ONE, p.getSize() - MultiIndex::ONE,
			m_funcTemp, MultiIndex( -1, -1 ) );
}

#include "gridfunction.h"
#include <cmath>
#include <cstdio>

#ifndef NULL
#define NULL 0
#endif

GridFunction::GridFunction()
	: m_pGridValues( NULL ), m_miSize( 0, 0 )
{
}

GridFunction::GridFunction( const MultiIndex& size )
	: m_miSize( size )
{
	allocMem();
}

GridFunction::GridFunction( index_t dimX, index_t dimY )
	: m_miSize( dimX, dimY )
{
	allocMem();
}

GridFunction::~GridFunction()
{
	freeMem();
}

void GridFunction::allocMem()
{
	index_t totalSize = m_miSize.x * m_miSize.y;
	if( totalSize > 0 )
	{
		m_pGridValues = new real[ totalSize ];
	}
}

void GridFunction::freeMem()
{
	if( m_pGridValues != NULL )
	{
		delete[] m_pGridValues;
	}
}

const MultiIndex& GridFunction::getSize() const
{
	return m_miSize;
}

real& GridFunction::operator()( index_t i, index_t j )
{
	return m_pGridValues[ j * m_miSize.x + i ];
}

const real& GridFunction::operator()( index_t i, index_t j ) const
{
	return m_pGridValues[ j * m_miSize.x + i ];
}

void GridFunction::set( const MultiIndex& begin,
		const MultiIndex& end,
		real value )
{
	for( index_t j = begin.y; j < end.y; ++j )
	{
		for( index_t i = begin.x; i < end.x; ++i )
		{
			m_pGridValues[ j * m_miSize.x + i ] = value;
		}
	}
}

void GridFunction::scale( const MultiIndex& begin,
		const MultiIndex& end,
		real factor )
{
	for( index_t j = begin.y; j < end.y; ++j )
	{
		for( index_t i = begin.x; i < end.x; ++i )
		{
			m_pGridValues[ j * m_miSize.x + i ] *= factor;
		}
	}
}

void GridFunction::offsetScale( const MultiIndex& begin,
		const MultiIndex& end,
		const MultiIndex& offset,
		real factor )
{
	// negative offset = evil = not allowed!
	assert( offset.x >= 0 && offset.x >= 0 );
	// check for out-of-bound indices
	assert( end.x + offset.x <= m_miSize.x || end.y + offset.y <= m_miSize.y );

	index_t linearOffset = offset.y * m_miSize.x + offset.x;
	for( index_t j = begin.y; j < end.y; ++j )
	{
		for( index_t i = begin.x; i < end.x; ++i )
		{
			index_t idx = j * m_miSize.x + i;
			m_pGridValues[ idx ] =
				factor * m_pGridValues[ idx + linearOffset ];
		}
	}
}

void GridFunction::copy( const MultiIndex& begin,
		const MultiIndex& end,
		const GridFunction& source )
{
	for( index_t j = begin.y; j < end.y; ++j )
	{
		for( index_t i = begin.x; i < end.x; ++i )
		{
			m_pGridValues[ j * m_miSize.x + i ] =
				source.m_pGridValues[ j * source.m_miSize.x + i ];
		}
	}
}

void GridFunction::scaledCopy( const MultiIndex& begin,
		const MultiIndex& end,
		const GridFunction& source,
		real factor )
{
	for( index_t j = begin.y; j < end.y; ++j )
	{
		for( index_t i = begin.x; i < end.x; ++i )
		{
			m_pGridValues[ j * m_miSize.x + i ] =
				factor * source.m_pGridValues[ j * source.m_miSize.x + i ];
		}
	}
}

void GridFunction::scaledOffsetCopy( const MultiIndex& begin,
		const MultiIndex& end,
		const GridFunction& source,
		const MultiIndex& offset,
		real factor )
{
	index_t linearOffset = offset.y * source.m_miSize.x + offset.x;
	for( index_t j = begin.y; j < end.y; ++j )
	{
		for( index_t i = begin.x; i < end.x; ++i )
		{
			m_pGridValues[ j * m_miSize.x + i ] =
				factor * source.m_pGridValues[
					j * source.m_miSize.x + i + linearOffset ];
		}
	}
}

void GridFunction::affineOffsetCopy( const MultiIndex& begin,
		const MultiIndex& end,
		const GridFunction& source,
		const MultiIndex& offset,
		real factor,
		real constant )
{
	index_t linearOffset = offset.y * source.m_miSize.x + offset.x;
	for( index_t j = begin.y; j < end.y; ++j )
	{
		for( index_t i = begin.x; i < end.x; ++i )
		{
			m_pGridValues[ j * m_miSize.x + i ] = constant +
				factor * source.m_pGridValues[
					j * source.m_miSize.x + i + linearOffset];
		}
	}
}

void GridFunction::addOffsetCopy( const MultiIndex& begin,
		const MultiIndex& end,
		const GridFunction& source,
		const MultiIndex& offset )
{
	index_t linearOffset = offset.y * source.m_miSize.x + offset.x;
	for( index_t j = begin.y; j < end.y; ++j )
	{
		for( index_t i = begin.x; i < end.x; ++i )
		{
			m_pGridValues[ j * m_miSize.x + i ] +=
				source.m_pGridValues[ j * source.m_miSize.x + i + linearOffset ];
		}
	}
}

void GridFunction::addScaledCopy( const MultiIndex& begin,
		const MultiIndex& end,
		const GridFunction& source,
		real factor )
{
	for( index_t j = begin.y; j < end.y; ++j )
	{
		for( index_t i = begin.x; i < end.x; ++i )
		{
			m_pGridValues[ j * m_miSize.x + i ] +=
				factor * source.m_pGridValues[ j * source.m_miSize.x + i ];
		}
	}
}

real GridFunction::getMaxValue( const MultiIndex& begin,
		const MultiIndex& end ) const
{
	assert( end.x - begin.x > 0 && end.y - begin.y > 0 );

	real max = m_pGridValues[ begin.y * m_miSize.x + begin.x ];
	for( index_t j = begin.y; j < end.y; ++j )
	{
		for( index_t i = begin.x; i < end.x; ++i )
		{
			real val = std::abs( m_pGridValues[ j * m_miSize.x + i ] );
			if( val > max )
			{
				max = val;
			}
		}
	}

	return max;
}

void GridFunction::multiply( const MultiIndex& begin,
		const MultiIndex& end,
		const GridFunction& source )
{
	for( index_t j = begin.y; j < end.y; ++j )
	{
		for( index_t i = begin.x; i < end.x; ++i )
		{
			m_pGridValues[ j * m_miSize.x + i ] *=
				source.m_pGridValues[ j * source.m_miSize.x + i ];
		}
	}
}

void GridFunction::abs( const MultiIndex& begin,
		const MultiIndex& end )
{
	for( index_t j = begin.y; j < end.y; ++j )
	{
		for( index_t i = begin.x; i < end.x; ++i )
		{
			index_t idx = j * m_miSize.x + i;
			if( m_pGridValues[ idx ] < 0.0 )
			{
				m_pGridValues[ idx ] = -m_pGridValues[ idx ];
			}
		}
	}
}

void GridFunction::print( const MultiIndex& begin, const MultiIndex& end ) const
{
	for( index_t j = begin.y; j < end.y; ++j )
	{
		for( index_t i = begin.x; i < end.x; ++i )
		{
			printf( "%12.4e", m_pGridValues[ j * m_miSize.x + i ] );
		}
		printf( "\n" );
	}
}

#include "stencil.h"

#ifdef DEBUG
#define MY_ASSERT(x) x
#else
#define MY_ASSERT(x)
#endif

Stencil Stencil::dx( real h )
{
	real _2h = 2 * h;

	Stencil s;
	s.m_Values[ 0 ][ 0 ] =  0.0;
	s.m_Values[ 0 ][ 1 ] =  0.0;
	s.m_Values[ 0 ][ 2 ] =  0.0;

	s.m_Values[ 1 ][ 0 ] = -1.0 / _2h;
	s.m_Values[ 1 ][ 1 ] =  0.0;
	s.m_Values[ 1 ][ 2 ] =  1.0 / _2h;

	s.m_Values[ 2 ][ 0 ] =  0.0;
	s.m_Values[ 2 ][ 1 ] =  0.0;
	s.m_Values[ 2 ][ 2 ] =  0.0;

	return s;
}

Stencil Stencil::dy( real h )
{
	real _2h = 2 * h;

	Stencil s;
	s.m_Values[ 0 ][ 0 ] =  0.0;
	s.m_Values[ 0 ][ 1 ] = -1.0 / _2h;
	s.m_Values[ 0 ][ 2 ] =  0.0;

	s.m_Values[ 1 ][ 0 ] =  0.0;
	s.m_Values[ 1 ][ 1 ] =  0.0;
	s.m_Values[ 1 ][ 2 ] =  0.0;

	s.m_Values[ 2 ][ 0 ] =  0.0;
	s.m_Values[ 2 ][ 1 ] =  1.0 / _2h;
	s.m_Values[ 2 ][ 2 ] =  0.0;

	return s;
}

Stencil Stencil::dx_bw( real h )
{
	Stencil s;
	s.m_Values[ 0 ][ 0 ] =  0.0;
	s.m_Values[ 0 ][ 1 ] =  0.0;
	s.m_Values[ 0 ][ 2 ] =  0.0;

	s.m_Values[ 1 ][ 0 ] = -1.0 / h;
	s.m_Values[ 1 ][ 1 ] =  1.0 / h;
	s.m_Values[ 1 ][ 2 ] =  0.0;

	s.m_Values[ 2 ][ 0 ] =  0.0;
	s.m_Values[ 2 ][ 1 ] =  0.0;
	s.m_Values[ 2 ][ 2 ] =  0.0;

	return s;
}

Stencil Stencil::dy_bw( real h )
{
	Stencil s;
	s.m_Values[ 0 ][ 0 ] =  0.0;
	s.m_Values[ 0 ][ 1 ] = -1.0 / h;
	s.m_Values[ 0 ][ 2 ] =  0.0;

	s.m_Values[ 1 ][ 0 ] =  0.0;
	s.m_Values[ 1 ][ 1 ] =  1.0 / h;
	s.m_Values[ 1 ][ 2 ] =  0.0;

	s.m_Values[ 2 ][ 0 ] =  0.0;
	s.m_Values[ 2 ][ 1 ] =  0.0;
	s.m_Values[ 2 ][ 2 ] =  0.0;

	return s;
}

Stencil Stencil::dx_fw( real h )
{
	Stencil s;
	s.m_Values[ 0 ][ 0 ] =  0.0;
	s.m_Values[ 0 ][ 1 ] =  0.0;
	s.m_Values[ 0 ][ 2 ] =  0.0;

	s.m_Values[ 1 ][ 0 ] =  0.0;
	s.m_Values[ 1 ][ 1 ] = -1.0 / h;
	s.m_Values[ 1 ][ 2 ] =  1.0 / h;

	s.m_Values[ 2 ][ 0 ] =  0.0;
	s.m_Values[ 2 ][ 1 ] =  0.0;
	s.m_Values[ 2 ][ 2 ] =  0.0;

	return s;
}

Stencil Stencil::dy_fw( real h )
{
	Stencil s;
	s.m_Values[ 0 ][ 0 ] =  0.0;
	s.m_Values[ 0 ][ 1 ] =  0.0;
	s.m_Values[ 0 ][ 2 ] =  0.0;

	s.m_Values[ 1 ][ 0 ] =  0.0;
	s.m_Values[ 1 ][ 1 ] = -1.0 / h;
	s.m_Values[ 1 ][ 2 ] =  0.0;

	s.m_Values[ 2 ][ 0 ] =  0.0;
	s.m_Values[ 2 ][ 1 ] =  1.0 / h;
	s.m_Values[ 2 ][ 2 ] =  0.0;

	return s;
}

Stencil Stencil::dxx( real h )
{
	real h2 = h * h;

	Stencil s;
	s.m_Values[ 0 ][ 0 ] =  0.0;
	s.m_Values[ 0 ][ 1 ] =  0.0;
	s.m_Values[ 0 ][ 2 ] =  0.0;

	s.m_Values[ 1 ][ 0 ] =  1.0 / h2;
	s.m_Values[ 1 ][ 1 ] = -2.0 / h2;
	s.m_Values[ 1 ][ 2 ] =  1.0 / h2;

	s.m_Values[ 2 ][ 0 ] =  0.0;
	s.m_Values[ 2 ][ 1 ] =  0.0;
	s.m_Values[ 2 ][ 2 ] =  0.0;

	return s;
}

Stencil Stencil::dyy( real h )
{
	real h2 = h * h;

	Stencil s;
	s.m_Values[ 0 ][ 0 ] =  0.0;
	s.m_Values[ 0 ][ 1 ] =  1.0 / h2;
	s.m_Values[ 0 ][ 2 ] =  0.0;

	s.m_Values[ 1 ][ 0 ] =  0.0;
	s.m_Values[ 1 ][ 1 ] = -2.0 / h2;
	s.m_Values[ 1 ][ 2 ] =  0.0;

	s.m_Values[ 2 ][ 0 ] =  0.0;
	s.m_Values[ 2 ][ 1 ] =  1.0 / h2;
	s.m_Values[ 2 ][ 2 ] =  0.0;

	return s;
}

Stencil Stencil::lerpx( real h )
{
	Stencil s;
	s.m_Values[ 0 ][ 0 ] =  0.0;
	s.m_Values[ 0 ][ 1 ] =  0.0;
	s.m_Values[ 0 ][ 2 ] =  0.0;

	s.m_Values[ 1 ][ 0 ] =  0.0;
	s.m_Values[ 1 ][ 1 ] =  h;
	s.m_Values[ 1 ][ 2 ] =  ( 1.0 - h );

	s.m_Values[ 2 ][ 0 ] =  0.0;
	s.m_Values[ 2 ][ 1 ] =  0.0;
	s.m_Values[ 2 ][ 2 ] =  0.0;

	return s;
}

Stencil Stencil::lerpy( real h )
{
	Stencil s;
	s.m_Values[ 0 ][ 0 ] =  0.0;
	s.m_Values[ 0 ][ 1 ] =  0.0;
	s.m_Values[ 0 ][ 2 ] =  0.0;

	s.m_Values[ 1 ][ 0 ] =  0.0;
	s.m_Values[ 1 ][ 1 ] =  h;
	s.m_Values[ 1 ][ 2 ] =  0.0;

	s.m_Values[ 2 ][ 0 ] =  0.0;
	s.m_Values[ 2 ][ 1 ] =  ( 1.0 - h );
	s.m_Values[ 2 ][ 2 ] =  0.0;

	return s;
}

Stencil Stencil::sor( const Point& h )
{
	real dx2 = h.x * h.x;
	real dy2 = h.y * h.y;

	Stencil s;
	s.m_Values[ 0 ][ 0 ] =  0.0;
	s.m_Values[ 0 ][ 1 ] =  1.0 / dy2;
	s.m_Values[ 0 ][ 2 ] =  0.0;

	s.m_Values[ 1 ][ 0 ] =  1.0 / dx2;
	s.m_Values[ 1 ][ 1 ] =  0.0;
	s.m_Values[ 1 ][ 2 ] =  1.0 / dx2;

	s.m_Values[ 2 ][ 0 ] =  0.0;
	s.m_Values[ 2 ][ 1 ] =  1.0 / dy2;
	s.m_Values[ 2 ][ 2 ] =  0.0;

	return s;
}

void Stencil::apply( const GridFunction& source,
		const MultiIndex& readBegin,
		const MultiIndex& readEnd,
		GridFunction& image,
		const MultiIndex& writeBegin,
		const MultiIndex& writeEnd ) const
{
	// maximum size of the region that the stencil can be applied to
	MultiIndex readSize( readEnd.x - readBegin.x,
			readEnd.y - readBegin.y );

	// check if read size equals write size
	assert( readSize.x == writeEnd.x - writeBegin.x &&
			readSize.y == writeEnd.y - writeBegin.y );

	// check if non-zero entries require out-of-bounds accesses
	MY_ASSERT(
	for( index_t j = 0; j < readSize.y; ++j )
	{
		index_t j_abs = readBegin.y + j - SIZE / 2;
		bool j_out = ( j_abs < 0 || j_abs >= source.getSize().y );
		for( index_t i = 0; i < readSize.x; ++i )
		{
			index_t i_abs = readBegin.x + i - SIZE / 2;
			bool i_out = ( i_abs < 0 || i_abs >= source.getSize().x );
			assert( !i_out || !j_out || m_Values[ j ][ i ] != 0.0 );
		}
	}
	);

	/*
	// hard-coding stencil application inside double loop for SIZE = 3 to avoid
	// additional loops
	assert( SIZE == 3 );
	*/

	for( index_t j = 0; j < readSize.y; ++j )
	{
		for( index_t i = 0; i < readSize.x; ++i )
		{
			/*
			image( writeBegin.x + i, writeBegin.y + j ) =
				m_Values[ 0 ][ 0 ] * source( readBegin.x + i - 1, readBegin.y + j - 1 ) +
				m_Values[ 0 ][ 1 ] * source( readBegin.x + i    , readBegin.y + j - 1 ) +
				m_Values[ 0 ][ 2 ] * source( readBegin.x + i + 1, readBegin.y + j - 1 ) +

				m_Values[ 1 ][ 0 ] * source( readBegin.x + i - 1, readBegin.y + j     ) +
				m_Values[ 1 ][ 1 ] * source( readBegin.x + i    , readBegin.y + j     ) +
				m_Values[ 1 ][ 2 ] * source( readBegin.x + i + 1, readBegin.y + j     ) +

				m_Values[ 2 ][ 0 ] * source( readBegin.x + i - 1, readBegin.y + j + 1 ) +
				m_Values[ 2 ][ 1 ] * source( readBegin.x + i    , readBegin.y + j + 1 ) +
				m_Values[ 2 ][ 2 ] * source( readBegin.x + i + 1, readBegin.y + j + 1 );
			*/

			real sum = 0.0;
			for( index_t k = 0; k < SIZE; ++k )
			{
				for( index_t l = 0; l < SIZE; ++l )
				{
					if( m_Values[ k ][ l ] != 0.0 )
					{
						sum += m_Values[ k ][ l ] *	source(
								readBegin.x + l - SIZE / 2,
								readBegin.y + k - SIZE / 2 );
					}
				}
			}
			image( writeBegin.x + i, writeBegin.y + j ) = sum;
		}
	}
}

StencilSet::StencilSet( const Point& h )
	: h( h ),
	dx( Stencil::dx( h.x ) ),
	dy( Stencil::dy( h.y ) ),
	dx_bw( Stencil::dx_bw( h.x ) ),
	dy_bw( Stencil::dy_bw( h.y ) ),
	dx_fw( Stencil::dx_fw( h.x ) ),
	dy_fw( Stencil::dy_fw( h.y ) ),
	dxx( Stencil::dxx( h.x ) ),
	dyy( Stencil::dyy( h.y ) ),
	lerpx( Stencil::lerpx( 0.5 ) ),
	lerpy( Stencil::lerpy( 0.5 ) )
{
}

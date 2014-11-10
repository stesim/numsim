#ifndef STENCIL_H
#define STENCIL_H

#include "typedef.h"
#include "gridfunction.h"

class Stencil
{
public:
	static Stencil dx( real h );
	static Stencil dy( real h );
	static Stencil dx_bw( real h );
	static Stencil dy_bw( real h );
	static Stencil dx_fw( real h );
	static Stencil dy_fw( real h );
	static Stencil dxx( real h );
	static Stencil dyy( real h );
	static Stencil lerpx( real h );
	static Stencil lerpy( real h );
	static Stencil sor( const Point& h );

	void apply( const GridFunction& source,
			const MultiIndex& readBegin,
			const MultiIndex& readEnd,
			GridFunction& image,
			const MultiIndex& writeBegin,
			const MultiIndex& writeEnd ) const;

private:
	static const index_t SIZE = 3;
	real m_Values[ SIZE ][ SIZE ];
};

struct StencilSet
{
public:
	StencilSet( const Point& h );

public:
	Point   h;

	Stencil dx;
	Stencil dy;
	Stencil dx_bw;
	Stencil dy_bw;
	Stencil dx_fw;
	Stencil dy_fw;
	Stencil dxx;
	Stencil dyy;
	Stencil lerpx;
	Stencil lerpy;
};

#endif

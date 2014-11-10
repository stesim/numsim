#ifndef GRIDFUNCTION_H
#define GRIDFUNCTION_H

#include "typedef.h"

class GridFunction
{
public:
	GridFunction();
	GridFunction( const MultiIndex& size );
	GridFunction( index_t dimX, index_t dimY );

	~GridFunction();

	const MultiIndex& getSize() const;

	real& operator ()( index_t i, index_t j );
	const real& operator ()( index_t i, index_t j ) const;

	void set( const MultiIndex& begin,
			const MultiIndex& end,
			real value );

	void scale( const MultiIndex& begin,
			const MultiIndex& end,
			real factor );

	void offsetScale( const MultiIndex& begin,
			const MultiIndex& end,
			const MultiIndex& offset,
			real factor );

	void copy( const MultiIndex& begin,
			const MultiIndex& end,
			const GridFunction& source );

	void scaledCopy( const MultiIndex& begin,
			const MultiIndex& end,
			const GridFunction& source,
			real factor );

	void scaledOffsetCopy( const MultiIndex& begin,
			const MultiIndex& end,
			const GridFunction& source,
			const MultiIndex& offset,
			real factor );

	void affineOffsetCopy( const MultiIndex& begin,
			const MultiIndex& end,
			const GridFunction& source,
			const MultiIndex& offset,
			real factor,
			real constant );

	void addOffsetCopy( const MultiIndex& begin,
			const MultiIndex& end,
			const GridFunction& source,
			const MultiIndex& offset );

	void addScaledCopy( const MultiIndex& begin,
			const MultiIndex& end,
			const GridFunction& source,
			real factor );

	void multiply( const MultiIndex& begin,
			const MultiIndex& end,
			const GridFunction& source );

	void abs( const MultiIndex& begin,
			const MultiIndex& end );

	real getMaxValue( const MultiIndex& begin,
			const MultiIndex& end ) const;

private:
	inline void allocMem();
	inline void freeMem();

private:
	real*        m_pGridValues;
	MultiIndex   m_miSize;
};

#endif

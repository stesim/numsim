#include "geometry.h"
#include <metamath/metamath.h>

using namespace mm;

Geometry::Geometry( const Point& domainSize, const MultiIndex& gridSize )
	: m_DomainSize( domainSize ),
	m_GridSize( gridSize ),
	m_Raster( gridSize ),
	m_iShapeCounter( 1 ),
	m_MaskU( gridSize + MultiIndex( 3, 2 ) ),
	m_MaskV( gridSize + MultiIndex( 2, 3 ) ),
	m_MaskP( gridSize + MultiIndex( 2, 2 ) )
{
	RasterCell cell;
	cell.type = CellType::Fluid;
	cell.shapeId = 0;

	m_Raster = constant( cell );

	m_MaskU = constant( false );
	m_MaskV = constant( false );
	m_MaskP = constant( false );
	set( m_MaskU, MultiIndex( 2, 1 ), m_MaskU.size() - MultiIndex( 2, 1 ), constant( true ) );
	set( m_MaskV, MultiIndex( 1, 2 ), m_MaskV.size() - MultiIndex( 1, 2 ), constant( true ) );
	set( m_MaskP, MultiIndex::ONE, m_MaskP.size() - MultiIndex::ONE, constant( true ) );
}

void Geometry::generateSampleGeometry()
{
	MultiIndex gridSize = m_GridSize;

	MultiIndex squarePos = MultiIndex::ZERO;//MultiIndex( gridSize.x * 3 / 4, gridSize.y / 4 );
	MultiIndex squareSize = MultiIndex( gridSize.x / 2, gridSize.y / 2 );

	m_SettersU = std::vector<BoundarySetter>(
			2 *  ( gridSize.x + 3 ) + 2 * gridSize.y
			+ 2 * ( squareSize.x - 1 + squareSize.y ) );
	m_SettersV = std::vector<BoundarySetter>(
			2 *  ( gridSize.y + 3 ) + 2 * gridSize.x
			+ 2 * ( squareSize.x + squareSize.y - 1 ) );
	m_SettersP = std::vector<BoundarySetter>(
			2 * gridSize.x + 2 * gridSize.y
			+ 2 * ( squareSize.x + squareSize.y ) );

	BoundarySetter setter;

	// u boundaries
	for( index_t i = 1; i <= gridSize.y; ++i )
	{
		// left
		setter.index = MultiIndex( 1, i );
		setter.sourceIndex = MultiIndex( -1, -1 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = 0.0;

		m_SettersU.push_back( setter );

		// right
		setter.index = MultiIndex( gridSize.x + 1, i );
		setter.sourceIndex = MultiIndex( -1, -1 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = 0.0;

		m_SettersU.push_back( setter );
	}
	for( index_t i = 0; i < gridSize.x + 3; ++i )
	{
		// top
		setter.index = MultiIndex( i, gridSize.y + 1 );
		setter.sourceIndex = setter.index - MultiIndex( 0, 1 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 2.0;
		setter.factor = -1.0;

		m_SettersU.push_back( setter );

		// bottom
		setter.index = MultiIndex( i, 0 );
		setter.sourceIndex = setter.index + MultiIndex( 0, 1 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = -1.0;

		m_SettersU.push_back( setter );
	}
	// square
	for( index_t i = 0; i < squareSize.x - 1; ++i )
	{
		setter.index = MultiIndex( squarePos.x + 2 + i, squarePos.y + 1 );
		setter.sourceIndex = setter.index - MultiIndex( 0, 1 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = -1.0;

		m_SettersU.push_back( setter );

		setter.index = MultiIndex( squarePos.x + 2 + i,
				squarePos.y + 1 + squareSize.y - 1 );
		setter.sourceIndex = setter.index + MultiIndex( 0, 1 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = -1.0;

		m_SettersU.push_back( setter );
	}
	for( index_t i = 0; i < squareSize.y; ++i )
	{
		setter.index = MultiIndex( squarePos.x + 1, squarePos.y + 1 + i );
		setter.sourceIndex = MultiIndex( -1, -1 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = 0.0;

		m_SettersU.push_back( setter );

		setter.index = MultiIndex( squarePos.x + 1 + squareSize.x,
				squarePos.y + 1 + i );
		setter.sourceIndex = MultiIndex( -1, -1 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = 0.0;

		m_SettersU.push_back( setter );
	}

	// v boundaries
	for( index_t i = 0; i < gridSize.y + 3; ++i )
	{
		// left
		setter.index = MultiIndex( 0, i );
		setter.sourceIndex = setter.index + MultiIndex( 1, 0 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = -1.0;

		m_SettersV.push_back( setter );

		// right
		setter.index = MultiIndex( gridSize.x + 1, i );
		setter.sourceIndex = setter.index - MultiIndex( 1, 0 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = -1.0;

		m_SettersV.push_back( setter );
	}
	for( index_t i = 1; i <= gridSize.x; ++i )
	{
		// top
		setter.index = MultiIndex( i, gridSize.y + 1 );
		setter.sourceIndex = MultiIndex( -1, -1 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = 0.0;

		m_SettersV.push_back( setter );

		// bottom
		setter.index = MultiIndex( i, 1 );
		setter.sourceIndex = MultiIndex( -1, -1 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = 0.0;

		m_SettersV.push_back( setter );
	}
	// square
	for( index_t i = 0; i < squareSize.y - 1; ++i )
	{
		setter.index = MultiIndex( squarePos.x + 1, squarePos.y + 2 + i );
		setter.sourceIndex = setter.index - MultiIndex( 1, 0 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = -1.0;

		m_SettersV.push_back( setter );

		setter.index = MultiIndex( squarePos.x + 1 + squareSize.x - 1,
				squarePos.y + 2 + i );
		setter.sourceIndex = setter.index + MultiIndex( 1, 0 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = -1.0;

		m_SettersV.push_back( setter );
	}
	for( index_t i = 0; i < squareSize.x; ++i )
	{
		setter.index = MultiIndex( squarePos.x + 1 + i, squarePos.y + 1 );
		setter.sourceIndex = MultiIndex( -1, -1 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = 0.0;

		m_SettersV.push_back( setter );

		setter.index = MultiIndex( squarePos.x + 1 + i,
				squarePos.y + 1 + squareSize.y );
		setter.sourceIndex = MultiIndex( -1, -1 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = 0.0;

		m_SettersV.push_back( setter );
	}

	// p boundaries
	for( index_t i = 1; i <= gridSize.x; ++i )
	{
		// top
		setter.index = MultiIndex( i, gridSize.y + 1 );
		setter.sourceIndex = setter.index - MultiIndex( 0, 1 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = 1.0;

		m_SettersP.push_back( setter );

		// bottom
		setter.index = MultiIndex( i, 0 );
		setter.sourceIndex = setter.index + MultiIndex( 0, 1 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = 1.0;

		m_SettersP.push_back( setter );
	}
	for( index_t i = 1; i <= gridSize.y; ++i )
	{
		// left
		setter.index = MultiIndex( 0, i );
		setter.sourceIndex = setter.index + MultiIndex( 1, 0 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = 1.0;

		m_SettersP.push_back( setter );

		// right
		setter.index = MultiIndex( gridSize.x + 1, i );
		setter.sourceIndex = setter.index - MultiIndex( 1, 0 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = 1.0;

		m_SettersP.push_back( setter );
	}
	// square
	for( index_t i = 0; i < squareSize.x - 2; ++i )
	{
		setter.index = squarePos + MultiIndex( 2 + i, 1 );
		setter.sourceIndex = setter.index - MultiIndex( 0, 1 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = 1.0;

		m_SettersP.push_back( setter );

		setter.index = squarePos + MultiIndex( 2 + i, 1 + squareSize.y - 1 );
		setter.sourceIndex = setter.index + MultiIndex( 0, 1 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = 1.0;

		m_SettersP.push_back( setter );
	}

	setter.index = squarePos + MultiIndex::ONE;
	setter.sourceIndex = setter.index - MultiIndex( 0, 1 );
	setter.sourceIndex2 = setter.index - MultiIndex( 1, 0 );
	setter.constant = 0.0;
	setter.factor = 0.5;

	m_SettersP.push_back( setter );

	setter.index = squarePos + MultiIndex( 1 + squareSize.x - 1, 1 );
	setter.sourceIndex = setter.index - MultiIndex( 0, 1 );
	setter.sourceIndex2 = setter.index + MultiIndex( 1, 0 );
	setter.constant = 0.0;
	setter.factor = 0.5;

	m_SettersP.push_back( setter );
	
	for( index_t i = 0; i < squareSize.y - 2; ++i )
	{
		setter.index = squarePos + MultiIndex( 1, 2 + i );
		setter.sourceIndex = setter.index - MultiIndex( 1, 0 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = 1.0;

		m_SettersP.push_back( setter );

		setter.index = squarePos + MultiIndex( 1 + squareSize.x - 1, 2 + i );
		setter.sourceIndex = setter.index + MultiIndex( 1, 0 );
		setter.sourceIndex2 = MultiIndex( -1, -1 );
		setter.constant = 0.0;
		setter.factor = 1.0;

		m_SettersP.push_back( setter );
	}

	setter.index = squarePos + MultiIndex( 1, 1 + squareSize.y - 1 );
	setter.sourceIndex = setter.index - MultiIndex( 1, 0 );
	setter.sourceIndex2 = setter.index + MultiIndex( 0, 1 );
	setter.constant = 0.0;
	setter.factor = 0.5;

	m_SettersP.push_back( setter );

	setter.index = squarePos + squareSize;
	setter.sourceIndex = setter.index + MultiIndex( 1, 0 );
	setter.sourceIndex2 = setter.index + MultiIndex( 0, 1 );
	setter.constant = 0.0;
	setter.factor = 0.5;

	m_SettersP.push_back( setter );

	// mask functions
	set( m_MaskU, MultiIndex::ZERO, MultiIndex( 2, gridSize.y + 2 ),
			constant( false ) );
	set( m_MaskU, MultiIndex( gridSize.x + 1, 0 ),
			gridSize + MultiIndex( 3, 2 ), constant( false ) );
	set( m_MaskU, MultiIndex( 2, 0 ), MultiIndex( gridSize.x + 1, 1 ),
			constant( false ) );
	set( m_MaskU, MultiIndex( 2, gridSize.y + 1 ),
			gridSize + MultiIndex( 1, 2 ), constant( false ) );

	set( m_MaskU, squarePos + MultiIndex::ONE,
			squarePos + squareSize + MultiIndex( 2, 1 ), constant( false ) );

	set( m_MaskV, MultiIndex::ZERO, MultiIndex( gridSize.x + 2, 2 ),
			constant( false ) );
	set( m_MaskV, MultiIndex( 0, gridSize.y + 1 ),
			gridSize + MultiIndex( 2, 3 ), constant( false ) );
	set( m_MaskV, MultiIndex( 0, 2 ), MultiIndex( 1, gridSize.y + 1 ),
			constant( false ) );
	set( m_MaskV, MultiIndex( gridSize.x + 1, 2 ),
			gridSize + MultiIndex( 2, 1 ), constant( false ) );

	set( m_MaskV, squarePos + MultiIndex::ONE,
			squarePos + squareSize + MultiIndex( 1, 2 ), constant( false ) );

	set( m_MaskP, MultiIndex::ZERO, MultiIndex( 1, gridSize.y + 2 ),
			constant( false ) );
	set( m_MaskP, MultiIndex( 0, gridSize.y + 1 ),
			gridSize + MultiIndex( 2, 2 ), constant( false ) );
	set( m_MaskP, MultiIndex( 0, 1 ), MultiIndex( 1, gridSize.y + 1 ),
			constant( false ) );
	set( m_MaskP, MultiIndex( gridSize.x + 1, 1 ),
			gridSize + MultiIndex( 2, 1 ), constant( false ) );

	set( m_MaskP, squarePos + MultiIndex::ONE,
			squarePos + squareSize + MultiIndex::ONE, constant( false ) );
}

void Geometry::addPolygon( const Polygon& polygon )
{
	index_t id = m_iShapeCounter++;

	if( polygon.m_Lines.empty() )
	{
		return;
	}

	MultiIndex min;
	MultiIndex max;

	rasterizeLine( polygon.m_StartPoint,
			polygon.m_Lines[ 0 ].endPoint,
			polygon.m_Lines[ 0 ].type,
			min,
			max,
			id,
			polygon.m_Lines[ 0 ].value );

	for( std::size_t i = 1; i < polygon.m_Lines.size(); ++i )
	{
		rasterizeLine( polygon.m_Lines[ i - 1 ].endPoint,
				polygon.m_Lines[ i ].endPoint,
				polygon.m_Lines[ i ].type,
				min,
				max,
				id,
				polygon.m_Lines[ i ].value );
	}

	// naive scanline algorithm, assumes no intersecting lines
	for( index_t j = min.y; j <= max.y; ++j )
	{
		index_t imin = max.x;
		index_t imax = min.x;
		for( index_t i = min.x; i <= max.x; ++i )
		{
			if( m_Raster( i, j ).type == CellType::Boundary
					&& m_Raster( i, j ).shapeId == id )
			{
				imin = i;
				break;
			}
		}
		for( index_t i = max.x; i >= min.x; --i )
		{
			if( m_Raster( i, j ).type == CellType::Boundary
					&& m_Raster( i, j ).shapeId == id )
			{
				imax = i;
				break;
			}
		}

		for( index_t i = imin + 1; i < imax; ++i )
		{
			m_Raster( i, j ).type = CellType::Obstacle;
		}
	}
}

void Geometry::addRectangle( const Point& start, const Point& end,
		BoundaryType boundary, const Point& value )
{
	Polygon poly;
	poly.setStartPoint( start );
	poly.addLine( Point( end.x, start.y ), boundary, value );
	poly.addLine( end, boundary, value );
	poly.addLine( Point( start.x, end.y ), boundary, value );
	poly.addLine( start, boundary, value );

	addPolygon( poly );
}

void Geometry::addRectangle( const Point& start, const Point& end,
		BoundaryType boundary )
{
	addRectangle( start, end, boundary, Point::ZERO );
}

void Geometry::setLeftBoundaryCondition( BoundaryType condition,
		const Point& value )
{
	BoundarySetter setter;
	switch( condition )
	{
		case BoundaryType::Noslip:
		{
			for( index_t i = 1; i <= m_GridSize.y + 1; ++i )
			{
				setter.index = MultiIndex( 0, i );
				setter.sourceIndex = setter.index + MultiIndex( 1, 0 );
				setter.sourceIndex2 = MultiIndex( -1, -1 );
				setter.constant = 0.0;
				setter.factor = -1.0;

				m_SettersV.push_back( setter );
			}
			for( index_t i = 0; i < m_GridSize.y; ++i )
			{
				setter.index = MultiIndex( 1, 1 + i );
				setter.sourceIndex = MultiIndex( -1, -1 );
				setter.constant = 0.0;

				m_SettersU.push_back( setter );

				setter.index = MultiIndex( 0, 1 + i );
				setter.sourceIndex = setter.index + MultiIndex( 1, 0 );
				setter.constant = 0.0;
				setter.factor = 1.0;

				m_SettersP.push_back( setter );
			}
			break;
		}
		case BoundaryType::Inflow:
		{
			for( index_t i = 1; i <= m_GridSize.y + 1; ++i )
			{
				setter.index = MultiIndex( 0, i );
				setter.sourceIndex = setter.index + MultiIndex( 1, 0 );
				setter.sourceIndex2 = MultiIndex( -1, -1 );
				setter.constant = 2.0 * value.y;
				setter.factor = -1.0;

				m_SettersV.push_back( setter );
			}
			for( index_t i = 0; i < m_GridSize.y; ++i )
			{
				setter.index = MultiIndex( 1, 1 + i );
				setter.sourceIndex = MultiIndex( -1, -1 );
				setter.constant = value.x;

				m_SettersU.push_back( setter );

				setter.index = MultiIndex( 0, 1 + i );
				setter.sourceIndex = setter.index + MultiIndex( 1, 0 );
				setter.constant = 0.0;
				setter.factor = 1.0;

				m_SettersP.push_back( setter );
			}
			break;
		}
		default:
		{
			break;
		}
	};
}

void Geometry::setRightBoundaryCondition( BoundaryType condition,
		const Point& value )
{
	BoundarySetter setter;
	switch( condition )
	{
		case BoundaryType::Noslip:
		{
			for( index_t i = 1; i <= m_GridSize.y + 1; ++i )
			{
				setter.index = MultiIndex( m_GridSize.x + 1, i );
				setter.sourceIndex = setter.index - MultiIndex( 1, 0 );
				setter.sourceIndex2 = MultiIndex( -1, -1 );
				setter.constant = 0.0;
				setter.factor = -1.0;

				m_SettersV.push_back( setter );
			}
			for( index_t i = 0; i < m_GridSize.y; ++i )
			{
				setter.index = MultiIndex( m_GridSize.x + 1, 1 + i );
				setter.sourceIndex = MultiIndex( -1, -1 );
				setter.constant = 0.0;

				m_SettersU.push_back( setter );

				setter.index = MultiIndex( m_GridSize.x + 1, 1 + i );
				setter.sourceIndex = setter.index - MultiIndex( 1, 0 );
				setter.constant = 0.0;
				setter.factor = 1.0;

				m_SettersP.push_back( setter );
			}
			break;
		}
		case BoundaryType::Inflow:
		{
			for( index_t i = 1; i <= m_GridSize.y + 1; ++i )
			{
				setter.index = MultiIndex( m_GridSize.x + 1, i );
				setter.sourceIndex = setter.index - MultiIndex( 1, 0 );
				setter.sourceIndex2 = MultiIndex( -1, -1 );
				setter.constant = 2.0 * value.y;
				setter.factor = -1.0;

				m_SettersV.push_back( setter );
			}
			for( index_t i = 0; i < m_GridSize.y; ++i )
			{
				setter.index = MultiIndex( m_GridSize.x + 1, 1 + i );
				setter.sourceIndex = MultiIndex( -1, -1 );
				setter.constant = value.x;

				m_SettersU.push_back( setter );

				setter.index = MultiIndex( m_GridSize.x + 1, 1 + i );
				setter.sourceIndex = setter.index - MultiIndex( 1, 0 );
				setter.constant = 0.0;
				setter.factor = 1.0;

				m_SettersP.push_back( setter );
			}
			break;
		}
		default:
		{
			break;
		}
	};
}

void Geometry::setBottomBoundaryCondition( BoundaryType condition,
		const Point& value )
{
	BoundarySetter setter;
	switch( condition )
	{
		case BoundaryType::Noslip:
		{
			for( index_t i = 1; i <= m_GridSize.x + 1; ++i )
			{
				setter.index = MultiIndex( i, 0 );
				setter.sourceIndex = setter.index + MultiIndex( 0, 1 );
				setter.sourceIndex2 = MultiIndex( -1, -1 );
				setter.constant = 0.0;
				setter.factor = -1.0;

				m_SettersU.push_back( setter );
			}
			for( index_t i = 0; i < m_GridSize.x; ++i )
			{
				setter.index = MultiIndex( 1 + i, 1 );
				setter.sourceIndex = MultiIndex( -1, -1 );
				setter.constant = 0.0;

				m_SettersV.push_back( setter );

				setter.index = MultiIndex( 1 + i, 0 );
				setter.sourceIndex = setter.index + MultiIndex( 0, 1 );
				setter.constant = 0.0;
				setter.factor = 1.0;

				m_SettersP.push_back( setter );
			}
			break;
		}
		case BoundaryType::Inflow:
		{
			for( index_t i = 1; i <= m_GridSize.x + 1; ++i )
			{
				setter.index = MultiIndex( i, 0 );
				setter.sourceIndex = setter.index + MultiIndex( 0, 1 );
				setter.sourceIndex2 = MultiIndex( -1, -1 );
				setter.constant = 2.0 * value.x;
				setter.factor = -1.0;

				m_SettersU.push_back( setter );
			}
			for( index_t i = 0; i < m_GridSize.x; ++i )
			{
				setter.index = MultiIndex( 1 + i, 1 );
				setter.sourceIndex = MultiIndex( -1, -1 );
				setter.constant = value.y;

				m_SettersV.push_back( setter );

				setter.index = MultiIndex( 1 + i, 0 );
				setter.sourceIndex = setter.index + MultiIndex( 0, 1 );
				setter.constant = 0.0;
				setter.factor = 1.0;

				m_SettersP.push_back( setter );
			}
			break;
		}
		default:
		{
			break;
		}
	};
}

void Geometry::setTopBoundaryCondition( BoundaryType condition,
		const Point& value )
{
	BoundarySetter setter;
	switch( condition )
	{
		case BoundaryType::Noslip:
		{
			for( index_t i = 1; i <= m_GridSize.x + 1; ++i )
			{
				setter.index = MultiIndex( i, m_GridSize.y + 1 );
				setter.sourceIndex = setter.index - MultiIndex( 0, 1 );
				setter.sourceIndex2 = MultiIndex( -1, -1 );
				setter.constant = 0.0;
				setter.factor = -1.0;

				m_SettersU.push_back( setter );
			}
			for( index_t i = 0; i < m_GridSize.x; ++i )
			{
				setter.index = MultiIndex( 1 + i, m_GridSize.y + 1 );
				setter.sourceIndex = MultiIndex( -1, -1 );
				setter.constant = 0.0;

				m_SettersV.push_back( setter );

				setter.index = MultiIndex( 1 + i, m_GridSize.y + 1 );
				setter.sourceIndex = setter.index - MultiIndex( 0, 1 );
				setter.constant = 0.0;
				setter.factor = 1.0;

				m_SettersP.push_back( setter );
			}
			break;
		}
		case BoundaryType::Inflow:
		{
			for( index_t i = 1; i <= m_GridSize.x + 1; ++i )
			{
				setter.index = MultiIndex( i, m_GridSize.y + 1 );
				setter.sourceIndex = setter.index - MultiIndex( 0, 1 );
				setter.sourceIndex2 = MultiIndex( -1, -1 );
				setter.constant = 2.0 * value.x;
				setter.factor = -1.0;

				m_SettersU.push_back( setter );
			}
			for( index_t i = 0; i < m_GridSize.x; ++i )
			{
				setter.index = MultiIndex( 1 + i, m_GridSize.y + 1 );
				setter.sourceIndex = MultiIndex( -1, -1 );
				setter.constant = value.y;

				m_SettersV.push_back( setter );

				setter.index = MultiIndex( 1 + i, m_GridSize.y + 1 );
				setter.sourceIndex = setter.index - MultiIndex( 0, 1 );
				setter.constant = 0.0;
				setter.factor = 1.0;

				m_SettersP.push_back( setter );
			}
			break;
		}
		default:
		{
			break;
		}
	};
}

void Geometry::bake()
{
	/*
	m_SettersU.clear();
	m_SettersV.clear();
	m_SettersP.clear();

	set( m_MaskU, MultiIndex( 2, 1 ), m_MaskU.size() - MultiIndex( 2, 1 ), constant( true ) );
	set( m_MaskV, MultiIndex( 1, 2 ), m_MaskV.size() - MultiIndex( 1, 2 ), constant( true ) );
	set( m_MaskP, MultiIndex::ONE, m_MaskP.size() - MultiIndex::ONE, constant( true ) );
	*/

	std::size_t numBoundaryCells = 0;

	// determine boundary cell orientation
	for( index_t j = 0; j < m_Raster.size().y; ++j )
	{
		for( index_t i = 0; i < m_Raster.size().x; ++i )
		{
			if( m_Raster( i, j ).type == CellType::Boundary )
			{
				m_Raster( i, j ).orientation = static_cast<Orientation>(
						( j > 0 && m_Raster( i, j - 1 ).type == CellType::Fluid )
						| ( i < m_Raster.size().x - 1 && m_Raster( i + 1, j ).type == CellType::Fluid ) << 1
						| ( j < m_Raster.size().y - 1 && m_Raster( i, j + 1 ).type == CellType::Fluid ) << 2
						| ( i > 0 && m_Raster( i - 1, j ).type == CellType::Fluid ) << 3 );
				++numBoundaryCells;
			}
			// mark preliminary obstacle nodes in masks
			if( m_Raster( i, j ).type != CellType::Fluid )
			{
				m_MaskU( i + 2, j + 1 ) = false;
				m_MaskV( i + 1, j + 2 ) = false;
				m_MaskP( i + 1, j + 1 ) = false;
			}
		}
	}

	m_SettersU.reserve( numBoundaryCells );
	m_SettersV.reserve( numBoundaryCells );
	m_SettersP.reserve( numBoundaryCells );

	static const MultiIndex mi10( 1, 0 );
	static const MultiIndex mi01( 0, 1 );
	static const MultiIndex mi12( 1, 2 );
	static const MultiIndex mi21( 2, 1 );
	static const MultiIndex miNeg( -1, -1 );

	for( index_t j = 0; j < m_Raster.size().y; ++j )
	{
		for( index_t i = 0; i < m_Raster.size().x; ++i )
		{
			const RasterCell& cell = m_Raster( i, j );
			if( cell.type != CellType::Boundary )
			{
				continue;
			}

			MultiIndex idx( i, j );
			BoundarySetter setter;

			switch( cell.boundary )
			{
				case BoundaryType::Noslip:
				{
					switch( cell.orientation )
					{
						case Orientation::North:
						case Orientation::South:
						{
							int dir = ( cell.orientation == Orientation::North ) ? 1 : -1;

							setter.index = idx + MultiIndex::ONE;
							setter.sourceIndex =
								setter.index + MultiIndex( 0, dir );
							setter.sourceIndex2 = miNeg;
							setter.constant = 0.0;
							setter.factor = -1.0;

							m_SettersU.push_back( setter );

							setter.index = idx + MultiIndex( 1, 1 + ( cell.orientation == Orientation::North ) );
							setter.sourceIndex = miNeg;
							setter.constant = 0.0;

							m_SettersV.push_back( setter );

							setter.index = idx + MultiIndex::ONE;
							setter.sourceIndex =
								setter.index + MultiIndex( 0, dir );
							setter.constant = 0.0;
							setter.factor = 1.0;

							m_SettersP.push_back( setter );
							break;
						}
						case Orientation::East:
						case Orientation::West:
						{
							int dir = ( cell.orientation == Orientation::East ) ? 1 : -1;

							setter.index = idx + MultiIndex::ONE;
							setter.sourceIndex =
								setter.index + MultiIndex( dir, 0 );
							setter.sourceIndex2 = miNeg;
							setter.constant = 0.0;
							setter.factor = -1.0;

							m_SettersV.push_back( setter );

							setter.index = idx + MultiIndex( 1 + ( cell.orientation == Orientation::East ), 1 );
							setter.sourceIndex = miNeg;
							setter.constant = 0.0;

							m_SettersU.push_back( setter );

							setter.index = idx + MultiIndex::ONE;
							setter.sourceIndex =
								setter.index + MultiIndex( dir, 0 );
							setter.constant = 0.0;
							setter.factor = 1.0;

							m_SettersP.push_back( setter );
							break;
						}
						case Orientation::NorthEast:
						case Orientation::SouthEast:
						case Orientation::SouthWest:
						case Orientation::NorthWest:
						{
							int dirH =
								( cell.orientation == Orientation::NorthEast
								  || cell.orientation == Orientation::SouthEast )
								? 1 : -1;
							int dirV =
								( cell.orientation == Orientation::NorthEast
								  || cell.orientation == Orientation::NorthWest )
								? 1 : -1;
							int posH = ( dirH > 0 );
							int posV = ( dirV > 0 );

							if( posH )
							{
								setter.index = idx + MultiIndex::ONE;
								setter.sourceIndex =
									setter.index + MultiIndex( 0, dirV );
								setter.sourceIndex2 = miNeg;
								setter.constant = 0.0;
								setter.factor = -1.0;

								m_SettersU.push_back( setter );
							}

							setter.index = idx + MultiIndex( 1 + posH, 1 );
							setter.sourceIndex = miNeg;
							setter.sourceIndex2 = miNeg;
							setter.constant = 0.0;

							m_SettersU.push_back( setter );

							if( posV )
							{
								setter.index = idx + MultiIndex::ONE;
								setter.sourceIndex =
									setter.index + MultiIndex( dirH, 0 );
								setter.sourceIndex2 = miNeg;
								setter.constant = 0.0;
								setter.factor = -1.0;

								m_SettersV.push_back( setter );
							}

							setter.index = idx + MultiIndex( 1, 1 + posV );
							setter.sourceIndex = miNeg;
							setter.sourceIndex2 = miNeg;
							setter.constant = 0.0;

							m_SettersV.push_back( setter );

							setter.index = idx + MultiIndex::ONE;
							setter.sourceIndex =
								setter.index + MultiIndex( dirH, 0 );
							setter.sourceIndex2 =
								setter.index + MultiIndex( 0, dirV );
							setter.constant = 0.0;
							setter.factor = 0.5;

							m_SettersP.push_back( setter );
							break;
						}
					}
					break;
				}
				case BoundaryType::Inflow:
				{
					switch( cell.orientation )
					{
						case Orientation::North:
						case Orientation::South:
						{
							int dir = ( cell.orientation == Orientation::North ) ? 1 : -1;

							setter.index = idx + MultiIndex::ONE;
							setter.sourceIndex =
								setter.index + MultiIndex( 0, dir );
							setter.sourceIndex2 = miNeg;
							setter.constant = 2.0 * cell.value.x;
							setter.factor = -1.0;

							//printf( "--- %f\n", setter.constant );

							m_SettersU.push_back( setter );

							setter.index = idx + MultiIndex( 1, 1 + ( cell.orientation == Orientation::North ) );
							setter.sourceIndex = miNeg;
							setter.constant = cell.value.y;

							m_SettersV.push_back( setter );

							setter.index = idx + MultiIndex::ONE;
							setter.sourceIndex =
								setter.index + MultiIndex( 0, dir );
							setter.constant = 0.0;
							setter.factor = 1.0;

							m_SettersP.push_back( setter );
							break;
						}
						case Orientation::East:
						case Orientation::West:
						{
							int dir = ( cell.orientation == Orientation::East ) ? 1 : -1;

							setter.index = idx + MultiIndex::ONE;
							setter.sourceIndex =
								setter.index + MultiIndex( dir, 0 );
							setter.sourceIndex2 = miNeg;
							setter.constant = 2.0 * cell.value.y;
							setter.factor = -1.0;

							m_SettersV.push_back( setter );

							setter.index = idx + MultiIndex( 1 + ( cell.orientation == Orientation::East ), 1 );
							setter.sourceIndex = miNeg;
							setter.constant = cell.value.x;

							m_SettersU.push_back( setter );

							setter.index = idx + MultiIndex::ONE;
							setter.sourceIndex =
								setter.index + MultiIndex( dir, 0 );
							setter.constant = 0.0;
							setter.factor = 1.0;

							m_SettersP.push_back( setter );
							break;
						}
						case Orientation::NorthEast:
						case Orientation::SouthEast:
						case Orientation::SouthWest:
						case Orientation::NorthWest:
						{
							int dirH =
								( cell.orientation == Orientation::NorthEast
								  || cell.orientation == Orientation::SouthEast )
								? 1 : -1;
							int dirV =
								( cell.orientation == Orientation::NorthEast
								  || cell.orientation == Orientation::NorthWest )
								? 1 : -1;
							int posH = ( dirH > 0 );
							int posV = ( dirV > 0 );

							if( posH )
							{
								setter.index = idx + MultiIndex::ONE;
								setter.sourceIndex =
									setter.index + MultiIndex( 0, dirV );
								setter.sourceIndex2 = miNeg;
								setter.constant = 2.0 * cell.value.x;
								setter.factor = -1.0;

								m_SettersU.push_back( setter );
							}

							setter.index = idx + MultiIndex( 1 + posH, 1 );
							setter.sourceIndex = miNeg;
							setter.sourceIndex2 = miNeg;
							setter.constant = cell.value.x;

							m_SettersU.push_back( setter );

							if( posV )
							{
								setter.index = idx + MultiIndex::ONE;
								setter.sourceIndex =
									setter.index + MultiIndex( dirH, 0 );
								setter.sourceIndex2 = miNeg;
								setter.constant = 2.0 * cell.value.y;
								setter.factor = -1.0;

								m_SettersV.push_back( setter );
							}

							setter.index = idx + MultiIndex( 1, 1 + posV );
							setter.sourceIndex = miNeg;
							setter.sourceIndex2 = miNeg;
							setter.constant = cell.value.y;

							m_SettersV.push_back( setter );

							setter.index = idx + MultiIndex::ONE;
							setter.sourceIndex =
								setter.index + MultiIndex( dirH, 0 );
							setter.sourceIndex2 =
								setter.index + MultiIndex( 0, dirV );
							setter.constant = 0.0;
							setter.factor = 0.5;

							m_SettersP.push_back( setter );
							break;
						}
					}
					break;
				}
				default:
				{
					break;
				}
			}
		}
	}

	// mark all additional boundary nodes in masks
	for( std::size_t i = 0; i < m_SettersU.size(); ++i )
	{
		m_MaskU( m_SettersU[ i ].index.x, m_SettersU[ i ].index.y ) = false;
	}
	for( std::size_t i = 0; i < m_SettersV.size(); ++i )
	{
		m_MaskV( m_SettersV[ i ].index.x, m_SettersV[ i ].index.y ) = false;
	}
	for( std::size_t i = 0; i < m_SettersP.size(); ++i )
	{
		m_MaskP( m_SettersP[ i ].index.x, m_SettersP[ i ].index.y ) = false;
	}

	for( index_t j = 0; j < m_Raster.size().y; ++j )
	{
		for( index_t i = 0; i < m_Raster.size().x; ++i )
		{
			printf( "%d ", m_Raster( i, j ).type == CellType::Obstacle );
		}
		printf( "\n" );
	}
}

void Geometry::applyVelocityBoundary( GridFunction& u, GridFunction& v ) const
{
	applySetters( m_SettersU, u );
	applySetters( m_SettersV, v );
}

void Geometry::applyPressureBoundary( GridFunction& p ) const
{
	applySetters( m_SettersP, p );
}

const MaskFunction& Geometry::getComputationMaskU() const
{
	return m_MaskU;
}

const MaskFunction& Geometry::getComputationMaskV() const
{
	return m_MaskV;
}

const MaskFunction& Geometry::getComputationMaskP() const
{
	return m_MaskP;
}

void Geometry::rasterizeLine( const Point& start, const Point& end,
		BoundaryType boundary, MultiIndex& minIdx, MultiIndex& maxIdx,
		index_t id, const Point& value )
{
	int x1 = std::round( start.x / m_DomainSize.x * ( m_GridSize.x - 1 ) );
	int y1 = std::round( start.y / m_DomainSize.y * ( m_GridSize.y - 1 ) );
	int x2 = std::round( end.x / m_DomainSize.x * ( m_GridSize.x - 1 ) );
	int y2 = std::round( end.y / m_DomainSize.y * ( m_GridSize.y - 1 ) );

	minIdx.x = std::min( minIdx.x, std::min( x1, x2 ) );
	minIdx.y = std::min( minIdx.y, std::min( y1, y2 ) );
	maxIdx.x = std::max( maxIdx.x, std::max( x1, x2 ) );
	maxIdx.y = std::max( maxIdx.y, std::max( y1, y2 ) );

	// rasterize line using Bresenham's algorithm
	int x, y, dx, dy, dx1, dy1, px, py, xe, ye, i;
	dx = x2 - x1;
	dy = y2 - y1;
	dx1 = abs( dx );
	dy1 = abs( dy );
	px = 2 * dy1 - dx1;
	py = 2 * dx1 - dy1;
	if( dy1 <= dx1 )
	{
		if( dx >= 0 )
		{
			x = x1;
			y = y1;
			xe = x2;
		}
		else
		{
			x = x2;
			y = y2;
			xe = x1;
		}

		m_Raster( x, y ).type = CellType::Boundary;
		m_Raster( x, y ).boundary = boundary;
		m_Raster( x, y ).shapeId = id;
		m_Raster( x, y ).value = value;

		for( i = 0; x < xe; i++ )
		{
			x = x + 1;
			if( px < 0 )
			{
				px = px + 2 * dy1;
			}
			else
			{
				if( ( dx < 0 && dy < 0 ) || ( dx > 0 && dy > 0 ) )
				{
					y = y + 1;
				}
				else
				{
					y = y - 1;
				}
				px = px + 2 * ( dy1 - dx1 );
			}

			m_Raster( x, y ).type = CellType::Boundary;
			m_Raster( x, y ).boundary = boundary;
			m_Raster( x, y ).shapeId = id;
			m_Raster( x, y ).value = value;
		}
	}
	else
	{
		if( dy >= 0 )
		{
			x = x1;
			y = y1;
			ye = y2;
		}
		else
		{
			x = x2;
			y = y2;
			ye = y1;
		}

		m_Raster( x, y ).type = CellType::Boundary;
		m_Raster( x, y ).boundary = boundary;
		m_Raster( x, y ).shapeId = id;
		m_Raster( x, y ).value = value;

		for( i = 0; y < ye; i++ )
		{
			y = y + 1;
			if( py <= 0 )
			{
				py = py + 2 * dx1;
			}
			else
			{
				if( ( dx < 0 && dy < 0 ) || ( dx > 0 && dy > 0 ) )
				{
					x = x + 1;
				}
				else
				{
					x = x - 1;
				}
				py = py + 2 * ( dx1 - dy1 );
			}

			m_Raster( x, y ).type = CellType::Boundary;
			m_Raster( x, y ).boundary = boundary;
			m_Raster( x, y ).shapeId = id;
			m_Raster( x, y ).value = value;
		}
	}
}

void Geometry::applySetters( const std::vector<BoundarySetter>& setters,
		GridFunction& fun )
{
	for( std::size_t i = 0; i < setters.size(); ++i )
	{
		const BoundarySetter& setter = setters[ i ];
		real& val = fun( setter.index.x, setter.index.y );
		val = setter.constant;
		if( setter.sourceIndex.x >= 0 )
		{
			val += setter.factor *
				fun( setter.sourceIndex.x, setter.sourceIndex.y );
		}
		if( setter.sourceIndex2.x >= 0 )
		{
			val += setter.factor *
				fun( setter.sourceIndex2.x, setter.sourceIndex2.y );
		}
	}
}


Geometry::Polygon::Line::Line( const Point& pt, BoundaryType t,
		const Point& val )
	: endPoint( pt ), type( t ), value( val )
{
}

void Geometry::Polygon::setStartPoint( const Point& point )
{
	m_StartPoint = point;
}

void Geometry::Polygon::addLine( const Point& endPoint, BoundaryType boundary )
{
	m_Lines.push_back( Line( endPoint, boundary, Point::ZERO ) );
}

void Geometry::Polygon::addLine( const Point& endPoint, BoundaryType boundary,
		const Point& boundaryValue )
{
	m_Lines.push_back( Line( endPoint, boundary, boundaryValue ) );
}

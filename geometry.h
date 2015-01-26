#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "typedef.h"
#include <vector>

class Geometry
{
public:
	enum class BoundaryType
	{
		Noslip = 0,
		Slip,
		Inflow,
		Outflow
	};

	class Polygon;

private:
	enum class CellType
	{
		Fluid = 0,
		Boundary,
		Obstacle
	};

	enum class Orientation
	{
		North = 1,
		NorthEast = 3,
		East = 2,
		SouthEast = 6,
		South = 4,
		SouthWest = 12,
		West = 8,
		NorthWest = 9
	};

	struct RasterCell
	{
		CellType type;
		BoundaryType boundary;
		Point value;
		index_t shapeId;
		Orientation orientation;
	};

	typedef mm::Function<RasterCell> RasterFunction;

	struct BoundarySetter
	{
		MultiIndex index;
		MultiIndex sourceIndex;
		MultiIndex sourceIndex2;
		real constant;
		real factor;
	};

public:
	Geometry( const Point& domainSize, const MultiIndex& gridSize );

	void generateSampleGeometry();

	void addPolygon( const Polygon& polygon );

	void addRectangle( const Point& start, const Point& end,
			BoundaryType boundary, const Point& value );

	void addRectangle( const Point& start, const Point& end,
			BoundaryType boundary );

	void setLeftBoundaryCondition( BoundaryType condition, const Point& value = Point::ZERO );

	void setRightBoundaryCondition( BoundaryType condition, const Point& value = Point::ZERO );

	void setTopBoundaryCondition( BoundaryType condition, const Point& value = Point::ZERO );

	void setBottomBoundaryCondition( BoundaryType condition, const Point& value = Point::ZERO );

	void bake();

	void applyVelocityBoundary( GridFunction& u, GridFunction& v ) const;

	void applyPressureBoundary( GridFunction& p ) const;

	const MaskFunction& getComputationMaskU() const;

	const MaskFunction& getComputationMaskV() const;

	const MaskFunction& getComputationMaskP() const;

private:
	void rasterizeLine( const Point& start, const Point& end,
			BoundaryType boundary, MultiIndex& minIdx, MultiIndex& maxIdx,
			index_t id, const Point& value );

	void setBoundaryCondition( bool horizontal, bool positiveDir,
			BoundaryType condition, const Point& value );

	static void applySetters( const std::vector<BoundarySetter>& setters,
			GridFunction& fun );

private:
	Point m_DomainSize;
	MultiIndex m_GridSize;

	RasterFunction m_Raster;
	index_t m_iShapeCounter;

	std::vector<BoundarySetter> m_SettersU;
	std::vector<BoundarySetter> m_SettersV;
	std::vector<BoundarySetter> m_SettersP;
	MaskFunction m_MaskU;
	MaskFunction m_MaskV;
	MaskFunction m_MaskP;
};

class Geometry::Polygon
{
friend class Geometry;

private:
	struct Line
	{
		Point endPoint;
		BoundaryType type;
		Point value;

		Line( const Point& pt, BoundaryType t, const Point& val );
	};

public:
	void setStartPoint( const Point& point );

	void addLine( const Point& endPoint, BoundaryType boundary );

	void addLine( const Point& endPoint, BoundaryType boundary,
			const Point& boundaryValue );

private:
	Point m_StartPoint;
	std::vector<Line> m_Lines;
};

#endif

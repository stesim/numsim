#include "computation.h"
#include <algorithm>

Computation::Computation( const Params& params )
	: m_Params( params ),
	h( params.domainSize.x / params.gridSize.x,
			params.domainSize.y / params.gridSize.y ),
	u( m_Params.gridSize.x + 1, m_Params.gridSize.y + 2 ),
	v( m_Params.gridSize.x + 2, m_Params.gridSize.y + 1 ),
	f( m_Params.gridSize.x + 1, m_Params.gridSize.y + 2 ),
	g( m_Params.gridSize.x + 2, m_Params.gridSize.y + 1 ),
	p( m_Params.gridSize.x + 2, m_Params.gridSize.y + 2 ),
	rhs( m_Params.gridSize.x, m_Params.gridSize.y ),
	m_funcTemp1( m_Params.gridSize.x + 2, m_Params.gridSize.y + 2 ),
	m_funcTemp2( m_Params.gridSize.x + 2, m_Params.gridSize.y + 2 ),
	m_funcTemp3( m_Params.gridSize.x + 2, m_Params.gridSize.y + 2 ),
	m_stDxBw( Stencil::dx_bw( h.x ) ),
	m_stDyBw( Stencil::dy_bw( h.y ) ),
	m_stDxx( Stencil::dxx( h.x ) ),
	m_stDyy( Stencil::dyy( h.y ) ),
	m_stInterpx( Stencil::lerpx( 0.5 ) ),
	m_stInterpy( Stencil::lerpy( 0.5 ) )
{
	u.set( MultiIndex::ZERO, u.getSize(), params.initialVelocity.x );
	v.set( MultiIndex::ZERO, v.getSize(), params.initialVelocity.y );
	p.set( MultiIndex::ZERO, p.getSize(), params.initialPressure );
}

real Computation::computeTimeStep()
{
	real hx2 = h.x * h.x;
	real hy2 = h.y * h.y;

	real uMax = u.getMaxValue( MultiIndex( 0, 0 ), u.getSize() );
	real vMax = v.getMaxValue( MultiIndex( 0, 0 ), v.getSize() );

	return m_Params.tau * std::min( m_Params.Re / 2.0 * hx2 * hy2 /
			( hx2 + hy2 ), std::min( h.x / uMax, h.y / vMax ) );
}

void Computation::computeNewVelocities( real dt )
{
	// ------------------------- compute u -------------------------

	MultiIndex uWriteEnd = u.getSize() - MultiIndex::ONE;
	// u = dp/dx
	m_stDxBw.apply( p, MultiIndex( 2, 1 ), p.getSize() - MultiIndex::ONE,
			u, MultiIndex::ONE, uWriteEnd );
	// u = -dt * dp/dx
	u.scale( MultiIndex::ONE, uWriteEnd, -dt );
	// u = f - dt * dp/dx
	u.addScaledCopy( MultiIndex::ONE, uWriteEnd, f, 1.0 );

	// ------------------------- compute v -------------------------

	MultiIndex vWriteEnd = v.getSize() - MultiIndex::ONE;
	// v = dp/dy
	m_stDyBw.apply( p, MultiIndex( 1, 2 ), p.getSize() - MultiIndex::ONE,
			v, MultiIndex::ONE, vWriteEnd );
	// v = -dt * dp/dy
	v.scale( MultiIndex::ONE, vWriteEnd, -dt );
	// v = g - dt * dp/dy
	v.addScaledCopy( MultiIndex::ONE, vWriteEnd, g, 1.0 );
}

void Computation::computeMomentumEquations( real dt )
{
	{
	// ========================= compute f =========================
	MultiIndex uReadEnd = u.getSize() - MultiIndex::ONE;

	// ------------ compute 1/Re*(d^2u/dx^2+d^2u/dy^2)  ------------
	// f = d^2u/dx^2
	m_stDxx.apply( u, MultiIndex::ONE, uReadEnd,
			f, MultiIndex::ONE, uReadEnd );
	// tmp1 = d^2u/dy^2
	m_stDyy.apply( u, MultiIndex::ONE, uReadEnd,
			m_funcTemp1, MultiIndex::ONE, uReadEnd );
	// f = f + tmp1 = d^2u/dx^2 + d^2u/dy^2
	f.addScaledCopy( MultiIndex::ONE, uReadEnd, m_funcTemp1, 1.0f );
	// f = f * 1 / Re = 1 / Re * ( d^2u/dx^2 + d^2u/dy^2 )
	f.scale( MultiIndex::ONE, uReadEnd, 1.0 / m_Params.Re );

	// --------------------- compute d(u^2)/dx ---------------------
	MultiIndex tmpBegin = MultiIndex( 0, 1 );
	// tmp1 = (u_(i+1) + u_i)/2
	m_stInterpx.apply( u, tmpBegin, uReadEnd,
			m_funcTemp1, tmpBegin, uReadEnd );
	// tmp2 = ((u_(i+1) + u_i)/2)^2 = tmp1 * tmp1
	m_funcTemp2.copy( tmpBegin, uReadEnd, m_funcTemp1 );
	m_funcTemp2.multiply( tmpBegin, uReadEnd, m_funcTemp2 );
	// tmp3 = d/dx tmp2
	m_stDxBw.apply( m_funcTemp2, MultiIndex::ONE, uReadEnd,
			m_funcTemp3, MultiIndex::ONE, uReadEnd );

	// tmp1 = |u_(i+1) + u_i|/2 = |tmp1|
	m_funcTemp1.abs( tmpBegin, uReadEnd );
	// tmp2 = (u_i-u_(i+1))/2 = -du/dx * dx
	m_stDxBw.apply( u, MultiIndex::ONE, u.getSize() - tmpBegin,
			m_funcTemp2, tmpBegin, uReadEnd );
	m_funcTemp2.scale( tmpBegin, uReadEnd, -h.x );
	// tmp1 = tmp1 * tmp2 = |u_(i+1) + u_i|/2 * (u_i-u_(i+1))/2
	m_funcTemp1.multiply( tmpBegin, uReadEnd, m_funcTemp2 );
	// tmp2 = d/dx tmp1
	m_stDxBw.apply( m_funcTemp1, MultiIndex::ONE, uReadEnd,
			m_funcTemp2, MultiIndex::ONE, uReadEnd );
	// tmp3 = tmp3 + alpha * tmp2 = d(u^2)/dx
	m_funcTemp3.addScaledCopy( MultiIndex::ONE, uReadEnd, m_funcTemp2,
			m_Params.alpha );
	// f = f + tmp3 = 1 / Re * ( d^2u/dx^2 + d^2u/dy^2 ) - d(u^2)/dx
	f.addScaledCopy( MultiIndex::ONE, uReadEnd, m_funcTemp3, -1.0 );

	// --------------------- compute d(u*v)/dy ---------------------
	tmpBegin = MultiIndex( 1, 0 );
	// tmp1 = (v_i + v_(i+1))/2
	m_stInterpx.apply( v, tmpBegin, uReadEnd,
			m_funcTemp1, tmpBegin, uReadEnd );
	// tmp2 = (u_j + u_(j+1))/2
	m_stInterpy.apply( u, tmpBegin, uReadEnd,
			m_funcTemp2, tmpBegin, uReadEnd );
	// tmp2 = tmp2 * tmp1 = (v_i + v_(i+1))/2 * (u_j + u_(j+1))/2
	m_funcTemp2.multiply( tmpBegin, uReadEnd, m_funcTemp1 );
	// tmp3 = d/dy tmp2
	m_stDyBw.apply( m_funcTemp2, MultiIndex::ONE, uReadEnd,
			m_funcTemp3, MultiIndex::ONE, uReadEnd );

	// tmp1 = |v_i + v_(i+1)|/2 = |tmp1|
	m_funcTemp1.abs( tmpBegin, uReadEnd );
	// tmp2 = (u_j-u_(j+1))/2 = -du/dy * dy
	m_stDyBw.apply( u, MultiIndex::ONE, u.getSize() - tmpBegin,
			m_funcTemp2, tmpBegin, uReadEnd );
	m_funcTemp2.scale( tmpBegin, uReadEnd, -h.y );
	// tmp1 = tmp1 * tmp2 = |v_i + v_(i+1)|/2 * (u_j-u_(j+1))/2
	m_funcTemp1.multiply( tmpBegin, uReadEnd, m_funcTemp2 );
	// tmp2 = d/dy tmp1
	m_stDyBw.apply( m_funcTemp1, MultiIndex::ONE, uReadEnd,
			m_funcTemp2, MultiIndex::ONE, uReadEnd );
	// tmp3 = tmp3 + alpha * tmp2 = d(u*v)/dy
	m_funcTemp3.addScaledCopy( MultiIndex::ONE, uReadEnd, m_funcTemp2,
			m_Params.alpha );
	// f = f + tmp3 = 1 / Re * ( d^2u/dx^2 + d^2u/dy^2 ) - d(u^2)/dx - d(u*v)/dy
	f.addScaledCopy( MultiIndex::ONE, uReadEnd, m_funcTemp3, -1.0 );

	// ---------------------- compute final f ----------------------
	f.scale( MultiIndex::ONE, uReadEnd, dt );
	f.addScaledCopy( MultiIndex::ONE, uReadEnd, u, 1.0 );
	}
	{
	// ========================= compute g =========================
	MultiIndex vReadEnd = v.getSize() - MultiIndex::ONE;
	// ------------ compute 1/Re*(d^2v/dx^2+d^2v/dy^2)  ------------
	// g = d^2v/dx^2
	m_stDxx.apply( v, MultiIndex::ONE, vReadEnd,
			g, MultiIndex::ONE, vReadEnd );
	// tmp1 = d^2v/dy^2
	m_stDyy.apply( v, MultiIndex::ONE, vReadEnd,
			m_funcTemp1, MultiIndex::ONE, vReadEnd );
	// g = g + tmp1 = d^2v/dx^2 + d^2v/dy^2
	g.addScaledCopy( MultiIndex::ONE, vReadEnd, m_funcTemp1, 1.0f );
	// g = g * 1 / Re = 1 / Re * ( d^2u/dx^2 + d^2u/dy^2 )
	g.scale( MultiIndex::ONE, vReadEnd, 1.0 / m_Params.Re );

	// --------------------- compute d(v^2)/dy ---------------------
	MultiIndex tmpBegin = MultiIndex( 1, 0 );
	// tmp1 = (v_(j+1) + v_j)/2
	m_stInterpy.apply( v, tmpBegin, vReadEnd,
			m_funcTemp1, tmpBegin, vReadEnd );
	// tmp2 = ((v_(j+1) + v_j)/2)^2 = tmp1 * tmp1
	m_funcTemp2.copy( tmpBegin, vReadEnd, m_funcTemp1 );
	m_funcTemp2.multiply( tmpBegin, vReadEnd, m_funcTemp2 );
	// tmp3 = d/dy tmp2
	m_stDyBw.apply( m_funcTemp2, MultiIndex::ONE, vReadEnd,
			m_funcTemp3, MultiIndex::ONE, vReadEnd );

	// tmp1 = |v_(j+1) + v_j|/2 = |tmp1|
	m_funcTemp1.abs( tmpBegin, vReadEnd );
	// tmp2 = (v_j-v_(j+1))/2 = -dv/dy * dy
	m_stDyBw.apply( v, MultiIndex::ONE, v.getSize() - tmpBegin,
			m_funcTemp2, tmpBegin, vReadEnd );
	m_funcTemp2.scale( tmpBegin, vReadEnd, -h.y );
	// tmp1 = tmp1 * tmp2 = |v_(j+1) + v_j|/2 * (v_j-v_(j+1))/2
	m_funcTemp1.multiply( tmpBegin, vReadEnd, m_funcTemp2 );
	// tmp2 = d/dy tmp1
	m_stDyBw.apply( m_funcTemp1, MultiIndex::ONE, vReadEnd,
			m_funcTemp2, MultiIndex::ONE, vReadEnd );
	// tmp3 = tmp3 + alpha * tmp2 = d(v^2)/dy
	m_funcTemp3.addScaledCopy( MultiIndex::ONE, vReadEnd, m_funcTemp2,
			m_Params.alpha );
	// g = g + tmp3 = 1 / Re * ( d^2v/dy^2 + d^2v/dx^2 ) - d(v^2)/dy
	g.addScaledCopy( MultiIndex::ONE, vReadEnd, m_funcTemp3, -1.0 );

	// --------------------- compute d(u*v)/dx ---------------------
	tmpBegin = MultiIndex( 0, 1 );
	// tmp1 = (u_j + u_(j+1))/2
	m_stInterpy.apply( u, tmpBegin, vReadEnd,
			m_funcTemp1, tmpBegin, vReadEnd );
	// tmp2 = (v_i + v_(i+1))/2
	m_stInterpx.apply( v, tmpBegin, vReadEnd,
			m_funcTemp2, tmpBegin, vReadEnd );
	// tmp2 = tmp2 * tmp1 = (u_j + u_(j+1))/2 * (v_i + v_(i+1))/2
	m_funcTemp2.multiply( tmpBegin, vReadEnd, m_funcTemp1 );
	// tmp3 = d/dx tmp2
	m_stDxBw.apply( m_funcTemp2, MultiIndex::ONE, vReadEnd,
			m_funcTemp3, MultiIndex::ONE, vReadEnd );

	// tmp1 = |u_j + u_(j+1)|/2 = |tmp1|
	m_funcTemp1.abs( tmpBegin, vReadEnd );
	// tmp2 = (v_i-v_(i+1))/2 = -dv/dx * dx
	m_stDxBw.apply( v, MultiIndex::ONE, v.getSize() - tmpBegin,
			m_funcTemp2, tmpBegin, vReadEnd );
	m_funcTemp2.scale( tmpBegin, vReadEnd, -h.x );
	// tmp1 = tmp1 * tmp2 = |u_j + u_(j+1)|/2 * (v_i-v_(i+1))/2
	m_funcTemp1.multiply( tmpBegin, vReadEnd, m_funcTemp2 );
	// tmp2 = d/dx tmp1
	m_stDxBw.apply( m_funcTemp1, MultiIndex::ONE, vReadEnd,
			m_funcTemp2, MultiIndex::ONE, vReadEnd );
	// tmp3 = tmp3 + alpha * tmp2 = d(v*u)/dx
	m_funcTemp3.addScaledCopy( MultiIndex::ONE, vReadEnd, m_funcTemp2,
			m_Params.alpha );
	// g = g + tmp3 = 1 / Re * ( d^2v/dy^2 + d^2v/dx^2 ) - d(v^2)/dy - d(v*u)/dx
	g.addScaledCopy( MultiIndex::ONE, vReadEnd, m_funcTemp3, -1.0 );

	// ---------------------- compute final g ----------------------
	g.scale( MultiIndex::ONE, vReadEnd, dt );
	g.addScaledCopy( MultiIndex::ONE, vReadEnd, v, 1.0 );
	}
}

void Computation::computeRightHandSide( real dt )
{
	// rhs = df/dx
	m_stDxBw.apply( f, MultiIndex::ONE, f.getSize() - MultiIndex( 0, 1 ),
			rhs, MultiIndex::ZERO, rhs.getSize() );
	// tmp1 = dg/dy
	m_stDyBw.apply( g, MultiIndex::ONE, g.getSize() - MultiIndex( 1, 0 ),
			m_funcTemp1, MultiIndex::ZERO, rhs.getSize() );
	// rhs = rhs + tmp1 = (df/dx + dg/dy)
	rhs.addScaledCopy( MultiIndex::ZERO, rhs.getSize(), m_funcTemp1, 1.0 );
	// rhs = 1/dt * rhs = 1/dt * (df/dx + dg/dy)
	rhs.scale( MultiIndex::ZERO, rhs.getSize(), 1 / dt );
}

void Computation::setVelocityBoundary( GridFunction& u, GridFunction& v )
{
	static const real FREE_FLOW_VELOCITY = 1;

	// left
	u.set( MultiIndex::ZERO, MultiIndex( 1, u.getSize().y ), 0.0 );
	// right
	u.set( MultiIndex( u.getSize().x - 1, 1 ), u.getSize(), 0.0 );
	// bottom
	u.scaledOffsetCopy( MultiIndex( 1, 0 ), MultiIndex( u.getSize().x - 1, 1 ),
			u, MultiIndex( 0, 1 ), -1.0 );
	// top
	u.affineOffsetCopy( MultiIndex( u.getSize().y - 1, 1 ),
			u.getSize() - MultiIndex( 1, 0 ), u, MultiIndex( 0, -1 ), -1.0,
			2.0 * FREE_FLOW_VELOCITY );

	// left
	v.scaledOffsetCopy( MultiIndex( 0, 1 ), MultiIndex( 1, v.getSize().y - 1 ),
			v, MultiIndex( 1, 0 ), -1.0 );
	// right
	v.scaledOffsetCopy( MultiIndex( v.getSize().x - 1, 1 ),
			v.getSize() - MultiIndex( 0, 1 ), v, MultiIndex( -1, 0 ), -1.0 );
	// bottom
	v.set( MultiIndex::ZERO, MultiIndex( v.getSize().x, 1 ), 0.0 );
	// top
	v.set( MultiIndex( 0, v.getSize().y - 1 ), v.getSize(), 0.0 );
}

void Computation::setPressureBoundary()
{
	// left
	p.scaledOffsetCopy( MultiIndex( 0, 1 ), MultiIndex( 1, p.getSize().y - 1 ),
			p, MultiIndex( 1, 0 ), 1.0 );
	// right
	p.scaledOffsetCopy( MultiIndex( p.getSize().x - 1, 1 ),
			p.getSize() - MultiIndex( 0, 1 ), p, MultiIndex( -1, 0 ), 1.0 );
	// bottom
	p.scaledOffsetCopy( MultiIndex( 1, 0 ), MultiIndex( p.getSize().x - 1, 1 ),
			p, MultiIndex( 0, 1 ), 1.0 );
	// top
	p.scaledOffsetCopy( MultiIndex( 1, p.getSize().y - 1 ),
			p.getSize() - MultiIndex( 1, 0 ), p, MultiIndex( 0, -1 ), 1.0 );
}

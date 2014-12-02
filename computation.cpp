#include "computation.h"
#include <algorithm>
#include "metamath/metamath.h"
#include "metamath/mmutils.h"
#include "communication.h"

using namespace mm;

real Computation::computeTimeStep( const GridFunction& u,
		const GridFunction& v,
		const Point& h,
		real Re,
		real tau )
{
	real hx2 = h.x * h.x;
	real hy2 = h.y * h.y;

	real uMax = max( abs( u ), MultiIndex::ZERO, u.size() );
	real vMax = max( abs( v ), MultiIndex::ZERO, v.size() );

	real dt_local = tau * std::min( Re / 2.0 * hx2 * hy2 / ( hx2 + hy2 ),
			std::min( h.x / uMax, h.y / vMax ) );

	return Communication::min( dt_local );
}

void Computation::computeNewVelocities( GridFunction& u,
		GridFunction& v,
		const GridFunction& f,
		const GridFunction& g,
		const GridFunction& p,
		const Point& h,
		real dt )
{
	// ------------------------- compute u -------------------------
	set( u, MultiIndex::ONE, u.size() - MultiIndex::ONE,
			f - dt % utils::diffX_Bw( p, h.x ) );

	// ------------------------- compute v -------------------------
	set( v, MultiIndex::ONE, v.size() - MultiIndex::ONE,
			g - dt % utils::diffY_Bw( p, h.y ) );
}

void Computation::computeMomentumEquations( GridFunction& f,
		GridFunction& g,
		const GridFunction& u,
		const GridFunction& v,
		const Point& h,
		real dt,
		real Re,
		real alpha )
{
	// ========================= compute f =========================
	// --------------------- compute d(u^2)/dx ---------------------
	auto u2x = utils::diffX_Bw( sqr( 0.5 % ( u + eval<+1,0>( u ) ) ), h.x )
		+ alpha % utils::diffX_Bw( ( 0.5 % abs( u + eval<+1,0>( u ) ) ) * ( 0.5 % ( u - eval<+1,0>( u ) ) ), h.x );

	// --------------------- compute d(u*v)/dy ---------------------
	auto uvy = utils::diffY_Bw( ( 0.5 % ( eval<-1,+1>( v ) + eval<0,+1>( v ) ) ) * ( 0.5 % ( u + eval<0,+1>( u ) ) ), h.y )
		 + alpha % utils::diffY_Bw( ( 0.5 % abs( eval<-1,+1>( v ) + eval<0,+1>( v ) ) ) * ( 0.5 % ( u - eval<0,+1>( u ) ) ), h.y );

	// ---------------------- compute final f ----------------------
	set( f, MultiIndex::ONE, f.size() - MultiIndex::ONE,
			u + dt % ( 1.0 / Re % utils::diffXX_YY( u, h )
				- u2x - uvy ) );

	// ========================= compute g =========================
	// --------------------- compute d(v^2)/dy ---------------------
	auto v2y = utils::diffY_Bw( sqr( 0.5 % ( v + eval<0,+1>( v ) ) ), h.y )
		+ alpha % utils::diffY_Bw( ( 0.5 % abs( v + eval<0,+1>( v ) ) ) * ( 0.5 % ( v - eval<0,+1>( v ) ) ), h.y );

	// --------------------- compute d(u*v)/dx ---------------------
	auto uvx = utils::diffX_Bw( ( 0.5 % ( eval<+1,-1>( u ) + eval<+1,0>( u ) ) ) * ( 0.5 % ( v + eval<+1,0>( v ) ) ), h.x )
		 + alpha % utils::diffX_Bw( ( 0.5 % abs( eval<+1,-1>( u ) + eval<+1,0>( u ) ) ) * ( 0.5 % ( v - eval<+1,0>( v ) ) ), h.x );

	// ---------------------- compute final g ----------------------
	set( g, MultiIndex::ONE, g.size() - MultiIndex::ONE,
			v + dt % ( 1.0 / Re % utils::diffXX_YY( v, h )
				- v2y - uvx ) );
}

void Computation::computeRightHandSide( GridFunction& rhs,
		const GridFunction& f,
		const GridFunction& g,
		const Point& h,
		real dt )
{
	// das ist Backward stencil weil in der Formel es so steht, ist aber trotzdem
	// eine *zentrale* Differenz, da wir entsprechend staggered sind beim Druck

	// rhs = 1/dt * (df/dx + dg/dy)
	rhs = 1.0 / dt % ( utils::diffX_Fw( eval<+1,+1>( f ), h.x )
			+ utils::diffY_Fw( eval<+1,+1>( g ), h.y ) );
}

void Computation::setVelocityBoundary( GridFunction& u, GridFunction& v,
		bool left, bool top, bool right, bool bottom )
{
	static const real FREE_FLOW_VELOCITY = 1;

	// left
	if( left )
	{
		set( u, MultiIndex( 1, 0 ), MultiIndex( 2, u.size().y - 1 ),
				constant( 0.0 ) );
		set( v, MultiIndex::ZERO, MultiIndex( 1, v.size().y ),
				-eval<+1,0>( v ) );
	}
	// right
	if( right )
	{
		set( u, MultiIndex( u.size().x - 2, 1 ), u.size() - MultiIndex::ONE,
				constant( 0.0 ) );
		set( v, MultiIndex( v.size().x - 1, 0 ), v.size(),
				-eval<-1,0>( v ) );
	}
	// bottom
	if( bottom )
	{
		set( u, MultiIndex::ZERO, MultiIndex( u.size().x, 1 ),
				-eval<0,+1>( u ) );
		set( v, MultiIndex::ONE, MultiIndex( v.size().x - 1, 2 ),
				constant( 0.0 ) );
	}
	// top
	if( top )
	{
		set( u, MultiIndex( 0, u.size().y - 1 ), u.size(),
				constant( 2.0 * FREE_FLOW_VELOCITY ) - eval<0,-1>( u ) );
		set( v, MultiIndex( 1, v.size().y - 2 ), v.size() - MultiIndex::ONE,
				constant( 0.0 ) );
	}
}

void Computation::setPressureBoundary( GridFunction& p,
		bool left, bool top, bool right, bool bottom )
{
	// left
	if( left )
	{
		set( p, MultiIndex( 0, 1 ), MultiIndex( 1, p.size().y - 1 ),
				eval<+1,0>( p ) );
	}
	// right
	if( right )
	{
		set( p, MultiIndex( p.size().x - 1, 1 ), p.size() - MultiIndex( 0, 1 ),
				eval<-1,0>( p ) );
	}
	// bottom
	if( bottom )
	{
		set( p, MultiIndex( 1, 0 ), MultiIndex( p.size().x - 1, 1 ),
				eval<0,+1>( p ) );
	}
	// top
	if( top )
	{
		set( p, MultiIndex( 1, p.size().y - 1 ), p.size() - MultiIndex( 1, 0 ),
				eval<0,-1>( p ) );
	}
}

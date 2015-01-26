#include "computation.h"
#include <algorithm>
#include "metamath/metamath.h"
#include "metamath/mmutils.h"
#include "communication.h"
#include <cmath>

#include <iostream>

using namespace mm;

real Computation::computeTimeStep( const GridFunction& u,
		const GridFunction& v,
		const Point& h,
		real Re,
		real Pr,
		real tau )
{
	static const real EPS = 1.0 / ( 1 << 20 );

	real hx2 = h.x * h.x;
	real hy2 = h.y * h.y;

	real uMax = std::max( max( abs( u ), MultiIndex::ZERO, u.size() ), EPS );
	real vMax = std::max( max( abs( v ), MultiIndex::ZERO, v.size() ), EPS );
	real tLim = 0.5 * Re * Pr / ( 1.0 / hx2 + 1 / hy2 );

	real dt_local = tau * std::min( tLim, std::min( h.x / uMax, h.y / vMax ) );

	return Communication::min( dt_local );
}

void Computation::computeNewVelocities( GridFunction& u,
		GridFunction& v,
		const GridFunction& f,
		const GridFunction& g,
		const GridFunction& p,
		const MaskFunction& uMask,
		const MaskFunction& vMask,
		const Point& h,
		real dt )
{
	// ------------------------- compute u -------------------------
	utils::setMasked( u, MultiIndex::ONE, u.size() - MultiIndex::ONE, uMask,
			f - dt * utils::diffX_Bw( p, h.x ) );

	// ------------------------- compute v -------------------------
	utils::setMasked( v, MultiIndex::ONE, v.size() - MultiIndex::ONE, vMask,
			g - dt * utils::diffY_Bw( p, h.y ) );
}

void Computation::computeMomentumEquations( GridFunction& f,
		GridFunction& g,
		const GridFunction& u,
		const GridFunction& v,
		const GridFunction& T,
		const Point& vf,
		const MaskFunction& uMask,
		const MaskFunction& vMask,
		const Point& h,
		real dt,
		real Re,
		real alpha,
		real beta )
{
	// ========================= compute f =========================
	// --------------------- compute d(u^2)/dx ---------------------
	auto u2x = utils::diffX_Bw( sqr( 0.5 * ( u + eval<+1,0>( u ) ) ), h.x )
		+ alpha * utils::diffX_Bw( ( 0.5 * abs( u + eval<+1,0>( u ) ) ) * ( 0.5 * ( u - eval<+1,0>( u ) ) ), h.x );

	// --------------------- compute d(u*v)/dy ---------------------
	auto uvy = utils::diffY_Bw( ( 0.5 * ( eval<-1,+1>( v ) + eval<0,+1>( v ) ) ) * ( 0.5 * ( u + eval<0,+1>( u ) ) ), h.y )
		 + alpha * utils::diffY_Bw( ( 0.5 * abs( eval<-1,+1>( v ) + eval<0,+1>( v ) ) ) * ( 0.5 * ( u - eval<0,+1>( u ) ) ), h.y );

	auto Tvfx = ( constant( 1.0 ) - beta * ( eval<-1,0>( T ) + T ) ) * vf.x;

	// ---------------------- compute final f ----------------------
	utils::setMasked( f, MultiIndex::ONE, f.size() - MultiIndex::ONE, uMask,
			u + dt * ( 1.0 / Re * utils::diffXX_YY( u, h )
				- u2x - uvy + Tvfx ) );

	// ========================= compute g =========================
	// --------------------- compute d(v^2)/dy ---------------------
	auto v2y = utils::diffY_Bw( sqr( 0.5 * ( v + eval<0,+1>( v ) ) ), h.y )
		+ alpha * utils::diffY_Bw( ( 0.5 * abs( v + eval<0,+1>( v ) ) ) * ( 0.5 * ( v - eval<0,+1>( v ) ) ), h.y );

	// --------------------- compute d(u*v)/dx ---------------------
	auto uvx = utils::diffX_Bw( ( 0.5 * ( eval<+1,-1>( u ) + eval<+1,0>( u ) ) ) * ( 0.5 * ( v + eval<+1,0>( v ) ) ), h.x )
		 + alpha * utils::diffX_Bw( ( 0.5 * abs( eval<+1,-1>( u ) + eval<+1,0>( u ) ) ) * ( 0.5 * ( v - eval<+1,0>( v ) ) ), h.x );

	auto Tvfy = ( constant( 1.0 ) - beta * ( eval<0,-1>( T ) + T ) ) * vf.y;

	// ---------------------- compute final g ----------------------
	utils::setMasked( g, MultiIndex::ONE, g.size() - MultiIndex::ONE, vMask,
			v + dt * ( 1.0 / Re * utils::diffXX_YY( v, h )
				- v2y - uvx + Tvfy ) );
}

void Computation::computeTemperatureEquations( GridFunction& T,
		const GridFunction& TT,
		const GridFunction& u,
		const GridFunction& v,
		const MaskFunction& pMask,
		const Point& h,
		real dt,
		real Re,
		real Pr,
		real gamma )
{
	auto uu = eval<+1,0>( u );
	auto vv = eval<0,+1>( v );

	// ========================= compute T =========================
	// --------------------- compute d(uT)/dx ---------------------
	auto uTx = utils::diffX_Bw( uu * 0.5 * ( eval<+1,0>( TT ) + TT )
		+ gamma * ( abs( uu ) * 0.5 * ( TT - eval<+1,0>( TT ) ) ), h.x );

	// --------------------- compute d(vT)/dy ---------------------
	auto vTy = utils::diffY_Bw( vv * 0.5 * ( eval<0,+1>( TT ) + TT )
		+ gamma * ( abs( vv ) * 0.5 * ( TT - eval<0,+1>( TT ) ) ), h.y );

	// ---------------------- compute final f ----------------------
	utils::setMasked( T, MultiIndex::ONE, T.size() - MultiIndex::ONE, pMask,
			TT + dt * ( 1.0 / ( Re * Pr ) * utils::diffXX_YY( TT, h )
				- uTx - vTy ) );
}

void Computation::computeRightHandSide( GridFunction& rhs,
		const GridFunction& f,
		const GridFunction& g,
		const MaskFunction& pMask,
		const Point& h,
		real dt )
{
	// das ist Backward stencil weil in der Formel es so steht, ist aber trotzdem
	// eine *zentrale* Differenz, da wir entsprechend staggered sind beim Druck

	// rhs = 1/dt * (df/dx + dg/dy)
	utils::setMasked( rhs, eval<+1,+1>( pMask ),
			1.0 / dt * ( utils::diffX_Fw( eval<+1,+1>( f ), h.x )
				+ utils::diffY_Fw( eval<+1,+1>( g ), h.y ) ) );
}

void Computation::computeVelocityPotential( GridFunction& psi,
		const GridFunction& u,
		const GridFunction& v,
		const Point& h )
{
	// use d(psi)/dx=-v instead of d(psi)/dy=u for more efficient memory access patterns
	//  => psi_ij = psi_(i-1)j - v_ij*dx, psi_0j = <boundary value>

	// "HACK": use psi as source and image function to accumulate values
	set( psi, MultiIndex( 1, 0 ), psi.size(),
			eval<-1,0>( psi ) - h.x * eval<0,+1>( v ) );
}

void Computation::computeVorticity( GridFunction& zeta,
		const GridFunction& u,
		const GridFunction& v,
		const Point& h )
{
	zeta = utils::diffY_Fw( eval<+2,+1>( u ), h.y ) - utils::diffX_Fw( eval<+1,+2>( v ), h.x );
}

void Computation::setVelocityBoundary( GridFunction& u, GridFunction& v,
		bool left, bool top, bool right, bool bottom )
{
	static const real FREE_FLOW_VELOCITY = 1;

	// left
	if( left )
	{
		set( u, MultiIndex( 1, 1 ), MultiIndex( 2, u.size().y - 1 ),
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

#ifndef PARAMS_H
#define PARAMS_H

#include "typedef.h"

struct Params
{
	Point      domainSize;
	MultiIndex gridSize;
	real       T;
	real       dt;
	real       tau;
	real       deltaVec;
	index_t    maxIter;
	real       eps;
	real       omega;
	real       alpha;
	real       Re;
	Point      initialVelocity;
	real       initialPressure;
};

#endif

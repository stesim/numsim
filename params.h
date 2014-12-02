#ifndef PARAMS_H
#define PARAMS_H

#include "typedef.h"
#include <iostream>

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
  
  void print( std::ostream& s ) const // man hätte auch den operator<< für Params überladen können...
  {
    s
      << "  Domain size:           "            << domainSize.x       << " x " << domainSize.y << std::endl
      << "  Grid size:             "            << gridSize.x         << " x " << gridSize.y   << std::endl
      << "  T:                     "            << T                  <<                          std::endl
      << "  Re:                    "            << Re                 <<                          std::endl
      << "  dt:                    "            << dt                 <<                          std::endl
      << "  tau:                   "            << tau                <<                          std::endl
      << "  Output time step:      "            << deltaVec           <<                          std::endl
      << "  Max SOR iterations:    "            << maxIter            <<                          std::endl
      << "  SOR error tolerance:   "            << eps                <<                          std::endl
      << "  Relaxation factor:     "            << omega              <<                          std::endl
      << "  Upwinding factor:      "            << alpha              <<                          std::endl
      << "  Initial velocity:      ("           << initialVelocity.x  << "," << initialVelocity.y << ")" << std::endl
      << "  Initial pressure:      "            << initialPressure    << std::endl;
  }
};

#endif

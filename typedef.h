#ifndef TYPEDEF_H
#define TYPEDEF_H

#include "metamath/mmfunction.h"

/*! Creates a type name for RealType */
typedef double real;

/*! Creates a type name for IndexType */
typedef int index_t;

/*! Creates a type name for MultiIndexType */
typedef mm::Tuple<index_t, 2> MultiIndex;

/*! Creates a type name for PointType */
typedef mm::Tuple<real, 2> Point;

typedef mm::Function<real, 2> GridFunction;

typedef mm::Function<bool, 2> MaskFunction;

#endif

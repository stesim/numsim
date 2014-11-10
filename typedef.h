//! Typedefs for the apllication 
/*!
 * @author diehlpk
 * @date 2012
 */

#ifndef TYPEDEF
#define TYPEDEF

#include "template.h"

/*! Creates a type name for RealType */
typedef double real;

/*! Creates a type name for IndexType */
typedef int index_t;

/*! Creates a type name for MultiIndexType */
typedef Tuple2<index_t> MultiIndex;

/*! Creates a type name for PointType */
typedef Tuple2<real> Point;

#endif

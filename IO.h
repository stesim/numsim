//! The class implements the IO 
/*!
 * @author diehlpk
 * @date 2012
 */

#ifndef IO_H_
#define IO_H_

#include <iostream>
#include <fstream>
#include "typedef.h"
#include "params.h"

class IO
{
public:
  //! Method writes the GridFunctions u,v,p in the vtk data format to the hard disk. 
  //! The files are named in the following convention: field_(instance)_(step).vts
  /*!
   * \param griddimension The dimension of the gridfunctions.
   * \param u The gridfunction u.
   * \param v The gridfunction v.
   * \param p The gridfunction p.
   * \param delta The mesh width in x direction and y direction
   * \param params The input parameters.
   * \param instance The number of instance.
   * \param step The number of the timestep.
   */
  static void writeVTKFile (const MultiIndex & griddimension,
		     const GridFunction & u, const GridFunction & v,
		     const GridFunction & p, const Point & delta, const Params& params,
			 int instance, int step);

  //! Method writes the GridFunctions u,v,p in a text data format to the hard disk. 
  //! The files are named in the following convention: field-(instance)-(step).txt
  /*!
   * \param griddimension The dimension of the gridfunctions.
   * \param u The gridfunction u.
   * \param v The gridfunction v.
   * \param p The gridfunction p.
   * \param delta The mesh width in x direction and y direction
   * \param instance The number of instance.
   * \param step The number of the timestep.
   */
  static void writeRawOutput( const GridFunction& u, const GridFunction& v,
		  const GridFunction& p, const Point& delta, const Point& offset,
		  index_t step, const MultiIndex& rank );

  /*!
   * Methods reads the simulation parameters from the specified input file.
   * 
   * @param filename The name of the file with the simulations paremters
   */
  static Params readInputFile (const char *filename);


  //! Method interpolates the velocity for u in the staggered grid.
  /*!
   * \param x Value of the x coordinate.
   * \param y Value of the y coordinate.
   * \param u The gridfunction u.
   * \param delta The mesh width in x direction and y direction.
   */
  static real interpolateVelocityU (real x, real y, const GridFunction & u,
				 const Point & delta);

  //! Method interpolates the velocity for v in the staggered grid.
  /*!
   * \param x Value of the x coordinate.
   * \param y Value of the y coordinate.
   * \param u The gridfunction v.
   * \param delta The mesh width in x direction and y direction.
   */
  static real interpolateVelocityV (real x, real y, const GridFunction & v,
				 const Point & delta);
};

#endif

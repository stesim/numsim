#include "IO.h"
#include <cstring>
#include <string>
#include <cstdlib>

using std::size_t;

Params
IO::readInputFile (const char *filename)
{
  Params p;

  std::ifstream stream( filename );

  std::string line;
  while( std::getline( stream, line ) )
  {
	  size_t beginName = line.find_first_not_of( " \t", 0 );
	  if( beginName >= std::string::npos )
	  {
		  continue;
	  }
	  size_t endName = line.find_first_of( " \t=", beginName );
	  if( endName >= std::string::npos )
	  {
		  continue;
	  }
	  size_t beginValue = line.find_first_not_of( " \t=", endName );
	  if( beginValue >= std::string::npos )
	  {
		  continue;
	  }
	  size_t endValue = line.find_first_of( " \t\n\0", beginValue );

	  std::string name = line.substr( beginName, endName - beginName );
	  std::string value = line.substr( beginValue, endValue - beginValue );

#define SET_PARAM_DOUBLE(_member, _name) if( name == _name ) {\
	  p._member = strtod( value.c_str(), NULL ); }
#define SET_PARAM_INT(_member, _name) if( name == _name ) {\
	  p._member = strtol( value.c_str(), NULL, 10 ); }

	  SET_PARAM_DOUBLE(domainSize.x, "xLength")
	  else SET_PARAM_DOUBLE(domainSize.y, "yLength")
	  else SET_PARAM_INT(gridSize.x, "iMax")
	  else SET_PARAM_INT(gridSize.y, "jMax")
	  else SET_PARAM_DOUBLE(T, "tEnd")
	  else SET_PARAM_DOUBLE(dt, "deltaT")
	  else SET_PARAM_DOUBLE(tau, "tau")
	  else SET_PARAM_DOUBLE(deltaVec, "deltaVec")
	  else SET_PARAM_INT(maxIter, "iterMax")
	  else SET_PARAM_DOUBLE(eps, "eps")
	  else SET_PARAM_DOUBLE(omega, "omg")
	  else SET_PARAM_DOUBLE(alpha, "alpha")
	  else SET_PARAM_DOUBLE(Re, "re")
	  else SET_PARAM_DOUBLE(initialVelocity.x, "ui")
	  else SET_PARAM_DOUBLE(initialVelocity.y, "vi")
	  else SET_PARAM_DOUBLE(initialPressure, "pi")

#undef SET_PARAM_DOUBLE
#undef SET_PARAM_INT
  }

  return p;
}


#define Element(field,ic) ((field)((ic)[0],(ic)[1]))

real
  IO::interpolateVelocityU (real x, real y, const GridFunction & u,
			    const Point & delta)
{

  real deltaX = delta[0];
  real deltaY = delta[1];

  MultiIndex index;

  // Computation of u(x,y)
  index[0] = ((int) (x / deltaX)) + 1;
  index[1] = ((int) ((y + (deltaY / 2)) / deltaY)) + 1;

  // The coordinates of the cell corners

  real x1 = (index[0] - 1) * deltaX;
  real x2 = index[0] * deltaX;
  real y1 = ((index[1] - 1) - 0.5) * deltaY;
  real y2 = (index[1] - 0.5) * deltaY;

  MultiIndex offset;

  offset[0] = index[0] - 1;
  offset[1] = index[1] - 1;

  real u1 = Element (u, offset);	// datafields->u->getField ()[i - 1][j - 1];

  offset[0] = index[0];
  offset[1] = index[1] - 1;

  real u2 = Element (u, offset);	//datafields->u->getField ()[i][j - 1];

  offset[0] = index[0] - 1;
  offset[1] = index[1];

  real u3 = Element (u, offset);	//datafields->u->getField ()[i - 1][j];
  real u4 = Element (u, index);

  real
    uInterploated =
    (1.0 / (deltaX * deltaY)) * ((x2 - x) * (y2 - y) *
				 u1 + (x - x1) * (y2 -
						  y) *
				 u2 + (x2 - x) * (y -
						  y1) *
				 u3 + (x - x1) * (y - y1) * u4);

  return uInterploated;
}


real
  IO::interpolateVelocityV (real x, real y, const GridFunction & v,
			    const Point & delta)
{
  real deltaX = delta[0];
  real deltaY = delta[1];

  // Computation of v(x,y)
  MultiIndex index;
  index[0] = ((int) ((x + (deltaX / 2)) / deltaX)) + 1;
  index[1] = ((int) (y / deltaY)) + 1;

  // The coordinates of the cell corners

  real x1 = ((index[0] - 1) - 0.5) * deltaX;
  real x2 = (index[0] - 0.5) * deltaX;
  real y1 = (index[1] - 1) * deltaY;
  real y2 = index[1] * deltaY;

  MultiIndex offset;

  offset[0] = index[0] - 1;
  offset[1] = index[1] - 1;

  real v1 = Element (v, offset);	//datafields->v->getField ()[i - 1][j - 1];

  offset[0] = index[0];
  offset[1] = index[1] - 1;

  real v2 = Element (v, offset);	//datafields->v->getField ()[i][j - 1];

  offset[0] = index[0] - 1;
  offset[1] = index[1];

  real v3 = Element (v, offset);	//datafields->v->getField ()[i - 1][j];


  real v4 = Element (v, index);	//datafields->v->getField ()[i][j];

  real
    vInterpolated =
    (1.0 / (deltaX * deltaY)) * ((x2 - x) * (y2 - y) *
				 v1 + (x - x1) * (y2 -
						  y) *
				 v2 + (x2 - x) * (y -
						  y1) *
				 v3 + (x - x1) * (y - y1) * v4);
  return vInterpolated;
}

void
IO::writeVTKFile (const MultiIndex & griddimension, const GridFunction & u,
		  const GridFunction & v, const GridFunction & p,
		  const Point & delta, int step)
{
  real deltaX = delta[0];
  real deltaY = delta[1];

  index_t iMax = griddimension[0] - 1;
  index_t jMax = griddimension[1] - 1;

  char numstr[21];
  sprintf (numstr, "%d", step);
  std::string filename;
  filename.append ("./");
  filename.append ("field_");
  filename.append (numstr);
  filename.append (".vts");

  std::filebuf fb;
  fb.open (filename.c_str (), std::ios::out);
  std::ostream os (&fb);

  os << "<?xml version=\"1.0\"?>" << std::endl
    << "<VTKFile type=\"StructuredGrid\">" << std::endl
    << "<StructuredGrid WholeExtent=\""
    << "0" << " " << (iMax - 1) << " "
    << "0" << " " << (jMax - 1) << " "
    << "0" << " " << "0" << " "
    << "\" GhostLevel=\"" << "1" << "\">" << std::endl
    << "<Piece Extent=\""
    << "0" << " " << (iMax - 1) << " "
    << "0" << " " << (jMax - 1) << " "
    << "0" << " " << "0" << " "
    << "\">" << std::endl
    << "<Points>" << std::endl
    <<
    "<DataArray type=\"Float64\" format=\"ascii\" NumberOfComponents=\"3\"> "
    << std::endl;
  for (int i = 0; i < iMax; ++i)
    {
      for (int j = 0; j < jMax; ++j)
	{
	  os << std::scientific << i * deltaX << " " << j *
	    deltaY << " " << 0.0 << std::endl;
	}
    }
  os << "</DataArray>" << std::endl
    << "</Points>" << std::endl
    << "<PointData Vectors=\"field\"  Scalars=\"P\">"
    << std::endl <<
    "<DataArray Name=\"field\" NumberOfComponents=\"3\" type=\"Float64\" >" <<
    std::endl;
  for (int i = 0; i < iMax; ++i)
    {
      real x = i * deltaX;

      for (int j = 0; j < jMax; ++j)
	{
	  real y = j * deltaY;

	  os << std::scientific << interpolateVelocityU (x, y, u,
							 delta) << " " <<
	    interpolateVelocityV (x, y, v, delta) << " " << 0. << std::endl;
	}

    }
  os << "</DataArray>" << std::endl
    << "<DataArray type=\"Float64\" Name=\"P\" format=\"ascii\">" <<
    std::endl;
  for (int i = 0; i <= iMax; ++i)
    {
      for (int j = 0; j <= jMax; ++j)
	{
	  os << std::scientific << p(i,j) << " ";

	}
      os << std::endl;

    }

  os << "</DataArray>" << std::endl
    << "</PointData>" << std::endl
    << "</Piece>" << std::endl
    << "</StructuredGrid>" << std::endl << "</VTKFile>" << std::endl;
  fb.close ();
}

void IO::writeRawOutput( const MultiIndex& griddimension,
		const GridFunction& u, const GridFunction& v,
		const GridFunction& p, const Point& delta, index_t step )
{
	char filename[ 256 ];
	sprintf( filename, "field-uvp-%d.txt", step );

	std::ofstream os( filename );

	os << griddimension.x << ' ' << griddimension.y << std::endl;

	for( int j = 0; j < griddimension.y; ++j )
	{
		real y = j * delta.y + delta.y / 2.0;
		for( int i = 0; i < griddimension.x; ++i )
		{
			real x = i * delta.x + delta.x / 2.0;
			real interpU = ( u( i + 1, j + 1 ) + u( i, j + 1 ) ) / 2.0;
			real interpV = ( v( i + 1, j + 1 ) + v( i + 1, j ) ) / 2.0;
			os << std::scientific << x << ' ' << y << ' ' << interpU << ' '
				<< interpV << ' ' << p( i + 1, j + 1 ) << std::endl;
		}
	}

	os.close();
}

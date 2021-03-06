#include "IO.h"
#include <cstring>
#include <string>
#include <cstdlib>
#include <cmath>

using std::size_t;


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
		  const Point & delta, const Params& params, int instance, int step)
{
  real deltaX = delta[0];
  real deltaY = delta[1];

  index_t iMax = griddimension[0] - 1;
  index_t jMax = griddimension[1] - 1;

  char numstr[32];
  sprintf (numstr, "%d_%d", instance, step);
  std::string filename;
  filename.append ("./");
  filename.append ("field_");
  filename.append (numstr);
  filename.append (".vts");

  std::filebuf fb;
  fb.open (filename.c_str (), std::ios::out);
  std::ostream os (&fb);
  if(os.fail())
  {
    std::cout << "Fehler beim \"Offnen von " << filename;
    exit(1);
  }

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
     //inputfile as comment:
     << "<!-- " ;
  os << params;
  os << "--> " << std::endl;
  os << "<DataArray type=\"Float64\" format=\"ascii\" NumberOfComponents=\"3\"> "
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
    << "<PointData Vectors=\"u\"  Scalars=\"P\">"
    << std::endl <<
    "<DataArray Name=\"u\" NumberOfComponents=\"3\" type=\"Float64\" >" <<
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

void IO::writeRawOutput( const GridFunction& u, const GridFunction& v,
		const GridFunction& p, const GridFunction& T, const GridFunction& psi,
		const GridFunction& zeta, const Point& delta, const Point& offset,
		index_t instance, index_t step, const MultiIndex& rank )
{
	char filename[ 256 ];
	MultiIndex griddimension;

	sprintf( filename, "field-uvpT-instance-%03d-step-%06d-rank-%03d-%03d.txt",
			instance, step, rank.x, rank.y );
	std::ofstream os( filename );

	griddimension = p.size() - MultiIndex( 2, 2 );

	os << griddimension.x << ' ' << griddimension.y << std::endl;
	for( int j = 0; j < griddimension.y; ++j )
	{
		real y = offset.y + j * delta.y + delta.y / 2.0;
		for( int i = 0; i < griddimension.x; ++i )
		{
			real x = offset.x + i * delta.x + delta.x / 2.0;
			real interpU = ( u( i + 1, j + 1 ) + u( i + 2, j + 1 ) ) / 2.0;
			real interpV = ( v( i + 1, j + 1 ) + v( i + 1, j + 2 ) ) / 2.0;
			os << std::scientific << x << ' ' << y << ' ' << interpU << ' '
				<< interpV << ' ' << p( i + 1, j + 1 ) << ' ' << T( i + 1, j + 1 ) << std::endl;
		}
	}
	os.close();

	sprintf( filename, "field-psi-instance-%03d-step-%06d-rank-%03d-%03d.txt",
			instance, step, rank.x, rank.y );
	std::ofstream os2( filename );

	griddimension = psi.size();

	os2 << griddimension.x << ' ' << griddimension.y << std::endl;
	for( int j = 0; j < griddimension.y; ++j )
	{
		real y = offset.y + j * delta.y;
		for( int i = 0; i < griddimension.x; ++i )
		{
			real x = offset.x + i * delta.x;
			os2 << std::scientific << x << ' ' << y << ' ' << psi( i, j )
				<< std::endl;
		}
	}
	os2.close();

	sprintf( filename, "field-zeta-instance-%03d-step-%06d-rank-%03d-%03d.txt",
			instance, step, rank.x, rank.y );
	std::ofstream os3( filename );

	griddimension = zeta.size();

	os3 << griddimension.x << ' ' << griddimension.y << std::endl;
	for( int j = 0; j < griddimension.y; ++j )
	{
		real y = offset.y + ( j + 1 ) * delta.y;
		for( int i = 0; i < griddimension.x; ++i )
		{
			real x = offset.x + ( i + 1 ) * delta.x;
			os3 << std::scientific << x << ' ' << y << ' ' << zeta( i, j )
				<< std::endl;
		}
	}
	os3.close();
}

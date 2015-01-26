#include "params.h"

#include <fstream>
#include <string>
#include <cstring>
#include <iostream>

using std::size_t;

void Params::parseFile( const char *filename )
{
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
	_member = strtod( value.c_str(), NULL ); }
#define SET_PARAM_INT(_member, _name) if( name == _name ) {\
	_member = strtol( value.c_str(), NULL, 10 ); }

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
		else SET_PARAM_DOUBLE(gamma, "gamma")
		else SET_PARAM_DOUBLE(Re, "re")
		else SET_PARAM_DOUBLE(Pr, "pr")
		else SET_PARAM_DOUBLE(beta, "beta")
		else SET_PARAM_DOUBLE(initialVelocity.x, "ui")
		else SET_PARAM_DOUBLE(initialVelocity.y, "vi")
		else SET_PARAM_DOUBLE(initialPressure, "pi")
		else SET_PARAM_DOUBLE(volForce.x, "gx")
		else SET_PARAM_DOUBLE(volForce.y, "gy")

#undef SET_PARAM_DOUBLE
#undef SET_PARAM_INT
	}
}

int Params::parseCmdArgs( int argc, char* argv[] )
{
	int i = 0;
	while( i < argc )
	{
		if( argv[ i ][ 0 ] != '-' )
		{
			++i;
			continue;
		}

#define SET_PARAM_DOUBLE(_member, _name) \
		if( strcmp( name, _name ) == 0 && i + 1 < argc )\
		{\
			str2double( argv[ ++i ], _member );\
		}
#define SET_PARAM_INT(_member, _name) \
		if( strcmp( name, _name ) == 0 && i + 1 < argc )\
		{\
			str2int( argv[ ++i ], _member );\
		}

		char* name = &argv[ i ][ 1 ];
		if( strcmp( name, "paramfile" ) == 0 && i + 1 < argc )
		{
			parseFile( argv[ ++i ] );
		}
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
		else SET_PARAM_DOUBLE(gamma, "gamma")
		else SET_PARAM_DOUBLE(Re, "re")
		else SET_PARAM_DOUBLE(Pr, "pr")
		else SET_PARAM_DOUBLE(beta, "beta")
		else SET_PARAM_DOUBLE(initialVelocity.x, "ui")
		else SET_PARAM_DOUBLE(initialVelocity.y, "vi")
		else SET_PARAM_DOUBLE(initialPressure, "pi")
		else SET_PARAM_DOUBLE(volForce.x, "gx")
		else SET_PARAM_DOUBLE(volForce.y, "gy")
		else if( strcmp( name, "-" ) == 0 )
		{
			return i + 1;
		}

#undef SET_PARAM_DOUBLE
#undef SET_PARAM_INT

		++i;
	}

	return i;
}

void Params::str2double( const char* str, double& out )
{
	char* end;
	double res = strtod( str, &end );
	if( end > str )
	{
		out = res;
	}
}

void Params::str2int( const char* str, int& out )
{
	char* end;
	int res = strtol( str, &end, 10 );
	if( end > str )
	{
		out = res;
	}
}

std::ostream& operator <<( std::ostream& s, const Params& p )
{
	s
		<< "  Domain size:           "            << p.domainSize.x       << " x " << p.domainSize.y << std::endl
		<< "  Grid size:             "            << p.gridSize.x         << " x " << p.gridSize.y   << std::endl
		<< "  T:                     "            << p.T                  <<                          std::endl
		<< "  Re:                    "            << p.Re                 <<                          std::endl
		<< "  Pr:                    "            << p.Pr                 <<                          std::endl
		<< "  Beta:                  "            << p.beta               <<                          std::endl
		<< "  dt:                    "            << p.dt                 <<                          std::endl
		<< "  tau:                   "            << p.tau                <<                          std::endl
		<< "  Output time step:      "            << p.deltaVec           <<                          std::endl
		<< "  Max SOR iterations:    "            << p.maxIter            <<                          std::endl
		<< "  SOR error tolerance:   "            << p.eps                <<                          std::endl
		<< "  Relaxation factor:     "            << p.omega              <<                          std::endl
		<< "  Upwinding factor:      "            << p.alpha              <<                          std::endl
		<< "  Upwinding factor (T):  "            << p.gamma              <<                          std::endl
		<< "  Initial velocity:      ("           << p.initialVelocity.x  << "," << p.initialVelocity.y << ")" << std::endl
		<< "  Initial pressure:      "            << p.initialPressure    << std::endl
		<< "  Volumetric force:      ("           << p.volForce.x  << "," << p.volForce.y << ")" << std::endl;

	return s;
}

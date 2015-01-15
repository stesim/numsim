#ifndef PARAMS_H
#define PARAMS_H

#include "typedef.h"
#include <ostream>

struct Params
{
public:
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

	/*!
	* Method reads the simulation parameters from the specified input file.
	* 
	* @param filename The name of the file with the simulations parameters
	*/
	void parseFile( const char *filename );

	int parseCmdArgs( int argc, char* argv[] );

private:
	void str2double( const char* str, double& out );

	void str2int( const char* str, int& out );
};

std::ostream& operator <<( std::ostream& s, const Params& p );

#endif

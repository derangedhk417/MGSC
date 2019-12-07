// Author:  Adam J. Robinson
// Summary: This file contains definition for the integral that 
//          characterizes the interaction between nuclei and electrons,
//          as defined in the mathematical documentation.
// Performs the Electron Nuclei integral for the given 
// Gaussian wavefunction.

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_errno.h>

// Can't just include "gaussian.h", it will create a circular
// reference.
class GaussianWFN;

class ENIntegrator {

public:

	gsl_monte_vegas_state * state;

	double lowerBounds[3];
	double upperBounds[3];
	GaussianWFN * wavefn;

	int    maxCalls;
	double width;
	double epsilon;

	// Initializes an integral calculator.
	//     fn     = The wavefunction that the computation is being carried
	//              out for.
	//     width  = The number of half width half maxes to 
	//              limit the integral range to. Monte Carlo
	//              integration does not work for bounds at infinity.
	//              typically, 3 < width < 6
	//
	//     ncalls = The maximum number of calls to the integrand function
	//              that the integrator is permitted to make in order to 
	//              calculate the integral. Higher values are slower and
	//              more accurate.
	ENIntegrator(GaussianWFN * fn, double width, int ncalls);

	~ENIntegrator();

	// These are set temporarily before an integral is performed.
	double A1;
	double A2;
	double A3;

	double s1;
	double s2;
	double s3;

	double R1;
	double R2;
	double R3;

	// Given the three indices that define the integral being calculated,
	// sets the lowerBounds and upperBounds members of the class to 
	// values appropriate for the given width, nuclei location, gaussian
	// center and gaussian width.
	void setBounds(int iu, int w, int l);

	double Integrate(int iu, int w, int l);

private:

	static double integrand(double * x, size_t dim, void * params);

};
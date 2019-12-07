// Author:  Adam J. Robinson
// Summary: This file contains definitions for the Gaussian
//          trial wavefunction used by MGSC.

#include "integral_en.h"

// Stores all variables related to a Gaussian trial wavefunction and
// provides methods for calculating expectation values.
class GaussianWFN {
public:
	// These variables need to be assigned a value before the 
	// trial wavefunction can be used.
	int       m;  // Number of terms
	int       n;  // Number of inputs (3 * # of electrons)
	double ** A;  // Shape matrices (really vectors due to simplifications)
	double ** s;  // Shift matrices (really vectors)
	double  * C;  // Coefficients on each Gaussian term
	double ** R;  // List of nuclear coordinates
	int       Nu; // Number of nuclei
	double  * Q;  // List of nuclear charges


	ENIntegrator integrator;

	// Initializes a Gaussian trial wavefunction with a given number of 
	// Gaussian terms. The variables annotated above need to be given a value
	// before the wavefunction is actually useable.
	GaussianWFN(int nterms, int nelectrons);

	~GaussianWFN();

public:

	// This section contains the primary functionality of the class,
	// calculation of the expectation value of the Hamiltonian.

	// Gets the expectation value of the Hamiltonian based on the
	// current set of parameters defining the wavefunction.
	double getHamiltonianExpectation();

	// This calculates the value denoted 'B' in the documentation. 
	// The value is calculated such that |B|^2 multiplied by the 
	// integral over all space of |wavefunction|^2 will always be
	// unity.
	double getNormalizationConstant();

	// Given the appropriate normalization constant, returns the
	// expectation value of the kinetic energy for the given set 
	// of parameters.
	double getKineticExpectation(double B);

	// Given the appropriate normalization constant, returns the
	// expectation value of the potential energy for the given set 
	// of parameters.
	double getPotentialExpectation(double B);
};
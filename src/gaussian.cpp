// Author:  Adam J. Robinson
// Summary: This file contains implementation of the Gaussian trial 
//          wavefunction used by MGSC. Some of the code is delegated to other
//          files for readability.

#include <iostream>   // Used only in DEBUG mode for error messages
#include <mathimf.h>
#include <chrono>     // Used in profile mode 

// Comment out the undef statements to see corresponding messages.
#define DEBUG
#undef  DEBUG

#define PROFILE
#undef  PROFILE

// Stores all variables related to a Gaussian trial wavefunction and
// provides methods for calculating expectation values.
class GaussianTrialWavefunction {
public:
	// These variables need to be assigned a value before the 
	// trial wavefunction can be used.
	int      m;  // Number of terms
	int      n;  // Number of inputs (3 * # of electrons)
	float ** A;  // Shape matrices (really vectors due to simplifications)
	float ** s;  // Shift matrices (really vectors)
	float  * C;  // Coefficients on each Gaussian term
	float ** R;  // List of nuclear coordinates
	int      Nu; // Number of nuclei
	float  * Q;  // List of nuclear charges

	// Initializes a Gaussian trial wavefunction with a given number of 
	// Gaussian terms. The variables annotated above need to be given a value
	// before the wavefunction is actually useable.
	GaussianWavefunction(int nterms, int nelectrons) {
		// In a future version, this will be 3 times and number of electrons.
		// Right now it is just going to be 3, corresponding to a single 
		// electron.
		n = nelectrons;
		m = nterms;
	}


}
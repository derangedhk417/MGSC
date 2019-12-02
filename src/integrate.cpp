// This file defines the class that is used to perform multi-dimensional 
// integration numerically.

#include <iostream>
#include <mathimf.h>
#include <chrono>
#include <gsl/gsl_integration.h>

// Integrates the specified function over the specified range 
// of independent variable values using a triangular integral.
// The step size is adjusted dynamically according to the 
// slope of the integrand.
class DynamicIntegrator {
public:
	// These are the minimum and maximum increment of the 
	// variable of integration as a function of the total
	// range of integration. For example, if an integral
	// is on [0, 2] and min_step is 1e-2, the minimum step
	// size will be 2e-2.
	// This effectively creates an upper bound on the number
	// of steps in the integration, while allowing the number
	// of steps to be significantly lower in cases where the
	// integrand has a low slope.
	float min_step    = 0.05;
	float max_step    = 0.1;
	float min_slope   = 0.05; // The slope under which the max_step will be used.
	float max_slope   = 1.0;  // The slope above which the min_step will be used.
	float finite_step = 1e-7;    // The finite step to use for the derivative. 

	void *voidData = NULL;

	int     nVars;
	float * lowerBounds;
	float * upperBounds;
	float (*integrand)(float *, void *);

	// These are used when scaling the step size.
	float a;
	float b;

	// Arguments:
	// 1) The integrand, the first parameter should be the input coordinates in order.
	//    The second parameter should be a void pointer containing any additional info
	//    that the function needs. Call DynamicIntegrator.setVoidPointer to determine
	//    what information should be passed.
	// 2) An array of lower bounds for each variable.
	// 3) An array of upper bounds for each variable.
	// 4) The number of variables to integrate over.
	DynamicIntegrator(float (*integrand)(float *, void *), float * lower, float * upper, int n) {
		this->integrand = integrand;
		nVars           = n;
		lowerBounds     = lower;
		upperBounds     = upper;
	}

	// Used to set the additional data provided to the integrand.
	void setVoidPointer(void *p) {
		voidData = p;
	}

	// Gets the appropriate step for the integration.
	float getStep(float slope) {
		float abs_slope = fabsf(slope);

		if (abs_slope <= min_slope) {
			return max_step;
		} else if (abs_slope >= max_slope) {
			return min_step;
		} else {
			return a * abs_slope * abs_slope + b;
		}
	}

	// Returns the proper integral of the function, under the given conditions.
	float integrate() {
		a = (min_step - max_step) / ((max_slope*max_slope) - (min_slope*min_slope));
		b = max_step - a*min_slope*min_slope;

		float result = 0.0;

		float * location = new float[nVars];

		location[0] = lowerBounds[0];
		while (location[0] <= upperBounds[0]) {
			result += integral_step(location, 0);
		}

		return result;
	}

	// Recursively integrates each independent variable. Hard to explain,
	// read the function for a full understanding.
	float integral_step(float *location, int var) {
		if (var == nVars - 1) {
			float point       = integrand(location, voidData);
			location[var]    += finite_step;
			float incremented = integrand(location, voidData);
			location[var]    -= finite_step;
			float slope       = (incremented - point) / finite_step;

			float proper_step = getStep(slope);
			location[var] += proper_step;
			float proper_incremented_value = integrand(location, voidData);

			return (point + proper_incremented_value) * (proper_step / 2);
		} else {
			// First, we get the appropriate value of the integral at this point.

			float point = 0.0;

			// Get a value at this point.
			float save = location[var + 1];
			location[var + 1] = lowerBounds[var + 1];
			while (location[var + 1] < upperBounds[var + 1]) {
				point += integral_step(location, var + 1);
			}
			location[var + 1] = save;

			// Now we calculate the slope at this point.
			location[var] += finite_step;
			float incval = 0.0;
			location[var + 1] = lowerBounds[var + 1];
			while (location[var + 1] < upperBounds[var + 1]) {
				incval += integral_step(location, var + 1);
			}
			location[var + 1] = save;
			location[var]    -= finite_step;

			float slope = (incval - point) / finite_step;

			float proper_step = getStep(slope);

			// We now know how far to step. Take the step, evaluate
			// the integrand and determine the appropriate value to add.
			location[var] += proper_step;
			float proper_incremented_value = 0.0;
			location[var + 1] = lowerBounds[var + 1];
			while (location[var + 1] < upperBounds[var + 1]) {
				proper_incremented_value += integral_step(location, var + 1);
			}
			location[var + 1] = save;

			float step_result = (point + proper_incremented_value) * (proper_step / 2);

			return step_result;
		}
	}
};
#include <iostream>
#include <mathimf.h>
#include <chrono>

using namespace std;

const float pi = 3.141592653589793238;

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
	float min_step    = 1e-3;
	float max_step    = 1e-1;
	float min_slope   = 0.25; // The slope under which the max_step will be used.
	float max_slope   = 25.0; // The slope above which the min_step will be used.
	float finite_step = 1e-7; // The finite step to use for the derivative. 

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


// This is the same as the dynamic integrator, except that
// the step size is constant.
class StaticIntegrator {
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
	float step = 0.5e-2;

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
	StaticIntegrator(float (*integrand)(float *, void *), float * lower, float * upper, int n) {
		this->integrand = integrand;
		nVars           = n;
		lowerBounds     = lower;
		upperBounds     = upper;
	}

	// Used to set the additional data provided to the integrand.
	void setVoidPointer(void *p) {
		voidData = p;
	}

	// Returns the proper integral of the function, under the given conditions.
	float integrate() {
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
			
			location[var] += step;
			float proper_incremented_value = integrand(location, voidData);

			return (point + proper_incremented_value) * (step / 2);
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


			// We now know how far to step. Take the step, evaluate
			// the integrand and determine the appropriate value to add.
			location[var] += step;
			float proper_incremented_value = 0.0;
			location[var + 1] = lowerBounds[var + 1];
			while (location[var + 1] < upperBounds[var + 1]) {
				proper_incremented_value += integral_step(location, var + 1);
			}
			location[var + 1] = save;

			float step_result = (point + proper_incremented_value) * (step / 2);

			return step_result;
		}
	}
};


float xCubedIntegrand(float * input, void * unused) {
	return input[0] * input[0] * input[0];
}

float xCubedIntegrand2D(float * input, void * unused) {
	return input[0] * input[0] * input[0] + input[1] * input[1] * input[1];
}

float hydrogenGroundState(float * input, void * a_ptr) {
	float a = *((float *)a_ptr);
	return input[0]*expf(-a*input[0])*sinf(input[1]);
}

struct ENParams {
	float * A;
	float * s;
	float * R;
};

float ENIntegrand(float * input, void *EN_ptr) {
	float x1 = input[0];
	float x2 = input[1];
	float x3 = input[2];

	struct ENParams * ENP = (struct ENParams * )EN_ptr;
	float * A = ENP->A;
	float * s = ENP->s;
	float * R = ENP->R;

	float numerator = -(x1*x1)*A[0] + x1*s[0] - (x2*x2)*A[1] + x2*s[1]  -(x3*x3)*A[2] + x3*s[2];
	numerator = expf(numerator);

	float d1 = R[0] - x1;
	float d2 = R[1] - x2;
	float d3 = R[2] - x3;

	float denominator = d1*d1 + d2*d2 + d3*d3;
	denominator = fmax(sqrtf(denominator), 1e-7);


	return numerator / denominator;
}

int main(int argc, char ** argv) {
	// What follows is a series of tests of the integration method.
	// All of these integrals are known analytically, so this should
	// work.

	// ----------------------------------------------------
	// x cubed 1D
	// ----------------------------------------------------
	float correct_value = powrf(5.0, 4.0) * (1.0 / 4.0);
	float lower[] = {0.0};
	float upper[] = {5.0};
	DynamicIntegrator d1(xCubedIntegrand, lower, upper, 1);

	d1.min_step = 1e-3;
	d1.max_step = 1e-2;

	auto begin = chrono::high_resolution_clock::now();
	float result = d1.integrate();
	auto end      = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::nanoseconds>(end - begin).count();

	cout << "Cubic 1D" << endl;
	cout << "Analytical Value: " << correct_value << endl;
	cout << "Numerical Value:  " << result        << endl;
	cout << "Relative Error:   " << (result - correct_value) / correct_value << endl;
	cout << "Time:             " << duration / 1000.0 << "μs" << endl;
	cout << endl << endl;

	// ----------------------------------------------------
	// x cubed 2D
	// ----------------------------------------------------
    correct_value = 10 * correct_value;
	float lower2[] = {0.0, 0.0};
	float upper2[] = {5.0, 5.0};
	DynamicIntegrator d2(xCubedIntegrand2D, lower2, upper2, 2);

	d2.min_step = 0.004;
	d2.max_step = 0.01;

	begin    = chrono::high_resolution_clock::now();
	result   = d2.integrate();
	end      = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::nanoseconds>(end - begin).count();

	cout << "Cubic 2D" << endl;
	cout << "Analytical Value: " << correct_value << endl;
	cout << "Numerical Value:  " << result        << endl;
	cout << "Relative Error:   " << (result - correct_value) / correct_value << endl;
	cout << "Time:             " << duration / 1000.0 << "μs" << endl;
	cout << endl << endl;

	// ----------------------------------------------------
	// Hydrogen Ground State
	// ----------------------------------------------------
	float a       = 1.0;
    correct_value = (4 * pi) / (a * a);
	float lower3[] = {0.0, 0.0, 0.0};
	float upper3[] = {25.0, pi,  2*pi};
	DynamicIntegrator d3(hydrogenGroundState, lower3, upper3, 3);

	d3.min_step  = 0.004;
	d3.max_step  = 0.5; 
	d3.setVoidPointer(&a);

	begin    = chrono::high_resolution_clock::now();
	result   = d3.integrate();
	end      = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::nanoseconds>(end - begin).count();

	cout << "Hydrogen Wavefunction" << endl;
	cout << "Analytical Value: "    << correct_value << endl;
	cout << "Numerical Value:  "    << result        << endl;
	cout << "Relative Error:   "    << (result - correct_value) / correct_value << endl;
	cout << "Time:             "    << duration / 1000.0 << "μs" << endl;
	cout << endl << endl;

	// ----------------------------------------------------
	// Electron - Nucleus Integral
	// ----------------------------------------------------
	struct ENParams ENP;
	float A[] = {1.0, 1.0, 1.0};
	float s[] = {0.0, 0.0, 0.0};
	float R[] = {0.0, 0.0, 0.0};
	ENP.A     = A;
	ENP.s     = s;
	ENP.R     = R;

    correct_value = 1.0;

    cout << "Electron - Nucleus Integral"  << endl;

    float range = 5.0;
    float xlower = fmin(R[0], -s[0]) - range;
    float xupper = fmax(R[0], -s[0]) + range;

    float ylower = fmin(R[1], -s[1]) - range;
    float yupper = fmax(R[1], -s[1]) + range;

    float zlower = fmin(R[2], -s[2]) - range;
    float zupper = fmax(R[2], -s[2]) + range;

    cout << "x -> " << xlower << " - " << xupper << endl;
    cout << "y -> " << ylower << " - " << yupper << endl;
    cout << "z -> " << zlower << " - " << zupper << endl;

	float lower4[] = {xlower, ylower, zlower};
	float upper4[] = {xupper, yupper, zupper};
	DynamicIntegrator d4(ENIntegrand, lower4, upper4, 3);

	d4.min_step  = 0.005;
	d4.max_step  = 0.5;
	d4.min_slope = 0.05;
	d4.max_slope = 1.0;
	d4.setVoidPointer(&ENP);

	begin    = chrono::high_resolution_clock::now();
	result   = d4.integrate();
	end      = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::nanoseconds>(end - begin).count();

	
	cout << "Analytical Value: "           << correct_value << endl;
	cout << "Numerical Value:  "           << result        << endl;
	cout << "Relative Error:   "           << (result - correct_value) / correct_value << endl;
	cout << "Time:             "           << duration / 1000.0 / 1000.0 << "ms" << endl;
	cout << endl << endl;

	cout << "Beginning Convergence Test" << endl;

	float test_vals[] = {0.25, 0.1, 0.05, 0.01, 0.0075, 0.005};
	int   n           = 6;

	for (int i = 0; i < n; ++ i) {
		d4.min_step = test_vals[i];
		begin    = chrono::high_resolution_clock::now();
		result   = d4.integrate();
		end      = chrono::high_resolution_clock::now();
		duration = chrono::duration_cast<chrono::nanoseconds>(end - begin).count();
		cout << d4.min_step << " " << duration / 1000.0 / 1000.0 << " " << result << endl;
	}
	

	return 0;
}
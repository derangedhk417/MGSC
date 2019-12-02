#include <iostream>
#include <mathimf.h>

using namespace std;

const float pi = 3.141592653589793238;

// Given a rather large set of parameters, performs the specified
// electron - electron integral. 
class EEIntegrator {
public:
	float dx_min      = 1e-4;
	float dx_max      = 1e-1;
	float slope_min   = 0.1;
	float slope_max   = 25.0;
	float finite_step = 1e-8;

	int i;
	int j;

	float * Aw;
	float * Al;
	float * sw;
	float * sl;

	float r1;
	float r2;
	float theta1;
	float theta2;
	float phi1;
	float phi2;

	float r_upper_bound = 500.0;

	float a;
	float b;

	EEIntegrator(int i, int j) {
		this->i = i;
		this->j = j;
	}

	float Integrate() {
		a = (dx_min - dx_max) / ((slope_max*slope_max) - (slope_min*slope_min));
		b = dx_max - a*slope_min*slope_min;

		r1     = 0.0;
		r2     = 0.0;
		theta1 = 0.0;
		theta2 = 0.0;
		phi1   = 0.0;
		phi2   = 0.0;
		
		float total = 0.0;

		while (r1 <= r_upper_bound) {
			total += r1Step();
		}

		return 0.0;
	}

private:

	// Given the current values of the integrand and the shape
	// matrices, will calculate the value of the integrand.
	float integrand() {
		// Figure out what the Y and V should be.
		float Y1 = Aw[3 * i]     + Al[3 * i];
		float Y2 = Aw[3 * i + 1] + Al[3 * i + 1];
		float Y3 = Aw[3 * i + 2] + Al[3 * i + 2];

		float V1 = sw[3 * i]     + sl[3 * i];
		float V2 = sw[3 * i + 1] + sl[3 * i + 1];
		float V3 = sw[3 * i + 2] + sl[3 * i + 2];

		float Y1_2 = Aw[3 * j]     + Al[3 * j];
		float Y2_2 = Aw[3 * j + 1] + Al[3 * j + 1];
		float Y3_2 = Aw[3 * j + 2] + Al[3 * j + 2];

		float V1_2 = sw[3 * j]     + sl[3 * j];
		float V2_2 = sw[3 * j + 1] + sl[3 * j + 1];
		float V3_2 = sw[3 * j + 2] + sl[3 * j + 2];

		// Calculate P1.
		float sintheta1 = sinf(theta1);
		float costheta1 = cosf(theta1);
		float sinphi1   = sinf(phi1);
		float cosphi1   = cosf(phi1);

		float t1 = -r1*sintheta1*cosphi1*cosphi1*Y1 + cosphi1*V1 - r1*sintheta1*sinphi1*sinphi1*Y2 + sinphi1*V2;
		float t2 = -r1*costheta1*Y3 + V3;

		float P1 = expf(r1*(sintheta1*t1 + costheta1*t2));

		// Calculate P2.
		float sintheta2 = sinf(theta2);
		float costheta2 = cosf(theta2);
		float sinphi2   = sinf(phi2);
		float cosphi2   = cosf(phi2);

		float t1_2 = -r2*sintheta2*cosphi2*cosphi2*Y1_2 + cosphi2*V1_2 - r2*sintheta2*sinphi2*sinphi2*Y2_2 + sinphi2*V2_2;
		float t2_2 = -r2*costheta2*Y3_2 + V3_2;

		float P2 = expf(r2*(sintheta2*t1_2 + costheta2*t2_2));

		// Calculate the full integrand.
		float result = P1 * P2 * r1*r1 * r2*r2 * sintheta1 * sintheta2
		result      /= sqrtf(r1*r1 + r2*r2 - 2*r1*r2*cosf(theta1 - theta2))

		return result;
	}

	// Given the slope of the function at the current point,
	// will return the appropriate step in the independent
	// variable.
	float getStep(float slope) {
		float abs_slope = fabsf(slope);

		if (abs_slope <= slope_min) {
			return dx_max;
		} else if (abs_slope >= slope_max) {
			return dx_min;
		} else {
			return a * abs_slope * abs_slope + b;
		}
	}

	float integralStep(float **varlist, int nvars, int current) {
		// First, get the derivative.

		float currentval = 0.0;

		// Get a value at this point.
		float original_inner = *varlist[current + 1];
		while (*varlist[current + 1] < r_upper_bound) {
			currentval += integralStep(varlist, nvars, current + 1);
		}
		*varlist[current + 1] = original_inner;

		// increment r1 and get a value.
		*varlist[current] += finite_step;
		float incval = 0.0;
		original_r2 = r2;
		while (r2 < r_upper_bound) {
			incval += r2Step();
		}
		r2 = original_r2;
		r1 -= finite_step;

		float slope = (incval - currentval) / finite_step;

		float proper_step = getStep(slope);

		// We now know how far to step. Take the step, evaluate
		// the integrand and determine the appropriate value to add.
		r1 += proper_step;
		float proper_incremented_value = 0.0;
		original_r2 = r2;
		while (r2 < r_upper_bound) {
			proper_incremented_value += r2Step();
		}
		r2 = original_r2;

		float step_result = (currentval + proper_incremented_value) * (proper_step / 2);

		return step_result;
	}


	// This process is as follows:
	//     - Calculate the derivative and choose an appropriate step.
	//     - Evaluate the integrand both here and at that step.
	//     - Calculate the triangular integral.
	//     - Increment by the appropriate step.
	float r1Step() {
		// First, get the derivative.

		float currentval = 0.0;

		// Get a value at this point.
		float original_r2 = r2;
		while (r2 < r_upper_bound) {
			currentval += r2Step();
		}
		r2 = original_r2;

		// increment r1 and get a value.
		r1 += finite_step;
		float incval = 0.0;
		original_r2 = r2;
		while (r2 < r_upper_bound) {
			incval += r2Step();
		}
		r2 = original_r2;
		r1 -= finite_step;

		float slope = (incval - currentval) / finite_step;

		float proper_step = getStep(slope);

		// We now know how far to step. Take the step, evaluate
		// the integrand and determine the appropriate value to add.
		r1 += proper_step;
		float proper_incremented_value = 0.0;
		original_r2 = r2;
		while (r2 < r_upper_bound) {
			proper_incremented_value += r2Step();
		}
		r2 = original_r2;

		float step_result = (currentval + proper_incremented_value) * (proper_step / 2);

		return step_result;
	}

	float r2Step() {
		// First, get the derivative.

		float currentval = 0.0;

		// Get a value at this point.
		float original_theta1 = theta1;
		while (theta1 < pi) {
			currentval += theta1Step();
		}
		theta1 = original_theta1;

		// increment r2 and get a value.
		r2 += finite_step;
		float incval = 0.0;
		original_theta1 = theta1;
		while (theta1 < pi) {
			incval += theta1Step();
		}
		theta1 = original_theta1;
		r2 -= finite_step;

		float slope = (incval - currentval) / finite_step;

		float proper_step = getStep(slope);

		// We now know how far to step. Take the step, evaluate
		// the integrand and determine the appropriate value to add.
		r2 += proper_step;
		float proper_incremented_value = 0.0;
		original_theta1 = theta1;
		while (theta1 < pi) {
			proper_incremented_value += theta1Step();
		}
		theta1 = original_theta1;

		float step_result = (currentval + proper_incremented_value) * (proper_step / 2);

		return step_result;
	}

};



int main(int argc, char ** argv) {
	float A1[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	float s1[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	float A2[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	float s2[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

	return 0;
}
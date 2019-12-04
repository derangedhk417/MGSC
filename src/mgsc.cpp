#include <iostream>
#include <mathimf.h>
#include <chrono>
#include "gaussian.cpp"
#include <gsl/gsl_multimin.h>

using namespace std;

double ExpectationValue(const gsl_vector * x, void * params) {
	float _A1 = gsl_vector_get(x, 0);
	float _A2 = gsl_vector_get(x, 1);
	float _C1 = gsl_vector_get(x, 2);
	float _C2 = gsl_vector_get(x, 3);

	float A1[]  = {_A1, _A1, _A1};
	float A2[]  = {_A2, _A2, _A2};
	float ** A  = new float*[2];
	A[0]        = A1;
	A[1]        = A2;

	float s1[] = {0.0, 0.0, 0.0};
	float ** s = new float*[1];
	s[0]       = s1;

	float C[]  = {_C1, _C2};
	float R1[] = {0.0, 0.0, 0.0};
	float ** R = new float*[1];
	R[0]       = R1;

	float Q[] = {N_qe};
	
	GaussianWavefunction wavefn(2, 0.01, 256);
	wavefn.A  = A;
	wavefn.s  = s;
	wavefn.C  = C;
	wavefn.R  = R;
	wavefn.Nu = 1;
	wavefn.Q  = Q; 

	float err;

	return wavefn.getHamiltonianExpectation(&err);
}

void ExpectationValueGradient(const gsl_vector * x, void * params, gsl_vector * grad) {
	double x, y;
	double *p = (double *)params;

	x = gsl_vector_get(v, 0);
	y = gsl_vector_get(v, 1);

	float _A1 = gsl_vector_get(x, 0);
	float _A2 = gsl_vector_get(x, 1);
	float _C1 = gsl_vector_get(x, 2);
	float _C2 = gsl_vector_get(x, 3);
	double dx = 1e-4; 

	gsl_vector_set(df, 0, 2.0 * p[2] * (x - p[0]));
	gsl_vector_set(df, 1, 2.0 * p[3] * (y - p[1]));

	// Get the gradient with respect to all four variables.
	float dA1;
	float dA2;
	float dC1;
	float dC2;

	float currentValue = ExpectationValue(x, NULL);

	gsl_vector_set(x, 0, _A1 + dx);
	float offsetA1 = ExpectationValue(x, NULL);
	gsl_vector_set(x, 0, _A1);
	dA1 = (offsetA1 - currentValue) / dx;

	gsl_vector_set(x, 1, _A2 + dx);
	float offsetA2 = ExpectationValue(x, NULL);
	gsl_vector_set(x, 1, _A2);
	dA2 = (offsetA2 - currentValue) / dx;

	gsl_vector_set(x, 2, _C1 + dx);
	float offsetC1 = ExpectationValue(x, NULL);
	gsl_vector_set(x, 2, _C1);
	dC1 = (offsetC1 - currentValue) / dx;

	gsl_vector_set(x, 3, _C2 + dx);
	float offsetC2 = ExpectationValue(x, NULL);
	gsl_vector_set(x, 3, _C2);
	dC2 = (offsetC2 - currentValue) / dx;

	gsl_vector_set(x, 0, dA1);
	gsl_vector_set(x, 1, dA2);
	gsl_vector_set(x, 2, dC1);
	gsl_vector_set(x, 3, dC2);

}

void Combined(const gsl_vector * x, void * params, double * f, gsl_vector * grad) {
 	*f = ExpectationValue(x, NULL);
 	ExpectationValueGradient(x, NULL, grad);
}


int main(int argc, char ** argv) {
	size_t iter = 0;
	int status;

	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;

	/* Position of the minimum (1,2), scale factors
	10,20, height 30. */
	double params[] = {};

	gsl_vector * x;
	gsl_multimin_function_fdf F;

	F.n      = 4;
	F.f      = ExpectationValue;
	F.df     = ExpectationValueGradient;
	F.fdf    = Combined;
	F.params = params;

	x = gsl_vector_alloc(4);
	gsl_vector_set(x, 0, 3.93e6);
	gsl_vector_set(x, 1, 3.93e7);
	gsl_vector_set(x, 2, 1.0);
	gsl_vector_set(x, 3, 0.2);


	T = gsl_multimin_fdfminimizer_conjugate_fr;
	s = gsl_multimin_fdfminimizer_alloc(T, 4);

	gsl_multimin_fdfminimizer_set(s, &F, x, 0.01, 1e-4);

	do {
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(s);

		if (status) break;

		status = gsl_multimin_test_gradient(s->gradient, 1e-3);

		if (status == GSL_SUCCESS) {
			printf("Minimum found at:\n");
			printf(
				"%5d %.5f %.5f %.5f %.5f %10.5f\n", 
				iter,
				gsl_vector_get(s->x, 0),
				gsl_vector_get(s->x, 1),
				gsl_vector_get(s->x, 2),
				gsl_vector_get(s->x, 3),
				s->f
			);
		}

	} while (status == GSL_CONTINUE && iter < 100);

	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(x);

	return 0;
}
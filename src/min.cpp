#include <iostream>
#include <mathimf.h>
#include <chrono>
#include "gaussian.h"
#include "constants.h"
#include <gsl/gsl_multimin.h>

using namespace std;

double ExpecationValueFinal(const gsl_vector * x, double * err) {
	double _A1 = gsl_vector_get(x, 0);
	double _A2 = gsl_vector_get(x, 1);
	double _C1 = gsl_vector_get(x, 2);
	double _C2 = gsl_vector_get(x, 3);

	double A1[]  = {_A1, _A1, _A1};
	double A2[]  = {_A2, _A2, _A2};
	double ** A  = new double*[2];
	A[0]        = A1;
	A[1]        = A2;

	double s1[] = {0.0, 0.0, 0.0};
	double s2[] = {0.0, 0.0, 0.0};
	double ** s = new double*[2];
	s[0]       = s1;
	s[1]       = s2;

	double C[]  = {_C1, _C2};
	double R1[] = {0.0, 0.0, 0.0};
	double ** R = new double*[1];
	R[0]       = R1;

	double Q[] = {N_qe};
	
	GaussianWFN wavefn(2, 1);
	wavefn.A  = A;
	wavefn.s  = s;
	wavefn.C  = C;
	wavefn.R  = R;
	wavefn.Nu = 1;
	wavefn.Q  = Q; 

	double res = wavefn.getHamiltonianExpectation(err);
	return res;
}

double ExpectationValue(const gsl_vector * x, void * params) {
	double _A1 = gsl_vector_get(x, 0);
	double _A2 = gsl_vector_get(x, 1);
	double _C1 = gsl_vector_get(x, 2);
	double _C2 = gsl_vector_get(x, 3);

	double A1[]  = {_A1, _A1, _A1};
	double A2[]  = {_A2, _A2, _A2};
	double ** A  = new double*[2];
	A[0]        = A1;
	A[1]        = A2;

	double s1[] = {0.0, 0.0, 0.0};
	double s2[] = {0.0, 0.0, 0.0};
	double ** s = new double*[2];
	s[0]       = s1;
	s[1]       = s2;

	double C[]  = {_C1, _C2};
	double R1[] = {0.0, 0.0, 0.0};
	double ** R = new double*[1];
	R[0]       = R1;

	double Q[] = {N_qe};
	
	GaussianWFN wavefn(2, 1);
	wavefn.A  = A;
	wavefn.s  = s;
	wavefn.C  = C;
	wavefn.R  = R;
	wavefn.Nu = 1;
	wavefn.Q  = Q; 

	double err;
	double res = wavefn.getHamiltonianExpectation(&err);
	return res;
}

void ExpectationValueGradient(const gsl_vector * x, void * params, gsl_vector * grad) {
	double _A1 = gsl_vector_get(x, 0);
	double _A2 = gsl_vector_get(x, 1);
	double _C1 = gsl_vector_get(x, 2);
	double _C2 = gsl_vector_get(x, 3);
	double dx  = 1e-2; 
	double dx2 = 1e-2;

	// Get the gradient with respect to all four variables.
	double dA1;
	double dA2;
	double dC1;
	double dC2;

	double currentValue = ExpectationValue(x, NULL);

	gsl_vector_set((gsl_vector *)x, 0, _A1 + dx);
	double offsetA1 = ExpectationValue(x, NULL);
	gsl_vector_set((gsl_vector *)x, 0, _A1);
	dA1 = (offsetA1 - currentValue) / dx;

	gsl_vector_set((gsl_vector *)x, 1, _A2 + dx);
	double offsetA2 = ExpectationValue(x, NULL);
	gsl_vector_set((gsl_vector *)x, 1, _A2);
	dA2 = (offsetA2 - currentValue) / dx;

	gsl_vector_set((gsl_vector *)x, 2, _C1 + dx2);
	double offsetC1 = ExpectationValue(x, NULL);
	gsl_vector_set((gsl_vector *)x, 2, _C1);
	dC1 = (offsetC1 - currentValue) / dx2;

	gsl_vector_set((gsl_vector *)x, 3, _C2 + dx2);
	double offsetC2 = ExpectationValue(x, NULL);
	gsl_vector_set((gsl_vector *)x, 3, _C2);
	dC2 = (offsetC2 - currentValue) / dx2;

	gsl_vector_set(grad, 0, dA1);
	gsl_vector_set(grad, 1, dA2);
	gsl_vector_set(grad, 2, dC1);
	gsl_vector_set(grad, 3, dC2);
}

void Combined(const gsl_vector * x, void * params, double * f, gsl_vector * grad) {
 	*f = ExpectationValue(x, NULL);
 	ExpectationValueGradient(x, NULL, grad);
}

int main() {

	int N_A1 = 16;
	int N_A2 = 16;

	double A1_min = 1e5;
	double A1_max = 1e7;
	double A2_min = 1e5;
	double A2_max = 1e7;

	double best = 1.0;

	cout << "i, j, E, error" << endl;

	for (int i = 0; i < N_A1; ++i) {
		double A1_inc     = (A1_max - A1_min) / N_A1;
		double current_A1 = A1_min + i * A1_inc;

		for (int j = 0; j < N_A2; ++j) {
			double A2_inc     = (A2_max - A2_min) / N_A2;
			double current_A2 = A2_max - i * A2_inc;

			size_t iter = 0;
			int status;

			const gsl_multimin_fdfminimizer_type * T;
			gsl_multimin_fdfminimizer * s;

			double params[] = {1.0};

			gsl_vector * x;
			gsl_multimin_function_fdf F;

			F.n      = 4;
			F.f      = &ExpectationValue;
			F.df     = &ExpectationValueGradient;
			F.fdf    = &Combined;
			F.params = params;

			x = gsl_vector_alloc(4);
			gsl_vector_set(x, 0, current_A1);
			gsl_vector_set(x, 1, current_A2);
			gsl_vector_set(x, 2, 1.0);
			gsl_vector_set(x, 3, 1.0);

			T = gsl_multimin_fdfminimizer_conjugate_fr;
			s = gsl_multimin_fdfminimizer_alloc(T, 4);

			gsl_multimin_fdfminimizer_set(s, &F, x, 0.01, 1e-4);

			do {
				iter++;
				status = gsl_multimin_fdfminimizer_iterate(s);

				if (status) break;

				status = gsl_multimin_test_gradient(s->gradient, 1e-3);

			} while (status == GSL_CONTINUE && iter < 100);

			double err;
			double exp = ExpecationValueFinal(x, &err);

			if (exp < best) {
				best = exp;
			}

			cout << i << ", " << j << ", " << exp << ", " << err << endl;

			gsl_multimin_fdfminimizer_free(s);
			gsl_vector_free(x);
		}
	}

	printf("Best Result: %f eV\n", best);

	return 0;
}
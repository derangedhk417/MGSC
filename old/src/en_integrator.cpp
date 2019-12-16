#include <iostream>
#include <mathimf.h>
#include <chrono>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

class ElectronNucleiIntegrator {

public:
	float A1;
	float A2;
	float A3;
	float s1;
	float s2;
	float s3;
	float *R;

	int   size;
	float relerr;

	double x;
	double y;
	double z;

	gsl_integration_workspace * wkspace0;
	gsl_integration_workspace * wkspace1;
	gsl_integration_workspace * wkspace2;

	ElectronNucleiIntegrator(int size) {
		gsl_set_error_handler(ElectronNucleiIntegrator::handler);
		wkspace0 = gsl_integration_workspace_alloc(size);
		wkspace1 = gsl_integration_workspace_alloc(size);
		wkspace2 = gsl_integration_workspace_alloc(size);
		this->size = size;
	}

	~ElectronNucleiIntegrator() {
		gsl_integration_workspace_free(wkspace0);
		gsl_integration_workspace_free(wkspace1);
		gsl_integration_workspace_free(wkspace2);
	}

	// Integrates with the given relative error tolerance.
	float Integrate(double relerr, float * err) {
		this->relerr = relerr;

		gsl_function F0;

		double result;
		double error;

		F0.function = &ElectronNucleiIntegrator::integral0;
		F0.params   = this;

		gsl_integration_qagi(
			&F0, 0, relerr, size, wkspace0, &result, &error
		);

		*err = error;

		return (float)result;
	}

private:
	static double integral0(double x, void * params) {
		ElectronNucleiIntegrator * self = (ElectronNucleiIntegrator *)params;
		self->x = x;

		gsl_function F1;

		double result;
		double error;

		F1.function = &ElectronNucleiIntegrator::integral1;
		F1.params   = self;

		int err = gsl_integration_qagi(
			&F1, 0, self->relerr, 
			self->size, self->wkspace1, 
			&result, &error
		);
	}

	static double integral1(double y, void * params) {
		ElectronNucleiIntegrator * self = (ElectronNucleiIntegrator *)params;
		self->y = y;

		gsl_function F2;

		double result;
		double error;

		F2.function = &ElectronNucleiIntegrator::integrand;
		F2.params   = self;

		int err = gsl_integration_qagi(
			&F2, 0, self->relerr, 
			self->size, self->wkspace2, 
			&result, &error
		);
	}

	static double integrand(double z, void * params) {
		ElectronNucleiIntegrator * self = (ElectronNucleiIntegrator *)params;

		float x1 = self->x;
		float x2 = self->y;
		float x3 = z;

		float A1  = self->A1;
		float A2  = self->A2;
		float A3  = self->A3;
		float s1  = self->s1;
		float s2  = self->s2;
		float s3  = self->s3;
		float * R = self->R;

		float numerator = -(x1*x1)*A1 + x1*s1 - (x2*x2)*A2 + x2*s2 - (x3*x3)*A3 + x3*s3;
		numerator = expf(numerator);

		float d1 = R[0] - x1;
		float d2 = R[1] - x2;
		float d3 = R[2] - x3;

		float denominator = d1*d1 + d2*d2 + d3*d3;
		// denominator = fmax(sqrtf(denominator), 1e-7);
		denominator = sqrtf(denominator);

		return numerator / denominator;
	}

	static void handleError(int err) {
		if (err == GSL_EDIVERGE) {
			std::cout << "Integral Diverges" << std::endl;
		} else if (err == GSL_ESING) {
			std::cout << "Non-Integrable Singularity" << std::endl;
		} else {
			std::cout << "Unspecified Integration Error" << std::endl;
		}
	}

	static void handler(
		const char * reason,
		const char * file,
		int line,
		int gsl_errno
	) {
		// Do nothing.
		// Some of the individual integrals will appear to 
		// be divergent accoring to the algorithm. The full integral
		// in three dimensions is convergent, but the gsl library
		// can't handle it.
	}
};
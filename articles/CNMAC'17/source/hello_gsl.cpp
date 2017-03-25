#include <stdio.h>
#include <stdlib.h>
#include<gsl/gsl_multimin.h>

double fun(const gsl_vector *x, void *params)
{
	double *y = x->data;
	double *p = (double *)params;

	return (y[0] - 1) * (y[0] - 1) + (y[1] - 2.3) * (y[1] - 2.3);
}

void grad(const gsl_vector *x, void *params, gsl_vector *g)
{
	double *y = x->data;
	double *p = (double *)params;

	gsl_vector_set(g, 0, 2.0 * (y[0] - 1));
	gsl_vector_set(g, 1, 2.0 * (y[1] - 2.3));
}

/* Compute both f and g together. */
void fun_grad(const gsl_vector *x, void *params, double *f, gsl_vector *g)
{
	*f = fun(x, params);
	grad(x, params, g);
}

int main(int argc, char **argv)
{
	gsl_vector *x = gsl_vector_alloc(2);

	// init x
	gsl_vector_set(x, 0, 5.2);
	gsl_vector_set(x, 1, 3.6);

	// set function
	gsl_multimin_function_fdf fg;
	fg.n = 2;
	fg.f = fun;
	fg.df = grad;
	fg.fdf = fun_grad;
	fg.params = NULL;

	// optimization type and method
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;

	T = gsl_multimin_fdfminimizer_conjugate_fr;
	s = gsl_multimin_fdfminimizer_alloc(T, 2);

	gsl_multimin_fdfminimizer_set(s, &fg, x, 0.01, 1E-4);

	int iter = 0, status;
	do {
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(s);
		if (status) break;
		status = gsl_multimin_test_gradient(s->gradient, 1E-3);
		double fx = gsl_multimin_fdfminimizer_minimum(s);
		printf("[%2d] x = [%3.2g %3.2g], fx = %g\n",
			iter, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1), fx);
	} while (status == GSL_CONTINUE && iter < 100);

	// free memory
	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(x);

	return EXIT_SUCCESS;
}
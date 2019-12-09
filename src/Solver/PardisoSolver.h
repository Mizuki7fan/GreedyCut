#pragma once
#include "Solver.h"

extern "C" {
	/* PARDISO prototype. */
	void pardisoinit(void   *, int    *, int *, int *, double *, int *);
	void pardiso(void   *, int    *, int *, int *, int *, int *,
		double *, int    *, int *, int *, int *, int *,
		int *, double *, double *, int *, double *);
	void pardiso_chkmatrix(int *, int *, double *, int *, int *, int *);
	void pardiso_chkvec(int *, int *, double *, int *);
	void pardiso_printstats(int *, int *, double *, int *, int *, int *,
		double *, int *);
}

class PardisoSolver: public Solver
{
private:
    /* data */
public:
    PardisoSolver(/* args */);
    ~PardisoSolver();

	void pardiso_init();
	bool factorize();
	void pardiso_solver();
	void free_numerical_factorization_memory();
};
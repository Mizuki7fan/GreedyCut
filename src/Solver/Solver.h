#pragma once
#include <vector>

class Solver
{
private:
	/* data */
public:
	virtual ~Solver() {};
	virtual void pardiso_init() = 0;
	virtual bool factorize() = 0;
	virtual void pardiso_solver() {};
	virtual void free_numerical_factorization_memory() {};

	std::vector<double> result;

	std::vector<int> ia;
	std::vector<int> ja;
	std::vector<double> a;
	std::vector<double> rhs;

	int nnz;
	int num;

protected:
	//pardiso stuff
	/*
	1: real and structurally symmetric, supernode pivoting
	2: real and symmetric positive definite
	-2: real and symmetric indefinite, diagonal or Bunch-Kaufman pivoting
	11: real and nonsymmetric, complete supernode pivoting
	*/
	int mtype;       /* Matrix Type */

	// Remember if matrix is symmetric or not, to
	// decide whether to eliminate the non-upper-
	// diagonal entries from the input II,JJ,SS
	bool is_symmetric;
	bool is_upper_half;

	int nrhs;     /* Number of right hand sides. */
	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void* pt[64];
	/* Pardiso control parameters. */
	int iparm[64];
	double   dparm[64];
	int maxfct, mnum, phase, error, msglvl, solver;
	/* Number of processors. */
	int      num_procs;
	/* Auxiliary variables. */
	char* var;
	int i, k;
	double ddum;          /* Double dummy */
	int idum;         /* Integer dummy. */
};

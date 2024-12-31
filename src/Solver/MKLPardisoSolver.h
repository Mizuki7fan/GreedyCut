#pragma once
#include "Solver.h"
#include "mkl_pardiso.h"
#include "mkl_types.h"

class MKLPardisoSolver : public Solver {
public:
  MKLPardisoSolver();
  ~MKLPardisoSolver();

  void pardiso_init();
  bool factorize();
  void pardiso_solver();
  void free_numerical_factorization_memory();
};
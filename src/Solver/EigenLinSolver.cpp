#include "EigenLinSolver.h"

EigenLinSolver::EigenLinSolver() {}

EigenLinSolver::~EigenLinSolver() {}

void EigenLinSolver::pardiso_init() {}

bool EigenLinSolver::factorize() {
  update_coef();
  simplicialLDLT.compute(coefMtr);
  return true;
}

void EigenLinSolver::pardiso_solver() {
  Eigen::Map<Eigen::VectorXd> b(rhs.data(), rhs.size(), 1);

  Eigen::VectorXd res_ = simplicialLDLT.solve(b);
  result.resize(num);
  for (std::size_t i = 0; i < rhs.size(); i++) {
    result[i] = res_[i];
  }
}

void EigenLinSolver::free_numerical_factorization_memory() {}

void EigenLinSolver::update_coef() {
  typedef Eigen::Triplet<double> T;
  vector<T> tripletvec;
  for (std::size_t i = 0; i < num; i++) {
    int j = ia[i];
    tripletvec.emplace_back(i, i, a[j]);
    for (j = ia[i] + 1; j < ia[i + 1]; j++) {
      tripletvec.emplace_back(i, ja[j], a[j]);
      tripletvec.emplace_back(ja[j], i, a[j]);
    }
  }
  coefMtr.resize(num, num);
  coefMtr.setFromTriplets(tripletvec.begin(), tripletvec.end());
}

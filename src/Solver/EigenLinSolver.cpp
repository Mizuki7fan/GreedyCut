#include "EigenLinSolver.h"
#include <atomic>
#include <cstddef>

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
  for (std::size_t ii = 0; ii < rhs.size(); ii++)
    result[ii] = res_[ii];
}

void EigenLinSolver::free_numerical_factorization_memory() {}

void EigenLinSolver::update_coef() {
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletvec;
  for (int ii = 0; ii < num; ii++) {
    std::size_t j = ia[ii];
    tripletvec.emplace_back(ii, ii, a[j]);
    for (j = ia[ii] + 1; j < ia[ii + 1]; j++) {
      tripletvec.emplace_back(ii, static_cast<int>(ja[j]), a[j]);
      tripletvec.emplace_back(static_cast<int>(ja[j]), ii, a[j]);
    }
  }
  coefMtr.resize(num, num);
  coefMtr.setFromTriplets(tripletvec.begin(), tripletvec.end());
}

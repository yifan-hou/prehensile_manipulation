#include "solvehfvc.h"

#include <iostream>
#include <algorithm>
#include <Eigen/LU>
#include <vector>

#include "eiquadprog.hpp"

#define TOL 1e-9

using std::cout;
using std::endl;
using Eigen::VectorXd;
using Eigen::MatrixXd;

#define MatrixDecomposition Eigen::FullPivLU<MatrixXd>

bool solvehfvc(const MatrixXd &N,
  const MatrixXd &G, const VectorXd &b_G,
  const VectorXd &F,
  const int kDimActualized, const int kDimUnActualized,
  const int kNumSeeds, const int kPrintLevel,
  HFVC *action) {

  /* Size checking */
  const int kDimGeneralized = kDimActualized + kDimUnActualized;
  assert(N.cols() == kDimGeneralized);
  assert(G.rows() == b_G.rows());
  assert(G.cols() == kDimGeneralized);
  assert(F.rows() == kDimGeneralized);

  if (kPrintLevel >= 2) {
    cout << "Begin solving for velocity commands" << endl;
    cout << "  [1] Determine Possible Dimension of control" << endl;
  }

  MatrixXd NG(N.rows()+G.rows(), N.cols());
  NG << N, G;

  // matrix decomposition
  //  There are two things to care about Eigen::FullPivLU.
  //  Firstly, if the matrix to be decomposed is empty or all zero, the kernel()
  //    method will throw a run time error.
  //  Secondly, if the null space has zero dimension, the output is NOT a zero
  //    dimensional matrix, but a n x 1 vector with all zeros.
  MatrixDecomposition lu_decomp_N(N);
  int rank_N = lu_decomp_N.rank();
  MatrixXd basis_N; // columns of basis_N forms a basis of the null-space of N
  if (rank_N > 0)
      basis_N = lu_decomp_N.kernel();
  else
      basis_N = MatrixXd::Identity(kDimGeneralized, kDimGeneralized);

  MatrixDecomposition lu_decomp_NG(NG);
  int rank_NG = lu_decomp_NG.rank();
  assert(rank_NG > 0);

  int n_av = rank_NG - rank_N;
  int n_af = kDimActualized - n_av;
  assert(rank_N + kDimActualized >= kDimGeneralized);

  MatrixXd basis_c;
  if (rank_NG < kDimGeneralized) {
    // null_NG is not empty
    // so C_c is also not empty
    MatrixXd null_NG = lu_decomp_NG.kernel(); // columns of null_NG forms a basis
                                             // of the null-space of NG
    MatrixXd C_c(null_NG.cols()+kDimUnActualized, kDimGeneralized);
    C_c << null_NG.transpose(),
        MatrixXd::Identity(kDimUnActualized,kDimUnActualized),
        MatrixXd::Zero(kDimUnActualized,kDimActualized);
    MatrixDecomposition lu_decomp_C_c(C_c);
    basis_c = lu_decomp_C_c.kernel(); // columns of basis_C_c forms a basis
                                             // of the null-space of C_c
  } else if (kDimUnActualized > 0) {
    MatrixXd C_c(kDimUnActualized, kDimGeneralized);
    C_c << MatrixXd::Identity(kDimUnActualized,kDimUnActualized),
        MatrixXd::Zero(kDimUnActualized,kDimActualized);
    MatrixDecomposition lu_decomp_C_c(C_c);
    basis_c = lu_decomp_C_c.kernel(); // columns of basis_C_c forms a basis
  } else {
    basis_c = MatrixXd::Identity(kDimGeneralized, kDimGeneralized);
  }

  MatrixXd R_a(kDimActualized, kDimActualized);
  MatrixXd T(kDimGeneralized, kDimGeneralized);
  VectorXd w_av;
  if (n_av == 0) {
    if (kPrintLevel >= 2)
      cout << "  [2] n_av = 0, no need for velocity control" << endl;
    R_a = MatrixXd::Identity(kDimActualized, kDimActualized);
    T = MatrixXd::Identity(kDimGeneralized, kDimGeneralized);
    w_av = VectorXd(0);
  } else {
    if (kPrintLevel >= 2)
      cout << "  [2] Solving for Directions by PGD" << endl;
    assert(basis_c.norm() > 0.1);// this shouldn't happen
    int NIter   = 50;
    int n_c     = rank_NG - kDimUnActualized;
    MatrixXd BB = basis_c.transpose()*basis_c;
    MatrixXd NN = basis_N*basis_N.transpose();

    std::vector<MatrixXd> k_all;
    float cost_all[kNumSeeds] = {0};
    for (int seed = 0; seed < kNumSeeds; ++seed)  {
      MatrixXd k  = MatrixXd::Random(n_c, n_av); // initial solution
      MatrixXd bck = basis_c*k;
      for (int i = 0; i < bck.cols(); i++) {
          float bck_col_norm = bck.col(i).norm();
          k.col(i) /= bck_col_norm;
      }
      MatrixXd g(n_c, n_av);
      float costs = 0;
      for (int iter = 0; iter < NIter; ++iter) {
        // compute gradient
        g = MatrixXd::Zero(n_c, n_av);
        costs = 0;
        for (int i = 0; i < n_av; ++i) {
          for (int j = 0; j < n_av; ++j) {
              if (i == j) continue;
              float tempcost = (k.col(i).transpose()*BB*k.col(j)).norm();
              costs += tempcost*tempcost;
              g.col(i) += 2.0f*(k.col(i).transpose()*BB*k.col(j))(0)*BB*k.col(j);
          }
          g.col(i) -= 2.0f*basis_c.transpose()*NN*basis_c*k.col(i);
          costs -= k.col(i).transpose()*basis_c.transpose()*NN*basis_c*k.col(i);
        }
        // descent
        k -= 10.0f*g;
        // project
        bck = basis_c*k;
        for (int i = 0; i < bck.cols(); i++) {
            float bck_col_norm = bck.col(i).norm();
            k.col(i) /= bck_col_norm;
        }
        // cout << "     cost: " << costs << ", grad: " << g.norm() << endl;
      }
      cost_all[seed] = costs;
      k_all.push_back(k);
    }
    float *cost_best    = std::min_element(cost_all, cost_all + kNumSeeds);
    int min_id          = std::distance(cost_all, cost_best);
    MatrixXd k_best     = k_all[min_id];
    MatrixXd C_best     = (basis_c*k_best).transpose();

    // R_a = [null(C_best(:, kDimUnActualized+1:end))';
    //         C_best(:, kDimUnActualized+1:end)];
    // For this decomposition, the input C_best_actualized won't be empty
    // because it has kDimActualized cols; its output basis_C_best_actualized
    // also won't be empty as we are conditioned on n_av > 0
    MatrixXd C_best_actualized = C_best.rightCols(kDimActualized);
    MatrixDecomposition lu_decomp_C_best_actualized(C_best_actualized);
    MatrixXd basis_C_best_actualized;
    int rank_C_best_actualized = lu_decomp_C_best_actualized.rank();
    assert(rank_C_best_actualized > 0);
    basis_C_best_actualized = // columns of basis_C_best_actualized forms the
        lu_decomp_C_best_actualized.kernel(); // null space of C_best_actualized
    R_a = MatrixXd::Zero(kDimActualized, kDimActualized);
    R_a << basis_C_best_actualized.transpose(), C_best_actualized;
    T = MatrixXd::Zero(kDimGeneralized, kDimGeneralized);
    T.topLeftCorner(kDimUnActualized, kDimUnActualized) =
        MatrixXd::Identity(kDimUnActualized, kDimUnActualized);
    T.bottomRightCorner(kDimActualized,kDimActualized) = R_a;

    // b_NG = [zeros(size(N, 1), 1); b_G];
    VectorXd b_NG = VectorXd::Zero(N.rows() + b_G.rows());
    b_NG.tail(b_G.rows()) = b_G;

    // v_star = NG\b_NG;
    VectorXd v_star = NG.fullPivLu().solve(b_NG);
    // cout << "C_best: " << C_best.rows() << ", " << C_best.cols() << endl;
    // cout << C_best;
    w_av = C_best*v_star;
  }

  action->n_av   = n_av;
  action->n_af   = n_af;
  action->R_a    = R_a;
  action->w_av   = w_av;

  return true;
}

// Find a unit orthogonal basis of the column space of A
// A must have full column rank
Eigen::MatrixXd QRWrapper(const Eigen::MatrixXd &A, Eigen::HouseholderQR<Eigen::MatrixXd> *qr) {
  assert(A.cols() <= A.rows());
  qr->compute(A);
  Eigen::MatrixXd thinQ(Eigen::MatrixXd::Identity(A.rows(), A.cols()));
  return qr->householderQ() * thinQ;
}

//  There are three things to care about Eigen::FullPivLU.
//  Firstly, if the matrix to be decomposed is empty or all zero, the kernel()
//    method will throw a run time error.
//  Secondly, if the null space has zero dimension, the output is NOT a zero
//    dimensional matrix, but a n x 1 vector with all zeros.
//  Thirdly, if the input is almost all zero, the rank computation will be too sensitive.
//    e.g. the rank for [0, 0, 1e-15] will be 1, not zero.

// todo: implement manual row reduction to find basis, instead of using QR
// todo: test setThreshold() for decompositions
bool solvehfvc_new(const MatrixXd &N,
  const MatrixXd &G, const VectorXd &b_G,
  const VectorXd &F,
  const int kDimActualized, const int kDimUnActualized,
  const int kNumSeeds, const int kPrintLevel,
  HFVC *action) {
  /* Size checking */
  const int kDimGeneralized = kDimActualized + kDimUnActualized;
  assert(N.cols() == kDimGeneralized);
  assert(G.rows() == b_G.rows());
  assert(G.cols() == kDimGeneralized);
  assert(F.rows() == kDimGeneralized);
  assert(b_G.norm() > TOL);
  Eigen::HouseholderQR<Eigen::MatrixXd> qr;
  Eigen::FullPivLU<MatrixXd> lu;

  /**
   * Step one: gather constraints
   */
  // regularize N
  lu.compute(N.transpose());
  int rank_N = lu.rank();
  if (rank_N == kDimGeneralized) return false;
  Eigen::MatrixXd N_reg = QRWrapper(lu.image(N.transpose()), &qr);
  N_reg.transposeInPlace();

  MatrixXd NG(N_reg.rows()+G.rows(), N_reg.cols());
  NG << N_reg, G;
  VectorXd b_NG = VectorXd::Zero(N_reg.rows() + b_G.rows());
  b_NG.tail(b_G.rows()) = b_G;

  lu.compute(NG);
  int rank_NG = lu.rank();
  assert(rank_NG > 0);
  Eigen::MatrixXd null_space_NG = lu.kernel();
  // check if empty
  if (null_space_NG.norm() == 0) {
    null_space_NG = Eigen::MatrixXd(kDimGeneralized, 0);
  }

  // get a special solution
  Eigen::VectorXd v_star = lu.solve(b_NG);
  // check if no solution
  if ((NG*v_star - b_NG).norm() > TOL) return false;
  /**
   * Step two: Handle un-actuated DOF by adding their linear generators to the solution set
   */
  Eigen::MatrixXd unactuated_linear_generators =
      Eigen::MatrixXd::Identity(kDimGeneralized, kDimUnActualized);

  // add generators from N to reduce unnecessary dims in C
  Eigen::MatrixXd null_space_C(null_space_NG.rows(), null_space_NG.cols() + unactuated_linear_generators.cols());
  null_space_C << null_space_NG, unactuated_linear_generators;

  Eigen::MatrixXd null_space_C_proj_N = N_reg * null_space_C;

  lu.compute(null_space_C_proj_N.transpose());
  int rank_NCProj = lu.rank();
  // check if the added generators will break the goal
  assert(rank_NCProj > 0); // This is probably not going to happen
  assert(rank_NCProj <= kDimUnActualized); // projection all comes from unactuated_linear_generators
  if (rank_NCProj < kDimUnActualized) return false; // no vel action can achieve the goal
  if (rank_NCProj < N_reg.rows()) {
    Eigen::MatrixXd NCProj_comp = lu.kernel(); // might need QR here to get an orthogonal basis
    Eigen::MatrixXd temp(kDimGeneralized, null_space_C.cols() + NCProj_comp.cols());
    temp << null_space_C, N_reg.transpose()*NCProj_comp;
    null_space_C = temp;
  }

  /**
   * Step three: extract control constraint from the solution set
   */
  lu.compute(null_space_C.transpose());
  assert(lu.rank() < kDimGeneralized); // this shouldn't happen as long as b_G is not 0
  Eigen::MatrixXd C = lu.kernel().transpose();
  Eigen::VectorXd b_C = C*v_star;

  // Get the force-controlled directions
  // and get orthogonal
  assert(C.rows() > 0);
  Eigen::MatrixXd C_actualized = C.rightCols(kDimActualized);
  // lu.compute(C_actualized);
  Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(C_actualized.transpose());
  assert(cod.rank() == C.rows());
  Eigen::MatrixXd image = cod.matrixQ();
  C_actualized = image.leftCols(C.rows()).transpose();
  Eigen::MatrixXd force_control_directions = image.rightCols(image.cols() - C.rows()).transpose();

  // R_a is a unitary matrix, R_a^T * R_a = I
  action->R_a = MatrixXd::Zero(kDimActualized, kDimActualized);
  action->R_a << force_control_directions, C_actualized;
  action->n_av   = C.rows();
  action->n_af   = kDimActualized - C.rows();
  action->w_av   = b_C;
  action->C      = C;
  action->b_C    = b_C;

  /**
   * Check the solution
   */
  // std::cout << "Velocity Command:" << std::endl;
  // std::cout << "  C:\n" << C << std::endl;
  // std::cout << "  b_C:\n" << b_C << std::endl;

  // MatrixXd NC(N_reg.rows()+C.rows(), N_reg.cols());
  // NC << N_reg, C;
  // VectorXd b_NC = VectorXd::Zero(N_reg.rows() + b_C.rows());
  // b_NC.tail(b_C.rows()) = b_C;

  // lu.compute(NC);
  // Eigen::MatrixXd sol_homo = lu.kernel();
  // Eigen::VectorXd sol_sp = lu.solve(b_NC);
  // std::cout << "NC solution:" << std::endl;
  // std::cout << "  sol_sp:\n" << sol_sp << std::endl;
  // std::cout << "  sol_homo:\n" << sol_homo << std::endl;
  // std::cout << "  NG*sol_sp:\n" << NG*sol_sp << std::endl;
  // std::cout << "  NG*sol_homo:\n" << NG*sol_homo << std::endl;

  return true;
}
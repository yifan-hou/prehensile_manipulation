#include "polyhedron.h"
#include "setoper.h"
#include "cdd.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include <vector>

bool Poly::vertexEnumeration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, Eigen::MatrixXd *R) {
  /**
   * cddlib initialization
   */
  dd_PolyhedraPtr poly;
  dd_MatrixPtr A_, G_;
  dd_rowrange m;
  dd_colrange d;
  dd_ErrorType err;
  dd_set_global_constants();  /* First, this must be called to use cddlib. */

  /**
   * Eigen to cdd format conversion
   */
  Eigen::MatrixXd A_eigen(A.rows(), A.cols() + 1);
  A_eigen.leftCols<1>() = b;
  A_eigen.rightCols(A.cols()) = A;

  m = A_eigen.rows();
  d = A_eigen.cols();
  A_ = dd_CreateMatrix(m, d);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < d; ++j)
      dd_set_d(A_->matrix[i][j], A_eigen(i, j));
  /**
   * Call vertex enumeration from cddlib
   */
  A_->representation = dd_Inequality;
  poly=dd_DDMatrix2Poly(A_, &err);  /* compute the second (generator) representation */
  if (err!=dd_NoError) {
    dd_WriteErrorMessages(stdout,err);
    dd_free_global_constants();
    return false;
  }
  G_=dd_CopyGenerators(poly);

  // printf("\nInput is H-representation:\n");
  // dd_WriteMatrix(stdout,A_);  printf("\n");
  // dd_WriteMatrix(stdout,G_);

  /**
   * Read results to Eigen format
   */
  // get linearity
  int lin_num = set_card(G_->linset);
  std::vector<int> lin_elements;
  if (lin_num > 0){
    for (long elem=1;elem<=G_->linset[0];elem++) {
      if (set_member(elem,G_->linset)) lin_elements.push_back(int(elem));
    }
  }
  int Grow = G_->rowsize;
  int Gcol = G_->colsize;
  *R = Eigen::MatrixXd(Grow + lin_num, Gcol);
  for (int i = 0; i < Grow; ++i)
    for (int j = 0; j < Gcol; ++j)
      (*R)(i,j) = dd_get_d(G_->matrix[i][j]);
  for (int i = 0; i < lin_num; ++i) {
    R->middleRows<1>(Grow + i) = -R->middleRows<1>(lin_elements[i]-1);
  }

  /**
   * cddlib clean up
   */
  dd_FreeMatrix(A_);
  dd_FreeMatrix(G_);
  dd_FreePolyhedra(poly);

  if (err!=dd_NoError) {
    dd_WriteErrorMessages(stdout,err);
    dd_free_global_constants();
    return false;
  }
  return true;
}

bool Poly::facetEnumeration(const Eigen::MatrixXd &R, Eigen::MatrixXd *A, Eigen::VectorXd *b) {
  /**
   * cddlib initialization
   */
  dd_PolyhedraPtr poly;
  dd_MatrixPtr R_, A_;
  dd_rowrange m;
  dd_colrange d;
  dd_ErrorType err;
  dd_set_global_constants();  /* First, this must be called to use cddlib. */

  /**
   * Eigen to cdd format conversion
   */
  m = R.rows();
  d = R.cols();
  R_ = dd_CreateMatrix(m, d);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < d; ++j)
      dd_set_d(R_->matrix[i][j], R(i, j));
  /**
   * Call face enumeration from cddlib
   */
  R_->representation = dd_Generator;
  poly=dd_DDMatrix2Poly(R_, &err);  /* compute the second (Inequality) representation */
  if (err!=dd_NoError) {
    dd_WriteErrorMessages(stdout,err);
    dd_free_global_constants();
    return false;
  }
  A_=dd_CopyInequalities(poly);

  /**
   * Read results to Eigen format
   */
  int Arow = A_->rowsize;
  int Acol = A_->colsize;
  *A = Eigen::MatrixXd(Arow, Acol-1);
  *b = Eigen::VectorXd(Arow);
  for (int i = 0; i < Arow; ++i) {
    (*b)(i) = dd_get_d(A_->matrix[i][0]);
    for (int j = 1; j < Acol; ++j)
      (*A)(i,j-1) = dd_get_d(A_->matrix[i][j]);
  }

  /**
   * cddlib clean up
   */
  dd_FreeMatrix(R_);
  dd_FreeMatrix(A_);
  dd_FreePolyhedra(poly);

  if (err!=dd_NoError) {
    dd_WriteErrorMessages(stdout,err);
    dd_free_global_constants();
    return false;
  }
  return true;
}

bool Poly::intersection(const Eigen::MatrixXd &R1, const Eigen::MatrixXd &R2, Eigen::MatrixXd *R) {
  /*
    A = [facetEnumeration(R1); facetEnumeration(R2)];
    R = vertexEnumeration(A);
   */
  // get inequalities
  Eigen::MatrixXd A1, A2;
  Eigen::VectorXd b1, b2;
  facetEnumeration(R1, &A1, &b1);
  facetEnumeration(R2, &A2, &b2);
  // stitch
  Eigen::MatrixXd A(A1.rows() + A2.rows(), A1.cols());
  Eigen::VectorXd b(b1.size() + b2.size());
  A << A1, A2;
  b << b1, b2;
  // get generators back
  return vertexEnumeration(A, b, R);
}


#include "polyhedron.h"
#include "setoper.h"
#include "cdd.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include <vector>
#include <Eigen/LU>


#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/Qhull.h"

using orgQhull::Qhull;
using orgQhull::QhullPoint;
using orgQhull::QhullVertexList;
using orgQhull::QhullVertexListIterator;

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
  A_eigen.rightCols(A.cols()) = -A;

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
  // std::cout << "A_: " << std::endl;
  // dd_WriteMatrix(stdout,A_);

  /**
   * Read results to Eigen format
   */
  // get linearity number
  int lin_num = set_card(A_->linset);
  // read id of linearity rows
  std::vector<int> lin_elements;
  if (lin_num > 0){
    for (long elem=1;elem<=A_->linset[0];elem++) {
      if (set_member(elem,A_->linset)) lin_elements.push_back(int(elem));
    }
  }
  // read id of linearity rows
  int Arow = A_->rowsize;
  int Acol = A_->colsize;
  *A = Eigen::MatrixXd(Arow + lin_num, Acol-1);
  *b = Eigen::VectorXd(Arow + lin_num);
  for (int i = 0; i < Arow; ++i) {
    (*b)(i) = dd_get_d(A_->matrix[i][0]);
    for (int j = 1; j < Acol; ++j)
      (*A)(i,j-1) = -dd_get_d(A_->matrix[i][j]);
  }
  for (int i = 0; i < lin_num; ++i) {
    A->middleRows<1>(Arow + i) = -A->middleRows<1>(lin_elements[i]-1);
    b->segment<1>(Arow + i) = -b->segment<1>(lin_elements[i]-1);
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
  if ((R1.rows() == 0)||(R2.rows() == 0)) {
    *R = Eigen::MatrixXd(0, 0);
    return true;
  }
  if ((R1.cols() == 0)||(R2.cols() == 0)) {
    std::cerr << "[intersection] error: inputs has zero cols" << std::endl;
    return false;
  }

  // get inequalities
  Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

  Eigen::MatrixXd A1, A2;
  Eigen::VectorXd b1, b2;
  facetEnumeration(R1, &A1, &b1);
  facetEnumeration(R2, &A2, &b2);

  // stitch
  Eigen::MatrixXd A(A1.rows() + A2.rows(), A1.cols());
  Eigen::VectorXd b(b1.size() + b2.size());
  A << A1, A2;
  b << b1, b2;
  vertexEnumeration(A, b, R);
  // get generators back
  return vertexEnumeration(A, b, R);
}

bool Poly::coneIntersection(const Eigen::MatrixXd &C1, const Eigen::MatrixXd &C2, Eigen::MatrixXd *C) {
  if ((C1.rows() == 0) || (C2.rows() == 0)) {
    *C = Eigen::MatrixXd(0, 0);
    return true;
  }
  if ((C1.cols() == 0) || (C2.cols() == 0)) {
    std::cerr << "[coneIntersection] error: input has zero cols." << std::endl;
    return false;
  }
  Eigen::MatrixXd R1(C1.rows(), C1.cols() + 1);
  Eigen::MatrixXd R2(C2.rows(), C2.cols() + 1);
  R1 << Eigen::VectorXd::Zero(C1.rows()), C1;
  R2 << Eigen::VectorXd::Zero(C2.rows()), C2;

  Eigen::MatrixXd R;
  if (!intersection(R1, R2, &R)) return false;
  *C = R.rightCols(R.cols()-1);
  if (C->norm() > 1e-10)
    assert(R.leftCols<1>().norm() < 1e-10);
  else
    *C = Eigen::MatrixXd(0, 0);

  return true;
}

bool Poly::offsetPolytope(Eigen::MatrixXd *polytope, const Eigen::VectorXd &offset) {
  assert(offset.size() == polytope->cols());
  int num = polytope->rows();
  (*polytope) += Eigen::VectorXd::Ones(num) * offset.transpose();
  return true;
}

Eigen::MatrixXd Poly::convhull(const Eigen::MatrixXd &points) {

  int num_points = points.rows();
  int dim = points.cols();
  // check dimensions
  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(points.transpose());
  // lu_decomp.setThreshold(LU_THRESHOLD);

  int rank = lu_decomp.rank();

  Eigen::MatrixXd data;
  Eigen::MatrixXd basis;
  if (rank == dim) {
    data = points.transpose();
  } else {
    // lu_decomp.image(points): its columns form a basis of the column-space of points'
    basis = lu_decomp.image(points.transpose());
    for (int i = 0; i < basis.cols(); i++) {
        float bck_col_norm = basis.col(i).norm();
        basis.col(i) /= bck_col_norm;
    }
    data = (points * basis).transpose();
  }
  Qhull q("triangle", rank, num_points, data.data(), "");
  // Qhull q("", rank, num_points, data.data(), "");
  QhullVertexList vertices(q.beginVertex(), q.endVertex());
  int num_points_on_hull = vertices.count();
  Eigen::MatrixXd results(num_points_on_hull, rank);
  int row = 0;
  for (QhullVertexListIterator i = vertices; i.hasNext(); ) {
    double *data = i.next().point().coordinates();
    for (int j = 0; j < rank; ++j)
      results(row, j) = data[j];
    row ++;
  }
  if (rank == dim) return results;
  else return results*basis.transpose();
}

bool Poly::convhull(const std::vector<double> &vectors, int dim, int num,
    std::vector<double> *results) {
  Eigen::MatrixXd vec_mat(dim, num);
  vec_mat = Eigen::MatrixXd::Map(vectors.data(), vec_mat.rows(), vec_mat.cols()).transpose();
  Eigen::MatrixXd result_eigen = convhull(vec_mat);

  results->resize(result_eigen.rows()*dim);
  for (int i = 0; i < result_eigen.rows(); ++i)
    for (int j = 0; j < dim; ++j) (*results)[i*dim + j] = result_eigen(i, j);

  return true;
}

bool Poly::minkowskiSumOfVectors(const Eigen::MatrixXd &vectors, Eigen::MatrixXd *results) {
  /**
   * The algorithm:
   * 1. Initialize the set with a vector [0 p1]
   * 2. For each additional vector v, set <- [set, set + v] (double the size)
   *  2.1 For every KCONVHULL_ROUND round of step 2, call convhull to reduce redundancy
   * 3. When all vectors are added, call convhull one last time.
   */
  int num = vectors.rows();
  int dim = vectors.cols();
  if (num <= 1) {
    *results = vectors;
    return true;
  }
  std::vector<double> msum;
  // estimation of storage. Exponential growing only happens KCONVHULL_ROUND.
  msum.reserve(2 * dim * pow(2, KCONVHULL_ROUND) * num);
  // 1. Initialize
  msum.assign(dim, 0);
  for (int i = 0; i < dim; ++i) msum.push_back(vectors(0, i));

  std::vector<double> temp;
  for (int vid = 1; vid < num; ++vid) {
    /**
     * add the new vector vectors(vid, :) to msum
     */
    int msum_num_now = msum.size()/dim;
    std::cout << "msum_num_now: " << msum_num_now << std::endl;
    // self copy
    msum.insert(msum.end(), msum.begin(), msum.end());
    // add the new vector to the copy
    for (int i = 0; i < msum_num_now; ++i)
      for (int j = 0; j < dim; ++j)
        msum[(msum_num_now + i) * dim + j] += vectors(vid, j);
    /**
     * Call convhull to reduce dimension
     */
    if (vid % KCONVHULL_ROUND == 0) {

      std::cout << "calling convhull:\n ";
      convhull(msum, dim, 2*msum_num_now, &temp);
      std::cout << "call convhull, " << 2*msum_num_now << " to " << temp.size()/dim << std::endl;
      msum = temp;
    }
  }
  // Call Convhull one last time
  if ((num-1) % KCONVHULL_ROUND != 0) {
    convhull(msum, dim, msum.size()/dim, &temp);
    msum = temp;
  }

  *results = Eigen::MatrixXd::Map(msum.data(), dim, msum.size()/dim).transpose();
  return true;
}


bool Poly::minkowskiSum(const Eigen::MatrixXd &poly1, const Eigen::MatrixXd &poly2, Eigen::MatrixXd *results) {
  /**
   * 1. Convolution
   */
  int num1 = poly1.rows();
  int num2 = poly2.rows();
  int dim = poly1.cols();
  assert(poly2.cols() == dim);
  Eigen::MatrixXd convolution(num1*num2, dim);
  for (int i = 0; i < num1; i++) {
    /**
     * add a copy of poly2 to vertex i of poly1
     */
    convolution.block(num2*i, 0, num2, dim) = poly2 + Eigen::VectorXd::Ones(num2) * poly1.middleRows<1>(i);
  }
  /**
   * 2. Convex hull
   */
  *results = convolution;
  // *results = convhull(convolution);
  return true;
}

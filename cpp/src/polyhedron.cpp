#include "polyhedron.h"
#include "setoper.h"
#include "cdd.h"
#include "eiquadprog.hpp"
#include "RobotUtilities/utilities.h"

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



double Poly::angBTVec(const Eigen::VectorXd &x, const Eigen::VectorXd &b) {
  return acos(x.normalized().dot(b.normalized()));
}

double Poly::distP2Line(const Eigen::VectorXd &p, const Eigen::VectorXd &n) {
  double k = p.dot(n)/n.dot(n);
  double dist = (p - k*n).norm();
  return (k > 0) ? dist:-dist;
}

Eigen::VectorXd Poly::projectP2Line(const Eigen::VectorXd &p, const Eigen::VectorXd &n) {
  double k = p.dot(n)/n.dot(n);
  return k*n;
}

Eigen::VectorXd Poly::projectP2Hyperplane(const Eigen::VectorXd &p, const Eigen::VectorXd &a, double b) {
  double a_norm = a.norm();
  assert(a_norm > 1e-9);
  return p - (a.dot(p) - b)/a_norm/a_norm*a;
}

double Poly::distRay2ConeFromOutside(const Eigen::VectorXd &p, const Eigen::MatrixXd &A,
    const Eigen::MatrixXd &R) {
  /**
   * 1. Project to each facet, check if the projection is in the cone
   */
  Eigen::VectorXd proj_on_face;
  for (int r = 0; r < A.rows(); ++r) {
    if ((A.middleRows(r, 1)*p)(0) < 0) continue;
    proj_on_face = projectP2Hyperplane(p, A.middleRows(r, 1).transpose(), 0);
    if ((A*proj_on_face).maxCoeff() <= 1e-10) {
      return angBTVec(proj_on_face, p);
    }
  }
  // if not returned yet, the projection is not on any facet
  /**
   * 2. Project p to each generator
   */
  double min_dist = 999999.9;
  for (int r = 0; r < R.rows(); ++r) {
    double ang = angBTVec(R.middleRows(r, 1).transpose(), p);
    if (ang < min_dist) min_dist = ang;
  }
  assert(min_dist < 10);
  return min_dist;
}

double Poly::distP2Polyhedron(const Eigen::VectorXd &p, const Eigen::MatrixXd &A,
    const Eigen::VectorXd &b, const Eigen::VectorXd &x0) {
  // prepare the QP
  //  min 0.5 * x G0 x + g0 x
  //  s.t.
  //     CE^T x + ce0 = 0
  //     CI^T x + ci0 >= 0
  // Input:
  //  min ||x-p||^2
  //  s.t.
  //     Ax <= b
  // Reformulate as:
  //  min 0.5 * x I x -p^T x
  //  s.t.
  //    -A x + b >= 0
  int kDim = p.rows();
  // solve the QP
  Eigen::VectorXd x = x0;
  Eigen::MatrixXd G0 = Eigen::MatrixXd::Identity(kDim, kDim);
  Eigen::VectorXd g0 = -p.transpose();
  double cost = solve_quadprog(G0, g0,
      Eigen::MatrixXd(kDim, 0), Eigen::VectorXd(0), // no equality constraints
      -A.transpose(), b, x);
  // std::cout << "Debug: x: " << x.transpose() << std::endl;
  return (x - p).norm();
}

Eigen::MatrixXd Poly::hitAndRunSampleInPolytope(const Eigen::MatrixXd &A,
    const Eigen::VectorXd &b, const Eigen::VectorXd &x0, int N, int discard, int runup) {
  // https://www.mathworks.com/matlabcentral/fileexchange/34208-uniform-distribution-over-a-convex-polytope
  int dim = x0.rows();
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(N+runup+discard, dim);

  int n = 0; // num generated so far
  Eigen::VectorXd x = x0;
  Eigen::VectorXd M = Eigen::VectorXd::Zero(dim); // Incremental mean.
  Eigen::VectorXd u; // direction
  Eigen::VectorXd v, z, c; // temp
  while (n < N+runup+discard) {
      // test whether in runup or not
      if (n < runup) {
        // same as hitandrun
        u = Eigen::VectorXd::Random(dim).normalized();
      } else {
        // choose a previous point at random
        v = X.middleRows(RUT::randi(n), 1).transpose();
        // line sampling direction is from v to sample mean
        u = (v-M).normalized();
      }
      // proceed as in hit and run
      z = A*u;
      c = (b - A*x).cwiseQuotient(z);
      // std::cout << "\n(b-Ax) " << (b - A*x).transpose() << std::endl;
      // std::cout << "z " << z.transpose() << std::endl;
      // std::cout << "c " << c.transpose() << std::endl;

      // tmin = max(c(z<0));
      // tmax = min(c(z>0));
      double tmin = 99999999, tmax = - 99999999;
      for (int i = 0; i < z.rows(); ++i) {
        if (z(i) > 0) {
          if (c(i) < tmin) tmin = c(i);
        } else {
          if (c(i) > tmax) tmax = c(i);
        }
      }
      // Choose a random point on that line segment
      double distance = tmin+(tmax-tmin)*RUT::rand();
      // std::cout << "debug: tmin: " << tmin << ", tmax: " << tmax  << ", distance: " << distance << std::endl;
      // std::cout << "debug: x: " << x.transpose() << std::endl;
      // Eigen::VectorXd tempresult = A*x - b;
      // if (tempresult.maxCoeff() > 0) std::cout << "debug: INFEASIBLE!!!!!!!! x0 out of polytope!!" << std::endl;
      x = x + distance*u;
      X.middleRows(n,1) = x.transpose();
      n++;

      // Incremental mean and covariance updates
      M = M + (x - M)/n;     // sample mean
  }
  return X.bottomRows(N);
}


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

bool Poly::coneFacetEnumeration(const Eigen::MatrixXd &M, Eigen::MatrixXd *A) {
  Eigen::MatrixXd R(M.rows(), M.cols() + 1);
  R << Eigen::VectorXd::Zero(M.rows()), M;
  Eigen::VectorXd b;
  bool success = facetEnumeration(R, A, &b);
  assert(b.norm() < 1e-10);
  return success;
}

bool Poly::polytopeFacetEnumeration(const Eigen::MatrixXd &M, Eigen::MatrixXd *A,
    Eigen::VectorXd *b) {
  Eigen::MatrixXd R(M.rows(), M.cols() + 1);
  R << Eigen::VectorXd::Ones(M.rows()), M;
  bool success = facetEnumeration(R, A, b);
  return success;
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
  // std::cout << "[debug] A1: " << A1.rows() << " x " << A1.cols() << std::endl;
  // std::cout << "[debug] A2: " << A2.rows() << " x " << A2.cols() << std::endl;
  Eigen::MatrixXd A(A1.rows() + A2.rows(), A1.cols());
  Eigen::VectorXd b(b1.size() + b2.size());
  A << A1, A2;
  b << b1, b2;
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

bool Poly::polytopeIntersection(const Eigen::MatrixXd &P1, const Eigen::MatrixXd &P2, Eigen::MatrixXd *P) {
  if ((P1.rows() == 0) || (P2.rows() == 0)) {
    *P = Eigen::MatrixXd(0, 0);
    return true;
  }
  if ((P1.cols() == 0) || (P2.cols() == 0)) {
    std::cerr << "[polytopeIntersection] error: input has zero cols." << std::endl;
    assert(0);
    return false;
  }
  Eigen::MatrixXd R1(P1.rows(), P1.cols() + 1);
  Eigen::MatrixXd R2(P2.rows(), P2.cols() + 1);
  R1 << Eigen::VectorXd::Ones(P1.rows()), P1;
  R2 << Eigen::VectorXd::Ones(P2.rows()), P2;

  Eigen::MatrixXd R;
  if (!intersection(R1, R2, &R)) return false;
  *P = R.rightCols(R.cols()-1);
  if (P->norm() > 1e-10)
    assert((R.leftCols<1>() - Eigen::VectorXd::Ones(R.rows())).norm() < 1e-10);
  else
    *P = Eigen::MatrixXd(0, 0);

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
  Eigen::MatrixXd points_reduced = points;
  int rank = RUT::rowSpace(&points_reduced);
  Eigen::MatrixXd basis = points_reduced.topRows(rank).transpose();

  Eigen::MatrixXd data;
  if (rank == dim) {
    data = points.transpose();
  } else {
    data = (points * basis).transpose();
  }
  Qhull q("triangle", rank, num_points, data.data(), "");
  // Qhull q("", rank, num_points, data.data(), "");
  QhullVertexList vertices(q.beginVertex(), q.endVertex());
  int num_points_on_hull = vertices.count();
  Eigen::MatrixXd results(num_points_on_hull, rank);
  int row = 0;
  for (QhullVertexListIterator i = vertices; i.hasNext(); ) {
    double *result = i.next().point().coordinates();
    for (int j = 0; j < rank; ++j)
      results(row, j) = result[j];
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
    // std::cout << "msum_num_now: " << msum_num_now << std::endl;
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
      convhull(msum, dim, 2*msum_num_now, &temp);
      msum = temp;
      // std::cout << "call convhull, " << 2*msum_num_now << " to " << temp.size()/dim << std::endl;
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

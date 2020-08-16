#include "polyhedron.h"
#include <ppl.hh> // use PPL instead of CDD
#include <glpk.h> /* GNU GLPK linear/mixed integer solver */

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

namespace PPL = Parma_Polyhedra_Library;


// PPL utilities
void print_constraints(const PPL::Constraint_System& cs,
                       const std::string& intro, std::ostream& s) {
  if (!intro.empty())
    s << intro << "\n";
  PPL::Constraint_System::const_iterator i = cs.begin();
  PPL::Constraint_System::const_iterator cs_end = cs.end();
  bool printed_something = i != cs_end;
  while (i != cs_end) {
    using PPL::IO_Operators::operator<<;
    s << *i++;
    if (i != cs_end)
      s << ",\n";
  }
  s << (printed_something ? "." : "true.") << std::endl;
}

void print_constraints(const PPL::Polyhedron& ph,
                       const std::string& intro, std::ostream& s) {
  print_constraints(ph.constraints(), intro, s);
}

void print_generators(const PPL::Generator_System& gs,
                      const std::string& intro, std::ostream& s) {
  if (!intro.empty())
    s << intro << "\n";
  PPL::Generator_System::const_iterator i = gs.begin();
  PPL::Generator_System::const_iterator gs_end = gs.end();
  bool printed_something = i != gs_end;
  while (i != gs_end) {
    using PPL::IO_Operators::operator<<;
    s << *i++;
    if (i != gs_end)
      s << ",\n";
  }
  s << (printed_something ? "." : "false.") << std::endl;
}

void print_generators(const PPL::Polyhedron& ph,
                      const std::string& intro, std::ostream& s) {
  print_generators(ph.generators(), intro, s);
}


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

// Input:
//  min ||x-p||^2
//  s.t.
//     Ax <= b
double Poly::distP2Polyhedron(const Eigen::VectorXd &p, const Eigen::MatrixXd &A,
    const Eigen::VectorXd &b, const Eigen::VectorXd &x0) {
  // prepare the QP
  //  min 0.5 * x G0 x + g0 x
  //  s.t.
  //     CE^T x + ce0 = 0
  //     CI^T x + ci0 >= 0
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
  return (x - p).norm();
}

Eigen::MatrixXd Poly::hitAndRunSampleInPolytope(const Eigen::MatrixXd &A,
    const Eigen::VectorXd &b, const Eigen::VectorXd &x0, int N, int discard, int runup, double max_radius) {
  // https://www.mathworks.com/matlabcentral/fileexchange/34208-uniform-distribution-over-a-convex-polytope
  int dim = x0.rows();
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(N+runup+discard, dim);
  if ((A*x0 - b).maxCoeff() > 0) {
    std::cout << "[hitAndRunSampleInPolytope] x0 is outside of Ax<b" << std::endl;
    exit(1);
  }
  if ((max_radius > 0) && (x0.norm() > max_radius)) {
    std::cout << "[hitAndRunSampleInPolytope] x0 is outside of max_radius" << std::endl;
    exit(1);
  }

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

      // tmin = max(c(z<0));
      // tmax = min(c(z>0));
      double tmin = 9999999;
      double tmax = -9999999;
      if (max_radius > 0) {
        double a = u.dot(u);
        double bb = 2.0*x.dot(u);
        double c = x.dot(x) - max_radius*max_radius;
        double delta = sqrt(bb*bb - 4.0*a*c);
        tmin = (-bb + delta)/2/a;
        tmax = (-bb - delta)/2/a;
      }
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
      x = x + distance*u;
      if (max_radius > 0) {
        // std::cout << "x.norm(): " << x.norm() << ", max_radius: " << max_radius << std::endl;
        if(x.norm() >= max_radius) {
          std::cout << "[hitAndRunSampleInPolytope] assertion failed." << std::endl;
          getchar();
          exit(1);
        }
      }
      X.middleRows(n,1) = x.transpose();
      n++;

      // Incremental mean and covariance updates
      M = M + (x - M)/n;     // sample mean
  }
  return X.bottomRows(N);
}


bool Poly::vertexEnumeration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
    const Eigen::MatrixXd &Ae, const Eigen::VectorXd &be,  Eigen::MatrixXd *R) {
  /**
   * Convert the inputs to integers
   */
  Eigen::MatrixXi A_int = (A*PPL_MULTIPLIER_BEFORE_ROUNDING).cast<int>();
  Eigen::VectorXi b_int = (b*PPL_MULTIPLIER_BEFORE_ROUNDING).cast<int>();
  Eigen::MatrixXi Ae_int = (Ae*PPL_MULTIPLIER_BEFORE_ROUNDING).cast<int>();
  Eigen::VectorXi be_int = (be*PPL_MULTIPLIER_BEFORE_ROUNDING).cast<int>();
  /**
   * Eigen to PPL format conversion
   */
  PPL::Constraint_System cs;
  int dim;
  int dim1 = 0, dim2 = 0;

  int nrows = A_int.rows();
  if (nrows > 0) {
    dim1 = A_int.cols();
    for (int i = 0; i < nrows; ++i) {
      PPL::Linear_Expression e;
      for (unsigned j = dim1; j > 0; j--) {
        e += A_int(i,j-1) * PPL::Variable(j-1);
      }
      e -= b_int(i);
      cs.insert(e <= 0);
    }
  }
  nrows = Ae_int.rows();
  if (nrows > 0) {
    dim2 = Ae_int.cols();
    for (int i = 0; i < nrows; ++i) {
      PPL::Linear_Expression e;
      for (unsigned j = dim2; j > 0; j--) {
        e += Ae_int(i,j-1) * PPL::Variable(j-1);
      }
      e -= be_int(i);
      cs.insert(e == 0);
    }
  }
  if (dim1 != 0) {
    dim = dim1;
  } else if (dim2 != 0) {
    dim = dim2;
  } else {
    return true; // input is empty, nothing to do
  }
  PPL::C_Polyhedron ph(cs);

  /**
   * Call vertex enumeration from cddlib
   */
  ph.minimized_generators();
  // print_constraints(ph, "*** ph constraints ***", std::cout);
  // print_generators(ph, "*** ph generators ***", std::cout);
  /**
   * Read results to Eigen format
   */
  PPL::Generator_System gs = ph.generators();
  std::vector<double> R_vec;
  std::vector<int> t_vec;
  std::vector<double> R_row;
  std::vector<int> t_row;
  PPL::Generator_System::const_iterator ig = gs.begin();
  PPL::Generator_System::const_iterator gs_end = gs.end();

  while (ig != gs_end) {
    R_row.resize(dim);
    t_row.resize(1);
    t_row[0] = 0;
    if (ig->is_point()) {
      t_row[0] = 1;
      mpz_class divisor = ig->divisor();
      for (PPL::dimension_type j = ig->space_dimension(); j-- > 0; ) {
        mpq_class q(ig->coefficient(PPL::Variable(j)), divisor);
        R_row[j] = q.get_d();
      }
    } else if (ig->is_ray() || ig->is_line()) {
      mpz_class max = abs(ig->coefficient(PPL::Variable(0)));
      for (PPL::dimension_type j = ig->space_dimension(); j-- > 0; ) {
        if (cmp(max, abs(ig->coefficient(PPL::Variable(j)))) < 0)
          max = abs(ig->coefficient(PPL::Variable(j)));
      }
      for (PPL::dimension_type j = ig->space_dimension(); j-- > 0; ) {
        R_row[j] = mpq_class(ig->coefficient(PPL::Variable(j)), max).get_d();
      }
    } else {
      std::cout << "[vertexEnumeration] A Closure point! Not implemented yet, should be the same as point. Make sure you know why there is a Closure point\n";
      exit(1);
    }

    if (ig->is_line()) {
      R_row.resize(2*dim);
      for (int ii = 0; ii < dim; ++ii)
        R_row[dim + ii] = -R_row[ii];
      t_row.assign(2, 0);
    }

    for (int ii = 0; ii < R_row.size(); ++ii)
      R_vec.push_back(R_row[ii]);
    for (int ii = 0; ii < t_row.size(); ++ii)
      t_vec.push_back(t_row[ii]);
    ig++;
  }

  int num_generators = R_vec.size()/dim;
  Eigen::MatrixXd R_ = Eigen::MatrixXd::Map(R_vec.data(), dim, num_generators).transpose();
  Eigen::VectorXi t_ = Eigen::VectorXi::Map(t_vec.data(), t_vec.size());

  *R = Eigen::MatrixXd(R_.rows(), R_.cols() + 1);
  R->leftCols(1) = t_.cast<double>();
  R->rightCols(R_.cols()) = R_;
  return true;
}

bool Poly::facetEnumeration(const Eigen::MatrixXd &R_input, Eigen::MatrixXd *A, Eigen::VectorXd *b) {
  /**
   * Scale input matrix to make it numerically more stable
   */
  // Eigen::MatrixXd S = Eigen::MatrixXd::Identity(R_input.cols(), R_input.cols());
  // double inf_value = 0;
  // if (R_input.leftCols(1).minCoeff() > 0.1) {
  //   // no rays
  //   inf_value = std::numeric_limits<double>::infinity();
  // }
  // for (int i = 1; i < R_input.cols(); ++i) {
  //   double max = -inf_value;
  //   double min = inf_value;
  //   for (int j = 0; j < R_input.rows(); ++j) {
  //     if (R_input(j, i) > max) max = R_input(j, i);
  //     if (R_input(j, i) < min) min = R_input(j, i);
  //   }
  //   double range = max - min;
  //   if ( range > 1e-8) {
  //     S(i,i) = 1.0/range;
  //   }
  // }

  /**
   * Convert the inputs to integers
   */
  Eigen::MatrixXi R_int = (R_input*PPL_MULTIPLIER_BEFORE_ROUNDING).cast<int>();
  R_int.leftCols(1) = R_input.leftCols(1).cast<int>();

  /**
   * Convert input to PPL format
   */
  int nrows = R_int.rows();
  int dim = R_int.cols() - 1;
  PPL::Constraint_System cs;
  PPL::Generator_System gs;

  bool has_a_point = false;
  for (int i = 0; i < nrows; ++i) {
    PPL::Linear_Expression e;
    for (unsigned j = dim; j > 0; j--) {
      e += R_int(i,j) * PPL::Variable(j-1);
    }
    if (R_int(i,0) == 1) {
      // a vertex
      has_a_point = true;
      gs.insert(point(e, PPL_MULTIPLIER_BEFORE_ROUNDING));
    } else {
      gs.insert(ray(e));
    }
  }

  // Every non-empty generator system must have at least one point.
  if (nrows > 0 && !has_a_point) {
    gs.insert(PPL::point());
  }

  PPL::C_Polyhedron ph(gs);

  ph.minimized_constraints(); // V to H
  // print_generators(ph, "*** ph generators ***", std::cout);
  // print_constraints(ph, "*** ph constraints ***", std::cout);

  cs = ph.constraints();

  // convert constraints to Eigen Matrix format
  // from Ax + b >= 0
  // to   Ax < b
  std::vector<double> A_vec;
  std::vector<double> b_vec;
  std::vector<double> A_row;
  std::vector<double> b_row;
  PPL::Constraint_System::const_iterator ic = cs.begin();
  PPL::Constraint_System::const_iterator cs_end = cs.end();

  while (ic != cs_end) {
    A_row.resize(dim);
    b_row.resize(1);

    mpz_class b = ic->inhomogeneous_term();
    mpz_class max = abs(b);
    for (PPL::dimension_type j = ic->space_dimension(); j-- > 0; ) {
      if (cmp(max, abs(ic->coefficient(PPL::Variable(j)))) < 0) {
        max = abs(ic->coefficient(PPL::Variable(j)));
      }
    }

    for (PPL::dimension_type j = ic->space_dimension(); j-- > 0; ) {
      A_row[j] = - mpq_class(ic->coefficient(PPL::Variable(j)), max).get_d();
    }
    b_row[0] = mpq_class(b, max).get_d();

    if (ic->is_equality()) {
      A_row.resize(2*dim);
      for (int ii = 0; ii < dim; ++ii)
        A_row[dim + ii] = -A_row[ii];
      b_row.push_back(-b_row[0]);
    }

    for (int ii = 0; ii < A_row.size(); ++ii)
      A_vec.push_back(A_row[ii]);
    for (int ii = 0; ii < b_row.size(); ++ii)
      b_vec.push_back(b_row[ii]);
    ic++;
  }

  int num_generators = A_vec.size()/dim;
  *A = Eigen::MatrixXd::Map(A_vec.data(), dim, num_generators).transpose();
  *b = Eigen::VectorXd::Map(b_vec.data(), b_vec.size());

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
  if ((R1.rows() == 0)||(R2.rows() == 0)) {
    *R = Eigen::MatrixXd(0, 0);
    return true;
  }
  if ((R1.cols() == 0)||(R2.cols() == 0)) {
    std::cerr << "[intersection] error: inputs has zero cols" << std::endl;
    return false;
  }

  /**
   * Convert the inputs to integers
   */
  Eigen::MatrixXi R1_int = (R1*PPL_MULTIPLIER_BEFORE_ROUNDING).cast<int>();
  Eigen::MatrixXi R2_int = (R2*PPL_MULTIPLIER_BEFORE_ROUNDING).cast<int>();
  R1_int.leftCols(1) = R1.leftCols(1).cast<int>();
  R2_int.leftCols(1) = R2.leftCols(1).cast<int>();
  /**
   * Convert input to PPL format
   */
  PPL::Generator_System gs1;
  PPL::Generator_System gs2;

  int nrows = R1_int.rows();
  int dim = R1_int.cols() - 1;
  bool has_a_point = false;
  for (int i = 0; i < nrows; ++i) {
    PPL::Linear_Expression e;
    for (unsigned j = dim; j > 0; j--) {
      e += R1_int(i,j) * PPL::Variable(j-1);
    }
    if (R1_int(i,0) == 1) {
      // a vertex
      has_a_point = true;
      gs1.insert(point(e, PPL_MULTIPLIER_BEFORE_ROUNDING));
    } else {
      gs1.insert(ray(e));
    }
  }
  // Every non-empty generator system must have at least one point.
  if (nrows > 0 && !has_a_point) {
    gs1.insert(PPL::point());
  }

  nrows = R2_int.rows();
  dim = R2_int.cols() - 1;
  has_a_point = false;
  for (int i = 0; i < nrows; ++i) {
    PPL::Linear_Expression e;
    for (unsigned j = dim; j > 0; j--) {
      e += R2_int(i,j) * PPL::Variable(j-1);
    }
    if (R2_int(i,0) == 1) {
      // a vertex
      has_a_point = true;
      gs2.insert(point(e, PPL_MULTIPLIER_BEFORE_ROUNDING));
    } else {
      gs2.insert(ray(e));
    }
  }
  // Every non-empty generator system must have at least one point.
  if (nrows > 0 && !has_a_point) {
    gs2.insert(PPL::point());
  }


  PPL::C_Polyhedron ph1(gs1);
  PPL::C_Polyhedron ph2(gs2);
  ph1.minimized_constraints();
  ph2.minimized_constraints();

  /**
   * Intersection
   */
  // print_generators(ph1, "*** ph1 generators ***", std::cout);
  // print_generators(ph2, "*** ph2 generators ***", std::cout);

  // ph1.add_constraints(ph2.constraints());
  ph1.intersection_assign(ph2);
  // print_constraints(ph1, "*** ph constraints ***", std::cout);
  // print_generators(ph1, "*** ph generators ***", std::cout);


  /**
   * Read the results
   * Below are mostly copied from vertexEnumeration
   */
  ph1.minimized_generators();
  /**
   * Read results to Eigen format
   */
  PPL::Generator_System gs = ph1.generators();
  std::vector<double> R_vec;
  std::vector<int> t_vec;
  std::vector<double> R_row;
  std::vector<int> t_row;
  PPL::Generator_System::const_iterator ig = gs.begin();
  PPL::Generator_System::const_iterator gs_end = gs.end();

  while (ig != gs_end) {
    R_row.resize(dim);
    t_row.resize(1);
    t_row[0] = 0;
    if (ig->is_point()) {
      t_row[0] = 1;
      mpz_class divisor = ig->divisor();
      for (PPL::dimension_type j = ig->space_dimension(); j-- > 0; ) {
        mpq_class q(ig->coefficient(PPL::Variable(j)), divisor);
        R_row[j] = q.get_d();
      }
    } else if (ig->is_ray() || ig->is_line()) {
      mpz_class max = abs(ig->coefficient(PPL::Variable(0)));
      for (PPL::dimension_type j = ig->space_dimension(); j-- > 0; ) {
        if (cmp(max, abs(ig->coefficient(PPL::Variable(j)))) < 0)
          max = abs(ig->coefficient(PPL::Variable(j)));
      }
      for (PPL::dimension_type j = ig->space_dimension(); j-- > 0; ) {
        R_row[j] = mpq_class(ig->coefficient(PPL::Variable(j)), max).get_d();
      }
    } else {
      std::cout << "[Intersection] A Closure point! Not implemented yet, should be the same as point. Make sure you know why there is a Closure point\n";
      exit(1);
    }

    if (ig->is_line()) {
      R_row.resize(2*dim);
      for (int ii = 0; ii < dim; ++ii)
        R_row[dim + ii] = -R_row[ii];
      t_row.assign(2, 0);
    }

    for (int ii = 0; ii < R_row.size(); ++ii)
      R_vec.push_back(R_row[ii]);
    for (int ii = 0; ii < t_row.size(); ++ii)
      t_vec.push_back(t_row[ii]);
    ig++;
  }

  int num_generators = R_vec.size()/dim;
  Eigen::MatrixXd R_ = Eigen::MatrixXd::Map(R_vec.data(), dim, num_generators).transpose();
  Eigen::VectorXi t_ = Eigen::VectorXi::Map(t_vec.data(), t_vec.size());

  *R = Eigen::MatrixXd(R_.rows(), R_.cols() + 1);
  R->leftCols(1) = t_.cast<double>();
  R->rightCols(R_.cols()) = R_;

  return true;
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
  // get rid of the origin
  int id0 = -1;
  for (int i = 0; i < R.rows(); ++i) {
    if (int(R(i, 0)) == 1) {
      id0 = i;
      break;
    }
  }
  Eigen::MatrixXd R_;
  if (id0 >= 0) {
    R_ = Eigen::MatrixXd(R.rows() - 1, R.cols());
    R_.topRows(id0) = R.topRows(id0);
    R_.bottomRows(R.rows() - id0 - 1) = R.bottomRows(R.rows() - id0 - 1);
  }

  *C = R_.rightCols(R_.cols()-1);
  if (C->norm() > 1e-10) {
    assert(R_.leftCols<1>().norm() < 1e-10);
    C->rowwise().normalize();
  } else {
    *C = Eigen::MatrixXd(0, 0);
  }

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

// todo: make this a full lp. Make a separate function for
//  feasibility check
bool Poly::lpfeasibility(const Eigen::MatrixXd &A,
    const Eigen::VectorXd &b, Eigen::VectorXd *xs) {
  /* declare variables */
  glp_prob *lp;
  glp_smcp parm;
  glp_init_smcp(&parm);
  parm.presolve = GLP_OFF;
  parm.msg_lev = GLP_MSG_ERR; // error and warning only

  int *ia, *ja;
  double *ar;
  int rows = A.rows();
  int cols = A.cols();
  int size = rows * cols;
  ia = new int[size + 1000];
  ja = new int[size + 1000];
  ar = new double[size + 1000];
  double z, x1, x2;

  /* create problem */
  lp = glp_create_prob();
  glp_set_obj_dir(lp, GLP_MAX);

  /* fill problem */
  glp_add_rows(lp, rows);
  for (int r = 1; r <= rows; ++r) {
    glp_set_row_bnds(lp, r, GLP_UP, 0.0, b(r-1));
  }

  // set variable bounds and cost function
  glp_add_cols(lp, cols);
  for (int c = 1; c <= cols; ++c) {
    glp_set_col_bnds(lp, c, GLP_FR, 0.0, 0.0); // no boundary
    glp_set_obj_coef(lp, c, 0.0); // no cost function, so coef = 0
  }

  // fill in coefficient matrix
  int id = 0;
  for (int r = 1; r <= rows; ++r) {
    for (int c = 1; c <= cols; ++c) {
      id = (r-1)*cols + c;
      ia[id] = r, ja[id] = c, ar[id] = A(r-1, c-1);
    }
  }
  glp_load_matrix(lp, id, ia, ja, ar);
  /* solve problem */
  glp_simplex(lp, &parm);
  int result = glp_get_status(lp);

  /* housekeeping */
  glp_delete_prob(lp);
  glp_free_env();
  delete [] ia;
  delete [] ja;
  delete [] ar;

  if ((result == GLP_OPT) || (result == GLP_FEAS)) {
    // feasible
    // z = glp_get_obj_val(lp);
    for (int d = 0; d < cols; ++d) {
      (*xs)(d) = glp_get_col_prim(lp, d + 1);
    }
    // std::cout << "z: " << z << std::endl;
    // std::cout << "solution: " << xs.transpose() << std::endl;
    return true;
  } else {
    return false;
    // std::cout << "Infeasible." << std::endl;
  }

}
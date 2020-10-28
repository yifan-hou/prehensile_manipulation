#include "polyhedron.h"
#include <glpk.h> /* GNU GLPK linear/mixed integer solver */

#include "eiquadprog.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cmath>
#include <string.h>

#include <vector>
#include <Eigen/LU>

#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/Qhull.h"
#include "RobotUtilities/utilities.h"

using orgQhull::Qhull;
using orgQhull::QhullPoint;
using orgQhull::QhullVertexList;
using orgQhull::QhullVertexListIterator;

// pairs definition, used for trackign index during sorting
typedef std::pair<int,int> int_pair;
bool int_comparator_descending ( const int_pair& l, const int_pair& r) { return l.first > r.first; }

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
    const Eigen::VectorXd &b, const Eigen::VectorXd &x0, Eigen::VectorXd *x_closest) {
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
  if (x_closest != nullptr)
    *x_closest = x;
  return (x - p).norm();
}

// x must be feasible. For each P2, compute the closest point p2 to it, then form a half space constraint:
//    v'(x - p2) > 0
// where v = normalize(x - p2).
// This constraint can be written as Ax < b, where
//    A = -V',
//    b = -V'p2
double Poly::getAwayFromPolyhedrons(
    const std::vector<Eigen::MatrixXd> &A2, const std::vector<Eigen::VectorXd> &b2,
    const Eigen::MatrixXd &A1, const Eigen::VectorXd &b1,
    const Eigen::VectorXd &xl, const Eigen::VectorXd &xu,
    Eigen::VectorXd *x) {
  /**
   * Compute the closest point in each P2
   */
  int num_of_P2 = A2.size();
  int dim = x->rows();
  if (num_of_P2 > 0) {
    Eigen::MatrixXd A(num_of_P2, dim);
    Eigen::VectorXd b(num_of_P2);
    Eigen::VectorXd x_closest, v;
    for (int p2 = 0; p2 < A2.size(); ++p2) {
      double dist = distP2Polyhedron(*x, A2[p2], b2[p2], Eigen::VectorXd::Zero(dim), &x_closest);
      if (dist < 0) return -1; // input x must be feasible
      v = *x - x_closest;
      v.normalize();
      A.middleRows(p2, 1) = -v.transpose();
      b(p2) = -v.dot(x_closest);
    }
    /**
     * Solve for the largest inscribed sphere
     */
    return inscribedSphere(A, b, xl, xu, A1, b1, x);
  } else {
    // no P2. Return the center of P1
    assert(A1.rows() > 0); //P1 cannot also be empty
    Eigen::MatrixXd A = A1;
    Eigen::VectorXd b = b1;
    // normalize A, b
    for (int i = 0; i < A1.rows(); ++i) {
      double norm = A1.middleRows(i, 1).norm();
      A.middleRows(i, 1) /= norm;
      b(i) /= norm;
    }
    return inscribedSphere(A, b, xl, xu, Eigen::MatrixXd(0, dim), Eigen::VectorXd(0), x);
  }
}

Eigen::MatrixXd Poly::hitAndRunSampleInPolytope(const Eigen::MatrixXd &A,
    const Eigen::VectorXd &b, const Eigen::VectorXd &x0, int N, int discard, int runup, double max_radius) {
  // https://www.mathworks.com/matlabcentral/fileexchange/34208-uniform-distribution-over-a-convex-polytope
  int dim = x0.rows();
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(N+runup+discard, dim);
  if ((A*x0 - b).maxCoeff() > 1e-9) {
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
        // u = (v-M).normalized();
        u = v-M;
        if (u.norm() < 1e-9) {
          u = Eigen::VectorXd::Random(dim).normalized();
        } else {
          u.normalize();
        }
      }
      // proceed as in hit and run
      z = A*u;

      c = (b - A*x).cwiseQuotient(z);

      // tmin = max(c(z<0));
      // tmax = min(c(z>0));
      double tmin = -9999999;
      double tmax = 9999999;
      if (max_radius > 0) {
        double a = u.dot(u);
        double bb = 2.0*x.dot(u);
        double c = x.dot(x) - max_radius*max_radius;
        double delta = sqrt(bb*bb - 4.0*a*c);
        tmin = (-bb - delta)/2/a;
        tmax = (-bb + delta)/2/a;
      }
      for (int i = 0; i < z.rows(); ++i) {
        if (z(i) < 0) {
          if (c(i) > tmin) tmin = c(i);
        } else {
          if (c(i) < tmax) tmax = c(i);
        }
      }
      // Choose a random point on that line segment

      double distance = tmin+(tmax-tmin)*RUT::rand();
      x = x + distance*u;
      if (max_radius > 0) {
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

std::vector<Eigen::VectorXd> Poly::sampleInP1OutOfP2(const Eigen::MatrixXd &A1,
    const Eigen::VectorXd &b1, const std::vector<Eigen::MatrixXd> &A2,
    const std::vector<Eigen::VectorXd> &b2, const Eigen::VectorXd &x0,
    int N, double max_radius) {
  int dim = x0.rows();
  if ((A1*x0 - b1).maxCoeff() > 0) {
    std::cout << "[sampleInP1OutOfP2] x0 is outside of P1" << std::endl;
    exit(1);
  }
  if ((max_radius > 0) && (x0.norm() > max_radius)) {
    std::cout << "[sampleInP1OutOfP2] x0 is outside of max_radius" << std::endl;
    exit(1);
  }

  // Eigen::MatrixXd covar = 1000.0*Eigen::MatrixXd::Identity(dim,dim);
  // RUT::normal_random_variable sample { covar };
  Eigen::VectorXd vec_origin = x0.normalized();
  Eigen::VectorXd u, z, c; // temps
  std::vector<Eigen::VectorXd> X;
  std::vector<int_pair> P2_ranking;
  for (int i = 0; i < A2.size(); ++i) {
    P2_ranking.push_back(std::make_pair (0, i));
  }
  for (int s = 0; s < N; ++s) {
    // u = Eigen::VectorXd::Random(dim);
    // u = sample();
    while (true) {
      u = Eigen::VectorXd::Random(dim);
      if (u.norm() < 1.0) {
        break;
      }
    }
    // u = u - u.dot(vec_origin)*vec_origin; // (TODO) optional
    u.normalize();
    z = A1*u;
    c = (b1 - A1*x0).cwiseQuotient(z);

    /**
     * Compute the maximal motion along this direction within P1
     */
    // std::cout << "Direction: " << u.transpose() << std::endl;
    double tmin = -9999999;
    double tmax = 9999999;
    if (max_radius > 0) {
      double aa = u.dot(u);
      double bb = 2.0*x0.dot(u);
      double cc = x0.dot(x0) - max_radius*max_radius;
      double delta = sqrt(bb*bb - 4.0*aa*cc);
      tmin = (-bb - delta)/2/aa;
      tmax = (-bb + delta)/2/aa;
    }
    // tmin = max(c(z<0));
    // tmax = min(c(z>0));
    for (int i = 0; i < z.rows(); ++i) {
      if (z(i) < 0) {
        if (c(i) > tmin) tmin = c(i);
      } else {
        if (c(i) < tmax) tmax = c(i);
      }
    }
    // std::cout << "tmin: " << tmin << ", tmax: " << tmax << std::endl;
    std::vector<double> segments;
    segments.push_back(tmin);
    segments.push_back(tmax);
    for (int p = 0; p < A2.size(); ++p) {
      int id = P2_ranking[p].second;
      z = A2[id]*u;
      c = (b2[id] - A2[id]*x0).cwiseQuotient(z);
      /**
       * Compute the maximal motion along this direction within P2
       */
      tmin = -9999999;
      tmax = 9999999;
      // tmin = max(c(z<0));
      // tmax = min(c(z>0));
      for (int i = 0; i < z.rows(); ++i) {
        if (z(i) < 0) {
          if (c(i) > tmin) tmin = c(i);
        } else {
          if (c(i) < tmax) tmax = c(i);
        }
      }
      if (tmin > tmax) {
        // this line has no intersection with the polyhedron
        continue;
      }
      // std::cout << "round " << p << std::endl;
      // std::cout << "P2 tmin: " << tmin << ", tmax: " << tmax << std::endl;
      /**
       * Remove from the segment [tmin, tmax] any part that belongs to P2
       */
      for (int seg = 0; seg < segments.size()/2; ) {
        int i = 2*seg;
        int j = 2*seg + 1;
        double begin = segments[i];
        double end = segments[j];
        if (begin >= tmax) break;
        if (end <= tmin) {
          seg ++;
          continue;
        }
        if ((begin >= tmin) && (end <= tmax)) {
          segments.erase(segments.begin() + i, segments.begin() + i + 2);
          if (segments.size() == 0) break;
          continue;
        }
        if (begin >= tmin) {
          segments[i] = tmax;
          break;
        } else if (end <= tmax) {
          segments[j] = tmin;
          seg ++;
        } else {
          segments[j] = tmin;
          segments.insert(segments.begin() + j + 1, end);
          segments.insert(segments.begin() + j + 1, tmax);
          break;
        }
      } // end this P2
      // std::cout << "Segments: ";
      // for (int ii = 0; ii < segments.size(); ++ii) std::cout << segments[ii] << ", ";
      // std::cout << std::endl;
      if (segments.size() == 0) {
        // this line dies on the current P2. record it
        P2_ranking[p].first += 1;
        break;
      }
    } // end all P2
    // std::cout << "Finished processing all P2. Segments: ";
    // for (int ii = 0; ii < segments.size(); ++ii) std::cout << segments[ii] << ", ";
    // std::cout << std::endl;
    if (segments.size() == 0) {
      std::sort(P2_ranking.begin(), P2_ranking.end(), int_comparator_descending);
      continue;
    }
    // Pick the middle point of the longest segment
    double distance = 0;
    double max_length = 0;
    for (int seg = 0; seg < segments.size()/2; ++seg) {
      double seg_length = segments[2*seg + 1] - segments[2*seg];
      if (seg_length > max_length) {
        max_length = seg_length;
        distance = 0.5*(segments[2*seg + 1] + segments[2*seg]);
      }
    }
    // std::cout << "Distance: " << distance << std::endl;
    // std::cout << "x: " << (x0 + distance*u).transpose() << std::endl;
    // getchar();
    X.push_back(x0 + distance*u);
  } // end one direction
  // std::cout << "P2_ranking: " << std::endl;
  // for (int i = 0; i < P2_ranking.size(); ++i) {
  //   std::cout << "id: " << P2_ranking[i].second <<  ", #: " << P2_ranking[i].first << ", ";
  // }
  // std::cout << std::endl;
  return X;
}

bool Poly::constructPPLGeneratorsFromV(const Eigen::MatrixXd &R_input,
    PPL::Generator_System *gs) {
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

  bool has_a_point = false;
  for (int i = 0; i < nrows; ++i) {
    PPL::Linear_Expression e;
    for (unsigned j = dim; j > 0; j--) {
      e += R_int(i,j) * PPL::Variable(j-1);
    }
    if (R_int(i,0) == 1) {
      // a vertex
      has_a_point = true;
      gs->insert(point(e, PPL_MULTIPLIER_BEFORE_ROUNDING));
    } else if (R_int(i,0) == 2) {
      gs->insert(line(e));
    } else {
      gs->insert(ray(e));
    }
  }

  // Every non-empty generator system must have at least one point.
  if (nrows > 0 && !has_a_point) {
    gs->insert(PPL::point());
  }
  return true;
}

bool Poly::constructPPLConstraintsFromH(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
    const Eigen::MatrixXd &Ae, const Eigen::VectorXd &be, PPL::Constraint_System *cs) {
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
      cs->insert(e <= 0);
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
      cs->insert(e == 0);
    }
  }
  return true;
}

bool Poly::constructPPLPolyFromV(const Eigen::MatrixXd &R_input,
    PPL::C_Polyhedron *ph) {
  PPL::Generator_System gs;
  if (constructPPLGeneratorsFromV(R_input, &gs)) {
    ph->add_generators(gs);
    return true;
  } else {
    return false;
  }
}

bool Poly::constructPPLPolyFromH(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
    const Eigen::MatrixXd &Ae, const Eigen::VectorXd &be, PPL::C_Polyhedron *ph) {
  PPL::Constraint_System cs;
  if (constructPPLConstraintsFromH(A, b, Ae, be, &cs)) {
    // print_constraints(*ph, "*** ph constraints ***", std::cout);
    ph->add_constraints(cs);
    // print_constraints(*ph, "*** ph constraints ***", std::cout);
    return true;
  } else {
    return false;
  }
}

bool Poly::getVertexFromPPL(const PPL::C_Polyhedron &ph, Eigen::MatrixXd *R) {
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

  int dim = ph.space_dimension();
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

bool Poly::getFacetFromPPL(const PPL::C_Polyhedron &ph, Eigen::MatrixXd *A, Eigen::VectorXd *b) {
  // from Ax + b >= 0
  // to   Ax < b
  PPL::Constraint_System cs = ph.constraints();
  std::vector<double> A_vec;
  std::vector<double> b_vec;
  std::vector<double> A_row;
  std::vector<double> b_row;
  PPL::Constraint_System::const_iterator ic = cs.begin();
  PPL::Constraint_System::const_iterator cs_end = cs.end();

  int dim = ph.space_dimension();

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

bool Poly::vertexEnumeration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
    const Eigen::MatrixXd &Ae, const Eigen::VectorXd &be,  Eigen::MatrixXd *R) {
  int dim;
  if (A.cols() != 0) {
    dim = A.cols();
  } else if (Ae.cols() != 0) {
    dim = Ae.cols();
  } else {
    *R = Eigen::MatrixXd::Zero(0,0);
    return true; // input is empty, nothing to do
  }

  PPL::C_Polyhedron ph(dim, PPL::UNIVERSE); // start from the universe, then add constraints
  constructPPLPolyFromH(A, b, Ae, be, &ph);
  // print_constraints(ph, "*** ph constraints ***", std::cout);

  /**
   * Call vertex enumeration from cddlib
   */
  ph.minimized_generators();
  // print_constraints(ph, "*** ph constraints ***", std::cout);
  // print_generators(ph, "*** ph generators ***", std::cout);
  return getVertexFromPPL(ph, R);
}

bool Poly::vertexEnumeration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,  Eigen::MatrixXd *R) {
  Eigen::MatrixXd Ae(0, A.cols());
  Eigen::VectorXd be(0);
  return vertexEnumeration(A, b, Ae, be, R);
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

  int nrows = R_input.rows();
  int dim = R_input.cols() - 1;

  // std::cout << "Calling facetEnumeration." << std::endl;
  // std::cout << "R_input:\n" << R_input << std::endl;
  PPL::C_Polyhedron ph(dim, PPL::EMPTY); // start from empty, then add generators
  constructPPLPolyFromV(R_input, &ph);
  // print_generators(ph, "*** ph generators ***", std::cout);

  ph.minimized_constraints(); // V to H
  // print_generators(ph, "*** ph generators ***", std::cout);
  // print_constraints(ph, "*** ph constraints ***", std::cout);
  return getFacetFromPPL(ph, A, b);
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

  int dim = R1.cols() - 1;
  PPL::C_Polyhedron ph1(dim, PPL::EMPTY);
  PPL::C_Polyhedron ph2(dim, PPL::EMPTY);
  constructPPLPolyFromV(R1, &ph1);
  constructPPLPolyFromV(R2, &ph2);

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

  ph1.minimized_generators();

  /**
   * Read the results
   */

  return getVertexFromPPL(ph1, R);
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


// Solve:
//   Max:  r
//   S.t.  Vi' (x - Pi) >= r,   i = 1,...,n
//         A1 x <= b1
// Rewrite as
//   Max: r
//   s.t.  Ax + r <= b
//         A1 x <= b1
// Ax <= b describes the polytope. Rows of A must be normalized.
// xl, xu are additional bounds on x.
double Poly::inscribedSphere(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
    const Eigen::VectorXd &xl, const Eigen::VectorXd &xu,
    const Eigen::MatrixXd &A1, const Eigen::VectorXd &b1, Eigen::VectorXd *xc) {
  /**
   * construct the LP with variable [x, r]
   *  min: [0 -1] [x r]'
   *  s.t. [A 1] [x r]' <= b
   *       xl <= x <= xu
   *       r >= 0
   */
  // Constraints
  int dim = A.cols();
  int dim_ext = A.cols() + 1;
  Eigen::MatrixXd A_ext = Eigen::MatrixXd::Zero(A.rows() + A1.rows(), dim + 1);
  A_ext.block(0, 0, A.rows(), dim) = A;
  A_ext.block(0, dim, A.rows(), 1) = Eigen::VectorXd::Ones(A.rows());
  A_ext.block(A.rows(), 0, A1.rows(), dim) = A1;
  Eigen::VectorXd b_ext = Eigen::VectorXd::Zero(A.rows() + A1.rows());
  b_ext.head(A.rows()) = b;
  b_ext.tail(A1.rows()) = b1;
  // bounds
  Eigen::VectorXd xl_ext = Eigen::VectorXd(dim_ext) * nan("");
  Eigen::VectorXd xu_ext = Eigen::VectorXd(dim_ext) * nan("");
  if (xl.rows() > 0) {
    xl_ext.head(dim) = xl;
  }
  if (xu.rows() > 0) {
    xu_ext.head(dim) = xu;
  }
  xl_ext(dim) = 0; // r >= 0

  // costs
  Eigen::VectorXd C = Eigen::VectorXd::Zero(dim_ext);
  C(dim) = -1; // max r

  Eigen::MatrixXd Ae_ext(0, dim_ext);
  Eigen::VectorXd be_ext(0);
  Eigen::VectorXd xc_ext = Eigen::VectorXd::Zero(dim_ext);
  double optimal_cost;
  // std::cout << "xl:\n" << xl << std::endl;
  // std::cout << "xu:\n" << xu << std::endl;
  // std::cout << "A_ext:\n" << A_ext << std::endl;
  // std::cout << "b_ext:\n" << b_ext << std::endl;
  // std::cout << "xl_ext:\n" << xl_ext << std::endl;
  // std::cout << "xu_ext:\n" << xu_ext << std::endl;
  bool success = lp(C, A_ext, b_ext, Ae_ext, be_ext, xl_ext, xu_ext, &xc_ext, &optimal_cost);
  if (success) {
    *xc = xc_ext.head(dim);
    return -optimal_cost;
  } else {
    return -1;
  }
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
  // check singularity
  // todo: handle anti-podal situation
  Eigen::MatrixXd vectors_;
  Eigen::VectorXd v1 = vectors.middleRows(0,1).transpose();
  for (int i = 1; i < num; ++i) {
    Eigen::VectorXd v2 = vectors.middleRows(i,1).transpose();
    if ( v1.norm()*v2.norm() - v1.dot(v2) < 1e-7) {
      // v1, v2 are colinear
      // update v1
      v1 = v1+v2;
    } else {
      // colinearity ends
      vectors_ = Eigen::MatrixXd(num - i + 1, dim);
      vectors_.topRows(1) = v1.transpose();
      vectors_.bottomRows(num - i) = vectors.bottomRows(num - i);
      num = num - i + 1;
      break;
    }
  }
  if (vectors_.size() == 0) {
    // all colinear
    *results = v1.transpose();
    return true;
  }
  std::vector<double> msum;
  // estimation of storage. Exponential growing only happens KCONVHULL_ROUND.
  msum.reserve(2 * dim * pow(2, KCONVHULL_ROUND) * num);
  // 1. Initialize
  msum.assign(dim, 0);
  for (int i = 0; i < dim; ++i) msum.push_back(vectors_(0, i));

  std::vector<double> temp;
  for (int vid = 1; vid < num; ++vid) {
    /**
     * add the new vector vectors_(vid, :) to msum
     */
    int msum_num_now = msum.size()/dim;
    // std::cout << "msum_num_now: " << msum_num_now << std::endl;
    // self copy
    msum.insert(msum.end(), msum.begin(), msum.end());
    // add the new vector to the copy
    for (int i = 0; i < msum_num_now; ++i)
      for (int j = 0; j < dim; ++j)
        msum[(msum_num_now + i) * dim + j] += vectors_(vid, j);
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
    // double z = glp_get_obj_val(lp);
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

bool Poly::lp(const Eigen::VectorXd &C, const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
    const Eigen::MatrixXd &Ae, const Eigen::VectorXd &be,
    const Eigen::VectorXd &xl, const Eigen::VectorXd &xu, Eigen::VectorXd *xs, double *optimal_cost) {
  /* declare variables */
  if (xs->rows() <= 0) {
    std::cerr << "[lp] Error: xs is not initialized!!" << std::endl;
    exit(-1);
  }
  glp_prob *lp;
  glp_smcp parm;
  glp_init_smcp(&parm);
  parm.presolve = GLP_OFF;
  parm.msg_lev = GLP_MSG_ERR; // error and warning only
  int *ia, *ja;
  double *ar;
  int rows = A.rows();
  int cols = A.cols();
  int rows_e = Ae.rows();
  int cols_e = Ae.cols();
  assert(cols_e == cols);

  Eigen::VectorXd xu_expand = Eigen::VectorXd(cols) * nan("");
  Eigen::VectorXd xl_expand = Eigen::VectorXd(cols) * nan("");
  if (xu.rows() > 0) {
    xu_expand = xu;
  }
  if (xl.rows() > 0) {
    xl_expand = xl;
  }

  int size = rows * cols + rows_e * cols_e;
  ia = new int[size + 1000];
  ja = new int[size + 1000];
  ar = new double[size + 1000];

  /**
   * Create problem
   */
  lp = glp_create_prob();
  glp_set_obj_dir(lp, GLP_MIN); // minimization, not maximization

  /**
   * Fill problem
   */
  /* sign and right-hand-side of constraints */
  glp_add_rows(lp, rows + rows_e);
  for (int r = 1; r <= rows; ++r) {
    glp_set_row_bnds(lp, r, GLP_UP, 0.0, b(r-1)); // upper bound, <=
  }
  for (int r = 1; r <= rows_e; ++r) {
    glp_set_row_bnds(lp, rows+r, GLP_FX, be(r-1), be(r-1)); // equality =
  }

  /* cost function */
  glp_add_cols(lp, cols);
  /* variable bounds */
  for (int c = 1; c <= cols; ++c) {
    glp_set_obj_coef(lp, c, C(c-1)); // cost function
    if (isfinite(xl_expand(c-1))) {
      if (isfinite(xu_expand(c-1))) {
        glp_set_col_bnds(lp, c, GLP_DB, xl_expand(c-1), xu_expand(c-1)); // double bounded
      } else {
        glp_set_col_bnds(lp, c, GLP_LO, xl_expand(c-1), 0.0); // lower-bounded
      }
    } else {
      if (isfinite(xu_expand(c-1))) {
        glp_set_col_bnds(lp, c, GLP_UP, 0.0, xu_expand(c-1)); // upper-bounded
      } else {
        glp_set_col_bnds(lp, c, GLP_FR, 0.0, 0.0); // no boundary
      }
    }
  }
  /* fill in coefficient matrix */
  int id = 0;
  for (int r = 1; r <= rows; ++r) {
    for (int c = 1; c <= cols; ++c) {
      id = (r-1)*cols + c;
      ia[id] = r, ja[id] = c, ar[id] = A(r-1, c-1);
    }
  }
  for (int r = 1; r <= rows_e; ++r) {
    for (int c = 1; c <= cols_e; ++c) {
      id = (r + rows - 1)*cols + c;
      ia[id] = r + rows, ja[id] = c, ar[id] = Ae(r-1, c-1);
    }
  }
  glp_load_matrix(lp, id, ia, ja, ar);
  /**
   * solve problem
   */
  // glp_write_prob(lp, 0, "problem.txt");
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
    *optimal_cost = glp_get_obj_val(lp);
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

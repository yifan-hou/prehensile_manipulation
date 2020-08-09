#include <ppl.hh>

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include <Eigen/Dense>

namespace PPL = Parma_Polyhedra_Library;

using namespace PPL;


void
print_constraints(const Constraint_System& cs,
                  const std::string& intro, std::ostream& s) {
  if (!intro.empty())
    s << intro << "\n";
  Constraint_System::const_iterator i = cs.begin();
  Constraint_System::const_iterator cs_end = cs.end();
  bool printed_something = i != cs_end;
  while (i != cs_end) {
    using IO_Operators::operator<<;
    s << *i++;
    if (i != cs_end)
      s << ",\n";
  }
  s << (printed_something ? "." : "true.") << std::endl;
}

void
print_constraints(const Polyhedron& ph,
                  const std::string& intro, std::ostream& s) {
  print_constraints(ph.constraints(), intro, s);
}

void
print_generators(const Generator_System& gs,
                 const std::string& intro, std::ostream& s) {
  if (!intro.empty())
    s << intro << "\n";
  Generator_System::const_iterator i = gs.begin();
  Generator_System::const_iterator gs_end = gs.end();
  bool printed_something = i != gs_end;
  while (i != gs_end) {
    using IO_Operators::operator<<;
    s << *i++;
    if (i != gs_end)
      s << ",\n";
  }
  s << (printed_something ? "." : "false.") << std::endl;
}

void
print_generators(const Grid_Generator_System& gs,
                 const std::string& intro, std::ostream& s) {
  if (!intro.empty())
    s << intro << "\n";
  Grid_Generator_System::const_iterator i = gs.begin();
  Grid_Generator_System::const_iterator gs_end = gs.end();
  bool printed_something = i != gs_end;
  while (i != gs_end) {
    using IO_Operators::operator<<;
    s << *i++;
    if (i != gs_end)
      s << ",\n";
  }
  s << (printed_something ? "." : "false.") << std::endl;
}

void
print_generators(const Polyhedron& ph,
                 const std::string& intro, std::ostream& s) {
  print_generators(ph.generators(), intro, s);
}

void
print_generators(const Grid& gr,
                 const std::string& intro, std::ostream& s) {
  print_generators(gr.grid_generators(), intro, s);
}





int main(int argc, char *argv[]) {

// V to H
  Eigen::MatrixXi coefficients(4,4);
  coefficients << 1, 100000,  0, 0,
                  1, -100000, 0, 0,
                  1, 0,  100000, 1,
                  1, 0, -100000, 0;
  // coefficients << 0,      50000,        100000,  -45610,
  //          0,     -50000,        100000,  -51431,
  //          0,      47213,        94423,    48569,
  //          0,     -47213,        94423,    43067;
  // coefficients <<  0,   21350,   71180, -446067,
  //                 0,  -21350,   71180,    -503000,
  //                 0,   20160,   67210,     475000,
  //                 0,  -20160,   67210,    421200;


  int nrows = coefficients.rows();
  int dim = coefficients.cols() - 1;
  std::cout << "dim: " << dim << std::endl;
  PPL::Constraint_System cs;
  PPL::Generator_System gs;

  bool has_a_point = false;
  for (int i = 0; i < nrows; ++i) {
    PPL::Linear_Expression e;
    std::cout << "row " << i << " ";
    for (unsigned j = dim; j > 0; j--) {
      e += coefficients(i,j) * PPL::Variable(j-1);
      std::cout << j << ": " << coefficients(i,j) << ", ";
    }
    std::cout << std::endl;
    if (coefficients(i,0) == 1) {
      // a vertex
      has_a_point = true;
      gs.insert(point(e, 1));
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

  print_generators(ph, "*** ph generators ***", std::cout);
  print_constraints(ph, "*** ph constraints ***", std::cout);

  cs = ph.constraints();

  // convert constraints to Eigen Matrix format
  // from Ax + b >= 0
  // to   Ax < b
  std::cout << "Converting constraints to Eigen format:\n";
  std::vector<double> A_vec;
  std::vector<double> b_vec;
  std::vector<double> A_row;
  std::vector<double> b_row;
  Constraint_System::const_iterator ic = cs.begin();
  Constraint_System::const_iterator cs_end = cs.end();

  std::vector<mpz_class> g_coefficients;
  g_coefficients.resize(dim);
  while (ic != cs_end) {
    A_row.resize(dim);
    b_row.resize(1);

    mpz_class b = ic->inhomogeneous_term();
    mpz_class max = abs(b);
    for (dimension_type j = ic->space_dimension(); j-- > 0; ) {
      if (cmp(max, abs(ic->coefficient(Variable(j)))) < 0) {
        max = abs(ic->coefficient(Variable(j)));
      }
    }

    for (dimension_type j = ic->space_dimension(); j-- > 0; ) {
      A_row[j] = - mpq_class(ic->coefficient(Variable(j)), max).get_d();
    }
    b_row[0] = mpq_class(b, max).get_d();

    if (ic->is_equality()) {
      A_row.resize(2*dim);
      for (int ii = 0; ii < dim; ++ii)
        A_row[dim + ii] = -A_row[ii];
      b_row.push_back(-b_row[0]);
    }
    for (int ii = 0; ii < dim; ++ii) {
      std::cout << A_row[ii] << ", ";
    }
    std::cout << std::endl;

    for (int ii = 0; ii < A_row.size(); ++ii)
      A_vec.push_back(A_row[ii]);
    for (int ii = 0; ii < b_row.size(); ++ii)
      b_vec.push_back(b_row[ii]);
    ic++;
  }

  int num_generators = A_vec.size()/dim;
  Eigen::MatrixXd A_ = Eigen::MatrixXd::Map(A_vec.data(), dim, num_generators).transpose();
  Eigen::VectorXd b_ = Eigen::VectorXd::Map(b_vec.data(), b_vec.size());

  std::cout << "A_:\n" << A_ << std::endl;
  std::cout << "b_:\n" << b_.transpose() << std::endl;





  // H to V
  std::cout << "\n\n\n\nH to V\n\n";
  // Eigen::MatrixXi A(3, 3);
  // A << -10, 0, 199,
  //       0, -11, 0,
  //       0, 0, -120;
  //       // 1,  1, 1;
  // Eigen::VectorXi b(3,1);
  // b << 0, 0, 0;

  // nrows = A.rows();
  // dim = A.cols();
  // std::cout << "dim: " << dim << std::endl;

  // cs.clear();
  // for (int i = 0; i < nrows; ++i) {
  //   PPL::Linear_Expression e;
  //   std::cout << "row " << i << " ";
  //   for (unsigned j = dim; j > 0; j--) {
  //     e += A(i,j-1) * PPL::Variable(j-1);
  //     std::cout << j-1 << ": " << A(i,j-1) << ", ";
  //   }
  //   std::cout << std::endl;
  //   e -= b(i);
  //   cs.insert(e <= 0);
  // }

  // ph = PPL::C_Polyhedron(cs);
  // ph.minimized_generators(); // V to H

  // print_generators(ph, "*** generators: ***", std::cout);
  // print_constraints(ph, "*** constraints: ***", std::cout);

  gs = ph.generators();
  // convert generators to Eigen Matrix format
  std::cout << "Converting Generators to Eigen format:\n";
  std::vector<double> R_vec;
  std::vector<int> t_vec;
  std::vector<double> R_row;
  std::vector<int> t_row;
  Generator_System::const_iterator ig = gs.begin();
  Generator_System::const_iterator gs_end = gs.end();

  while (ig != gs_end) {
    R_row.resize(dim);
    t_row.resize(1);
    t_row[0] = 0;
    if (ig->is_point()) {
      std::cout << "Point.";
      t_row[0] = 1;
      mpz_class divisor = ig->divisor();
      for (dimension_type j = ig->space_dimension(); j-- > 0; ) {
        mpq_class q(ig->coefficient(Variable(j)), divisor);
        R_row[j] = q.get_d();
      }
    } else if (ig->is_ray() || ig->is_line()) {
      std::cout << "Ray or line.";
      // needs to normalize
      mpz_class max = abs(ig->coefficient(Variable(0)));
      for (dimension_type j = ig->space_dimension(); j-- > 0; ) {
        if (cmp(max, abs(ig->coefficient(Variable(j)))) < 0)
          max = abs(ig->coefficient(Variable(j)));
      }
      for (dimension_type j = ig->space_dimension(); j-- > 0; ) {
        R_row[j] = mpq_class(ig->coefficient(Variable(j)), max).get_d();
      }
    } else {
      std::cout << "A Closure point! Not implemented yet, should be the same as point. Make sure you know why there is a Closure point\n";
      exit(1);
    }

    if (ig->is_line()) {
      R_row.resize(2*dim);
      for (int ii = 0; ii < dim; ++ii)
        R_row[dim + ii] = -R_row[ii];
      t_row.assign(2, 0);
    }
    for (int ii = 0; ii < dim; ++ii) {
      std::cout << R_row[ii] << ", ";
    }
    std::cout << std::endl;

    for (int ii = 0; ii < R_row.size(); ++ii)
      R_vec.push_back(R_row[ii]);
    for (int ii = 0; ii < t_row.size(); ++ii)
      t_vec.push_back(t_row[ii]);
    ig++;
  }

  num_generators = R_vec.size()/dim;
  Eigen::MatrixXd R_ = Eigen::MatrixXd::Map(R_vec.data(), dim, num_generators).transpose();
  Eigen::VectorXi t_ = Eigen::VectorXi::Map(t_vec.data(), t_vec.size());

  std::cout << "R_:\n" << R_ << std::endl;
  std::cout << "t_:\n" << t_.transpose() << std::endl;
  return 0;
}
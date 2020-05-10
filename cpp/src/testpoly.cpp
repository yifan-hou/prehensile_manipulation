#include <iostream>
#include "polyhedron.h"

Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");


int main(int argc, char *argv[]) {
  // Eigen::MatrixXd A(4, 2);
  // A << -3, 0,
  //       0, -3,
  //       1,  0,
  //       0,  1;
  // Eigen::VectorXd b(4);
  // b << 7, 7, 1, 1;

  Eigen::MatrixXd A(1, 3);
  A << 0,  0, 1;
  Eigen::VectorXd b(1);
  b << 0;

  Eigen::MatrixXd R;
  // if (!Poly::vertexEnumeration(A, b, &R)) return -1;

  // std::cout << "Vertex Enumeration:\n";
  // std::cout << "A:\n" << A.format(CleanFmt) << std::endl;
  // std::cout << "b:\n" << b.format(CleanFmt) << std::endl;
  // std::cout << "R:\n" << R.format(CleanFmt) << std::endl;

  // if (!Poly::facetEnumeration(R, &A, &b)) return -1;

  // std::cout << "Facet Enumeration:\n";
  // std::cout << "R:\n" << R.format(CleanFmt) << std::endl;
  // std::cout << "A:\n" << A.format(CleanFmt) << std::endl;
  // std::cout << "b:\n" << b.format(CleanFmt) << std::endl;

  // if (!Poly::vertexEnumeration(A, b, &R)) return -1;

  // std::cout << "Vertex Enumeration:\n";
  // std::cout << "A:\n" << A.format(CleanFmt) << std::endl;
  // std::cout << "b:\n" << b.format(CleanFmt) << std::endl;
  // std::cout << "R:\n" << R.format(CleanFmt) << std::endl;

  Eigen::MatrixXd R1(2, 3);
  Eigen::MatrixXd R2(2, 3);
  R1 << 0, 1, 0,
        0, 0, 1;
  R2 << 0, 1, 1,
        0, -1, 1;
  if (!Poly::intersection(R1, R2, &R)) return -1;
  std::cout << "Intersection:\n";
  std::cout << "R1:\n" << R1.format(CleanFmt) << std::endl;
  std::cout << "R2:\n" << R2.format(CleanFmt) << std::endl;
  std::cout << "R:\n" << R.format(CleanFmt) << std::endl;


  R1 = Eigen::MatrixXd(4, 3);
  R2 = Eigen::MatrixXd(4, 3);
  R1 << 1, 0, 0,
        1, 0, 1,
        1, 1, 0,
        1, 1, 1,
  R2 << 1, 0.5, 0.5,
        1, -0.5, 0.5,
        1, 0.5, -0.5,
        1, -0.5, -0.5;
  if (!Poly::intersection(R1, R2, &R)) return -1;
  std::cout << "Intersection:\n";
  std::cout << "R1:\n" << R1.format(CleanFmt) << std::endl;
  std::cout << "R2:\n" << R2.format(CleanFmt) << std::endl;
  std::cout << "R:\n" << R.format(CleanFmt) << std::endl;

  return 0;
}
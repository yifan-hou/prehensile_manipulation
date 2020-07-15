#include <iostream>
#include <vector>
#include "polyhedron.h"

Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

#include <string>

int main(int argc, char *argv[]) {
  std::cout << "############################################\n";
  std::cout << "    Distance to cone\n";
  std::cout << "############################################\n";
  Eigen::VectorXd p(3);
  p << -1, 0, 2;
  Eigen::VectorXd plane_a(3);
  plane_a << 1, 0, 0;
  double plane_b = 2;
  Eigen::VectorXd p_proj = Poly::projectP2Hyperplane(p, plane_a, plane_b);
  std::cout << "p: " << p.transpose() << std::endl;
  std::cout << "a: " << plane_a.transpose() << std::endl;
  std::cout << "p_proj: " << p_proj.transpose() << std::endl;

  std::cout << std::endl;
  p << -1, -1, -1;
  Eigen::MatrixXd cone_A(3, 3), cone_R(3, 3);
  cone_A << -1, 0, 0,
             0, -1, 0,
             0,  0, -1;
  cone_R << 1, 0, 0,
            0, 1, 0,
            0, 0, 1;
  std::cout << "p: " << p.transpose() << std::endl;
  double angle = Poly::distRay2ConeFromOutside(p, cone_A, cone_R);
  std::cout << "angle: " << angle << std::endl;

  std::cout << "############################################\n";
  std::cout << "    Distance to polyhedron\n";
  std::cout << "############################################\n";
  Eigen::VectorXd p2(3);
  p2 << -1, -1, -1;
  Eigen::MatrixXd A2(3,3);
  A2 << 1, 0, 0,
        0, 1, 0,
        0, 0, 1;
  Eigen::VectorXd b2(3);
  b2 << 0, 0, 0;
  std::cout << "p: " << p2.transpose() << std::endl;
  double dist = Poly::distP2Polyhedron(p2, A2, b2, Eigen::VectorXd::Random(3));
  std::cout << "dist: " << dist << std::endl;

  
  return 1;

  std::cout << "############################################\n";
  std::cout << "    Basic vertex and facet enumeration\n";
  std::cout << "############################################\n";

  Eigen::MatrixXd A1;
  Eigen::VectorXd b1;
  Eigen::MatrixXd Rhigh(1, 3);
  Rhigh << 0, 1, 0;
  if (!Poly::facetEnumeration(Rhigh, &A1, &b1)) return -1;

  std::cout << "Facet Enumeration:\n";
  std::cout << "R:\n" << Rhigh.format(CleanFmt) << std::endl;
  std::cout << "A:\n" << A1.format(CleanFmt) << std::endl;
  std::cout << "b:\n" << b1.format(CleanFmt) << std::endl;

  Eigen::MatrixXd A(1, 3);
  A << 0,  0, 1;
  Eigen::VectorXd b(1);
  b << 0;

  Eigen::MatrixXd R;
  if (!Poly::vertexEnumeration(A, b, &R)) return -1;

  std::cout << "Vertex Enumeration:\n";
  std::cout << "A:\n" << A.format(CleanFmt) << std::endl;
  std::cout << "b:\n" << b.format(CleanFmt) << std::endl;
  std::cout << "R:\n" << R.format(CleanFmt) << std::endl;

  if (!Poly::facetEnumeration(R, &A, &b)) return -1;

  std::cout << "Facet Enumeration:\n";
  std::cout << "R:\n" << R.format(CleanFmt) << std::endl;
  std::cout << "A:\n" << A.format(CleanFmt) << std::endl;
  std::cout << "b:\n" << b.format(CleanFmt) << std::endl;

  if (!Poly::vertexEnumeration(A, b, &R)) return -1;

  std::cout << "Vertex Enumeration:\n";
  std::cout << "A:\n" << A.format(CleanFmt) << std::endl;
  std::cout << "b:\n" << b.format(CleanFmt) << std::endl;
  std::cout << "R:\n" << R.format(CleanFmt) << std::endl;

  std::cout << "############################################\n";
  std::cout << "    Redundant facet enumeration\n";
  std::cout << "############################################\n";
  Eigen::MatrixXd R_redundant(4, 3);
  R_redundant << 0, 1, 0,
                 0, 0, 1,
                 0, 1, 1,
                 0, 1, 2;

  if (!Poly::facetEnumeration(R_redundant, &A, &b)) return -1;

  std::cout << "Facet Enumeration:\n";
  std::cout << "R:\n" << R_redundant.format(CleanFmt) << std::endl;
  std::cout << "A:\n" << A.format(CleanFmt) << std::endl;
  std::cout << "b:\n" << b.format(CleanFmt) << std::endl;

  std::cout << "############################################\n";
  std::cout << "   limited vertex enumeration\n";
  std::cout << "############################################\n";
  Eigen::MatrixXd Ai(3, 2);
  Ai << 0, 1,
        1, 0,
       -1, -1;
  Eigen::VectorXd bi(3);
  bi << 0, 0, 0;
  if (!Poly::vertexEnumeration(Ai, bi, &R)) return -1;

  std::cout << "Vertex Enumeration:\n";
  std::cout << "A:\n" << Ai.format(CleanFmt) << std::endl;
  std::cout << "b:\n" << bi.format(CleanFmt) << std::endl;
  std::cout << "R:\n" << R.format(CleanFmt) << std::endl;

  Eigen::MatrixXd Aii(3, 3);
  Aii << 0, 0, 1,
       0, 1, 0,
       0, -1, -1;
  Eigen::VectorXd bii(3);
  bii << 0, 0, -1;
  if (!Poly::vertexEnumeration(Aii, bii, &R)) return -1;

  std::cout << "Vertex Enumeration:\n";
  std::cout << "A:\n" << Aii.format(CleanFmt) << std::endl;
  std::cout << "b:\n" << bii.format(CleanFmt) << std::endl;
  std::cout << "R:\n" << R.format(CleanFmt) << std::endl;


  std::cout << "############################################\n";
  std::cout << "    Double Description high complexity\n";
  std::cout << "############################################\n";

  Eigen::MatrixXd M = Eigen::MatrixXd::Random(32, 6) + Eigen::MatrixXd::Ones(32, 6);
  Eigen::MatrixXd R3(M.rows(), M.cols()+1);
  R3 << Eigen::VectorXd::Zero(M.rows()), M;
  if (!Poly::facetEnumeration(R3, &A, &b)) return -1;
  std::cout << "Facet Enumeration:\n";
  std::cout << "R:\n" << R3.format(CleanFmt) << std::endl;
  std::cout << "A:\n" << A.format(CleanFmt) << std::endl;
  std::cout << "b^T:\n" << b.transpose().format(CleanFmt) << std::endl;

  std::cout << "############################################\n";
  std::cout << "    intersection\n";
  std::cout << "############################################\n";

  Eigen::MatrixXd R1(2, 3);
  Eigen::MatrixXd R2(2, 3);
  // R1 << 0, 1, 0,
  //       0, 0, 1;
  // R2 << 0, 1, 1,
  //       0, -1, 1;
  R1 << 0, 1, 0,
        0, 0, 1;
  R2 << 0, -1, -1,
        0, -1, 0;
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

  std::cout << "############################################\n";
  std::cout << "    coneIntersection\n";
  std::cout << "############################################\n";
  Eigen::MatrixXd C1(2, 2);
  Eigen::MatrixXd C2(2, 2);
  C1 << 0, 1,
        1, 0;
  C2 << 1, 1,
        -1, 1;
  Eigen::MatrixXd C;
  if(!Poly::coneIntersection(C1, C2, &C)) return -1;
  std::cout << "Cone Intersection:\n";
  std::cout << "C1:\n" << C1.format(CleanFmt) << std::endl;
  std::cout << "C2:\n" << C2.format(CleanFmt) << std::endl;
  std::cout << "C:\n" << C.format(CleanFmt) << std::endl;

  Eigen::MatrixXd C3(2, 2);
  Eigen::MatrixXd C4(2, 2);
  C3 << 0, 1,
         1, 0;
  C4 << -1, 0,
        0, -1;
  if(!Poly::coneIntersection(C3, C4, &C)) return -1;
  std::cout << "Cone Intersection:\n";
  std::cout << "C1:\n" << C3.format(CleanFmt) << std::endl;
  std::cout << "C2:\n" << C4.format(CleanFmt) << std::endl;
  std::cout << "C:\n" << C.format(CleanFmt) << std::endl;

  Eigen::MatrixXd C5(1, 6);
  Eigen::MatrixXd C6(1, 6);
  C5 << 1, 0, 0, 0, 0, 0;
  C6 << 0, 0, 1, 0, 0, 0;
  if(!Poly::coneIntersection(C5, C6, &C)) return -1;
  std::cout << "Cone Intersection:\n";
  std::cout << "C1:\n" << C5.format(CleanFmt) << std::endl;
  std::cout << "C2:\n" << C6.format(CleanFmt) << std::endl;
  std::cout << "C:\n" << C.format(CleanFmt) << std::endl;

  std::cout << "############################################\n";
  std::cout << "    offsetPolytope\n";
  std::cout << "############################################\n";
  Eigen::MatrixXd P1(2, 6);
  P1 << 1, 1, 1, 1, 1, 1,
        2, 2, 2, 2, 2, 2;
  Eigen::VectorXd v1(6);
  v1 << 4, 4, 4, 4, 4, 4;
  std::cout << "P1: " << P1 << std::endl;
  std::cout << "v1: " << v1 << std::endl;
  Poly::offsetPolytope(&P1, v1);
  std::cout << "P1 + v1: " << P1 << std::endl;

  std::cout << "############################################\n";
  std::cout << "    convhull\n";
  std::cout << "############################################\n";
  Eigen::MatrixXd points1(6, 3);
  points1 << 0, 0, 0,
            1, 0, 0,
            2, 0, 0,
            0, 0, 1,
            0, 2, 0,
            0, 0, 2;
  std::cout << "points:\n" << points1 << std::endl;
  Eigen::MatrixXd hull = Poly::convhull(points1);
  std::cout << "hull: \n" << hull << std::endl;;

  std::cout << "############################################\n";
  std::cout << "    convhull degeneration case\n";
  std::cout << "############################################\n";
  Eigen::MatrixXd points2(6, 3);
  points2 << 0, 0, 0,
            1, 0, 0,
            2, 0, 0,
            0, 0, 0,
            0, 2, 0,
            0, 0, 0;
  std::cout << "points:\n" << points2 << std::endl;
  hull = Poly::convhull(points2);
  std::cout << "hull: \n" << hull << std::endl;;

  std::cout << "############################################\n";
  std::cout << "    convhull vector\n";
  std::cout << "############################################\n";
  std::vector<double> points2_vec;
  points2_vec.resize(points2.size());
  Eigen::MatrixXd::Map(&points2_vec[0], points2.cols(), points2.rows()) = points2.transpose();
  std::cout << "points:\n";
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 3; ++j)
      std::cout << points2_vec[i*3 + j] << " ";
    std::cout << std::endl;
  }
  std::vector<double> results2;
  Poly::convhull(points2_vec, 3, 6,  &results2);
  std::cout << "hull: \n";
  for (int i = 0; i < results2.size()/3; ++i) {
    for (int j = 0; j < 3; ++j)
      std::cout << results2[i*3 + j] << " ";
    std::cout << std::endl;
  }


  std::cout << "############################################\n";
  std::cout << "    minkowskiSumOfVectors\n";
  std::cout << "############################################\n";
  Eigen::MatrixXd min_vec1(2, 3);
  min_vec1 << 1, 0, 0,
              0, 1, 0;
  std::cout << "vectors:\n" << min_vec1 << std::endl;
  Eigen::MatrixXd results_min;
  Poly::minkowskiSumOfVectors(min_vec1,  &results_min);
  std::cout << "results:\n" << results_min << std::endl;

  std::cout << "############################################\n";
  std::cout << "    minkowskiSum of polytopes\n";
  std::cout << "############################################\n";

  Eigen::MatrixXd mink_poly1 = Eigen::MatrixXd::Random(9, 6);
  Eigen::MatrixXd mink_poly2 = Eigen::MatrixXd::Random(9, 6);
  Eigen::MatrixXd mink_poly3 = Eigen::MatrixXd::Random(9, 6);
  Eigen::MatrixXd mink_poly4 = Eigen::MatrixXd::Random(9, 6);
  Eigen::MatrixXd results_12, results_123, results_1234;

  std::cout << "Minkowski Sum of poly1 and poly2:\n";
  Poly::minkowskiSum(mink_poly1, mink_poly2,  &results_12);
  std::cout << "rows: " << results_12.rows() << std::endl;
  std::cout << "Minkowski Sum of poly12 and poly3:\n";
  Poly::minkowskiSum(results_12, mink_poly3,  &results_123);
  std::cout << "rows: " << results_123.rows() << std::endl;
  std::cout << "Minkowski Sum of poly123 and poly4:\n";
  Poly::minkowskiSum(results_123, mink_poly4,  &results_1234);
  std::cout << "rows: " << results_1234.rows() << std::endl;



  std::cout << "\nDone." << std::endl;
  return 0;
}
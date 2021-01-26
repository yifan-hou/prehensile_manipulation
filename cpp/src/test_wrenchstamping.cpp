#include <iostream>
#include <stdio.h>

#include "shared_grasping_jacobian.h"

using namespace Eigen;

int test2d();
int test3d();

int main(int argc, char *argv[]) {
  return test3d();
}

int test2d() {
  std::cout << "[main] testing 2d." << std::endl;
  /**
   * Problem description
   */
  double kFrictionE = 0.3;
  double kFrictionH = 0.7;
  double kObjWeight = 10;
  double kW = 0.0435; // object width
  double kH = 0.0435; // object height
  MatrixXd CP_H_e(2,2);
  MatrixXd CN_H_e(2,2);
  MatrixXd CP_H_h(2,2);
  MatrixXd CN_H_h(2,2);
  CP_H_e << kW/2, 0, // each row is a point
            -kW/2, 0;
  CN_H_e << 0, 1,
            0, 1;
  CP_H_h << kW/2, 0,
           -kW/2, 0;
  CN_H_h << 0, -1,
            0, -1;
  CP_H_e.transposeInPlace();
  CN_H_e.transposeInPlace();
  CP_H_h.transposeInPlace();
  CN_H_h.transposeInPlace();

  VectorXd CP_H_G(2), v_HG(2);
  CP_H_G << 0, -kH/2.0;
  v_HG << 0, -kH/2.0;

  /**
   * Geometrical processing
   */
  MatrixXd N_e, T_e, N_h, T_h, eCone_allFix, hCone_allFix;
  VectorXd F_G;
  sharedGraspingGeometryProcessing2d(kFrictionE, kFrictionH,
      CP_H_e, CN_H_e, CP_H_h, CN_H_h, N_e, T_e, N_h, T_h,
      eCone_allFix, hCone_allFix);

  std::cout << "results: \n";
  std::cout << "N_e:\n" << N_e << std::endl;
  std::cout << "T_e:\n" << T_e << std::endl;
  std::cout << "N_h:\n" << N_h << std::endl;
  std::cout << "T_h:\n" << T_h << std::endl;
  std::cout << "eCone_allFix:\n" << eCone_allFix << std::endl;
  std::cout << "hCone_allFix:\n" << hCone_allFix << std::endl;
  // std::cout << "F_G:\n" << F_G << std::endl;
  return 1;
  /**
   * Contact mode enumeration
   */

  /**
   * Wrench stamping
   */

}

int test3d() {
  std::cout << "[main] testing 3d." << std::endl;
  /**
   * Problem description
   */
  double kFrictionE = 0.25;
  double kFrictionH = 0.5;
  double kNumSlidingPlanes = 2;
  double kObjWeight = 10;
  double kW = 0.0435; // object width
  double kH = 0.0435; // object height
  MatrixXd CP_H_e(4,3);
  MatrixXd CN_H_e(4,3);
  MatrixXd CP_H_h(4,3);
  MatrixXd CN_H_h(4,3);
  CP_H_e << kW/2, kW/2, -kH, // each row is a point
            kW/2, -kW/2, -kH,
            -kW/2, kW/2, -kH,
            -kW/2, -kW/2, -kH;
  CN_H_e << 0, 0, 1,
            0, 0, 1,
            0, 0, 1,
            0, 0, 1;
  CP_H_h << kW/2, kW/2, 0,
            kW/2, -kW/2, 0,
            -kW/2, kW/2, 0,
            -kW/2, -kW/2, 0;
  CN_H_h << 0, 0, -1,
            0, 0, -1,
            0, 0, -1,
            0, 0, -1;
  CP_H_e.transposeInPlace();
  CN_H_e.transposeInPlace();
  CP_H_h.transposeInPlace();
  CN_H_h.transposeInPlace();

  VectorXd CP_H_G(3), v_HG(3);
  CP_H_G << 0, 0, kH/2.0;
  v_HG << 0, 0, -1;

  /**
   * Geometrical processing
   */
  MatrixXd N_e, T_e, N_h, T_h, eCone_allFix, hCone_allFix;
  VectorXd F_G;
  sharedGraspingGeometryProcessing3d(kFrictionE, kFrictionH, kNumSlidingPlanes,
      CP_H_e, CN_H_e, CP_H_h, CN_H_h, N_e, T_e, N_h, T_h,
      eCone_allFix, hCone_allFix);

  std::cout << "results: \n";
  std::cout << "N_e:\n" << N_e << std::endl;
  std::cout << "T_e:\n" << T_e << std::endl;
  std::cout << "N_h:\n" << N_h << std::endl;
  std::cout << "T_h:\n" << T_h << std::endl;
  std::cout << "eCone_allFix:\n" << eCone_allFix << std::endl;
  std::cout << "hCone_allFix:\n" << hCone_allFix << std::endl;
  // std::cout << "F_G:\n" << F_G << std::endl;
  /**
   * Contact mode enumeration
   */
  
  return 1;
  /**
   * Wrench stamping
   */

}
#include <iostream>
#include <fstream>
#include <stdio.h>

#include "shared_grasping_jacobian.h" // 2d jacobian computation doesn't have wrapper yet
#include "wrench_space_analysis.h"

using namespace Eigen;

#define PI 3.1415926

int test3d();

int main(int argc, char *argv[]) {
  return test3d();
}

int test3d() {
  std::cout << "[main] testing 3d." << std::endl;

  /**
   * 1. call wrench stamping on the nominal problem
   */

  /**
   * Problem description
   */
  double kFrictionE = 0.8;
  double kFrictionH = 0.8;
  double kNumSlidingPlanes = 2;
  double kContactForce = 15;
  double kObjWeight = 10;
  double kL1 = 0.0435; // object width
  double kL2 = 0.0435; // object height
  double kF1 = 0.02; // finger 1 distance to corner
  double kF2 = 0.02; // finger 2 distance to corner
  double kCharacteristicLength = 1;
  MatrixXd CP_H_e(2, 3);
  MatrixXd CN_H_e(2, 3);
  MatrixXd CP_H_h(2, 3);
  MatrixXd CN_H_h(2, 3);
  int kNumContactsE = 2;
  int kNumContactsH = 2;
  double theta = 45*PI/180;
  CP_H_h << -kF1*sin(theta), 0, -kF1*cos(theta),
            kF2*cos(theta), 0, -kF2*sin(theta);
  CP_H_e << -kL1*sin(theta)+kL2*cos(theta), 0.0103, -kL1*cos(theta)-kL2*sin(theta),
            -kL1*sin(theta)+kL2*cos(theta), -0.0103, -kL1*cos(theta)-kL2*sin(theta);
  CN_H_h << cos(theta), 0, -sin(theta),
            -sin(theta), 0, -cos(theta);
  CN_H_e << 0, 0, 1,
            0, 0, 1;
  CP_H_e.transposeInPlace();
  CN_H_e.transposeInPlace();
  CP_H_h.transposeInPlace();
  CN_H_h.transposeInPlace();

  VectorXd CP_H_G(3), v_HG(3);
  CP_H_G << 0, 0, kL1/2.0;
  v_HG << 0, 0, -1;

  /**
   * Goal
   */
  MatrixXd G(1, 12);
  G << 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
  VectorXd b_G(1);
  b_G << 0.1;

  MatrixXi e_cs_modes_goal(1, kNumContactsE);
  e_cs_modes_goal << 0, 0;
  std::vector<MatrixXi> e_ss_modes_goal(1);
  e_ss_modes_goal[0] = MatrixXi::Zero(1, kNumContactsE*kNumSlidingPlanes);

  MatrixXi h_cs_modes_goal(1, kNumContactsH);
  h_cs_modes_goal << 0, 0;
  std::vector<MatrixXi> h_ss_modes_goal(1);
  h_ss_modes_goal[0] = MatrixXi::Zero(1, kNumContactsH*kNumSlidingPlanes);

  WrenchSpaceAnalysis wsa;

  wsa.updateConstants(kFrictionE, kFrictionH, kNumSlidingPlanes, kContactForce,
      kObjWeight, kCharacteristicLength);

  wsa.updateContactGeometry(kNumContactsE, kNumContactsH,
      CP_H_e, CN_H_e, CP_H_h, CN_H_h, CP_H_G, v_HG);

  wsa.computeContactModes();

  HFVC action;
  auto [g_margin, c_margin] = wsa.wrenchStampingWrapper(G, b_G,
    e_cs_modes_goal, e_ss_modes_goal, h_cs_modes_goal, h_ss_modes_goal, action);

  std::cout << "Geometrical stability margin: " << g_margin << std::endl;
  std::cout << "Control stability margin: " << c_margin << std::endl;

  auto [g_margin1, c_margin1] = wsa.wrenchStampingEvaluationWrapper(G, b_G,
    e_cs_modes_goal, e_ss_modes_goal, h_cs_modes_goal, h_ss_modes_goal, action);
  std::cout << "Geometrical stability margin: " << g_margin1 << std::endl;
  std::cout << "Control stability margin: " << c_margin1 << std::endl;

  // return 0;
  std::vector<double> kFrictionE_range = {0.1, 1.2};
  std::vector<double> kFrictionH_range = {0.1, 1.2};
  std::vector<double> kL1_range = {0.02, 0.1};
  std::vector<double> kL2_range = {0.02, 0.1};

  kNumSlidingPlanes = 4;
  e_ss_modes_goal[0] = MatrixXi::Zero(1, kNumContactsE*kNumSlidingPlanes);
  h_ss_modes_goal[0] = MatrixXi::Zero(1, kNumContactsH*kNumSlidingPlanes);

  // std::cout << "CP_H_e:\n" << CP_H_e << std::endl;
  // std::cout << "CN_H_e:\n" << CN_H_e << std::endl;
  // std::cout << "CP_H_h:\n" << CP_H_h << std::endl;
  // std::cout << "CN_H_h:\n" << CN_H_h << std::endl;
  // std::cout << "CP_H_G:\n" << CP_H_G << std::endl;
  // std::cout << "v_HG:\n" << v_HG << std::endl;
  /**
   * Evaluation 1: friction variations
   */
  std::ofstream fp_test1("/workspace/Git/prehensile_manipulation/data/test1_friction.txt");
  if (!fp_test1.is_open()) {
    std::cout << "Cannot open file.\n";
    return false;
  }
  int N = 10;
  fp_test1 << N << " " << kFrictionE << " " << kFrictionH << std::endl;
  fp_test1 << kFrictionE_range[0] << " " << kFrictionE_range[1] << std::endl;
  fp_test1 << kFrictionH_range[0] << " " << kFrictionH_range[1] << std::endl;
  for (int r1 = 0; r1 < N; ++r1) {
    double ratio1 = double(r1)/double(N);
    double frictionE = kFrictionE_range[0]*ratio1 + kFrictionE_range[1]*(1-ratio1);
    for (int r2 = 0; r2 < N; ++r2) {
      double ratio2 = double(r2)/double(N);
      double frictionH = kFrictionH_range[0]*ratio2 + kFrictionH_range[1]*(1-ratio2);

      // // debug
      // frictionE = kFrictionE;
      // frictionH = kFrictionH;
      wsa.updateConstants(frictionE, frictionH, kNumSlidingPlanes, kContactForce,
          kObjWeight, kCharacteristicLength);
      wsa.updateContactGeometry(kNumContactsE, kNumContactsH,
          CP_H_e, CN_H_e, CP_H_h, CN_H_h, CP_H_G, v_HG);
      wsa.computeContactModes();
      auto [g_margin, c_margin] = wsa.wrenchStampingEvaluationWrapper(G, b_G,
        e_cs_modes_goal, e_ss_modes_goal, h_cs_modes_goal, h_ss_modes_goal, action);

      if ((g_margin > 0) && (c_margin > 0)) {
        fp_test1 << "1 ";
      } else {
        fp_test1 << "0 ";
      }
      std::cout << "testing friction " << r1 << ", " << r2 << ", " << frictionE << ", " << frictionH << ", (" << g_margin << ", " << c_margin << ")" << std::endl;
    }
    fp_test1 << "\n";
  }
  fp_test1.close();


  /**
   * Evaluation 2: geometrical variations
   */
  std::ofstream fp_test2("/workspace/Git/prehensile_manipulation/data/test2_geometrical.txt");
  if (!fp_test2.is_open()) {
    std::cout << "Cannot open file.\n";
    return false;
  }
  N = 10;
  fp_test2 << N << " " << kL1 << " " << kL2  << std::endl;
  fp_test2 << kL1_range[0] << " " << kL1_range[1] << std::endl;
  fp_test2 << kL2_range[0] << " " << kL2_range[1] << std::endl;
  for (int r1 = 0; r1 < N; ++r1) {
    double ratio1 = double(r1)/double(N);
    double L1 = kL1_range[0]*ratio1 + kL1_range[1]*(1-ratio1);
    for (int r2 = 0; r2 < N; ++r2) {
      double ratio2 = double(r2)/double(N);
      double L2 = kL2_range[0]*ratio2 + kL2_range[1]*(1-ratio2);
      // kFrictionE = 0.25;
      // kFrictionH = 0.5;
      // kL1 = 0.0435; // object width
      // L2 = 0.0435; // object height
      // double theta = 45*PI/180;
      CP_H_h(0, 0) = -kF1*sin(theta);
      CP_H_h(1, 0) = 0;
      CP_H_h(2, 0) = -kF1*cos(theta);
      CP_H_h(0, 1) = kF2*cos(theta);
      CP_H_h(1, 1) = 0;
      CP_H_h(2, 1) = -kF2*sin(theta);

      CP_H_e(0, 0) = -L1*sin(theta)+L2*cos(theta);
      CP_H_e(1, 0) = 0.0103;
      CP_H_e(2, 0) = -L1*cos(theta)-L2*sin(theta);
      CP_H_e(0, 1) = -L1*sin(theta)+L2*cos(theta);
      CP_H_e(1, 1) = -0.0103;
      CP_H_e(2, 1) = -L1*cos(theta)-L2*sin(theta);

      CN_H_h(0, 0) = cos(theta);
      CN_H_h(1, 0) = 0;
      CN_H_h(2, 0) = -sin(theta);
      CN_H_h(0, 1) = -sin(theta);
      CN_H_h(1, 1) = 0;
      CN_H_h(2, 1) = -cos(theta);
      CP_H_G(2) = L1/2.0;

      // std::cout << "CP_H_e:\n" << CP_H_e << std::endl;
      // std::cout << "CN_H_e:\n" << CN_H_e << std::endl;
      // std::cout << "CP_H_h:\n" << CP_H_h << std::endl;
      // std::cout << "CN_H_h:\n" << CN_H_h << std::endl;
      // std::cout << "CP_H_G:\n" << CP_H_G << std::endl;
      // std::cout << "v_HG:\n" << v_HG << std::endl;

      wsa.updateConstants(kFrictionE, kFrictionH, kNumSlidingPlanes, kContactForce,
          kObjWeight, kCharacteristicLength);
      wsa.updateContactGeometry(kNumContactsE, kNumContactsH,
          CP_H_e, CN_H_e, CP_H_h, CN_H_h, CP_H_G, v_HG);
      wsa.computeContactModes();
      auto [g_margin, c_margin] = wsa.wrenchStampingEvaluationWrapper(G, b_G,
        e_cs_modes_goal, e_ss_modes_goal, h_cs_modes_goal, h_ss_modes_goal, action);

      if ((g_margin > 0) && (c_margin > 0)) {
        fp_test2 << "1 ";
      } else {
        fp_test2 << "0 ";
      }
      std::cout << "testing geo " << r1 << ", " << r2;
      std::cout << ", (" << g_margin << ", " << c_margin << ")" << std::endl;
      fp_test2 << "\n";
    }
  }
  fp_test2.close();
  return 0;
}
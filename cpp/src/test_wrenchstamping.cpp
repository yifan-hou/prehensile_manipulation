#include <iostream>
#include <stdio.h>

#include <modus/enumerate_modes.hpp>

#include "shared_grasping_jacobian.h"
#include "wrench_space_analysis.h"

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
  double kContactForce = 15;
  double kObjWeight = 10;
  double kW = 0.0435; // object width
  double kH = 0.0435; // object height
  double kCharacteristicLength = 0.15;
  MatrixXd CP_H_e(4,3);
  MatrixXd CN_H_e(4,3);
  MatrixXd CP_H_h(4,3);
  MatrixXd CN_H_h(4,3);
  int kNumContactsE = 4;
  int kNumContactsH = 4;
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
   * Goal
   */
  MatrixXd G(1, 12);
  G << 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  VectorXd b_G(1);
  b_G << 0.1;

  MatrixXi e_cs_modes_goal(1, kNumContactsE);
  e_cs_modes_goal << 0, 0, 1, 1;
  std::vector<MatrixXi> e_ss_modes_goal(1);
  e_ss_modes_goal[0] = MatrixXi::Zero(1, kNumContactsE*kNumSlidingPlanes);

  MatrixXi h_cs_modes_goal(1, kNumContactsH);
  h_cs_modes_goal << 0, 0, 0, 0;
  std::vector<MatrixXi> h_ss_modes_goal(1);
  h_ss_modes_goal[0] = MatrixXi::Zero(1, kNumContactsH*kNumSlidingPlanes);


  WrenchSpaceAnalysis wsa;

  wsa.updateConstants(kFrictionE, kFrictionH, kNumSlidingPlanes, kContactForce,
      kObjWeight, kCharacteristicLength);

  wsa.updateContactGeometry(kNumContactsE, kNumContactsH,
      CP_H_e, CN_H_e, CP_H_h, CN_H_h, CP_H_G, v_HG);

  wsa.updateContactModes();

  auto [g_margin, c_margin] = wsa.wrenchStampingWrapper(G, b_G,
    e_cs_modes_goal, e_ss_modes_goal, h_cs_modes_goal, h_ss_modes_goal);

  // /**
  //  * Geometrical processing
  //  */
  // MatrixXd N_e, T_e, N_h, T_h, eCone_allFix, hCone_allFix;
  // sharedGraspingGeometryProcessing3d(kFrictionE, kFrictionH, kNumSlidingPlanes,
  //     CP_H_e, CN_H_e, CP_H_h, CN_H_h, N_e, T_e, N_h, T_h,
  //     eCone_allFix, hCone_allFix);

  // VectorXd F_G = getGravityVector3d(kObjWeight, CP_H_G, v_HG);

  // std::cout << "results: \n";
  // std::cout << "N_e:\n" << N_e << std::endl;
  // std::cout << "T_e:\n" << T_e << std::endl;
  // std::cout << "N_h:\n" << N_h << std::endl;
  // std::cout << "T_h:\n" << T_h << std::endl;
  // std::cout << "eCone_allFix:\n" << eCone_allFix << std::endl;
  // std::cout << "hCone_allFix:\n" << hCone_allFix << std::endl;
  // std::cout << "F_G:\n" << F_G << std::endl;

  // /**
  //  * Contact mode enumeration
  //  */
  // VectorXd b_e = VectorXd::Zero(N_e.rows());
  // IncidenceGraph* cs_graph = enumerate_cs_modes(N_e, b_e, 1e-8);
  // std::vector<std::string> svs = cs_graph->get_proper_sign_vectors();
  // MatrixXi e_cs_modes = MatrixXi::Zero(svs.size(), kNumContactsE);
  // int ss_modes_count = 0, cs_modes_count = 0;
  // std::vector<MatrixXi> e_ss_modes(svs.size());
  // for (int k = cs_graph->rank() - 1; k >= 0; k--) {
  //   for (Node* u : cs_graph->rank(k)) {
  //     // get cs modes
  //     ModeEnumerationOptions opt;
  //     opt.interior_point = u->interior_point;
  //     std::string sv = get_sign_vector(N_e * u->interior_point - b_e, 1e-8);
  //     assert(kNumContactsE == sv.length());
  //     std::cout << "cs sign " << sv << std::endl;

  //     // decode cs modes
  //     for (int j = 0; j < kNumContactsE; ++j) {
  //       if (sv[j] == '0')
  //         e_cs_modes(cs_modes_count, j) = 0;
  //       else
  //         e_cs_modes(cs_modes_count, j) = 1;
  //     }

  //     // get ss modes
  //     IncidenceGraph* ss_graph =
  //         enumerate_ss_modes(N_e, b_e, T_e, u->sign_vector, 1e-8, &opt);

  //     auto ss_modes = ss_graph->get_sign_vectors(u->sign_vector);
  //     std::cout << "ss count: " << ss_modes.size() << std::endl;
  //     MatrixXi e_ss_mode = MatrixXi::Zero(ss_modes.size(), kNumContactsE*kNumSlidingPlanes);
  //     for (int ss_id = 0; ss_id < ss_modes.size(); ss_id++) {
  //       std::string ss = ss_modes[ss_id];
  //       std::cout << ss << std::endl;
  //       int sliding_count = 0;
  //       for (int j = 0; j < kNumContactsE; ++j) {
  //         if (sv[j] == '0') {
  //           for (int jj = 0; jj < kNumSlidingPlanes; ++jj) {
  //             // std::cout << "ss[" << kNumContactsE+sliding_count*kNumSlidingPlanes + jj << "] = " << ss[kNumContactsE+sliding_count*kNumSlidingPlanes + jj] << std::endl;
  //             if (ss[kNumContactsE+sliding_count*kNumSlidingPlanes + jj] == '+')
  //               e_ss_mode(ss_id, j*kNumSlidingPlanes + jj) = 1;
  //             else if (ss[kNumContactsE+sliding_count*kNumSlidingPlanes + jj] == '-')
  //               e_ss_mode(ss_id, j*kNumSlidingPlanes + jj) = -1;
  //             else
  //               e_ss_mode(ss_id, j*kNumSlidingPlanes + jj) = 0;
  //           }
  //           sliding_count++;
  //         }
  //       }
  //     }
  //     // std::cout << "e_ss_mode: \n" << e_ss_mode << std::endl;
  //     std::cout << std::endl;
  //     e_ss_modes[cs_modes_count] = e_ss_mode;

  //     ss_modes_count += ss_modes.size();
  //     delete ss_graph;

  //     if (u->sign_vector == "00-------0-0") {
  //         int*a = nullptr;
  //         std::cout << *a << std::endl;
  //     }
  //     cs_modes_count ++;
  //   }
  // }
  // delete cs_graph;

  // MatrixXi h_cs_modes = MatrixXi::Zero(1, kNumContactsH);
  // std::vector<MatrixXi> h_ss_modes(1);
  // h_ss_modes[0] = MatrixXi::Zero(1, kNumContactsH*kNumSlidingPlanes);

  
  // /**
  //  * Wrench stamping
  //  */
  // MatrixXd J_e(N_e.rows() + T_e.rows(), N_e.cols());
  // J_e << N_e, T_e;
  // MatrixXd J_h(N_h.rows() + T_h.rows(), N_h.cols());
  // J_h << N_h, T_h;

  // WrenchSpaceAnalysis wsa;
  // HFVC action;
  // auto [g_margin, c_margin] = wsa.wrenchStamping(
  //     J_e, J_h, eCone_allFix, hCone_allFix, F_G, kContactForce,
  //     kFrictionE, kFrictionH, kCharacteristicLength, kNumSlidingPlanes,
  //     e_cs_modes, e_ss_modes, h_cs_modes, h_ss_modes,
  //     G, b_G, e_cs_modes_goal, e_ss_modes_goal, h_cs_modes_goal, h_ss_modes_goal,
  //     &action);
  std::cout << "Geometrical stability margin: " << g_margin << std::endl;
  std::cout << "Control stability margin: " << c_margin << std::endl;
  return 0;
}
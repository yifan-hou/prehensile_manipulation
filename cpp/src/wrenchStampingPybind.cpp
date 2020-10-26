#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <iostream>
#include <Eigen/Dense>
#include <memory>
#include "wrench_space_analysis.h"

namespace py = pybind11;

using namespace Eigen;

std::shared_ptr<WrenchSpaceAnalysis> wsa_ptr;

void wrenchStamping_wrapper(Eigen::MatrixXd Jac_e, Eigen::MatrixXd Jac_h,
    Eigen::MatrixXd eCone_allFix, Eigen::MatrixXd hCone_allFix,
    const Eigen::VectorXd &F_G,
    const double kContactForce, const double kFrictionE, const double kFrictionH,
    const double kCharacteristicLength, const int kNumSlidingPlanes,
    const Eigen::MatrixXi &e_cs_modes, const std::vector<Eigen::MatrixXi> &e_ss_modes,
    const Eigen::MatrixXi &h_cs_modes, const std::vector<Eigen::MatrixXi> &h_ss_modes,
    Eigen::MatrixXd G, const Eigen::VectorXd &b_G,
    const Eigen::MatrixXi &e_cs_modes_goal, const std::vector<Eigen::MatrixXi> &e_ss_modes_goal,
    const Eigen::MatrixXi &h_cs_modes_goal, const std::vector<Eigen::MatrixXi> &h_ss_modes_goal) {
  HFVC action;
  wsa_ptr->wrenchStamping(Jac_e, Jac_h, eCone_allFix, hCone_allFix,
      F_G, kContactForce, kFrictionE, kFrictionH, kCharacteristicLength, kNumSlidingPlanes,
      e_cs_modes, e_ss_modes, h_cs_modes, h_ss_modes, G, b_G,
      e_cs_modes_goal, e_ss_modes_goal, h_cs_modes_goal, h_ss_modes_goal, &action);
  std::cout << "[wrapper] Finished!\n";
}

Eigen::VectorXd wrenchStamping_2d_wrapper(Eigen::MatrixXd Jac_e, Eigen::MatrixXd Jac_h,
    Eigen::MatrixXd eCone_allFix, Eigen::MatrixXd hCone_allFix,
    const Eigen::VectorXd &F_G,
    const double kContactForce, const double kFrictionE, const double kFrictionH,
    const double kCharacteristicLength,
    Eigen::MatrixXd G, const Eigen::VectorXd &b_G,
    const Eigen::MatrixXi &e_modes, const Eigen::MatrixXi &h_modes,
    const Eigen::VectorXi &e_modes_goal, const Eigen::VectorXi &h_modes_goal) {
  HFVC action;
  double margin = wsa_ptr->wrenchStamping_2d(Jac_e, Jac_h, eCone_allFix, hCone_allFix, F_G,
      kContactForce, kFrictionE, kFrictionH, kCharacteristicLength,
      G, b_G, e_modes, h_modes, e_modes_goal, h_modes_goal, &action);
  // [margin, n_af, n_av, R_a (9), eta_af, w_av]
  Eigen::VectorXd results = Eigen::VectorXd::Zero(15);
  results(0) = margin;
  if (margin > 0) {
    results(1) = action.n_af;
    results(2) = action.n_av;
    results(3) = action.R_a(0,0);
    results(4) = action.R_a(0,1);
    results(5) = action.R_a(0,2);
    results(6) = action.R_a(1,0);
    results(7) = action.R_a(1,1);
    results(8) = action.R_a(1,2);
    results(9) = action.R_a(2,0);
    results(10) = action.R_a(2,1);
    results(11) = action.R_a(2,2);
    for (int i = 0; i < action.n_af; ++i) {
      results(12+i) = action.eta_af(i);
    }
    for (int i = 0; i < action.n_av; ++i) {
      results(12+action.n_af+i) = action.w_av(i);
    }
  }
  // std::cout << "[wrapper] Finished! results: \n" << results.transpose()
  //     << std::endl;
  return results;
}

PYBIND11_MODULE(wrenchStampingLib, m) {
  std::cout << "Initializing wrenchStampingLib" << std::endl;
  m.doc() = "pybind11 wrenchStampingLib plugin"; // optional module docstring

  std::cout << "Initializing wrenchSpaceAnalysis" << std::endl;
  wsa_ptr = std::shared_ptr<WrenchSpaceAnalysis>(new WrenchSpaceAnalysis());
  std::cout << "wrenchSpaceAnalysis initialized." << std::endl;
  m.def("wrenchSpaceAnalysis", &wrenchStamping_wrapper, "Test mode cleaning");
  m.def("wrenchSpaceAnalysis_2d", &wrenchStamping_2d_wrapper, "Test mode cleaning");
  std::cout << "wrenchStampingLib initialized" << std::endl;
}


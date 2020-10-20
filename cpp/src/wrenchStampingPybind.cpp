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
  wsa_ptr->wrenchStamping(Jac_e, Jac_h, eCone_allFix, hCone_allFix,
      F_G, kContactForce, kFrictionE, kFrictionH, kCharacteristicLength, kNumSlidingPlanes,
      e_cs_modes, e_ss_modes, h_cs_modes, h_ss_modes, G, b_G,
      e_cs_modes_goal, e_ss_modes_goal, h_cs_modes_goal, h_ss_modes_goal);
  std::cout << "[wrapper] Finished!\n";
}

double wrenchStamping_2d_wrapper(Eigen::MatrixXd Jac_e, Eigen::MatrixXd Jac_h,
    Eigen::MatrixXd eCone_allFix, Eigen::MatrixXd hCone_allFix,
    const Eigen::VectorXd &F_G,
    const double kContactForce, const double kFrictionE, const double kFrictionH,
    const double kCharacteristicLength,
    Eigen::MatrixXd G, const Eigen::VectorXd &b_G,
    const Eigen::MatrixXi &e_modes, const Eigen::MatrixXi &h_modes,
    const Eigen::VectorXi &e_modes_goal, const Eigen::VectorXi &h_modes_goal) {
  double margin = wsa_ptr->wrenchStamping_2d(Jac_e, Jac_h, eCone_allFix, hCone_allFix, F_G,
      kContactForce, kFrictionE, kFrictionH, kCharacteristicLength,
      G, b_G, e_modes, h_modes, e_modes_goal, h_modes_goal);
  std::cout << "[wrapper] Finished!\n";
  return margin;
}

PYBIND11_MODULE(wrenchStampingLib, m) {
    m.doc() = "pybind11 wrenchStampingLib plugin"; // optional module docstring

    wsa_ptr = std::shared_ptr<WrenchSpaceAnalysis>(new WrenchSpaceAnalysis());
    m.def("wrenchSpaceAnalysis", &wrenchStamping_wrapper, "Test mode cleaning");
    m.def("wrenchSpaceAnalysis_2d", &wrenchStamping_2d_wrapper, "Test mode cleaning");
}


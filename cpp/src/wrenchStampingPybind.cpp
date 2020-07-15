#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <iostream>
#include <Eigen/Dense>
#include "wrench_space_analysis.h"

namespace py = pybind11;

using namespace Eigen;

std::vector<MatrixXi> modeCleaning_test(const MatrixXi &cs_modes, const std::vector<MatrixXi> &ss_modes, int kNumSlidingPlanes) {

  int num_cs_modes = cs_modes.rows();
  int num_contacts = cs_modes.cols();
  int total_sliding_directions = ss_modes[0].cols();
  // assert(sss_modes->size() == 0);
  /**
   * First, we need to build a cs mode representation with {f = sticking, 0 = sliding, 1 = separation}
   */
  std::cout << "[debug] cs_modes: " << cs_modes.rows() << " rows, " << cs_modes.cols() << " cols." << std::endl;
  std::cout << "[debug] ss_modes[0]: " << ss_modes[0].rows() << " rows, " << ss_modes[0].cols() << " cols." << std::endl;
  std::vector<std::string> sss_modes_; // sticking, sliding, separation
  std::vector<std::vector<int>> s_modes_;

  std::string sss_mode;
  sss_mode.resize(num_contacts);
  for (int cs = 0; cs < num_cs_modes; ++cs) {
    Eigen::VectorXi cs_mode = cs_modes.middleRows(cs, 1).transpose();
    // std::cout << "[debug] cs_mode cs: " << cs_mode << std::endl;
    int num_ss_modes = ss_modes[cs].rows();
    // std::cout << "[debug] num_ss_modes: " << num_ss_modes << std::endl;
    for (int row = 0; row < num_ss_modes; ++row) {
      // process one contact mode
      bool redundant = false;
      for (int i = 0; i < num_contacts; ++i) {
        if (cs_mode(i) == 0) {
          // std::cout << "block: " << ss_modes[cs].block(row,i*kNumSlidingPlanes, 1, kNumSlidingPlanes) << std::endl;
          int ss_mode_contact_i_sum = ss_modes[cs].block(row,i*kNumSlidingPlanes, 1, kNumSlidingPlanes).cwiseAbs().sum();
          if (ss_mode_contact_i_sum == 0) {
            sss_mode[i] = 'f';
            // std::cout << "assign value to sss_mode[" << i << "]. sss_mode = " << sss_mode[i] << std::endl;
          }
          else if (ss_mode_contact_i_sum == kNumSlidingPlanes) {
            sss_mode[i] = '0';
            // std::cout << "assign value to sss_mode[" << i << "]. sss_mode = " << sss_mode[i] << std::endl;
          } else {
            // std::cout << "redundant\n";
            redundant = true;
            break;
          }
        } else {
          sss_mode[i] = '1';
          // std::cout << "assign value to sss_mode[" << i << "]. sss_mode = " << sss_mode[i] << std::endl;
        }
      }
      if (!redundant) {
        // std::cout << "Not redundant! sssmode: " << sss_mode << std::endl;
        // std::cout << "sssmode[0]: " << sss_mode[0] << " [1]: " << sss_mode[1] << " [2]: " << sss_mode[2] << " [3]: " << sss_mode[3] << std::endl;
        // getchar();
        // update sss modes
        auto findIter = std::find(sss_modes_.begin(), sss_modes_.end(), sss_mode);
        if (findIter != sss_modes_.end()) {
          // this sss mode is already in our library
          // just record this mode
          int id = std::distance(sss_modes_.begin(), findIter);
          for (int i = 0; i < total_sliding_directions; ++i) s_modes_[id].push_back(ss_modes[cs](row, i));
        } else {
          // this sss mode is new
          // std::cout << "[debug] sss_modes_ is growing, size: " << sss_modes_.size() << std::endl;
          std::cout << "[debug] sss_mode: " << sss_mode << std::endl;
          // getchar();
          sss_modes_.push_back(sss_mode);
          s_modes_.emplace_back();
          for (int i = 0; i < total_sliding_directions; ++i) s_modes_.back().push_back(ss_modes[cs](row, i));
        }
      }
      // end processing one mode in a cs mode
    } // end processing all modes in a cs mode
    // end processing one cs mode
  } // end processing all cs modes
  // clean ups: package the result into output format
  // *sss_modes = MatrixXi(sss_modes_.size(), num_contacts);
  MatrixXi sss_modes(sss_modes_.size(), num_contacts);
  std::vector<MatrixXi> s_modes;
  for (int i = 0; i < sss_modes_.size(); ++i) {
    // std::cout << "[debug] sss_modes_[i]:" << sss_modes_[i] << std::endl;
    // std::cout << "[debug] a" << std::endl;
    // getchar();
    for (int j = 0; j < num_contacts; ++j) {
      switch(sss_modes_[i].at(j)) {
        case 'f':
          // std::cout << "[debug] f" << std::endl;
          (sss_modes)(i, j) = -1;
          break;
        case '0':
          // std::cout << "[debug] 0" << std::endl;
          (sss_modes)(i, j) = 0;
          break;
        case '1':
          // std::cout << "[debug] 1" << std::endl;
          (sss_modes)(i, j) = 1;
          break;
        default:
          std::cerr << "[modeCleaning] wrong sign: " << sss_modes_[i].at(j) << std::endl;
          // return false;
      }
    }
  }
  for (int i = 0; i < s_modes_.size(); ++i) {
    // std::cout << "[debug] first row of s_modes_[i]:" << std::endl;
    // for (int j = 0; j < total_sliding_directions; ++j) std::cout << s_modes_[i][j] << " ";
    Eigen::MatrixXi modes_i_colmajor(total_sliding_directions, s_modes_[i].size()/total_sliding_directions);
    modes_i_colmajor = MatrixXi::Map(s_modes_[i].data(), modes_i_colmajor.rows(), modes_i_colmajor.cols());
    // std::cout << "[debug] modes_i_colmajor:\n" << modes_i_colmajor << std::endl;
    // getchar();
    s_modes.push_back(modes_i_colmajor.transpose());
  }
  return s_modes;
}

void wrenchSpaceAnalysis_wrapper(Eigen::MatrixXd Jac_e, Eigen::MatrixXd Jac_h,
    Eigen::MatrixXd eCone_allFix, Eigen::MatrixXd hCone_allFix,
    const Eigen::VectorXd &F_G,
    const double kContactForce, const double kFrictionE, const double kFrictionH,
    const double kCharacteristicLength, const int kNumSlidingPlanes,
    const Eigen::MatrixXi &e_cs_modes, const std::vector<Eigen::MatrixXi> &e_ss_modes,
    const Eigen::MatrixXi &h_cs_modes, const std::vector<Eigen::MatrixXi> &h_ss_modes,
    Eigen::MatrixXd G, const Eigen::VectorXd &b_G,
    const Eigen::MatrixXi &e_cs_modes_goal, const std::vector<Eigen::MatrixXi> &e_ss_modes_goal,
    const Eigen::MatrixXi &h_cs_modes_goal, const std::vector<Eigen::MatrixXi> &h_ss_modes_goal) {
  wrenchSpaceAnalysis(Jac_e, Jac_h, eCone_allFix, hCone_allFix,
      F_G, kContactForce, kFrictionE, kFrictionH, kCharacteristicLength, kNumSlidingPlanes,
      e_cs_modes, e_ss_modes, h_cs_modes, h_ss_modes, G, b_G,
      e_cs_modes_goal, e_ss_modes_goal, h_cs_modes_goal, h_ss_modes_goal);
  std::cout << "[wrapper] Finished!\n";
}

void wrenchSpaceAnalysis_2d_wrapper(Eigen::MatrixXd Jac_e, Eigen::MatrixXd Jac_h,
    Eigen::MatrixXd eCone_allFix, Eigen::MatrixXd hCone_allFix,
    const Eigen::VectorXd &F_G,
    const double kContactForce, const double kFrictionE, const double kFrictionH,
    const double kCharacteristicLength,
    Eigen::MatrixXd G, const Eigen::VectorXd &b_G,
    const Eigen::MatrixXi &e_modes, const Eigen::MatrixXi &h_modes,
    const Eigen::VectorXi &e_modes_goal, const Eigen::VectorXi &h_modes_goal) {
  wrenchSpaceAnalysis_2d(Jac_e, Jac_h, eCone_allFix, hCone_allFix, F_G,
      kContactForce, kFrictionE, kFrictionH, kCharacteristicLength,
      G, b_G, e_modes, h_modes, e_modes_goal, h_modes_goal);
  std::cout << "[wrapper] Finished!\n";
}

PYBIND11_MODULE(wrenchStampingLib, m) {
    m.doc() = "pybind11 wrenchStampingLib plugin"; // optional module docstring

    m.def("modeCleaning", &modeCleaning_test, "Test mode cleaning");
    m.def("wrenchSpaceAnalysis", &wrenchSpaceAnalysis_wrapper, "Test mode cleaning");
    m.def("wrenchSpaceAnalysis_2d", &wrenchSpaceAnalysis_2d_wrapper, "Test mode cleaning");
}


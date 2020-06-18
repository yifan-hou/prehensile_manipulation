#include "wrench_space_analysis.h"
#include "timer.h"
#include "solvehfvc.h"

#include <list>
#include <string>
#include <algorithm>
#include <iostream>

// #include "solvehfvc.h"
#include "polyhedron.h"

#define TOL 1e-8
#define PI 3.1415926

using namespace Eigen;

typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;

// Geometrical parameters:
//  Jac_e, Jac_h
//  eCone_allFix_r, hCone_allFix_r: each row denotes a generator
//  F_G: gravity vector
// Magnitude parameters:
//  kContactForce
//  kCharacteristicLength
// List of all velocity-consistent contact modes
//  e_cs_modes, h_cs_modes: number of cs modes x number of contacts
//  e_ss_modes, h_ss_modes: a vector, each element is a (number of modes x number of total sliding directions) matrix describing the stick-sliding modes for a particular cs modes.
// Optional arguments (use an empty matrix to denote an unused argument):
//  A, b: additional force constraint
//  G, b_G: goal velocity description
//  e_mode_goal, h_mode_goal: a particular goal mode
void wrenchSpaceAnalysis(MatrixXd Jac_e, MatrixXd Jac_h,
    MatrixXd eCone_allFix_r, MatrixXd hCone_allFix_r,
    const VectorXd &F_G, const double kContactForce,
    const double kCharacteristicLength, const int kNumSlidingPlanes,
    const MatrixXi &e_cs_modes, const std::vector<MatrixXi> &e_ss_modes,
    const MatrixXi &h_cs_modes, const std::vector<MatrixXi> &h_ss_modes,
    MatrixXd G, const VectorXd &b_G,
    const MatrixXi &e_cs_modes_goal, const std::vector<MatrixXi> &e_ss_modes_goal,
    const MatrixXi &h_cs_modes_goal, const std::vector<MatrixXi> &h_ss_modes_goal) {

  std::cout << "[wrenchSpaceAnalysis] Calling..\n";
  // std::cout << "G:\n" << G << std::endl;
  // std::cout << "b_G:\n" << b_G << std::endl;
  // getchar();
  // std::cout << "Jac_e: " << Jac_e.rows() << " x " << Jac_e.cols() << std::endl;
  // std::cout << "Jac_h: " << Jac_h.rows() << " x " << Jac_h.cols() << std::endl;
  // std::cout << "eCone_allFix_r: " << eCone_allFix_r.rows() << " x " << eCone_allFix_r.cols() << std::endl;
  // std::cout << "hCone_allFix_r: " << hCone_allFix_r.rows() << " x " << hCone_allFix_r.cols() << std::endl;
  // std::cout << "e_cs_modes: " << e_cs_modes.rows() << " x " << e_cs_modes.cols() << std::endl;
  // std::cout << "h_cs_modes: " << h_cs_modes.rows() << " x " << h_cs_modes.cols() << std::endl;
  // std::cout << "e_ss_modes: " << e_ss_modes.size() << std::endl;
  // std::cout << "h_ss_modes: " << h_ss_modes.size() << std::endl;
  // std::cout << "G: " << G.rows() << " x " << G.cols() << std::endl;
  // std::cout << "F_G: " << F_G.size() << std::endl;
  // std::cout << "b_G: " << b_G.size() << std::endl;
  // std::cout << "e_mode_goal: " << e_mode_goal.size() << std::endl;
  // std::cout << "h_mode_goal: " << h_mode_goal.size() << std::endl;
  // getchar();
  Timer timer;
  timer.tic();

  bool flag_given_goal_velocity = false;
  if (b_G.size() > 0) flag_given_goal_velocity = true;

  // scaling for generalized velocity
  // V = gvscale * V_scaled
  Vector6d vscale_vec, vscale_inv_vec;
  vscale_vec << 1., 1., 1., 1./kCharacteristicLength, 1./kCharacteristicLength, 1./kCharacteristicLength;
  vscale_inv_vec << 1., 1., 1., kCharacteristicLength, kCharacteristicLength, kCharacteristicLength;
  MatrixXd vscale = vscale_vec.asDiagonal();
  MatrixXd vscale_inv = vscale_inv_vec.asDiagonal();
  MatrixXd gvscale = MatrixXd::Zero(12, 12);
  gvscale.block<6, 6>(0, 0) = vscale;
  gvscale.block<6, 6>(6, 6) = vscale;
  // N_ * V_scaled = 0, N*V = 0,
  Jac_e = Jac_e * vscale;
  Jac_h = Jac_h * vscale;
  eCone_allFix_r = eCone_allFix_r * vscale;
  hCone_allFix_r = hCone_allFix_r * vscale;
  if (flag_given_goal_velocity) G = G*gvscale;

  // for crashing check
  MatrixXd cone_allFix_r;
  Poly::coneIntersection(eCone_allFix_r, hCone_allFix_r, &cone_allFix_r);

  Eigen::MatrixXi e_sss_modes, h_sss_modes;
  std::vector<Eigen::MatrixXi> e_s_modes, h_s_modes;
  int kNumEContacts = e_cs_modes.cols();
  int kNumHContacts = h_cs_modes.cols();

  if (!modeCleaning(e_cs_modes, e_ss_modes, kNumSlidingPlanes, &e_sss_modes, &e_s_modes)) {
    std::cerr << "[wrenchSpaceAnalysis] failed to call modeCleaning for e contacts." << std::endl;
    return;
  }
  if (!modeCleaning(h_cs_modes, h_ss_modes, kNumSlidingPlanes, &h_sss_modes, &h_s_modes)) {
    std::cerr << "[wrenchSpaceAnalysis] failed to call modeCleaning for h contacts." << std::endl;
    return;
  }

  Eigen::MatrixXi e_sss_modes_goal, h_sss_modes_goal;
  std::vector<Eigen::MatrixXi> e_s_modes_goal, h_s_modes_goal;

  if (e_cs_modes_goal.size() == 0 ) {
    e_sss_modes_goal = e_sss_modes;
    h_sss_modes_goal = h_sss_modes;
    e_s_modes_goal = e_s_modes;
    h_s_modes_goal = h_s_modes;
  } else {
    if (!modeCleaning(e_cs_modes_goal, e_ss_modes_goal, kNumSlidingPlanes, &e_sss_modes_goal, &e_s_modes_goal)) {
      std::cerr << "[wrenchSpaceAnalysis] failed to call modeCleaning for goal e contacts." << std::endl;
      return;
    }
    if (!modeCleaning(h_cs_modes_goal, h_ss_modes_goal, kNumSlidingPlanes, &h_sss_modes_goal, &h_s_modes_goal)) {
      std::cerr << "[wrenchSpaceAnalysis] failed to call modeCleaning for goal h contacts." << std::endl;
      return;
    }
  }

  std::cout << "timer: modeCleaning time = " << timer.toc() << "ms" << std::endl;

  Eigen::VectorXd F = Eigen::VectorXd::Zero(12);
  Eigen::MatrixXd N;
  int kDimActualized = 6;
  int kDimUnActualized = 6;
  int kNumSeeds = 5;
  int kPrintLevel = 0;
  HFVC action;

  Eigen::VectorXi e_sss_mode_goal, h_sss_mode_goal;
  Eigen::VectorXi e_sss_mode, h_sss_mode;

  Eigen::FullPivLU<MatrixXd> lu;
  lu.setThreshold(TOL);

  // Eigen::VectorXi e_s_mode, h_s_mode;

  std::cout << "##         Begin loop         ##\n";
  timer.tic();
  for (int e_sss_i_goal = 0; e_sss_i_goal < e_sss_modes_goal.rows(); ++e_sss_i_goal) {
    e_sss_mode_goal = e_sss_modes_goal.middleRows(e_sss_i_goal, 1).transpose();
    for (int h_sss_i_goal = 0; h_sss_i_goal < h_sss_modes_goal.rows(); ++h_sss_i_goal) {
      h_sss_mode_goal = h_sss_modes_goal.middleRows(h_sss_i_goal, 1).transpose();
      std::cout << "[WrenchStamping] goal id: e " << e_sss_i_goal << ", h " << h_sss_i_goal << std::endl;
      std::cout << "[WrenchStamping] goal mdoe: e " << e_sss_mode_goal.transpose() << ", h " << h_sss_mode_goal.transpose() << std::endl;

      std::cout << "[WrenchStamping] 1. HFVC" << std::endl;
      Eigen::MatrixXd N, Nu, T; // TODO: Nu and T seem useless here
      getConstraintOfTheMode(Jac_e, Jac_h,
          e_sss_mode_goal, h_sss_mode_goal,
          e_s_modes[0].middleRows(0, 1).transpose(), // just temporary
          h_s_modes[0].middleRows(0, 1).transpose(),
          &N, &Nu, &T);

      if (!solvehfvc_new(N, G, b_G, F, kDimActualized, kDimUnActualized, kNumSeeds,
        kPrintLevel, &action)) {
        std::cout << "[WrenchStamping]    HFVC has no solution." << std::endl;
        continue;
      }
      assert(action.n_af < kDimActualized);
      assert(action.n_af > 0);
      // make sure all velocity commands >= 0
      for (int i = 0; i < action.n_av; ++i) {
        if (action.b_C(i) < 0) {
          action.b_C(i) = - action.b_C(i);
          action.C.middleRows(i, 1) = - action.C.middleRows(i, 1);
          action.R_a.middleRows(i + action.n_af, 1) = - action.R_a.middleRows(i + action.n_af, 1);
        }
      }
      action.w_av = action.b_C;
      MatrixXd V_control_directions_r = -action.R_a.bottomRows(action.n_av);
      MatrixXd F_control_directions_r(2*action.n_af, kDimActualized);
      F_control_directions_r << action.R_a.topRows(action.n_af),
          -action.R_a.topRows(action.n_af);
      // std::cout << "[WrenchStamping]    Debug: action.C\n" << action.C << std::endl;
      // std::cout << "[WrenchStamping]    Debug: action.b_C\n" << action.b_C << std::endl;

      // Crashing check
      MatrixXd R;
      Poly::coneIntersection(cone_allFix_r, V_control_directions_r, &R);
      if (R.rows() > 0) {
        std::cout << "[WrenchStamping]    Crashing." << std::endl;
        continue;
      }

      std::cout << "[WrenchStamping] 2. Check Velocity Feasibility." << std::endl;
      // How to filter out modes:
      // 1. If NC degenerates, mark this mode as incompatible;
      // 2. If NC gives unique solution, record the cone of this mode
      // 3. If NC gives multiple solutions, record the cone of the all sticking mode.
      for (int e_sss_i = 0; e_sss_i < e_sss_modes.rows(); ++e_sss_i) {
        e_sss_mode = e_sss_modes.middleRows(e_sss_i, 1).transpose();
        for (int h_sss_i = 0; h_sss_i < h_sss_modes.rows(); ++h_sss_i) {
          h_sss_mode = h_sss_modes.middleRows(h_sss_i, 1).transpose();
          std::cout << "[WrenchStamping]    Checking id: " << e_sss_i << ", " << h_sss_i;
          std::cout << " (e: " << e_sss_mode.transpose() << ", h: " << h_sss_mode.transpose() << ")\t";
          std::cout << e_s_modes[e_sss_i].rows() << " s modes:";

          getConstraintOfTheMode(Jac_e, Jac_h,
              e_sss_mode, h_sss_mode,
              e_s_modes[0].middleRows(0, 1).transpose(), // just temporary
              h_s_modes[0].middleRows(0, 1).transpose(),
              &N, &Nu, &T);

          Eigen::MatrixXd NC(N.rows() + action.C.rows(), N.cols());
          NC << N, action.C;
          Eigen::VectorXd b_NC = Eigen::VectorXd::Zero(NC.rows());
          b_NC.tail(action.b_C.rows()) = action.b_C;
          lu.compute(NC);
          Eigen::VectorXd sol_NC = lu.solve(b_NC);
          bool a_solution_exists = (NC*sol_NC).isApprox(b_NC, 10.*TOL);
          if (a_solution_exists) {
            Eigen::MatrixXd null_NC = lu.kernel();
            if (null_NC.norm() > TOL) {
              std::cout << "Multiple Solutions, kernel Dim = " << null_NC.cols() << "contact dim: ";
              // Multiple solutions
              for (int i = 0; i < kNumEContacts; ++i) {
                lu.compute(T.middleRows(2*i, 2)*null_NC);
                std::cout << lu.rank() << ", ";
              }
              std::cout << std::endl;
            } else {
              std::cout << "Unique Solution. ";
              if ((Nu.rows() > 0) && ((Nu*sol_NC).minCoeff() < -TOL)) {
                std::cout << "Violates inequalities. " << std::endl;
                continue;
              }
              std::cout << std::endl;
              // find the cone of this Unique solution

            }
          } else {
            std::cout << "No Solution. " << std::endl;
            // no solution. Continue
          }
        }
      }
      // done with inner SSS loop

      // for (int e_s_i = 0; e_s_i < e_s_modes[e_sss_i_goal].rows(); ++e_s_i) {
      //   e_s_mode = e_s_modes[e_sss_i_goal].middleRows(e_s_i, 1).transpose();
      //   e_cone = getConeOfTheMode(eCone_allFix_r, e_sss_mode, e_s_mode, kNumSlidingPlanes);
      //   for (int h_s_i = 0; h_s_i < h_s_modes[h_sss_i_goal].rows(); ++h_s_i) {
      //     h_s_mode = h_s_modes[h_sss_i_goal].middleRows(h_s_i, 1).transpose();
      //     h_cone = getConeOfTheMode(hCone_allFix_r, h_sss_mode, h_s_mode, kNumSlidingPlanes);
      //     Poly::coneIntersection(e_cone, h_cone, &R);
      //     // std::cout << "coneIntersection for " << e_sss_mode.transpose() << ", at time " << timer.toc() << "ms" << std::endl;
      //     if(R.rows() == 0) empty ++;
      //     else non_empty ++;
      //   }
      // }
    }
  }

  // Eigen::MatrixXd N = getConstraintOfTheMode(Jac_e, Jac_h,
  //     e_sss_modes.middleRows(0, 1).transpose(), h_sss_modes.middleRows(0, 1).transpose(),
  //     e_s_modes[0].middleRows(0, 1).transpose(), h_s_modes[0].middleRows(0, 1).transpose(),
  //     kNumSlidingPlanes);

  // Eigen::VectorXd F = Eigen::VectorXd::Zero(12);
  // int kDimActualized = 6;
  // int kDimUnActualized = 6;
  // int kNumSeeds = 5;
  // int kPrintLevel = 0;
  // HFVC action;

  // std::cout << "N:\n" << N << std::endl;
  // std::cout << "G:\n" << G << std::endl;
  // std::cout << "b_G:\n" << b_G << std::endl;
  // std::cout << "F:\n" << F << std::endl;

  // timer.tic();
  // for (int i = 0; i < 100; ++i) {
  //   solvehfvc(N, G, b_G, F, kDimActualized, kDimUnActualized, kNumSeeds,
  //     kPrintLevel, &action);
  // }
  // std::cout << "timer: HFVC time x 100 = " << timer.toc() << "ms" << std::endl;
  // std::cout << "solution:\n";
  // std::cout << "n_av: " << action.n_av << std::endl;
  // std::cout << "n_af: " << action.n_af << std::endl;
  // std::cout << "R_a:\n" << action.R_a << std::endl;
  // std::cout << "w_av:\n" << action.w_av << std::endl;

  // std::cout << "\n\nNew Algorithm\n";
  // timer.tic();
  // for (int i = 0; i < 100; ++i) {
  //   solvehfvc_new(N, G, b_G, F, kDimActualized, kDimUnActualized, kNumSeeds,
  //     kPrintLevel, &action);
  // }
  // std::cout << "timer: HFVC new time x 100 = " << timer.toc() << "ms" << std::endl;
  return;

  //     drawWrenchSpace(cone_allFix, eCone_allFix_r, hCone_allFix_r, ...
  //             V_control_directions, F_control_directions, eh_cones, eh_modes, ...
  //             eh_cone_feasible_mode_count, feasibilities, goal_id);

  //     eh_cones_goal_m = eh_cones{m};
  //     feasibilities(m) = 0;
  //     eh_cones_other_feasible_m = eh_cones(feasibilities);

  //     fprintf("===============================\n");
  //     fprintf("===  Compute Force Action   ===\n");
  //     fprintf("===============================\n");

  //     % action selection
  //     force_basis = R_a_inv(:, 1:n_af);
  //     [force_action, shape_margin] = forceControl(force_basis, ...
  //             eh_cones_goal_m, eh_cones_other_feasible_m);

  //     if ~isempty(force_action)
  //         disp('Force command:');
  //         disp(force_action);
  //         solution.eh_mode = eh_mode_goal;
  //         solution.n_af = n_af;
  //         solution.n_av = n_av;
  //         solution.w_av = w_av;
  //         solution.eta_af = -kForceMagnitude*force_action; % the minus sign comes from force balance
  //         solution.margin = min(shape_margin, margins(m));

  //         R_a_inv = R_a^-1;
  //         Cf_inv = vscale_inv*R_a_inv; Cf_inv = Cf_inv(:, 1:n_af);
  //         Cv_inv = vscale*R_a_inv; Cv_inv = Cv_inv(:, end-n_av+1:end);
  //         solution.R_a_inv = [Cf_inv Cv_inv];
  //         solution.R_a = solution.R_a_inv^-1;

  //         solutions_count = solutions_count + 1;
  //         solutions{solutions_count} = solution;

  //         if goal_mode_is_given
  //             break;
  //         end
  //     else
  //         disp('Failure: no distinguishable force command.');
  //         if goal_mode_is_given
  //             return;
  //         else
  //             continue;
  //         end
  //     end
  // end

  // fprintf("###############################################\n");
  // fprintf("##                  Results                  ##\n");
  // fprintf("###############################################\n");
  // fprintf("Total number of feasible modes found: %d\n", solutions_count);
  // solution = [];
  // if solutions_count > 0
  //     if goal_mode_is_given
  //         solution = solutions(1);
  //     else
  //         solutions = solutions(1:solutions_count);
  //         margins = zeros(1, solutions_count);
  //         for i = 1:solutions_count
  //             mode_text = printModes(solutions{i}.eh_mode, false);
  //             fprintf("Mode %d: %s\nmargin: %f\n", i, mode_text, solutions{i}.margin);
  //             margins(i) = solutions{i}.margin;
  //         end
  //         disp('Best solution:');
  //         [~, best_solution_id] = max(margins);
  //         solution = solutions{best_solution_id};
  //     end

  //     printModes(solution.eh_mode);
  //     V_T = solution.R_a_inv*[zeros(solution.n_af,1); solution.w_av];
  //     F_T = solution.R_a_inv*[solution.eta_af; zeros(solution.n_av, 1)];
  //     disp('R_a:');
  //     disp(solution.R_a);
  //     disp('V_T:');
  //     disp(V_T);
  //     disp('F_T:');
  //     disp(F_T);
  // end

}

bool modeCleaning(const MatrixXi &cs_modes, const std::vector<MatrixXi> &ss_modes, int kNumSlidingPlanes,
    MatrixXi *sss_modes, std::vector<MatrixXi> *s_modes) {
  int num_cs_modes = cs_modes.rows();
  int num_contacts = cs_modes.cols();
  int total_sliding_directions = ss_modes[0].cols();
  assert(sss_modes->size() == 0);
  /**
   * First, we need to build a cs mode representation with {f = sticking, 0 = sliding, 1 = separation}
   */
  std::vector<std::string> sss_modes_; // sticking, sliding, separation
  std::vector<std::vector<int>> s_modes_;

  Eigen::VectorXi cs_mode;
  std::string sss_mode;
  sss_mode.resize(num_contacts);
  for (int cs = 0; cs < num_cs_modes; ++cs) {
    cs_mode = cs_modes.middleRows(cs, 1).transpose();
    std::cout << "  cs_mode: " << cs_mode.transpose() << std::endl;
    int num_ss_modes = ss_modes[cs].rows();
    for (int row = 0; row < num_ss_modes; ++row) {
      // process one contact mode
      bool redundant = false;
      for (int i = 0; i < num_contacts; ++i) {
        if (cs_mode(i) == 0) {
          int ss_mode_contact_i_sum = ss_modes[cs].block(row,i*kNumSlidingPlanes, 1, kNumSlidingPlanes).cwiseAbs().sum();
          if (ss_mode_contact_i_sum == 0)
            sss_mode[i] = 'f';
          else if (ss_mode_contact_i_sum == kNumSlidingPlanes) {
            sss_mode[i] = '0';
          } else {
            redundant = true;
            break;
          }
        } else {
          sss_mode[i] = '1';
        }
      }
      if (!redundant) {
        // update sss modes
        auto findIter = std::find(sss_modes_.begin(), sss_modes_.end(), sss_mode);
        if (findIter != sss_modes_.end()) {
          // this sss mode is already in our library
          // just record this mode
          int id = std::distance(sss_modes_.begin(), findIter);
          for (int i = 0; i < total_sliding_directions; ++i) s_modes_[id].push_back(ss_modes[cs](row, i));
        } else {
          // this sss mode is new
          std::cout << "sss_mode: " << sss_mode << std::endl;
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
  *sss_modes = MatrixXi(sss_modes_.size(), num_contacts);
  for (int i = 0; i < sss_modes_.size(); ++i) {
    for (int j = 0; j < num_contacts; ++j) {
      switch(sss_modes_[i].at(j)) {
        case 'f':
          (*sss_modes)(i, j) = -1;
          break;
        case '0':
          (*sss_modes)(i, j) = 0;
          break;
        case '1':
          (*sss_modes)(i, j) = 1;
          break;
        default:
          std::cerr << "[modeCleaning] wrong sign: " << sss_modes_[i].at(j) << std::endl;
          return false;
      }
    }
  }
  for (int i = 0; i < s_modes_.size(); ++i) {
    Eigen::MatrixXi modes_i_colmajor(total_sliding_directions, s_modes_[i].size()/total_sliding_directions);
    modes_i_colmajor = MatrixXi::Map(s_modes_[i].data(), modes_i_colmajor.rows(), modes_i_colmajor.cols());
    s_modes->push_back(modes_i_colmajor.transpose());
  }
  return true;
}

Eigen::MatrixXd getConeOfTheMode(const Eigen::MatrixXd &cone_allFix,
    const Eigen::VectorXi &sss_mode, const Eigen::VectorXi &s_mode, int kNumSlidingPlanes) {

  int num_contacts = sss_mode.size();
  Eigen::VectorXi s_mode01 = s_mode;
  s_mode01.array() = (s_mode01.array() + 1)/2; // -1/1 -> 0/1
  std::vector<double> generators; // concatenation of all generators
  for (int i = 0; i < num_contacts; ++i) {
    // -1: sticking   0: sliding   1: separation
    int id_start = (2*kNumSlidingPlanes+1)*i;
    MatrixXd generators_add(0, 0);
    if (sss_mode[i] == -1) {
      // record all 2*kNumSlidingPlanes generators
      generators_add = cone_allFix.middleRows(id_start, 2*kNumSlidingPlanes);
    } else if (sss_mode[i] == 0) {
      // find two corresponding generators
      std::vector<double> v;
      v.resize(2*6);

      // Compute the order from bit-wise code:
      int id_s = kNumSlidingPlanes*i;
      int mode_sum = s_mode01.segment(id_s, kNumSlidingPlanes).sum();
      int id_g;
      if (s_mode01(id_s) == 1)
          id_g = mode_sum - 1;
      else
          id_g = 2*kNumSlidingPlanes - mode_sum - 1;
      generators_add = cone_allFix.middleRows(id_start + id_g, 2);
    }

    // add generators
    if (generators_add.size() > 0) {
      int g_size = generators.size();
      generators.resize(g_size + generators_add.rows()*6);
      Eigen::MatrixXd::Map(&generators[g_size], 6, generators_add.rows()) = generators_add.transpose();
    }
  }
  // convert generators to matrix form
  int num_generators = generators.size()/6;
  return MatrixXd::Map(generators.data(), 6, num_generators).transpose();
}

// for now, only implemented the bilateral part for hybrid servoing
void getConstraintOfTheMode(
    const Eigen::MatrixXd &J_e_AF, const Eigen::MatrixXd &J_h_AF,
    const Eigen::VectorXi &sss_mode_e, const Eigen::VectorXi &sss_mode_h,
    const Eigen::VectorXi &s_mode_e, const Eigen::VectorXi &s_mode_h,
    Eigen::MatrixXd *N, Eigen::MatrixXd *Nu, Eigen::MatrixXd *T) {
  const int kNumEContacts = sss_mode_e.size();
  const int kNumHContacts = sss_mode_h.size();
  const int kDim = J_h_AF.cols();
  assert(J_e_AF.cols() == kDim);

  // count number of rows
  // -1: sticking   0: sliding   1: separation
  int NeRows = 0, NhRows = 0;
  int NueRows = 0, NuhRows = 0;
  for (int i = 0; i < kNumEContacts; ++i) {
    if (sss_mode_e[i] == -1) {
      NeRows += 3;
    } else if (sss_mode_e[i] == 0) {
      NeRows += 1;
    } else {
      NueRows += 1;
    }
  }
  for (int i = 0; i < kNumHContacts; ++i) {
    if (sss_mode_h[i] == -1) {
      NhRows += 3;
    } else if (sss_mode_h[i] == 0) {
      NhRows += 1;
    } else {
      NuhRows += 1;
    }
  }

  *N  = Eigen::MatrixXd::Zero(NeRows+NhRows, 2*kDim);
  *Nu = Eigen::MatrixXd::Zero(NueRows+NuhRows, 2*kDim);
  *T = Eigen::MatrixXd::Zero(2*(kNumEContacts + kNumHContacts), 2*kDim);

  int n_count = 0, nu_count = 0;
  for (int i = 0; i < kNumEContacts; ++i) {
    T->block(2*i, 0, 2, kDim) = J_e_AF.middleRows(kNumEContacts + i*2, 2);
    // -1: sticking   0: sliding   1: separation
    if (sss_mode_e[i] == 1) {
      // separation
      Nu->block(nu_count++, 0, 1, kDim) = J_e_AF.middleRows(i, 1);
      continue;
    }
    N->block(n_count++, 0, 1, kDim) = J_e_AF.middleRows(i, 1);
    if (sss_mode_e[i] == -1) {
      // sticking
      N->block(n_count, 0, 2, kDim) = J_e_AF.middleRows(kNumEContacts + i*2, 2);
      n_count += 2;
    }
  }
  for (int i = 0; i < kNumHContacts; ++i) {
    // -1: sticking   0: sliding   1: separation
    T->block(2*i, 0, 2, kDim) = -J_h_AF.middleRows(kNumHContacts + i*2, 2);
    T->block(2*i, kDim, 2, kDim) = J_h_AF.middleRows(kNumHContacts + i*2, 2);
    if (sss_mode_h[i] == 1) {
      // separation
      Nu->block(nu_count, 0, 1, kDim) = -J_h_AF.middleRows(i, 1);
      Nu->block(nu_count++, kDim, 1, kDim) = J_h_AF.middleRows(i, 1);
      continue;
    }
    N->block(n_count, 0, 1, kDim) = -J_h_AF.middleRows(i, 1);
    N->block(n_count++, kDim, 1, kDim) = J_h_AF.middleRows(i, 1);
    if (sss_mode_h[i] == -1) {
      // sticking
      N->block(n_count, 0, 2, kDim) = -J_h_AF.middleRows(kNumHContacts + i*2, 2);
      N->block(n_count, kDim, 2, kDim) = J_h_AF.middleRows(kNumHContacts + i*2, 2);
      n_count += 2;
    }
  }
}
#include "wrench_space_analysis.h"
#include "timer.h"
#include "solvehfvc.h"

#include <list>
#include <string>
#include <algorithm>
#include <iostream>

// #include "solvehfvc.h"
#include "polyhedron.h"

#define TOL 1e-7
#define PI 3.1415926

using namespace Eigen;

typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;

// Geometrical parameters:
//  Jac_e, Jac_h
//  eCone_allFix, hCone_allFix: each row denotes a generator
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
    MatrixXd eCone_allFix, MatrixXd hCone_allFix,
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
  // std::cout << "eCone_allFix: " << eCone_allFix.rows() << " x " << eCone_allFix.cols() << std::endl;
  // std::cout << "hCone_allFix: " << hCone_allFix.rows() << " x " << hCone_allFix.cols() << std::endl;
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
  eCone_allFix = eCone_allFix * vscale;
  hCone_allFix = hCone_allFix * vscale;
  if (flag_given_goal_velocity) G = G*gvscale;

  // for crashing check
  MatrixXd cone_allFix;
  Poly::coneIntersection(eCone_allFix, hCone_allFix, &cone_allFix);

  Eigen::MatrixXi e_sss_modes, h_sss_modes;
  std::vector<Eigen::MatrixXi> e_s_modes, h_s_modes;
  int num_e_contacts = e_cs_modes.cols();
  int num_h_contacts = h_cs_modes.cols();

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
  // Eigen::VectorXi e_s_mode, h_s_mode;

  std::cout << "##         Begin loop         ##\n";
  timer.tic();
  for (int e_sss_i_goal = 0; e_sss_i_goal < e_sss_modes_goal.rows(); ++e_sss_i_goal) {
    e_sss_mode_goal = e_sss_modes_goal.middleRows(e_sss_i_goal, 1).transpose();
    for (int h_sss_i_goal = 0; h_sss_i_goal < h_sss_modes_goal.rows(); ++h_sss_i_goal) {
      h_sss_mode_goal = h_sss_modes_goal.middleRows(h_sss_i_goal, 1).transpose();
      std::cout << "[WrenchStamping] goal id: e " << e_sss_i_goal << ", h " << h_sss_i_goal << std::endl;

      std::cout << "[WrenchStamping] 1. HFVC" << std::endl;
      N = getConstraintOfTheMode(Jac_e, Jac_h,
          e_sss_mode_goal, h_sss_mode_goal,
          e_s_modes[0].middleRows(0, 1).transpose(), // just temporary
          h_s_modes[0].middleRows(0, 1).transpose(),
          kNumSlidingPlanes);
      if (!solvehfvc_new(N, G, b_G, F, kDimActualized, kDimUnActualized, kNumSeeds,
        kPrintLevel, &action)) {
        std::cout << "[WrenchStamping]    HFVC has no solution." << std::endl;
        continue;
      }
      

      // if isempty(R_a) || any(isnan(w_av))
      //     disp('[Hybrid servoing] solvehfvc returns no solution.')
      //     n_av = []; n_af = [];R_a = [];R_a_inv=[];w_av=[];Cv=[]; b_C=[];
      //     return;
      // end

      // fprintf('[Hybrid Servoing] force dimension: %d\n', n_af);
      // fprintf('[Hybrid Servoing] velocity dimension: %d\n', n_av);

      // assert(n_af < 3);
      // assert(n_af > 0);

      // % make sure all velocity commands >= 0
      // R_id_flip = find(w_av < 0);
      // R_a(n_af + R_id_flip, :) = - R_a(n_af + R_id_flip, :);
      // w_av(R_id_flip) = - w_av(R_id_flip);

      // Cv = [zeros(n_av, 3), R_a(n_af + 1: end, :)];
      // b_C = w_av;
      // R_a_inv = R_a^-1;







      // for (int e_s_i = 0; e_s_i < e_s_modes[e_sss_i_goal].rows(); ++e_s_i) {
      //   e_s_mode = e_s_modes[e_sss_i_goal].middleRows(e_s_i, 1).transpose();
      //   e_cone = getConeOfTheMode(eCone_allFix, e_sss_mode, e_s_mode, kNumSlidingPlanes);
      //   for (int h_s_i = 0; h_s_i < h_s_modes[h_sss_i_goal].rows(); ++h_s_i) {
      //     h_s_mode = h_s_modes[h_sss_i_goal].middleRows(h_s_i, 1).transpose();
      //     h_cone = getConeOfTheMode(hCone_allFix, h_sss_mode, h_s_mode, kNumSlidingPlanes);
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

  // /**
  //  * Preparation
  //  */
  // const int NE_SSS_MODES = e_sss_modes.rows();
  // const int NH_SSS_MODES = h_sss_modes.rows();

  // std::cout << "##         Check CS modes F-feasibility         ##\n";
  // Eigen::VectorXi e_sss_mode, h_sss_mode;
  // Eigen::MatrixXd e_cone, h_cone;
  // Eigen::MatrixXd R;

  // // check every s mode
  // Eigen::VectorXi e_s_mode, h_s_mode;
  // int empty = 0;
  // int non_empty = 0;
  // for (int e_sss_i = 0; e_sss_i < NE_SSS_MODES; ++e_sss_i) {
  //   e_sss_mode = e_sss_modes.middleRows(e_sss_i, 1).transpose();
  //   for (int h_sss_i = 0; h_sss_i < NH_SSS_MODES; ++h_sss_i) {
  //     h_sss_mode = h_sss_modes.middleRows(h_sss_i, 1).transpose();
  //     for (int e_s_i = 0; e_s_i < e_s_modes[e_sss_i].rows(); ++e_s_i) {
  //       e_s_mode = e_s_modes[e_sss_i].middleRows(e_s_i, 1).transpose();
  //       e_cone = getConeOfTheMode(eCone_allFix, e_sss_mode, e_s_mode, kNumSlidingPlanes);
  //       // std::cout << "\n\n\n";
  //       // std::cout << "debug, cone_allFix: \n" << eCone_allFix << std::endl;
  //       // std::cout << "debug, sss_mode: \n" << e_sss_mode << std::endl;
  //       // std::cout << "debug, s_mode: \n" << e_s_mode << std::endl;
  //       // std::cout << "debug, result cone: \n" << e_cone << std::endl;
  //       // getchar();
  //       for (int h_s_i = 0; h_s_i < h_s_modes[h_sss_i].rows(); ++h_s_i) {
  //         h_s_mode = h_s_modes[h_sss_i].middleRows(h_s_i, 1).transpose();
  //         h_cone = getConeOfTheMode(hCone_allFix, h_sss_mode, h_s_mode, kNumSlidingPlanes);
  //         Poly::coneIntersection(e_cone, h_cone, &R);
  //         // std::cout << "coneIntersection for " << e_sss_mode.transpose() << ", at time " << timer.toc() << "ms" << std::endl;
  //         if(R.rows() == 0) empty ++;
  //         else non_empty ++;
  //       }
  //     }
  //   }
  // }

  // std::cout << "\nEmpty: " << empty << std::endl;
  // std::cout << "Non-empty: " << non_empty << std::endl;
  // std::cout << "timer: All time = " << timer.toc() << "ms" << std::endl;



  // %%
  // %% Begin to check each cone
  // %%
  // fprintf("###############################################\n");
  // fprintf("##            Velocity Filtering             ##\n");
  // fprintf("###############################################\n");

  // solutions = cell(eh_cone_feasible_mode_count, 1);
  // solutions_count = 0;

  // for m = 1:eh_cone_feasible_mode_count
  //     if goal_mode_is_given
  //         m = goal_id;
  //     else
  //         disp('***********');
  //         fprintf("Mode %d of %d\n", m, eh_cone_feasible_mode_count);
  //         printModes(eh_modes(:, m));
  //     end

  //     eh_mode_goal = eh_modes(:, m);
  //     fprintf("====================================\n");
  //     fprintf("= Hybrid Servoing & Crashing Check =\n");
  //     fprintf("====================================\n");
  //     N_all = Jacs{m};
  //     [n_av, n_af, R_a, R_a_inv, w_av, Cv, b_C] = hybridServoing(N_all, G, b_G);
  //     if isempty(n_av)
  //         disp("Failure: Hybrid Servoing returns no solution.")
  //         if goal_mode_is_given
  //             return;
  //         else
  //             continue;
  //         end
  //     end

  //     V_control_directions = -R_a_inv(:, end-n_av+1:end);
  //     F_control_directions = [R_a_inv(:, 1:n_af) -R_a_inv(:, 1:n_af)];
  //     intersection = coneIntersection(cone_allFix, V_control_directions);
  //     % Crashing check
  //     if ~isempty(intersection) && norm(intersection) > TOL
  //         disp('Failure: Crashing. The mode is not feasible.');
  //         if goal_mode_is_given
  //             return;
  //         else
  //             continue;
  //         end
  //     end

  //     fprintf("=======================\n");
  //     fprintf("=== Mode Filtering ===\n");
  //     fprintf("=======================\n");

  //     % How to filter out modes:
  //     % 1. If NC degenerates, mark this mode as incompatible;
  //     % 2. If nominal velocity under NC exists, and it cause the contact point to
  //     %    slide in a different direction, remove this mode.
  //     feasibilities = false(eh_cone_feasible_mode_count, 1);

  //     % cone of possible actuation forces
  //     %   force controlled direction can have both positive and negative forces

  //     feasible_mode_count = 0;
  //     flag_goal_infeasible = false;
  //     flag_goal_impossible = false;

  //     for n = 1:eh_cone_feasible_mode_count
  //         % filter out modes using velocity command
  //         N = Jacs{n};
  //         Nu = Jacus{n};

  //         N = rref(N);
  //         rank_N = rank(N);
  //         N = N(1:rank_N, :);

  //         % compute possible sliding directions
  //         % these computations are based on 'Criteria for Maintaining Desired Contacts for Quasi-Static Systems'
  //         Lambda_bar = [Cv; N];
  //         b_Lambda_bar = [b_C; zeros(size(N,1), 1)];

  //         compatible = false;
  //         if rank([Lambda_bar b_Lambda_bar], TOL) == rank(Lambda_bar, TOL)
  //             v_star = linsolve(Lambda_bar, b_Lambda_bar);
  //             if any(Nu*v_star < -TOL)
  //                 % this mode can not exist
  //                 % V-Impossible
  //                 if n == m
  //                     flag_goal_impossible = true;
  //                     break;
  //                 else
  //                     continue;
  //                 end
  //             end
  //             compatible = true;
  //         else
  //             if n == m
  //                 flag_goal_infeasible = true;
  //                 break;
  //             end
  //         end

  //         % figure(1);clf(1);hold on;
  //         % drawCone(eh_cones{n},'g', true);
  //         % drawCone(W_action,'k', true);
  //         if compatible
  //             feasible_mode_count = feasible_mode_count + 1;
  //             feasibilities(n) = true;
  //         end
  //     end

  //     if flag_goal_infeasible
  //         disp('Goal mode violates velocity equality constraints. Discard this mode.');
  //         if goal_mode_is_given
  //             return;
  //         else
  //             continue;
  //         end
  //     end
  //     if flag_goal_impossible
  //         disp('Goal mode violates velocity inequality constraints. Discard this mode.');
  //         if goal_mode_is_given
  //             return;
  //         else
  //             continue;
  //         end
  //     end
  //     disp(['Remaining feasible modes: ' num2str(feasible_mode_count)]);

  //     drawWrenchSpace(cone_allFix, eCone_allFix, hCone_allFix, ...
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
Eigen::MatrixXd getConstraintOfTheMode(
    const Eigen::MatrixXd &J_e_AF, const Eigen::MatrixXd &J_h_AF,
    const Eigen::VectorXi &sss_mode_e, const Eigen::VectorXi &sss_mode_h,
    const Eigen::VectorXi &s_mode_e, const Eigen::VectorXi &s_mode_h,
    int kNumSlidingPlanes) {
  const int kNumEContacts = sss_mode_e.size();
  const int kNumHContacts = sss_mode_h.size();
  const int kDim = J_h_AF.cols();
  assert(J_e_AF.cols() == kDim);

  // count number of rows
  // -1: sticking   0: sliding   1: separation
  int NeRows = 0, NhRows = 0;
  for (int i = 0; i < kNumEContacts; ++i) {
    if (sss_mode_e[i] == -1) {
      NeRows += 1 + kNumSlidingPlanes;
    } else if (sss_mode_e[i] == 0) {
      NeRows += 1;
    }
  }
  for (int i = 0; i < kNumHContacts; ++i) {
    if (sss_mode_h[i] == -1) {
      NhRows += 1 + kNumSlidingPlanes;
    } else if (sss_mode_h[i] == 0) {
      NhRows += 1;
    }
  }
  Eigen::MatrixXd N;
  N = Eigen::MatrixXd::Zero(NeRows+NhRows, 2*kDim);
  int n_count = 0;
  for (int i = 0; i < kNumEContacts; ++i) {
    // -1: sticking   0: sliding   1: separation
    if (sss_mode_e[i] == 1) {
      // separation
      continue;
    }
    N.block(n_count++, 0, 1, kDim) = J_e_AF.middleRows(i, 1);
    if (sss_mode_e[i] == -1) {
      // sticking
      N.block(n_count, 0, kNumSlidingPlanes, kDim) = J_e_AF.middleRows(kNumEContacts + i*kNumSlidingPlanes, kNumSlidingPlanes);
      n_count += kNumSlidingPlanes;
    }
  }
  for (int i = 0; i < kNumHContacts; ++i) {
    // -1: sticking   0: sliding   1: separation
    if (sss_mode_h[i] == 1) {
      // separation
      continue;
    }
    N.block(n_count, 0, 1, kDim) = -J_h_AF.middleRows(i, 1);
    N.block(n_count++, kDim, 1, kDim) = J_h_AF.middleRows(i, 1);
    if (sss_mode_h[i] == -1) {
      // sticking
      N.block(n_count, 0, kNumSlidingPlanes, kDim) = -J_h_AF.middleRows(kNumHContacts + i*kNumSlidingPlanes, kNumSlidingPlanes);
      N.block(n_count, kDim, kNumSlidingPlanes, kDim) = J_h_AF.middleRows(kNumHContacts + i*kNumSlidingPlanes, kNumSlidingPlanes);
      n_count += kNumSlidingPlanes;
    }
  }
  return N;
}
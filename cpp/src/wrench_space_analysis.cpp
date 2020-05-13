#include <list>
#include <string>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "solvehfvc.h"
#include "polyhedron.h"

#define TOL 1e-7
#define PI 3.1415926

using namespace Eigen;

typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;


bool modeCleaning(const MatrixXi &cs_modes, const std::vector<MatrixXi> &ss_modes, int kNumSlidingPlanes,
    MatrixXi *sss_modes, std::vector<MatrixXi> *s_modes);

// Geometrical parameters:
//  Jac_e, Jac_h
//  eCone_allFix_raw, hCone_allFix_raw: each row denotes a generator
//  F_G: gravity vector
// Magnitude parameters:
//  kContactForce
//  kObjWeight
//  kCharacteristicLength
// List of all velocity-consistent contact modes
//  e_cs_modes, h_cs_modes: number of cs modes x number of contacts
//  e_ss_modes, h_ss_modes: a vector, each element is a (number of modes x number of total sliding directions) matrix describing the stick-sliding modes for a particular cs modes.
// Optional arguments (use an empty matrix to denote an unused argument):
//  A, b: additional force constraint
//  G, b_G: goal velocity description
//  e_mode_goal, h_mode_goal: a particular goal mode
void wrenchSpaceAnalysis(const MatrixXd &Jac_e_raw, const MatrixXd &Jac_h_raw,
    const MatrixXd &eCone_allFix_raw, const MatrixXd &hCone_allFix_raw,
    const VectorXd &F_G,
    const double kContactForce, const double kObjWeight,
    const double kCharacteristicLength, const int kNumSlidingPlanes,
    const MatrixXi &e_cs_modes, const std::vector<MatrixXi> &e_ss_modes,
    const MatrixXi &h_cs_modes, const std::vector<MatrixXi> &h_ss_modes,
    const MatrixXd &G, const VectorXd &b_G,
    const VectorXi &e_mode_goal, const VectorXi &h_mode_goal) {

  Eigen::MatrixXi e_sss_modes, h_sss_modes;
  std::vector<Eigen::MatrixXi> e_s_modes, h_s_modes;

  if (!modeCleaning(e_cs_modes, e_ss_modes, kNumSlidingPlanes, &e_sss_modes, &e_s_modes)) {
    std::cerr << "[wrenchSpaceAnalysis] failed to call modeCleaning for e contacts." << std::endl;
    return;
  }
  if (!modeCleaning(h_cs_modes, h_ss_modes, kNumSlidingPlanes, &h_sss_modes, &h_s_modes)) {
    std::cerr << "[wrenchSpaceAnalysis] failed to call modeCleaning for h contacts." << std::endl;
    return;
  }

  bool flag_given_goal_mode = false;
  if (e_mode_goal.size() > 0) flag_given_goal_mode = true;
  bool flag_given_goal_velocity = false;
  if (b_G.size() > 0) flag_given_goal_velocity = true;

  /**
   * Preparation
   */
  const int NE_SSS_MODES = e_sss_modes.rows();
  const int NH_SSS_MODES = h_sss_modes.rows();
  // scaling for generalized velocity
  // V = gvscale * V_scaled
  Matrix6d vscale = Vector6d(1., 1., 1., 1./kCharacteristicLength, 1./kCharacteristicLength, 1./kCharacteristicLength).asDiagonal();
  Matrix6d vscale_inv = Vector6d(1., 1., 1., kCharacteristicLength, kCharacteristicLength, kCharacteristicLength).asDiagonal();
  MatrixXd gvscale(12, 12);
  gvscale.block<6, 6>(0, 0) = vscale;
  gvscale.block<6, 6>(6, 6) = vscale;
  // N_ * V_scaled = 0, N*V = 0,
  // -> N_ = N*gvscale
  MatrixXd Jac_e = Jac_e_raw * gvscale;
  MatrixXd Jac_h = Jac_h_raw * gvscale;
  MatrixXd eCone_allFix = eCone_allFix_raw * vscale;
  MatrixXd hCone_allFix = hCone_allFix_raw * vscale;
  if (flag_given_goal_velocity) G = G*gvscale;

  // for crashing check
  MatrixXd cone_allFix;
  Poly::coneIntersection(eCone_allFix, hCone_allFix, &cone_allFix);

  printf("##         Compute Stability Margins         ##\n");
  Eigen::VectorXi e_sss_mode, h_sss_mode;
  Eigen::VectorXi e_s_mode, h_s_mode;
  Eigen::MatrixXd e_cone, h_cone;
  for (int e_sss_i = 0; e_sss_i < NE_SSS_MODES; ++e_sss_i) {
    e_sss_mode = e_sss_modes.middleRows(e_sss_i, 1);
    for (int h_sss_i = 0; h_sss_i < NH_SSS_MODES; ++h_sss_i) {
      h_sss_mode = h_sss_modes.middleRows(h_sss_i, 1);
      for (int e_s_i = 0; e_s_i < e_s_modes[e_sss_i].size(); ++e_s_i) {
        e_s_mode = e_s_modes[e_sss_i].middleRows(e_s_i, 1);
        for (int h_s_i = 0; h_s_i < h_s_modes[h_sss_i].size(); ++h_s_i) {
          h_s_mode = h_s_modes[h_sss_i].middleRows(h_s_i, 1);
          e_cone = getConeOfTheMode(eCone_allFix, e_sss_mode, e_s_mode, kNumSlidingPlanes);
          h_cone = getConeOfTheMode(hCone_allFix, h_sss_mode, h_s_mode, kNumSlidingPlanes);
          /* construct two polytopes */
          // scale kContactForce
          e_cone = kContactForce * e_cone;
          h_cone = - kContactForce * h_cone;
        }
      }
    }
  }

  // EH cone intersection
  //   * Compute safety margins
  //   * get rid of infeasible cone
  //

  eh_modes = [];
  margins = [];

  eh_cones = cell(size(e_modes, 2)*size(h_modes, 2), 1);
  cone_e = cell(size(e_modes, 2)*size(h_modes, 2), 1);
  cone_h = cell(size(e_modes, 2)*size(h_modes, 2), 1);
  Jacs = cell(size(e_modes, 2)*size(h_modes, 2), 1);
  Jacus = cell(size(e_modes, 2)*size(h_modes, 2), 1);
  Jacues = cell(size(e_modes, 2)*size(h_modes, 2), 1);

  eh_cone_feasible_mode_count = 0;
  goal_id = 0;
  for i = 1:size(e_modes, 2)
      for j = 1:size(h_modes, 2)
          [Je_, Jh_] = getFrictionalJacobianFromContacts(e_modes(:, i), h_modes(:, j), eCone_allFix, hCone_allFix);

          % check force balance
          R = coneIntersection(Je_', Jh_');
          if isempty(R) || norm(R) == 0
              continue;
          end

          % compute cone stability margin
          if sameConeCheck(Je_', R)
              margin_ = coneStabilityMargin(Jh_', R);
          elseif sameConeCheck(Jh_', R)
              margin_ = coneStabilityMargin(Je_', R);
          else
              margin_e = coneStabilityMargin(Je_', R);
              margin_h = coneStabilityMargin(Jh_', R);
              margin_ = min(margin_e, margin_h);
          end
          if margin_ < TOL
              continue;
          end

          % compute velocity Jacobian
          [N, Nu] = getJacobianFromContacts(e_modes(:, i), h_modes(:, j), Jac_e, Jac_h);

          % figure(1);clf(1);hold on;
          % printModes([e_modes(:, i); h_modes(:, j)]);
          % fprintf('Margin: %f\n', margin_);
          % drawCone(Je_','g', true);
          % drawCone(Jh_','b', true);

          eh_modes = [eh_modes [e_modes(:, i); h_modes(:, j)]];
          margins = [margins margin_];

          eh_cone_feasible_mode_count = eh_cone_feasible_mode_count + 1;
          eh_cones{eh_cone_feasible_mode_count} = R;
          Jacs{eh_cone_feasible_mode_count} = N;
          Jacus{eh_cone_feasible_mode_count} = Nu;

          if goal_mode_is_given
              if (i == e_mode_goal) && (j == h_mode_goal)
                  goal_id = eh_cone_feasible_mode_count;
              end
          end
      end
  end

  disp(['Modes with Margin > 0: ' num2str(eh_cone_feasible_mode_count)]);

  if goal_mode_is_given && (goal_id == 0)
      disp("Failure: the Cone of the goal mode is empty.");
      return;
  end

  % trim
  eh_cones = eh_cones(1:eh_cone_feasible_mode_count);
  Jacs = Jacs(1:eh_cone_feasible_mode_count);
  Jacus = Jacus(1:eh_cone_feasible_mode_count);

  %%
  %% Begin to check each cone
  %%
  fprintf("###############################################\n");
  fprintf("##            Velocity Filtering             ##\n");
  fprintf("###############################################\n");

  solutions = cell(eh_cone_feasible_mode_count, 1);
  solutions_count = 0;

  for m = 1:eh_cone_feasible_mode_count
      if goal_mode_is_given
          m = goal_id;
      else
          disp('***********');
          fprintf("Mode %d of %d\n", m, eh_cone_feasible_mode_count);
          printModes(eh_modes(:, m));
      end

      eh_mode_goal = eh_modes(:, m);
      fprintf("====================================\n");
      fprintf("= Hybrid Servoing & Crashing Check =\n");
      fprintf("====================================\n");
      N_all = Jacs{m};
      [n_av, n_af, R_a, R_a_inv, w_av, Cv, b_C] = hybridServoing(N_all, G, b_G);
      if isempty(n_av)
          disp("Failure: Hybrid Servoing returns no solution.")
          if goal_mode_is_given
              return;
          else
              continue;
          end
      end

      V_control_directions = -R_a_inv(:, end-n_av+1:end);
      F_control_directions = [R_a_inv(:, 1:n_af) -R_a_inv(:, 1:n_af)];
      intersection = coneIntersection(cone_allFix, V_control_directions);
      % Crashing check
      if ~isempty(intersection) && norm(intersection) > TOL
          disp('Failure: Crashing. The mode is not feasible.');
          if goal_mode_is_given
              return;
          else
              continue;
          end
      end

      fprintf("=======================\n");
      fprintf("=== Mode Filtering ===\n");
      fprintf("=======================\n");

      % How to filter out modes:
      % 1. If NC degenerates, mark this mode as incompatible;
      % 2. If nominal velocity under NC exists, and it cause the contact point to
      %    slide in a different direction, remove this mode.
      feasibilities = false(eh_cone_feasible_mode_count, 1);

      % cone of possible actuation forces
      %   force controlled direction can have both positive and negative forces

      feasible_mode_count = 0;
      flag_goal_infeasible = false;
      flag_goal_impossible = false;

      for n = 1:eh_cone_feasible_mode_count
          % filter out modes using velocity command
          N = Jacs{n};
          Nu = Jacus{n};

          N = rref(N);
          rank_N = rank(N);
          N = N(1:rank_N, :);

          % compute possible sliding directions
          % these computations are based on 'Criteria for Maintaining Desired Contacts for Quasi-Static Systems'
          Lambda_bar = [Cv; N];
          b_Lambda_bar = [b_C; zeros(size(N,1), 1)];

          compatible = false;
          if rank([Lambda_bar b_Lambda_bar], TOL) == rank(Lambda_bar, TOL)
              v_star = linsolve(Lambda_bar, b_Lambda_bar);
              if any(Nu*v_star < -TOL)
                  % this mode can not exist
                  % V-Impossible
                  if n == m
                      flag_goal_impossible = true;
                      break;
                  else
                      continue;
                  end
              end
              compatible = true;
          else
              if n == m
                  flag_goal_infeasible = true;
                  break;
              end
          end

          % figure(1);clf(1);hold on;
          % drawCone(eh_cones{n},'g', true);
          % drawCone(W_action,'k', true);
          if compatible
              feasible_mode_count = feasible_mode_count + 1;
              feasibilities(n) = true;
          end
      end

      if flag_goal_infeasible
          disp('Goal mode violates velocity equality constraints. Discard this mode.');
          if goal_mode_is_given
              return;
          else
              continue;
          end
      end
      if flag_goal_impossible
          disp('Goal mode violates velocity inequality constraints. Discard this mode.');
          if goal_mode_is_given
              return;
          else
              continue;
          end
      end
      disp(['Remaining feasible modes: ' num2str(feasible_mode_count)]);

      drawWrenchSpace(cone_allFix, eCone_allFix, hCone_allFix, ...
              V_control_directions, F_control_directions, eh_cones, eh_modes, ...
              eh_cone_feasible_mode_count, feasibilities, goal_id);

      eh_cones_goal_m = eh_cones{m};
      feasibilities(m) = 0;
      eh_cones_other_feasible_m = eh_cones(feasibilities);

      fprintf("===============================\n");
      fprintf("===  Compute Force Action   ===\n");
      fprintf("===============================\n");

      % action selection
      force_basis = R_a_inv(:, 1:n_af);
      [force_action, shape_margin] = forceControl(force_basis, ...
              eh_cones_goal_m, eh_cones_other_feasible_m);

      if ~isempty(force_action)
          disp('Force command:');
          disp(force_action);
          solution.eh_mode = eh_mode_goal;
          solution.n_af = n_af;
          solution.n_av = n_av;
          solution.w_av = w_av;
          solution.eta_af = -kForceMagnitude*force_action; % the minus sign comes from force balance
          solution.margin = min(shape_margin, margins(m));

          R_a_inv = R_a^-1;
          Cf_inv = vscale_inv*R_a_inv; Cf_inv = Cf_inv(:, 1:n_af);
          Cv_inv = vscale*R_a_inv; Cv_inv = Cv_inv(:, end-n_av+1:end);
          solution.R_a_inv = [Cf_inv Cv_inv];
          solution.R_a = solution.R_a_inv^-1;

          solutions_count = solutions_count + 1;
          solutions{solutions_count} = solution;

          if goal_mode_is_given
              break;
          end
      else
          disp('Failure: no distinguishable force command.');
          if goal_mode_is_given
              return;
          else
              continue;
          end
      end
  end

  fprintf("###############################################\n");
  fprintf("##                  Results                  ##\n");
  fprintf("###############################################\n");
  fprintf("Total number of feasible modes found: %d\n", solutions_count);
  solution = [];
  if solutions_count > 0
      if goal_mode_is_given
          solution = solutions(1);
      else
          solutions = solutions(1:solutions_count);
          margins = zeros(1, solutions_count);
          for i = 1:solutions_count
              mode_text = printModes(solutions{i}.eh_mode, false);
              fprintf("Mode %d: %s\nmargin: %f\n", i, mode_text, solutions{i}.margin);
              margins(i) = solutions{i}.margin;
          end
          disp('Best solution:');
          [~, best_solution_id] = max(margins);
          solution = solutions{best_solution_id};
      end

      printModes(solution.eh_mode);
      V_T = solution.R_a_inv*[zeros(solution.n_af,1); solution.w_av];
      F_T = solution.R_a_inv*[solution.eta_af; zeros(solution.n_av, 1)];
      disp('R_a:');
      disp(solution.R_a);
      disp('V_T:');
      disp(V_T);
      disp('F_T:');
      disp(F_T);
  end

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
    s_modes.push_back(modes_i_colmajor.transpose());
  }
  return true;
}

Eigen::MatrixXd getConeOfTheMode(const Eigen::MatrixXd &cone_allFix,
    const Eigen::VectorXi &sss_mode, const Eigen::VectorXi &s_mode, int kNumSlidingPlanes) {

  int number_of_edges = 0;
  int num_contacts = sss_mode.cols();

  s_mode.array() = (s_mode.array() + 1)/2; // -1/1 -> 0/1
  std::vector<double> generators; // concatenation of all generators
  for (int i = 0; i < count; ++i) {
    // -1: sticking   0: sliding   1: separation
    int id_start = (2*kNumSlidingPlanes+1)*(i-1);
    MatrixXd generators_add(0, 0);
    if (sss_mode[i] == -1) {
      // record all 2*kNumSlidingPlanes generators
      generators_add = cone_allFix.middleRows(id_start, 2*kNumSlidingPlanes);
    } else if (sss_mode[i] == 0) {
      // find two corresponding generators
      std::vector<double> v;
      v.resize(2*6);

      // Compute the order from bit-wise code:
      int id_s = kNumSlidingPlanes*(i-1);
      int mode_sum = s_mode.segment(id_s, kNumSlidingPlanes).sum();
      int id_g;
      if (s_mode(id_s) == 1)
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


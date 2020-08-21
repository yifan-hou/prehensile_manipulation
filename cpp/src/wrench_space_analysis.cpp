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

double wrenchSpaceAnalysis_2d(MatrixXd Jac_e, MatrixXd Jac_h,
    MatrixXd eCone_allFix_r, MatrixXd hCone_allFix_r,
    const VectorXd &F_G,
    const double kContactForce, const double kFrictionE, const double kFrictionH,
    const double kCharacteristicLength,
    MatrixXd G, const VectorXd &b_G,
    const MatrixXi &e_modes, const MatrixXi &h_modes,
    const VectorXi &e_mode_goal, const VectorXi &h_mode_goal,
    int print_level) {
  if (print_level > 0)
    std::cout << "[wrenchSpaceAnalysis_2d] Calling..\n";
  // std::cout << "Jac_e: " << Jac_e << std::endl;
  // std::cout << "Jac_h: " << Jac_h << std::endl;
  // std::cout << "eCone_allFix_r: " << eCone_allFix_r << std::endl;
  // std::cout << "hCone_allFix_r: " << hCone_allFix_r << std::endl;
  // std::cout << "F_G: " << F_G << std::endl;
  // std::cout << "kContactForce: " << kContactForce << std::endl;
  // std::cout << "kFrictionE: " << kFrictionE << std::endl;
  // std::cout << "kFrictionH: " << kFrictionH << std::endl;
  // std::cout << "kCharacteristicLength: " << kCharacteristicLength << std::endl;
  // std::cout << "G: " << G << std::endl;
  // std::cout << "b_G: " << b_G << std::endl;
  // std::cout << "e_modes: " << e_modes << std::endl;
  // std::cout << "h_modes: " << h_modes << std::endl;
  // std::cout << "e_mode_goal: " << e_mode_goal << std::endl;
  // std::cout << "h_mode_goal: " << h_mode_goal << std::endl;
  // getchar();
  Timer timer;
  timer.tic();

  Eigen::FullPivLU<MatrixXd> lu; // for quickly solving linear system
  lu.setThreshold(TOL);

  double time_stats_initialization;
  double time_stats_force_margin;
  double time_stats_hybrid_servoing;
  double time_stats_crashing_check;
  double time_stats_velocity_loops;

  bool flag_given_goal_velocity = false;
  if (b_G.size() > 0) flag_given_goal_velocity = true;

  int kNumEContacts = e_modes.cols();
  int kNumHContacts = h_modes.cols();
  int kDim = Jac_e.cols();
  // scaling for generalized velocity
  Vector3d kv_vec, kv_inv_vec, kf_vec;
  double scale = 1./kCharacteristicLength;
  kv_vec << scale, scale, 1;
  kv_inv_vec << 1./scale, 1./scale, 1;
  kf_vec << 1, 1, scale;
  Matrix3d Kv = kv_vec.asDiagonal();
  Matrix3d Kf = kf_vec.asDiagonal();
  Matrix3d Kv_inv = kv_inv_vec.asDiagonal();

  Jac_e = Jac_e * Kv_inv;
  Jac_h = Jac_h * Kv_inv;
  eCone_allFix_r = eCone_allFix_r * Kf;
  hCone_allFix_r = hCone_allFix_r * Kf;
  if (flag_given_goal_velocity) {
    G.leftCols(kDim) = G.leftCols(kDim) * Kv_inv;
    G.rightCols(kDim) = G.rightCols(kDim) * Kv_inv;
  }

  bool flag_given_goal_mode = false;
  if (e_mode_goal.rows() > 0) flag_given_goal_mode = true;

  // for crashing check
  MatrixXd cone_allFix_r;
  Poly::coneIntersection(eCone_allFix_r, hCone_allFix_r, &cone_allFix_r);

  // divide the big matrices
  // Jn: nContacts x 6,  Nn: nContacts x 12,  Jacobian for the normals
  // Jt: 2*nContacts x 6,  Nn: 2*nContacts x 12, Jacobian for the tangentials
  Eigen::MatrixXd Jn(kNumEContacts + kNumHContacts, kDim);
  Jn << Jac_e.topRows(kNumEContacts), Jac_h.topRows(kNumHContacts);
  Eigen::MatrixXd Jt(kNumEContacts + kNumHContacts, kDim);
  Jt << Jac_e.bottomRows(kNumEContacts), Jac_h.bottomRows(kNumHContacts);
  Eigen::MatrixXd Nn(Jn.rows(), 2*kDim);
  Nn << Jac_e.topRows(kNumEContacts), Eigen::MatrixXd::Zero(kNumEContacts, kDim),
       -Jac_h.topRows(kNumHContacts), Jac_h.topRows(kNumHContacts);
  Eigen::MatrixXd Nt(Jt.rows(), 2*kDim);
  Nt << Jac_e.bottomRows(kNumEContacts), Eigen::MatrixXd::Zero(kNumEContacts, kDim),
        -Jac_h.bottomRows(kNumHContacts), Jac_h.bottomRows(kNumHContacts);

  time_stats_initialization = timer.toc();
  timer.tic();
  /**
   * Filter out modes with empty Wrench Cones
   * Compute geometrical stability margin
   */
  if (print_level > 0)
    std::cout << "[WrenchStamping] 1. Wrench cone filtering" << std::endl;
  std::vector<VectorXi> eh_modes;
  std::vector<double> margins;
  std::vector<MatrixXd> polytope_of_the_modes;
  std::vector<MatrixXd> N_of_the_modes;
  std::vector<MatrixXd> Nu_of_the_modes;
  int goal_id = -1;
  double geometrical_stability_margin = 0;

  VectorXi e_mode_ii, h_mode_jj;
  MatrixXd e_cone_ii, e_polytope_ii, h_cone_jj, h_polytope_jj;
  MatrixXd polytope_ij_sum;
  MatrixXd polytope_ij_minus;
  MatrixXd polytope_ij_sum_A;
  VectorXd polytope_ij_sum_b;
  MatrixXd N, Nu;
  for (int ii = 0; ii < e_modes.rows(); ++ii) {
    e_mode_ii = e_modes.middleRows(ii, 1).transpose();
    e_cone_ii = kContactForce * getConeOfTheMode_2d(eCone_allFix_r, e_mode_ii);
    Poly::minkowskiSumOfVectors(e_cone_ii, &e_polytope_ii);

    // gravity offset
    Poly::offsetPolytope(&e_polytope_ii, F_G);

    bool is_goal_e = true;
    if (flag_given_goal_mode) {
      for (int i = 0; i < kNumEContacts; ++i) {
        if (e_mode_ii[i] != e_mode_goal[i]) {
          is_goal_e = false;
          break;
        }
      }
    }

    for (int jj = 0; jj < h_modes.rows(); ++jj) {
      h_mode_jj = h_modes.middleRows(jj, 1).transpose();
      h_cone_jj = - kContactForce * getConeOfTheMode_2d(hCone_allFix_r, h_mode_jj);
      Poly::minkowskiSumOfVectors(h_cone_jj, &h_polytope_jj);

      if (print_level > 0)
        std::cout << "  ii:" << ii << ", jj:" << jj << ", e mode:" << e_mode_ii.transpose() << ", h mode:" << h_mode_jj.transpose();

      Poly::minkowskiSum(e_polytope_ii, h_polytope_jj, &polytope_ij_sum);
      Poly::polytopeFacetEnumeration(polytope_ij_sum, &polytope_ij_sum_A, &polytope_ij_sum_b);


      if (polytope_ij_sum_b.minCoeff() <= 1e-5 ) {
        if (print_level > 0) std::cout << " F-Infeasible." << std::endl;
        continue;
      }
      // The cone is non-empty.
      // compute the stability margin
      // margin = min(b./normByRow(A));
      double margin = 9999;
      for (int i = 0; i < polytope_ij_sum_b.rows(); ++i) {
        double A_row_norm = polytope_ij_sum_A.middleRows(i, 1).norm();
        assert(A_row_norm > 1e-7);
        double margin_new = polytope_ij_sum_b(i)/A_row_norm;
        if (margin_new < margin) margin = margin_new;
      }
      if (print_level > 0) std::cout << " margin: " << margin;

      // compute the polytope of the mode (intersection, not minkowski sum)
      if (!Poly::polytopeIntersection(e_polytope_ii, -h_polytope_jj, &polytope_ij_minus)) {
        std::cerr << "Intersection is not found!!\n";
        return -1;
      }

      // check if this is a goal mode
      if (flag_given_goal_mode && is_goal_e) {
        bool is_goal_h = true;
        for (int i = 0; i < kNumHContacts; ++i) {
          if (h_mode_jj[i] != h_mode_goal[i]) {
            is_goal_h = false;
            break;
          }
        }
        if (is_goal_h) {
          if (print_level > 0) std::cout << " (Goal)";
          goal_id = margins.size();
          geometrical_stability_margin = margin;
        }
      }
      if (print_level > 0) std::cout << std::endl;

      // store the results
      margins.push_back(margin);
      // polytope_of_the_modes.push_back(polytope_ij_sum);
      polytope_of_the_modes.push_back(polytope_ij_minus);

      getConstraintOfTheMode_2d(Jac_e, Jac_h, e_mode_ii, h_mode_jj, &N, &Nu);

      N_of_the_modes.push_back(N);
      Nu_of_the_modes.push_back(Nu);

      // for debug purpose
      VectorXi eh_mode(kNumEContacts + kNumHContacts);
      eh_mode << e_mode_ii, h_mode_jj;
      eh_modes.push_back(eh_mode);
    }
  }
  if (goal_id < 0) {
    std::cout << "[WrenchStamping]    Goal mode does not have force-balance." << std::endl;
    return -1;
  }

  time_stats_force_margin = timer.toc();
  timer.tic();

  if (print_level > 0)
    std::cout << "[WrenchStamping] 2. HFVC" << std::endl;
  int kDimActualized = 3;
  int kDimUnActualized = 3;
  HFVC action;
  if (!solvehfvc_newer(N_of_the_modes[goal_id], G, b_G, kDimActualized, kDimUnActualized, &action)) {
    std::cout << "[WrenchStamping]    HFVC has no solution." << std::endl;
    return -1;
  }
  assert(action.n_af < kDimActualized); // shouldn't be all force
  assert(action.n_af > 0); // shouldn't be all velocity
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
  MatrixXd F_control_directions_r = action.R_a.topRows(action.n_af);

  time_stats_hybrid_servoing = timer.toc();
  timer.tic();

  // Crashing check
  MatrixXd R;
  Poly::coneIntersection(cone_allFix_r, V_control_directions_r, &R); // this line has errors sometimes
  if (R.rows() > 0) {
    std::cout << "[WrenchStamping]    Crashing." << std::endl;
    return -1;
  }

  time_stats_crashing_check = timer.toc();
  timer.tic();

  if (print_level > 0)  std::cout << "[WrenchStamping] 3. Begin Velocity loop." << std::endl;
  std::vector<int> feasible_ids;
  for (int id = 0; id < margins.size(); ++id) {
    if (print_level > 0) {
      std::cout << "[WrenchStamping]    Checking id: " << id;
      std::cout << " : " << eh_modes[id].transpose();
      if (id == goal_id) std::cout << " (Goal)";
    }

    MatrixXd N = N_of_the_modes[id];
    MatrixXd Nu = Nu_of_the_modes[id]; // Nu*v >= 0

    // check velocity constraints
    Eigen::MatrixXd NC(N.rows() + action.C.rows(), N.cols());
    NC << N, action.C;
    Eigen::VectorXd b_NC = Eigen::VectorXd::Zero(NC.rows());
    b_NC.tail(action.b_C.rows()) = action.b_C;
    Eigen::VectorXd b_Nu = Eigen::VectorXd(Nu.rows());
    if(!Poly::vertexEnumeration(-Nu, b_Nu, NC, b_NC, &R)) {
      std::cerr << "Error: vertexEnumeration returns error." << std::endl;
      return -1;
    }
    if (R.rows() == 0) {
      // no solution
      if (id == goal_id) {
        std::cout << " Goal is infeasible. Return." << std::endl;
        return -1;
      }
      if (print_level > 0) std::cout << " infeasible." << std::endl;
      continue; // otherwise, just discard this mode
    }
    // this mode is feasible. store its polytope
    feasible_ids.push_back(id);
    if (print_level > 0) std::cout << " stored." << std::endl;
  }

  if (print_level > 0)
    std::cout << "[WrenchStamping] 4. Compute force control and control-stability-margin." << std::endl;
  std::vector<MatrixXd> polytopes_projection_r;
  std::vector<MatrixXd> pps_A; // A x <= b
  std::vector<VectorXd> pps_b; // A x <= b
  MatrixXd pp_goal_A;
  VectorXd pp_goal_b;
  VectorXd pp_goal_point; // an inner point
  MatrixXd polytope_projection, pp_A;
  VectorXd pp_b;
  if (print_level > 0)
    std::cout << "[WrenchStamping]  4.1 Compute the projection of the polytopes." << std::endl;
  for (int c = 0; c < feasible_ids.size(); ++c) {
    if (print_level > 0)
      std::cout << "[WrenchStamping]    Polytope " << c << ": " << eh_modes[feasible_ids[c]].transpose() << ": ";
    // projection
    polytope_projection = polytope_of_the_modes[feasible_ids[c]] * F_control_directions_r.transpose();
    // get the inequality representations for the cone projections
    if(!Poly::polytopeFacetEnumeration(polytope_projection, &pp_A, &pp_b)) {
      std::cerr << "Error: polytopeFacetEnumeration return error." << std::endl;
      std::cout << "[debug] polytope_projection: " << polytope_projection << std::endl;
      exit(1);
    }
    if (print_level > 0) std::cout << polytope_projection.rows() << " vertices.";

    // std::cout << "F_control_directions_r:\n" << F_control_directions_r << std::endl;
    // std::cout << "polytope_of_the_modes[feasible_ids[c]]:\n" << polytope_of_the_modes[feasible_ids[c]] << std::endl;
    // std::cout << " A rows: " << pp_A.rows() << std::endl;
    // getchar();
    // return -1;

    if (goal_id == feasible_ids[c]) {
      if (print_level > 0) std::cout << " (id: Goal) ";
      pp_goal_A = pp_A;
      pp_goal_b = pp_b;
      // std::cout <<  "pp_goal_A:\n" << pp_goal_A << std::endl;
      // std::cout <<  "pp_goal_b:\n" << pp_goal_b << std::endl;
      pp_goal_point = (MatrixXd::Ones(1, polytope_projection.rows()) * polytope_projection).transpose() / polytope_projection.rows();
      assert(("Assertion fail: goal polytope projection is empty", polytope_projection.norm() > 1e-5));
      assert(("Assertion fail: goal polytope projection degenerates", polytope_projection.rows() >= action.n_af)); // ideally we should check its rank
    } else {
      polytopes_projection_r.push_back(polytope_projection);
      pps_A.push_back(pp_A);
      pps_b.push_back(pp_b);
      if (print_level > 0) std::cout << " (id:" << pps_A.size() - 1 << ") ";
    }
    if (print_level > 0) std::cout << " A rows: " << pp_A.rows() << std::endl;
    // std::cout << "polytope_projection:\n" << polytope_projection << std::endl;
    // std::cout << "polytope_of_the_modes[feasible_ids[c]]:\n" << polytope_of_the_modes[feasible_ids[c]] << std::endl;
  }
  if (print_level > 0)
    std::cout << "[WrenchStamping]  4.2 Sample wrenches and find feasible ones." << std::endl;
  // sample wrenches in the projection of the goal polytope
  int ns = 0;
  if (action.n_af == 1) ns = 5;
  else if (action.n_af == 2) ns = 20;
  else ns = 50;

  // int discard = 10;
  // int runup = 10;
  // MatrixXd wrench_samples = Poly::hitAndRunSampleInPolytope(pp_goal_A, pp_goal_b,
  //     pp_goal_point, ns, discard, runup);
  std::vector<VectorXd> wrench_samples = Poly::sampleInP1OutOfP2(
      pp_goal_A, pp_goal_b,
      pps_A, pps_b, pp_goal_point, ns);

  // std::cout << "wrench_samples: \n" << wrench_samples << std::endl;
  // std::cout << "wrench_sample_norms: \n" << wrench_samples.rowwise().norm() << std::endl;
  VectorXd wrench_sample, wrench_best;
  double control_stability_margin = -1;
  for (int s = 0; s < ns; ++s) { // make sure they sum to one
    // create the sample
    wrench_sample = wrench_samples[s];
    if (print_level > 0)
      std::cout << "[WrenchStamping]    Sample #" << s << ": ";
    // check if the sample is within any other polytopes
    bool infeasible_sample = false;
    for (int i = 0; i < pps_A.size(); ++i) {
      VectorXd Ax = pps_A[i]*wrench_sample - pps_b[i];
      if (Ax.maxCoeff() <= 0) {
        // the sample is within another polytope
        if (print_level > 0) std::cout << i << "th polytope infeasible." << std::endl;
        infeasible_sample = true;
        break;
      }
    }
    if (infeasible_sample) continue;

    // compute its distance to all other projections
    double min_dist = 999999.9;
    for (int i = 0; i < pps_A.size(); ++i) {
      double dist = Poly::distP2Polyhedron(wrench_sample, pps_A[i], pps_b[i], Eigen::VectorXd::Random(kDim));
      if (dist < min_dist) min_dist = dist;
    }
    if (min_dist > control_stability_margin) {
      control_stability_margin = min_dist;
      wrench_best = wrench_sample;
    }
    if (print_level > 0) std::cout << "min_dist:" << min_dist << std::endl;
  }
  if (control_stability_margin < 0) {
    std::cout << "[WrenchStamping] 4. Force control has no solution." << std::endl;
    return -1;
  }
  if (print_level > 0)
    std::cout << "[WrenchStamping] 4. Force control: " << wrench_best.transpose() << std::endl;
  // now we have the control stability margin
  action.eta_af = -wrench_best; // the minus sign comes from force balance
  // before scaling
  // R_scaled^T*eta_scaled = Kf*f
  // Eigen::VectorXd f0 = action.R_a.topRows(action.n_af).transpose() * action.eta_af;
  // Eigen::VectorXd F_T = Kf.inverse() * f0;

  // scale back so the control constraints work on normal units
  Eigen::VectorXd kv2_vec(6);
  kv2_vec << kv_vec, kv_vec;
  Eigen::MatrixXd Kv2 = kv2_vec.asDiagonal();
  action.R_a.topRows(action.n_af) *= Kf;
  action.R_a.bottomRows(action.n_av) *= Kv;
  action.C *= Kv2;

  // print the results
  Eigen::MatrixXd R_a_inv = action.R_a.inverse();
  Eigen::VectorXd V = Eigen::VectorXd::Zero(kDimActualized);
  V.tail(action.n_av) = action.w_av;
  Eigen::VectorXd V_T = R_a_inv*V;

  Eigen::VectorXd F = Eigen::VectorXd::Zero(kDimActualized);
  F.head(action.n_af) = action.eta_af;
  Eigen::VectorXd F_T = R_a_inv*F;

  if (print_level > 0) {
    // Eigen::MatrixXd N = N_of_the_modes[goal_id] * Kv2;
    // Eigen::MatrixXd NC(N.rows() + action.C.rows(), N.cols());
    // NC << N, action.C;
    // Eigen::VectorXd b_NC = Eigen::VectorXd::Zero(NC.rows());
    // b_NC.tail(action.b_C.rows()) = action.b_C;
    // assert(NC.norm() > 10*TOL); //otherwise lu won't be accurate
    // lu.compute(NC);
    // Eigen::VectorXd sol_NC = lu.solve(b_NC);
    // bool a_solution_exists = (NC*sol_NC).isApprox(b_NC, 10.*TOL);
    // if (!a_solution_exists) {
    //     // no solution
    //     std::cout << "No Solution. " << std::endl;
    //     continue;
    // }
    // Eigen::MatrixXd null_NC = lu.kernel();
    // std::cout << "null_NC: " << null_NC << std::endl;
    // std::cout << "Gv-bg: " << G*sol_NC - b_G << std::endl;

    std::cout << " 5. Results:" << std::endl;
    std::cout << "   geometrical_stability_margin: " << geometrical_stability_margin << std::endl;
    std::cout << "   control_stability_margin: " << control_stability_margin << std::endl;
    std::cout << "   R_a:\n" << action.R_a << std::endl;
    std::cout << "   w_av:\n" << action.w_av << std::endl;
    std::cout << "   eta_af:\n" << action.eta_af << std::endl;
    std::cout << "   V_T:" << V_T.transpose() << std::endl;
    std::cout << "   F_T:" << F_T.transpose() << std::endl;

    time_stats_velocity_loops = timer.toc();
    std::cout << "Timing statistics:\n";
    std::cout << "  time_stats_initialization: " << time_stats_initialization << " ms\n";
    std::cout << "  time_stats_force_margin: " << time_stats_force_margin << " ms\n";
    std::cout << "  time_stats_hybrid_servoing: " << time_stats_hybrid_servoing << " ms\n";
    std::cout << "  time_stats_crashing_check: " << time_stats_crashing_check << " ms\n";
    std::cout << "  time_stats_velocity_loops: " << time_stats_velocity_loops << " ms\n";
    std::cout << "  Total: " << time_stats_initialization + time_stats_force_margin
        + time_stats_hybrid_servoing + time_stats_crashing_check
        + time_stats_velocity_loops << " ms\n";
  }
  return std::min(control_stability_margin, geometrical_stability_margin);
}


void wrenchSpaceAnalysis(MatrixXd Jac_e, MatrixXd Jac_h,
    MatrixXd eCone_allFix_r, MatrixXd hCone_allFix_r,
    const VectorXd &F_G, const double kContactForce,
    const double kFrictionE, const double kFrictionH,
    const double kCharacteristicLength, const int kNumSlidingPlanes,
    const MatrixXi &e_cs_modes, const std::vector<MatrixXi> &e_ss_modes,
    const MatrixXi &h_cs_modes, const std::vector<MatrixXi> &h_ss_modes,
    MatrixXd G, const VectorXd &b_G,
    const MatrixXi &e_cs_modes_goal, const std::vector<MatrixXi> &e_ss_modes_goal,
    const MatrixXi &h_cs_modes_goal, const std::vector<MatrixXi> &h_ss_modes_goal,
    int print_level) {

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
  double time_stats_initialization;
  double time_stats_hybrid_servoing;
  double time_stats_velocity_filtering;
  double time_stats_projection;
  double time_stats_force_control;

  bool flag_given_goal_velocity = false;
  if (b_G.size() > 0) flag_given_goal_velocity = true;

  int kNumEContacts = e_cs_modes.cols();
  int kNumHContacts = h_cs_modes.cols();
  int kDim = Jac_e.cols();
  // scaling for generalized velocity
  Vector6d kv_vec, kv_inv_vec, kf_vec;
  double scale = 1./kCharacteristicLength;
  kv_vec << scale, scale, scale, 1, 1, 1;
  kv_inv_vec << 1./scale, 1./scale, 1./scale, 1, 1, 1;
  kf_vec << 1, 1, 1, scale, scale, scale;
  MatrixXd Kv = kv_vec.asDiagonal();
  MatrixXd Kf = kf_vec.asDiagonal();
  MatrixXd Kv_inv = kv_inv_vec.asDiagonal();

  Jac_e = Jac_e * Kv_inv;
  Jac_h = Jac_h * Kv_inv;
  eCone_allFix_r = eCone_allFix_r * Kf;
  hCone_allFix_r = hCone_allFix_r * Kf;
  if (flag_given_goal_velocity) {
    G.leftCols(kDim) = G.leftCols(kDim) * Kv_inv;
    G.rightCols(kDim) = G.rightCols(kDim) * Kv_inv;
  }

  // for crashing check
  // Instead of
  //    Poly::coneIntersection(eCone_allFix_r, hCone_allFix_r, &cone_allFix_r);
  // save the H-representation of the cone of the all-fixed mode
  MatrixXd Ae_allFix;
  MatrixXd Ah_allFix;
  Poly::coneFacetEnumeration(eCone_allFix_r, &Ae_allFix);
  Poly::coneFacetEnumeration(hCone_allFix_r, &Ah_allFix);
  MatrixXd A_allFix(Ae_allFix.rows() + Ah_allFix.rows(), Ae_allFix.cols());
  A_allFix << Ae_allFix, Ah_allFix;
  // do minimized_constraints here?



  Eigen::MatrixXi e_sss_modes, h_sss_modes;
  std::vector<Eigen::MatrixXi> e_s_modes, h_s_modes;

  // divide the big matrices
  // Jn: nContacts x 6,  Nn: nContacts x 12,  Jacobian for the normals
  // Jt: 2*nContacts x 6,  Nn: 2*nContacts x 12, Jacobian for the tangentials
  Eigen::MatrixXd Jn(kNumEContacts + kNumHContacts, kDim);
  Jn << Jac_e.topRows(kNumEContacts), Jac_h.topRows(kNumHContacts);
  Eigen::MatrixXd Jt(2*kNumEContacts + 2*kNumHContacts, kDim);
  Jt << Jac_e.bottomRows(2*kNumEContacts), Jac_h.bottomRows(2*kNumHContacts);
  Eigen::MatrixXd Nn(Jn.rows(), 2*kDim);
  Nn << Jac_e.topRows(kNumEContacts), Eigen::MatrixXd::Zero(kNumEContacts, kDim),
       -Jac_h.topRows(kNumHContacts), Jac_h.topRows(kNumHContacts);
  Eigen::MatrixXd Nt(Jt.rows(), 2*kDim);
  Nt << Jac_e.bottomRows(2*kNumEContacts), Eigen::MatrixXd::Zero(2*kNumEContacts, kDim),
        -Jac_h.bottomRows(2*kNumHContacts), Jac_h.bottomRows(2*kNumHContacts);

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

  std::vector<Eigen::MatrixXd> sample_grids;
  sample_grids.push_back(Eigen::MatrixXd(2, 1));
  sample_grids.push_back(Eigen::MatrixXd(4, 2));
  sample_grids.push_back(Eigen::MatrixXd(8, 3));
  sample_grids[0] << 1,
                    -1;
  sample_grids[1] << 1, 1,
                   1, -1,
                   -1, 1,
                   -1, -1;
  sample_grids[2] << 1, 1, 1,
                   1, 1, -1,
                   1, -1, 1,
                   1, -1, -1,
                   -1, 1, 1,
                   -1, 1, -1,
                   -1, -1, 1,
                   -1, -1, -1;

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
  std::cout << "e_sss_modes_goal: \n" << e_sss_modes_goal << std::endl;
  std::cout << "timer: modeCleaning time = " << timer.toc() << "ms" << std::endl;

  // Eigen::VectorXd F = Eigen::VectorXd::Zero(12);
  Eigen::MatrixXd N;
  int kDimActualized = 6;
  int kDimUnActualized = 6;
  HFVC action;

  Eigen::VectorXi e_sss_mode_goal, h_sss_mode_goal;
  Eigen::VectorXi e_sss_mode, h_sss_mode;

  Eigen::FullPivLU<MatrixXd> lu; // for quickly solving linear system
  lu.setThreshold(TOL);

  time_stats_initialization = timer.toc();
  timer.tic();

  /*******************************************************************
   *      First half: Velocity Filtering
   */

  std::cout << "##         Begin loop         ##\n";
  timer.tic();
  for (int e_sss_i_goal = 0; e_sss_i_goal < e_sss_modes_goal.rows(); ++e_sss_i_goal) {
    e_sss_mode_goal = e_sss_modes_goal.middleRows(e_sss_i_goal, 1).transpose();
    for (int h_sss_i_goal = 0; h_sss_i_goal < h_sss_modes_goal.rows(); ++h_sss_i_goal) {
      h_sss_mode_goal = h_sss_modes_goal.middleRows(h_sss_i_goal, 1).transpose();
      std::cout << "[WrenchStamping] goal id: e " << e_sss_i_goal << ", h " << h_sss_i_goal << std::endl;
      std::cout << "[WrenchStamping] goal mode: e " << e_sss_mode_goal.transpose() << ", h " << h_sss_mode_goal.transpose() << std::endl;
      int goal_id_e = findIdInModes(e_sss_mode_goal, e_sss_modes);
      int goal_id = -1;
      assert(goal_id_e >= 0);

      std::vector<Eigen::MatrixXd> e_cones_VFeasible;
      std::vector<Eigen::VectorXi> e_modes_VFeasible;
      std::vector<Eigen::VectorXi> h_modes_VFeasible;

      std::cout << "[WrenchStamping] 1. HFVC" << std::endl;
      Eigen::MatrixXd N, Nu;
      getConstraintOfTheMode(Jac_e, Jac_h,
          e_sss_mode_goal, h_sss_mode_goal,
          &N, &Nu);

      if (!solvehfvc_newer(N, G, b_G, kDimActualized, kDimUnActualized, &action)) {
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
      MatrixXd F_control_directions_r = action.R_a.topRows(action.n_af);
      MatrixXd R_a_inv = action.R_a.inverse();

      // Crashing check

      time_stats_hybrid_servoing = timer.toc();
      timer.tic();

      // this might be much faster:
      //    MatrixXd R = A_allFix * V_control_directions_r.transpose();
      MatrixXd A_V_cone;
      Poly::coneFacetEnumeration(V_control_directions_r, &A_V_cone);
      MatrixXd A_AF_V(A_V_cone.rows() + A_allFix.rows(), A_allFix.cols());
      A_AF_V << A_V_cone, A_allFix;
      VectorXd xs = VectorXd::Zero(A_V_cone.cols());
      bool is_feasible = Poly::lpfeasibility(A_AF_V, Eigen::VectorXd::Zero(A_AF_V.rows()), &xs);
      if (is_feasible && xs.norm() > TOL) {
        std::cout << "[WrenchStamping]    Crashing." << std::endl;
        continue;
      } else {
        std::cout << "  No crashing." << std::endl;
      }

      /**
       * Debug HS
       *
       */

      // Eigen::MatrixXd N0C(N.rows() + action.C.rows(), N.cols());
      // N0C << N, action.C;
      // Eigen::VectorXd b_N0C = Eigen::VectorXd::Zero(N0C.rows());
      // b_N0C.tail(action.b_C.rows()) = action.b_C;
      // assert(N0C.norm() > 10*TOL); //otherwise lu won't be accurate
      // lu.compute(N0C);
      // Eigen::VectorXd sol_N0C = lu.solve(b_N0C);
      // bool a_solution_exists = (N0C*sol_N0C).isApprox(b_N0C, 10.*TOL);
      // if (!a_solution_exists) {
      //     // no solution
      //     std::cout << "No Solution. " << std::endl;
      //     continue;
      // }
      // Eigen::MatrixXd null_N0C = lu.kernel();
      // std::cout << "sol_N0C: " << sol_N0C.transpose() << std::endl;
      // std::cout << "null_N0C: " << null_N0C << std::endl;
      // std::cout << "Gv-bg: " << G*sol_N0C - b_G << std::endl;

      // Eigen::VectorXd kv2_vec(12);
      // kv2_vec << kv_vec, kv_vec;
      // Eigen::MatrixXd Kv2 = kv2_vec.asDiagonal();
      // action.R_a.topRows(action.n_af) *= Kf;
      // action.R_a.bottomRows(action.n_av) *= Kv;
      // action.C *= Kv2;

      // Eigen::MatrixXd R_a_inv = action.R_a.inverse();
      // Eigen::VectorXd V = Eigen::VectorXd::Zero(kDimActualized);
      // V.tail(action.n_av) = action.w_av;
      // Eigen::VectorXd V_T = R_a_inv*V;
      // std::cout << "   V_T:" << V_T.transpose() << std::endl;
      // // Eigen::VectorXd F_T = action.R_a.topRows(action.n_af).transpose() * action.eta_af;
      // // std::cout << "   F_T:" << F_T.transpose() << std::endl;
      // return;

      // time_stats_hybrid_servoing = timer.toc();
      // timer.tic();

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

          getConstraintOfTheMode(Jac_e, Jac_h, e_sss_mode, h_sss_mode, &N, &Nu);

          Eigen::MatrixXd NC(N.rows() + action.C.rows(), N.cols());
          NC << N, action.C;
          Eigen::VectorXd b_NC = Eigen::VectorXd::Zero(NC.rows());
          b_NC.tail(action.b_C.rows()) = action.b_C;
          assert(NC.norm() > 10*TOL); //otherwise lu won't be accurate
          lu.compute(NC);
          Eigen::VectorXd sol_NC = lu.solve(b_NC);
          bool a_solution_exists = (NC*sol_NC).isApprox(b_NC, 10.*TOL);
          if (!a_solution_exists) {
              // no solution
              std::cout << "No Solution. " << std::endl;
              continue;
          }
          Eigen::MatrixXd null_NC = lu.kernel();
          bool has_penetration = false;
          if (Nu.rows() > 0) {
            // check inequalities
            Eigen::MatrixXd contact_normal_proj = Nu*null_NC; // this is a linear space, not cone
            Eigen::VectorXd contact_normal_sol = Nu*sol_NC;
            // std::cout << "debug: contact_normal_proj = " << contact_normal_proj << std::endl;
            // std::cout << "debug: contact_normal_sol = " << contact_normal_sol << std::endl;
            for (int i = 0; i < Nu.rows(); ++i) {
              if (contact_normal_proj.middleRows(i, 1).norm() < TOL) {
                if (contact_normal_sol(i) < -TOL) {
                  has_penetration = true;
                  break;
                }
              }
            }
          }
          if (has_penetration) {
            std::cout << "Violates inequalities. " << std::endl;
            continue;
          }

          // prepare to get the cone
          Eigen::MatrixXd e_cone_base = getConeOfTheMode(eCone_allFix_r, e_sss_mode, kNumSlidingPlanes);
          // First, decide where to get velocity samples.
          // If samples from contact, fill up g_sampled directly.
          // If samples from kernel, fill in vel_samples_in_kernel
          Eigen::MatrixXd vel_samples_in_kernel = Eigen::MatrixXd(0, 2*kDim);
          Eigen::MatrixXd g_sampled = Eigen::MatrixXd(0, kDim);
          if (null_NC.norm() < TOL) {
            // unique solution
            std::cout << "Unique Solution. " << std::endl;
            // find the cone of this Unique solution
            vel_samples_in_kernel = sol_NC.transpose();
          } else {
            // Multiple solutions
            // sample sliding velocities
            int dim_null_NC = null_NC.cols();
            std::cout << "Multiple Solutions, kernel Dim = " << dim_null_NC << ", contact dim = ";
            // check dimension of contact tangential projections for sliding contacts
            std::vector<int> contact_sliding_DOFs;
            std::vector<int> contact_ids;
            std::vector<Eigen::VectorXd> contact_vels;
            std::vector<Eigen::MatrixXd> contact_kernels;
            for (int i = 0; i < kNumEContacts; ++i) {
              if (e_sss_mode(i) == 0) {
                Eigen::MatrixXd contact_tangent_proj = Nt.middleRows(2*i, 2)*null_NC;
                if (contact_tangent_proj.norm() < 10.0*TOL) {
                  // rank = 0
                  continue;
                } else {
                  lu.compute(contact_tangent_proj);
                  int rank = lu.rank();
                  std::cout << rank << " ";
                  Eigen::VectorXd contact_tangent_sol = Nt.middleRows(2*i, 2)*sol_NC;
                  contact_sliding_DOFs.push_back(rank);
                  contact_ids.push_back(i);
                  contact_vels.push_back(contact_tangent_sol);
                  contact_kernels.push_back(lu.image(contact_tangent_proj));
                }
              }
            }
            std::cout << std::endl;
            // sample velocities
            //  1. no sliding: g_sampled = empty.
            //  2. has sliding:
            //    a. contact dim all zeros (x)
            //    b. contact dim = 1: sample from contact dim
            //    c. contact dim > 1: sample from kernel
            if (contact_sliding_DOFs.size() == 0) {
              // no sliding, do nothing here
            } else if (contact_sliding_DOFs.size() == 1) {
              // sample from contact, 1d or 2d
              assert(contact_sliding_DOFs[0] > 0);
              assert(contact_sliding_DOFs[0] <= 2);
              // fill g_sampled directly
              double mag = contact_vels[0].norm();
              if (mag < 100*TOL) mag = 1;
              Eigen::MatrixXd sample_grid = sample_grids[contact_sliding_DOFs[0]-1] * mag;
              Eigen::MatrixXd vel_samples_on_contact =
                  Eigen::MatrixXd::Ones(sample_grid.rows(), 1) * contact_vels[0].transpose()
                  + sample_grid * contact_kernels[0].transpose();
              // normalize, get cones
              double friction = (contact_ids[0] < kNumEContacts)? kFrictionE:kFrictionH;
              g_sampled = getSlidingGeneratorsFromOneContact(vel_samples_on_contact,
                  Jt.middleRows(2*contact_ids[0], 2), Jn.middleRows(contact_ids[0], 1), friction);
            } else {
              // sample from kernel
              double mag = sol_NC.norm();
              if (mag < 100*TOL) mag = 1;
              assert(dim_null_NC <= 3);
              Eigen::MatrixXd sample_grid = sample_grids[dim_null_NC-1] * mag;
              vel_samples_in_kernel =
                  Eigen::MatrixXd::Ones(sample_grid.rows(), 1) * sol_NC.transpose()
                  + sample_grid * null_NC.transpose();
            }
          }

          if ((g_sampled.rows() == 0) && (vel_samples_in_kernel.rows() != 0)) {
            // samples are drawn from kernel, in vel_samples_in_kernel
            // use it to fill up g_sampled

            // project kernel to every sliding contacts
            std::vector<Eigen::MatrixXd> g_sampled_each_contact;
            std::vector<int> g_sampled_rows;
            int g_sampled_rows_total = 0;
            for (int i = 0; i < kNumEContacts; ++i) {
              if (e_sss_mode(i) == 0) {
                Eigen::MatrixXd contact_tangent_proj = Nt.middleRows(2*i, 2)*vel_samples_in_kernel.transpose();
                if (contact_tangent_proj.norm() < 10*TOL) continue;
                double friction = kFrictionE;
                g_sampled_each_contact.push_back(getSlidingGeneratorsFromOneContact(contact_tangent_proj.transpose(),
                    Jt.middleRows(2*i, 2), Jn.middleRows(i, 1), friction));
                g_sampled_rows.push_back(g_sampled_each_contact.back().rows());
                g_sampled_rows_total += g_sampled_rows.back();
              }
            }

            g_sampled = Eigen::MatrixXd(g_sampled_rows_total, kDim);
            int g_sampled_rows_acc = 0;
            for (int i = 0; i < g_sampled_each_contact.size(); ++i) {
              g_sampled.middleRows(g_sampled_rows_acc, g_sampled_rows[i]) = g_sampled_each_contact[i];
              g_sampled_rows_acc += g_sampled_rows[i];
            }
          }
          // save the cone(s)
          Eigen::MatrixXd e_cone(e_cone_base.rows() + g_sampled.rows(), kDim);
          e_cone << e_cone_base, g_sampled;
          e_cones_VFeasible.push_back(e_cone);
          e_modes_VFeasible.push_back(e_sss_mode);
          h_modes_VFeasible.push_back(h_sss_mode);
          if (e_sss_i == goal_id_e) goal_id = e_cones_VFeasible.size() - 1;
        }
      }
      // done with inner SSS loop
      assert(goal_id >= 0);
      time_stats_velocity_filtering = timer.toc();
      timer.tic();

      /*******************************************************************
       *      Second half: Force Filtering
       */
      // projection cones of the modes onto force-controlled subspace
      std::cout << "[WrenchStamping] 3. Compute force control and control-stability-margin." << std::endl;
      // std::vector<Eigen::MatrixXd> cones_projection_r;
      std::vector<Eigen::MatrixXd> cp_A; // A x <= 0
      Eigen::MatrixXd cone_projection_goal_r;
      Eigen::MatrixXd cp_goal_A;
      Eigen::MatrixXd cone_of_the_mode, cone_projection, cp_A_temp;
      std::cout << "[WrenchStamping]  3.1 Compute cone of the modes." << std::endl;
      for (int c = 0; c < e_cones_VFeasible.size(); ++c) {
        std::cout << "[WrenchStamping]    Cone " << c << ": " << e_modes_VFeasible[c].transpose() << ": ";
        // compute cone of the modes
        PPL::C_Polyhedron ph1(6, PPL::EMPTY);
        Eigen::MatrixXd R1(e_cones_VFeasible[c].rows(), e_cones_VFeasible[c].cols() + 1);
        R1 << Eigen::VectorXd::Zero(e_cones_VFeasible[c].rows()), e_cones_VFeasible[c];
        Poly::constructPPLPolyFromV(R1, &ph1);

        PPL::C_Polyhedron ph2(6, PPL::EMPTY);
        Eigen::MatrixXd R2(hCone_allFix_r.rows(), hCone_allFix_r.cols() + 1);
        R2 << Eigen::VectorXd::Zero(hCone_allFix_r.rows()), hCone_allFix_r;
        Poly::constructPPLPolyFromV(R2, &ph2);

        ph1.minimized_constraints();
        ph2.minimized_constraints();
        ph1.intersection_assign(ph2);
        ph1.minimized_generators();

        // check intersection results
        MatrixXd R_mode;
        Poly::getVertexFromPPL(ph1, &R_mode);
        if (!((R_mode.rows() > 0) && (R_mode.rightCols(R_mode.cols()-1).norm() > TOL))) {
          std::cout << "Empty cone." << std::endl;
          assert(c != goal_id);
          continue;
        }
        /**
         * Projection to force controlled subspace
         */
        // cylindrificate the intersection in the V directions
        MatrixXd R_V_lines(V_control_directions_r.rows(), V_control_directions_r.cols() + 1);
        R_V_lines << 2*VectorXd::Ones(V_control_directions_r.rows()), V_control_directions_r;
        PPL::Generator_System gs_V;
        Poly::constructPPLGeneratorsFromV(R_V_lines, &gs_V);
        ph1.add_generators(gs_V);
        ph1.minimized_constraints();

        // project to force controlled subspace
        MatrixXd c_A_temp;
        VectorXd c_b_temp;
        Poly::getFacetFromPPL(ph1, &c_A_temp, &c_b_temp);
        assert(c_b_temp.norm() < TOL);
        c_A_temp = c_A_temp * R_a_inv;
        cp_A_temp = c_A_temp.leftCols(action.n_af);
        cp_A_temp.rowwise().normalize();

        // // construct the force controlled subspace
        // MatrixXd R_F_lines(F_control_directions_r.rows(), F_control_directions_r.cols() + 1);
        // R_F_lines << 2*VectorXd::Ones(F_control_directions_r.rows()), F_control_directions_r;
        // PPL::Generator_System gs_F;
        // Poly::constructPPLGeneratorsFromV(R_F_lines, &gs_F);
        // PPL::C_Polyhedron ph_F(6, PPL::EMPTY);
        // ph_F.add_generators(gs_F);

        // // Projection by intersection
        // ph1.intersection_assign(ph_F);


        // clean up R_mode
        // 1. get rid of point origin
        int id0 = -1;
        for (int i = 0; i < R_mode.rows(); ++i) {
          if (int(R_mode(i, 0)) == 1) {
            id0 = i;
            break;
          }
        }
        Eigen::MatrixXd R_;
        if (id0 >= 0) {
          R_ = Eigen::MatrixXd(R_mode.rows() - 1, R_mode.cols());
          R_.topRows(id0) = R_mode.topRows(id0);
          R_.bottomRows(R_mode.rows() - id0 - 1) = R_mode.bottomRows(R_mode.rows() - id0 - 1);
        }
        assert(R_.leftCols<1>().norm() < 1e-10); // make sure the all the rest are rays
        // 2. get rid of the first column
        cone_of_the_mode = R_.rightCols(R_.cols()-1);
        cone_of_the_mode.rowwise().normalize();

        // std::cout << cone_projection.rows() << " generators.";
        if (goal_id == c) {
          std::cout << " (id: Goal)" << std::endl;
          // generators projection
          cone_projection_goal_r = cone_of_the_mode * F_control_directions_r.transpose();
          cone_projection_goal_r.rowwise().normalize();
          cp_goal_A = cp_A_temp;
        } else {
          // cones_projection_r.push_back(cone_projection);
          cp_A.push_back(cp_A_temp);
          std::cout << " (id:" << cp_A.size() - 1 << ") ";;
          std::cout << " A rows: " << cp_A_temp.rows() << std::endl;
        }
      }
      assert(cone_projection_goal_r.norm() > 1e-5);
      assert(cone_projection_goal_r.rows() >= action.n_af); // ideally we should check its rank

      // get rid of cones that are not adjacent to goal cone
      std::vector<MatrixXd> cp_A_selected;
      std::cout << "reducing irrelevant cones. Total number of cones before: " << cp_A.size() << std::endl;
      for (int i = 0; i < cp_A.size(); ++i) {
        VectorXd xs = VectorXd::Zero(cp_goal_A.cols());
        MatrixXd cp_all_A(cp_goal_A.rows() + cp_A[i].rows(), cp_goal_A.cols());
        cp_all_A << cp_goal_A, cp_A[i];
        bool is_feasible = Poly::lpfeasibility(cp_all_A, -1e-7*Eigen::VectorXd::Ones(cp_all_A.rows()), &xs);
        std::cout << "is_feasible: " << is_feasible << ", xs: " << xs.transpose() << std::endl;
        if (is_feasible && xs.norm() > TOL) {
          cp_A_selected.push_back(cp_A[i]);
        } else {
        }
      }
      cp_A = cp_A_selected;
      std::cout << "Total number of cones after: " << cp_A.size() << std::endl;

      time_stats_projection = timer.toc();
      timer.tic();

      std::cout << "[WrenchStamping]  3.2 Sample wrenches and find feasible ones." << std::endl;
      // sample wrenches in the projection of the goal cone
      int ng = cone_projection_goal_r.rows();
      int ns = 0;
      if (action.n_af == 1) ns = 1;
      else if (action.n_af == 2) ns = 10;
      else if (action.n_af == 3) ns = 20;
      else if (action.n_af == 4) ns = 50;
      else if (action.n_af == 5) ns = 200;
      std::vector<VectorXd> cp_b;
      for (int i = 0; i < cp_A.size(); ++i) cp_b.push_back(VectorXd::Zero(cp_A[i].rows()));
      Eigen::VectorXd x0 = cone_projection_goal_r.colwise().mean().transpose();
      std::cout << "\ndebug: x0: " << x0.transpose() << "\n\n";
      double max_radius = 2;
      std::vector<VectorXd> wrench_samples = Poly::sampleInP1OutOfP2(
          cp_goal_A, VectorXd::Zero(cp_goal_A.rows()),
          cp_A, cp_b, x0, ns, max_radius);

      std::cout << "Found " << wrench_samples.size() << " feasible solutions." << std::endl;
      double control_stability_margin = -1;
      VectorXd wrench_best;
      for (int s = 0; s < wrench_samples.size(); ++s) {
        std::cout << "[WrenchStamping]    Sample #" << s << ": ";
        // check if the sample is within any other cones
        bool infeasible_sample = false;
        for (int i = 0; i < cp_A.size(); ++i) {
          Eigen::VectorXd cp_b = cp_A[i]*wrench_samples[s];
          if (cp_b.maxCoeff() <= 0) {
            std::cout << i << "th cone infeasible." << std::endl;
            infeasible_sample = true;
            break;
          }
        }
        if (infeasible_sample) continue;

        // compute its distance to all other cone projections
        double min_dist = 999999.9;
        for (int i = 0; i < cp_A.size(); ++i) {
          double dist = Poly::distP2Polyhedron(wrench_samples[s], cp_A[i],
              Eigen::VectorXd::Zero(cp_A[i].rows()), Eigen::VectorXd::Zero(wrench_samples[s].rows()));
          if (dist < min_dist) min_dist = dist;
        }
        if (min_dist > control_stability_margin) {
          control_stability_margin = min_dist;
          wrench_best = wrench_samples[s];
        }
        std::cout << "min_dist:" << min_dist << std::endl;
      }
      if (control_stability_margin < 0) {
        std::cout << "[WrenchStamping] 3. Force control has no solution." << std::endl;
        continue;
      }
      double time_stats_force_control = timer.toc();
      std::cout << "[WrenchStamping] 3. Force control: " << wrench_best.transpose() << std::endl;
      action.eta_af = -kContactForce*wrench_best; // the minus sign comes from force balance
      // now we have the control stability margin
      Eigen::VectorXd kv2_vec(12);
      kv2_vec << kv_vec, kv_vec;
      Eigen::MatrixXd Kv2 = kv2_vec.asDiagonal();
      action.R_a.topRows(action.n_af) *= Kf;
      action.R_a.bottomRows(action.n_av) *= Kv;
      action.C *= Kv2;
      // Eigen::MatrixXd R_a_inv = action.R_a.inverse();
      // Eigen::VectorXd V = Eigen::VectorXd::Zero(kDimActualized);
      // V.tail(action.n_av) = action.w_av;
      // Eigen::VectorXd V_T = R_a_inv*V;
      // Eigen::VectorXd F_T = action.R_a.topRows(action.n_af).transpose() * action.eta_af;
      if (print_level > 0) {
        std::cout << " 4. Results:" << std::endl;
        // std::cout << "   geometrical_stability_margin: " << geometrical_stability_margin << std::endl;
        std::cout << "   control_stability_margin: " << control_stability_margin << std::endl;
        std::cout << "   R_a:\n" << action.R_a << std::endl;
        std::cout << "   w_av:\n" << action.w_av << std::endl;
        std::cout << "   eta_af:\n" << action.eta_af << std::endl;
        // std::cout << "   V_T:" << V_T.transpose() << std::endl;
        // std::cout << "   F_T:" << F_T.transpose() << std::endl;

        std::cout << "Timing statistics:\n";
        std::cout << "  time_stats_initialization: " << time_stats_initialization << " ms\n";
        std::cout << "  time_stats_hybrid_servoing: " << time_stats_hybrid_servoing << " ms\n";
        std::cout << "  time_stats_velocity_filtering: " << time_stats_velocity_filtering << " ms\n";
        std::cout << "  time_stats_projection: " << time_stats_projection << " ms\n";
        std::cout << "  time_stats_force_control: " << time_stats_force_control << " ms\n";
        std::cout << "  Total: " << time_stats_initialization + time_stats_hybrid_servoing + time_stats_velocity_filtering
            + time_stats_projection + time_stats_force_control << " ms\n";
      }
    }
  }
  return;
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
    // std::cout << "  cs_mode: " << cs_mode.transpose() << std::endl;
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
          // std::cout << "sss_mode: " << sss_mode << std::endl;
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


Eigen::MatrixXd getConeOfTheMode_2d(const Eigen::MatrixXd &cone_allFix,
    const Eigen::VectorXi &modes) {
  int num_contacts = modes.size();
  std::vector<double> generators; // concatenation of all generators
  for (int i = 0; i < num_contacts; ++i) {
    // 0:separation 1:fixed 2/3: sliding
    int id_start = 2*i;
    MatrixXd generators_add(0, 0);
    if (modes[i] == 1) {
      // record all 2 generators
      generators_add = cone_allFix.middleRows(id_start, 2);
    } else if (modes[i] == 2) {
      // Right sliding, left edge of friction cone
      generators_add = cone_allFix.middleRows(id_start, 1);
    } else if (modes[i] == 3) {
      // Left sliding, right edge of friction cone
      generators_add = cone_allFix.middleRows(id_start + 1, 1);
    }
    // add generators
    if (generators_add.size() > 0) {
      int g_size = generators.size();
      generators.resize(g_size + generators_add.rows()*3);
      Eigen::MatrixXd::Map(&generators[g_size], 3, generators_add.rows()) = generators_add.transpose();
    }
  }
  // convert generators to matrix form
  int num_generators = generators.size()/3;
  return MatrixXd::Map(generators.data(), 3, num_generators).transpose();
}

// Eigen::MatrixXd getConeOfTheMode(const Eigen::MatrixXd &cone_allFix,
//     const Eigen::VectorXi &sss_mode, const Eigen::VectorXi &s_mode, int kNumSlidingPlanes) {
Eigen::MatrixXd getConeOfTheMode(const Eigen::MatrixXd &cone_allFix,
    const Eigen::VectorXi &sss_mode, int kNumSlidingPlanes) {

  int num_contacts = sss_mode.size();
  // Eigen::VectorXi s_mode01 = s_mode;
  // s_mode01.array() = (s_mode01.array() + 1)/2; // -1/1 -> 0/1
  std::vector<double> generators; // concatenation of all generators
  for (int i = 0; i < num_contacts; ++i) {
    // -1: sticking   0: sliding   1: separation
    int id_start = (2*kNumSlidingPlanes+1)*i;
    MatrixXd generators_add(0, 0);
    if (sss_mode[i] == -1) {
      // record all 2*kNumSlidingPlanes generators
      generators_add = cone_allFix.middleRows(id_start, 2*kNumSlidingPlanes);
    } else if (sss_mode[i] == 0) {
      // // find two corresponding generators
      // std::vector<double> v;
      // v.resize(2*6);

      // // Compute the order from bit-wise code:
      // int id_s = kNumSlidingPlanes*i;
      // int mode_sum = s_mode01.segment(id_s, kNumSlidingPlanes).sum();
      // int id_g;
      // if (s_mode01(id_s) == 1)
      //     id_g = mode_sum - 1;
      // else
      //     id_g = 2*kNumSlidingPlanes - mode_sum - 1;
      // generators_add = cone_allFix.middleRows(id_start + id_g, 2);
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

void getConstraintOfTheMode_2d(
    const Eigen::MatrixXd &J_e_AF, const Eigen::MatrixXd &J_h_AF,
    const Eigen::VectorXi &mode_e, const Eigen::VectorXi &mode_h,
    Eigen::MatrixXd *N, Eigen::MatrixXd *Nu) {
  const int kNumEContacts = mode_e.size();
  const int kNumHContacts = mode_h.size();
  const int kDim = J_h_AF.cols();
  assert(J_e_AF.cols() == kDim);

  // count number of rows
  // 0:separation 1:fixed 2/3: sliding
  int NRows = 0;
  int NuRows = 0;
  for (int i = 0; i < kNumEContacts; ++i) {
    if (mode_e[i] == 1) {
      NRows += 2;
    } else if (mode_e[i] == 0) {
      NuRows += 1;
    } else {
      NRows += 1;
      NuRows += 1;
    }
  }
  for (int i = 0; i < kNumHContacts; ++i) {
    if (mode_h[i] == 1) {
      NRows += 2;
    } else if (mode_h[i] == 0) {
      NuRows += 1;
    } else {
      NRows += 1;
      NuRows += 1;
    }
  }

  *N  = Eigen::MatrixXd::Zero(NRows, 2*kDim);
  *Nu = Eigen::MatrixXd::Zero(NuRows, 2*kDim);

  int n_count = 0, nu_count = 0;
  for (int i = 0; i < kNumEContacts; ++i) {
    // 0: separation, 1: sticking, 2:right sliding 3: left sliding
    if (mode_e[i] == 0) {
      // separation
      Nu->block(nu_count, 0, 1, kDim) = J_e_AF.middleRows(i, 1);
      nu_count ++;
    } else if (mode_e[i] == 1) {
      // sticking
      N->block(n_count, 0, 1, kDim) = J_e_AF.middleRows(i, 1);
      N->block(n_count+1, 0, 1, kDim) = J_e_AF.middleRows(kNumEContacts + i, 1);
      n_count += 2;
    } else if (mode_e[i] == 2) {
      // right sliding
      N->block(n_count, 0, 1, kDim) = J_e_AF.middleRows(i, 1);
      n_count ++;
      Nu->block(nu_count, 0, 1, kDim) = -J_e_AF.middleRows(kNumEContacts + i, 1);
      nu_count ++;
    } else {
      // left sliding
      N->block(n_count, 0, 1, kDim) = J_e_AF.middleRows(i, 1);
      n_count ++;
      Nu->block(nu_count, 0, 1, kDim) = J_e_AF.middleRows(kNumEContacts + i, 1);
      nu_count ++;
    }
  }
  for (int i = 0; i < kNumHContacts; ++i) {
    // 0: separation, 1: sticking, 2:right sliding 3: left sliding
    if (mode_h[i] == 0) {
      // separation
      Nu->block(nu_count, 0, 1, kDim) = -J_h_AF.middleRows(i, 1);
      Nu->block(nu_count, kDim, 1, kDim) = J_h_AF.middleRows(i, 1);
      nu_count ++;
    } else if (mode_h[i] == 1) {
      // sticking
      N->block(n_count, 0, 1, kDim) = -J_h_AF.middleRows(i, 1);
      N->block(n_count, kDim, 1, kDim) = J_h_AF.middleRows(i, 1);
      N->block(n_count+1, 0, 1, kDim) = -J_h_AF.middleRows(kNumHContacts + i, 1);
      N->block(n_count+1, kDim, 1, kDim) = J_h_AF.middleRows(kNumHContacts + i, 1);
      n_count += 2;
    } else if (mode_h[i] == 2) {
      // right sliding
      N->block(n_count, 0, 1, kDim) = -J_h_AF.middleRows(i, 1);
      N->block(n_count, kDim, 1, kDim) = J_h_AF.middleRows(i, 1);
      n_count ++;
      Nu->block(nu_count, 0, 1, kDim) = J_h_AF.middleRows(kNumHContacts + i, 1);
      Nu->block(nu_count, kDim, 1, kDim) = -J_h_AF.middleRows(kNumHContacts + i, 1);
      nu_count ++;
    } else {
      // left sliding
      N->block(n_count, 0, 1, kDim) = -J_h_AF.middleRows(i, 1);
      N->block(n_count, kDim, 1, kDim) = J_h_AF.middleRows(i, 1);
      n_count ++;
      Nu->block(nu_count, 0, 1, kDim) = -J_h_AF.middleRows(kNumHContacts + i, 1);
      Nu->block(nu_count, kDim, 1, kDim) = J_h_AF.middleRows(kNumHContacts + i, 1);
      nu_count ++;
    }
  }
}


// for now, only implemented the bilateral part for hybrid servoing
void getConstraintOfTheMode(
    const Eigen::MatrixXd &J_e_AF, const Eigen::MatrixXd &J_h_AF,
    const Eigen::VectorXi &sss_mode_e, const Eigen::VectorXi &sss_mode_h,
    Eigen::MatrixXd *N, Eigen::MatrixXd *Nu) {
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

  int n_count = 0, nu_count = 0;
  for (int i = 0; i < kNumEContacts; ++i) {
    // -1: sticking   0: sliding   1: separation
    if (sss_mode_e[i] == 1) {
      // separation
      Nu->block(nu_count, 0, 1, kDim) = J_e_AF.middleRows(i, 1);
      nu_count ++;
      continue;
    }
    N->block(n_count, 0, 1, kDim) = J_e_AF.middleRows(i, 1);
    n_count ++;
    if (sss_mode_e[i] == -1) {
      // sticking
      N->block(n_count, 0, 2, kDim) = J_e_AF.middleRows(kNumEContacts + i*2, 2);
      n_count += 2;
    }
  }
  for (int i = 0; i < kNumHContacts; ++i) {
    // -1: sticking   0: sliding   1: separation
    if (sss_mode_h[i] == 1) {
      // separation
      Nu->block(nu_count, 0, 1, kDim) = -J_h_AF.middleRows(i, 1);
      Nu->block(nu_count, kDim, 1, kDim) = J_h_AF.middleRows(i, 1);
      nu_count ++;
      continue;
    }
    N->block(n_count, 0, 1, kDim) = -J_h_AF.middleRows(i, 1);
    N->block(n_count, kDim, 1, kDim) = J_h_AF.middleRows(i, 1);
    n_count++;
    if (sss_mode_h[i] == -1) {
      // sticking
      N->block(n_count, 0, 2, kDim) = -J_h_AF.middleRows(kNumHContacts + i*2, 2);
      N->block(n_count, kDim, 2, kDim) = J_h_AF.middleRows(kNumHContacts + i*2, 2);
      n_count += 2;
    }
  }
}

Eigen::MatrixXd getSlidingGeneratorsFromOneContact(const Eigen::MatrixXd &vel_samples_on_contact,
    const Eigen::MatrixXd &Jt, const Eigen::MatrixXd &Jn, double friction) {
    int kDim = Jt.cols();
    if (vel_samples_on_contact.norm() < 100*TOL) {
        return Eigen::MatrixXd(0, kDim);
    }
    Eigen::MatrixXd g_sampled(vel_samples_on_contact.rows(), kDim);
    for (int i = 0; i < vel_samples_on_contact.rows(); ++i) {
      g_sampled.middleRows(i, 1) = (friction * vel_samples_on_contact.middleRows(i, 1).normalized() * Jt + Jn)/(std::sqrt(1+friction*friction));
  }
  return g_sampled;
}

int findIdInModes(const Eigen::VectorXi &target_mode, const Eigen::MatrixXi &modes) {
  assert(modes.cols() == target_mode.size());
  int dim = modes.cols();
  for (int row = 0; row < modes.rows(); ++row) {
    bool found = true;
    for (int i = 0; i < dim; ++i) {
      if (target_mode(i) != modes(row, i)) {
        found = false;
        break;
      }
    }
    if (found) return row;
  }
  return -1;
}


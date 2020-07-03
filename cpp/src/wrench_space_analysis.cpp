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
    const double kFrictionE, const double kFrictionH,
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
  std::cout << "Jac_e: " << Jac_e.rows() << " x " << Jac_e.cols() << std::endl;
  std::cout << "Jac_h: " << Jac_h.rows() << " x " << Jac_h.cols() << std::endl;
  std::cout << "eCone_allFix_r: " << eCone_allFix_r.rows() << " x " << eCone_allFix_r.cols() << std::endl;
  std::cout << "hCone_allFix_r: " << hCone_allFix_r.rows() << " x " << hCone_allFix_r.cols() << std::endl;
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

  int kNumEContacts = e_cs_modes.cols();
  int kNumHContacts = h_cs_modes.cols();
  int kDim = Jac_e.cols();
  // scaling for generalized velocity
  // Scaling is essentially changing units. For example, if we change m to mm, then
  // Translational velocity: 1 -> 1000 (m/s -> mm/s)
  // Rotational velocity: 1 -> 1 (rad/s -> rad/s)
  // Force: 1 -> 1 (N -> N)
  // Torque: 1 -> 1000 (N*m -> N*mm)
  // Define Kv = diag(1000,1000,1000,1,1,1), Kf=diag(1,1,1,1000,1000,1000)
  // Then we have
  //    V_scaled = Kv*V, f_scaled = Kf*f
  //    N*V=0 -> N*Kv_inv*V_scaled = 0, so  N_scaled = N*Kv_inv
  //    J'*tau=f -> J'*tau = Kf_inv*f_scaled -> (J*Kf)'*tau = f_scaled, so  J_scaled = J*Kf
  // For hfvc computed in scaled space, we have
  //    Rv_scaled*V_scaled = omega_scaled
  //    Rf_scaled*f_scaled = eta_scaled
  // We can transform these back to constraints on v, f:
  //    Rv_scaled*Kv*V = omega_scaled
  //    Rf_scaled*Kf*f = eta_scaled
  // So, in our output,
  //    Rv = Rv_scaled*Kv,   omega = omega_scaled
  //    Rf = Rf_scaled*Kf,   eta = eta_scaled
  // And R = [Rv; Rf]
  Vector6d kv_vec, kv_inv_vec, kf_vec;
  double scale = 1./kCharacteristicLength;
  kv_vec << scale, scale, scale, 1, 1, 1;
  kv_inv_vec << 1./scale, 1./scale, 1./scale, 1, 1, 1;
  kf_vec << 1, 1, 1, scale, scale, scale;
  MatrixXd Kv = kv_vec.asDiagonal();
  MatrixXd Kv_inv = kv_inv_vec.asDiagonal();
  MatrixXd Kf = kf_vec.asDiagonal();

  Jac_e = Jac_e * Kv_inv;
  Jac_h = Jac_h * Kv_inv;
  eCone_allFix_r = eCone_allFix_r * Kf;
  hCone_allFix_r = hCone_allFix_r * Kf;
  if (flag_given_goal_velocity) {
    G.leftCols(kDim) = G.leftCols(kDim) * Kv_inv;
    G.rightCols(kDim) = G.rightCols(kDim) * Kv_inv;
  }

  // for crashing check
  MatrixXd cone_allFix_r;
  Poly::coneIntersection(eCone_allFix_r, hCone_allFix_r, &cone_allFix_r);

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

  std::vector<Eigen::MatrixXd> sample_grids;
  sample_grids.push_back(Eigen::MatrixXd(2, 1));
  sample_grids.push_back(Eigen::MatrixXd(4, 2));
  sample_grids.push_back(Eigen::MatrixXd(8, 3));
  sample_grids[0] << 1, -1;
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

  Eigen::FullPivLU<MatrixXd> lu; // for quickly solving linear system
  lu.setThreshold(TOL);

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
      std::cout << "[WrenchStamping] goal mdoe: e " << e_sss_mode_goal.transpose() << ", h " << h_sss_mode_goal.transpose() << std::endl;
      int goal_id_e = findIdInModes(e_sss_mode_goal, e_sss_modes);
      int goal_id = -1;
      assert(goal_id_e >= 0);

      std::vector<Eigen::MatrixXd> e_cones_VFeasible;
      std::vector<Eigen::VectorXi> e_modes_VFeasible;
      std::vector<Eigen::VectorXi> h_modes_VFeasible;

      std::cout << "[WrenchStamping] 1. HFVC" << std::endl;
      Eigen::MatrixXd N, Nu; // TODO: Nu and T seem useless here
      getConstraintOfTheMode(Jac_e, Jac_h,
          e_sss_mode_goal, h_sss_mode_goal,
          &N, &Nu);

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
      MatrixXd F_control_directions_r = action.R_a.topRows(action.n_af);
      // std::cout << "[WrenchStamping]    Debug: action.C\n" << action.C << std::endl;
      // std::cout << "[WrenchStamping]    Debug: action.b_C\n" << action.b_C << std::endl;

      // Crashing check
      MatrixXd R;
      std::cout << "[debug] cone_allFix_r: " << cone_allFix_r.rows() << " x " << cone_allFix_r.cols() << std::endl;
      std::cout << "[debug] V_control_directions_r: " << V_control_directions_r.rows() << " x " << V_control_directions_r.cols() << std::endl;

      Poly::coneIntersection(cone_allFix_r, V_control_directions_r, &R); // this line has errors sometimes
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
              &N, &Nu);

          Eigen::MatrixXd NC(N.rows() + action.C.rows(), N.cols());
          NC << N, action.C;
          Eigen::VectorXd b_NC = Eigen::VectorXd::Zero(NC.rows());
          b_NC.tail(action.b_C.rows()) = action.b_C;
          assert(NC.norm() > 10*TOL); //otherwise lu won't be accurate
          lu.compute(NC);
          Eigen::VectorXd sol_NC = lu.solve(b_NC);
          bool a_solution_exists = (NC*sol_NC).isApprox(b_NC, 10.*TOL);
          if (a_solution_exists) {
            Eigen::MatrixXd null_NC = lu.kernel();
            bool has_penetration = false;
            if (Nu.rows() > 0) {
              // check inequalities
              Eigen::MatrixXd contact_normal_proj = Nu*null_NC;
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
                    assert(false); // this should not happen
                    std::cout << 0 << " ";
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

            if (g_sampled.rows() == 0) {
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
          } else {
              // no solution
              std::cout << "No Solution. " << std::endl;
          }// end if a solution exists
        }
      }
      // done with inner SSS loop
      assert(goal_id >= 0);
      /*******************************************************************
       *      Second half: Force Filtering
       */
      // projection cones of the modes onto force-controlled subspace
      std::cout << "[WrenchStamping] 3. Compute force control and control-stability-margin." << std::endl;
      std::vector<Eigen::MatrixXd> cones_projection_r;
      std::vector<Eigen::MatrixXd> cp_A; // A x <= 0
      Eigen::MatrixXd cone_projection_goal_r;
      // cones_projection_r.reserve(e_cones_VFeasible.size());
      // cp_A.resize(e_cones_VFeasible.size());
      Eigen::MatrixXd cone_of_the_mode, cone_projection, cp_A_temp;
      std::cout << "[WrenchStamping]  3.1 Compute cone of the modes." << std::endl;
      for (int c = 0; c < e_cones_VFeasible.size(); ++c) {
        std::cout << "[WrenchStamping]    Check cone " << c << ": ";
        // compute cone of the modes
        Poly::coneIntersection(e_cones_VFeasible[c], hCone_allFix_r, &cone_of_the_mode);
        if (cone_of_the_mode.rows() == 0) {
          std::cout << "Empty cone." << std::endl;
          assert(c != goal_id);
          continue;
        }
        // projection
        cone_projection = cone_of_the_mode * F_control_directions_r.transpose();
        cone_projection.rowwise().normalize();
        // get the inequality representations for the cone projections
        if(!Poly::coneFacetEnumeration(cone_projection, &cp_A_temp)) {
          std::cerr << "coneFacetEnumeration return error." << std::endl;
          std::cout << "[debug] cone_projection: " << cone_projection.rows() << " x " << cone_projection.cols() << std::endl;
          exit(1);
        }

        std::cout << cone_projection.rows() << " generators.";
        if (goal_id == c) {
          std::cout << " (Goal)" << std::endl;
          cone_projection_goal_r = cone_projection;
        } else {
          std::cout << " A rows: " << cp_A_temp.rows() << std::endl;
          cones_projection_r.push_back(cone_projection);
          cp_A.push_back(cp_A_temp);
        }
      }
      assert(cone_projection_goal_r.norm() > 1e-5);
      assert(cone_projection_goal_r.rows() >= action.n_af); // ideally we should check its rank

      std::cout << "[WrenchStamping]  3.2 Sample wrenches and find feasible ones." << std::endl;
      // sample wrenches in the projection of the goal cone
      int ng = cone_projection_goal_r.rows();
      int ns = 0;
      if (action.n_af == 1) ns = 3;
      else if (action.n_af == 2) ns = 10;
      else if (action.n_af == 3) ns = 20;
      else if (action.n_af == 4) ns = 50;
      else if (action.n_af == 5) ns = 100;
      Eigen::MatrixXd rand_weights = Eigen::MatrixXd::Random(ns, ng) + Eigen::MatrixXd::Ones(ns, ng);
      Eigen::VectorXd rand_weights_row_sum = rand_weights.rowwise().sum();
      Eigen::VectorXd wrench_sample, wrench_best;
      double control_stability_margin = -1;
      for (int s = 0; s < ns; ++s) { // make sure they sum to one
        // create the sample
        rand_weights.middleRows(s, 1) /= rand_weights_row_sum(s);
        wrench_sample = (rand_weights.middleRows(s, 1) * cone_projection_goal_r).normalized().transpose();
        std::cout << "[WrenchStamping]    Sample #" << s << ": " << wrench_sample.transpose() << ", ";
        // check if the sample is within any other cones
        bool infeasible_sample = false;
        for (int i = 0; i < cp_A.size(); ++i) {
          Eigen::VectorXd cp_b = cp_A[i]*wrench_sample;
          if (cp_b.maxCoeff() <= 0) {
            infeasible_sample = true;
            break;
          }
        }
        if (infeasible_sample) {
          std::cout << "infeasible." << std::endl;
          continue;
        }
        // compute its distance to all other cone projections
        double min_ang_dist = 999999.9;
        for (int i = 0; i < cp_A.size(); ++i) {
          double ang = Poly::distRay2ConeFromOutside(wrench_sample, cp_A[i], cones_projection_r[i]);
          if (ang < min_ang_dist) min_ang_dist = ang;
        }
        if (min_ang_dist > control_stability_margin) {
          control_stability_margin = min_ang_dist;
          wrench_best = wrench_sample;
        }
        std::cout << "min_ang_dist:" << min_ang_dist << std::endl;
      }
      if (control_stability_margin < 0) {
        std::cout << "[WrenchStamping] 3. Force control has no solution." << std::endl;
        continue;
      }
      std::cout << "[WrenchStamping] 3. Force control: " << wrench_best.transpose() << std::endl;
      // now we have the control stability margin
      action.eta_af = -kContactForce*wrench_best; // the minus sign comes from force balance
      // scale back so the control constraints work on normal units
      action.R_a.topRows(action.n_af) = action.R_a.topRows(action.n_af) * Kf;
      action.R_a.bottomRows(action.n_av) = action.R_a.bottomRows(action.n_av) * Kv;
      // end of handling a goal mode
      //     V_T = solution.R_a_inv*[zeros(solution.n_af,1); solution.w_av];
      //     F_T = solution.R_a_inv*[solution.eta_af; zeros(solution.n_av, 1)];
    }
  }

  return;
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


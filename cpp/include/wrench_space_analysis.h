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
//    R_scaled*V_scaled = omega_scaled
//    R_scaled^T*eta_scaled = f_scaled
// To retrieve the control, transform these back to constraints on v, f:
//    R_scaled*Kv*V = omega_scaled
//    R_scaled^T*eta_scaled = Kf*f -> (R_scaled*Kf_inv)^T * eta_scaled = f
//        ->  (R_scaled*Kf_inv*K)^T * eta_scaled/K = f
// Note R_scaled*Kf_inv*K = R_scaled*Kv, so the two agrees
// So, in our output,
//    R = R_scaled*Kv
//    omega = omega_scaled
//    eta = eta_scaled/K

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "solvehfvc.h"



class WrenchSpaceAnalysis {
public:
    WrenchSpaceAnalysis();
    ~WrenchSpaceAnalysis();
    // Geometrical parameters:
    //  Jac_e, Jac_h
    //  eCone_allFix_r, hCone_allFix_r: each row denotes a generator
    //  F_G: gravity vector
    // Magnitude parameters:
    //  kContactForce
    //  kCharacteristicLength
    // List of all velocity-consistent contact modes
    //  e_cs_modes, h_cs_modes: number of cs modes x number of contacts
    //  e_ss_modes, h_ss_modes: a std vector, each element is a (number of modes x number of total sliding directions) matrix describing the stick-sliding modes for a particular cs modes.
    // Optional arguments (use an empty matrix to denote an unused argument):
    //  A, b: additional force constraint
    //  G, b_G: goal velocity description: Gv = b_G
    //  e_mode_goal, h_mode_goal: a particular goal mode
    double wrenchStamping_2d(Eigen::MatrixXd Jac_e, Eigen::MatrixXd Jac_h,
        Eigen::MatrixXd eCone_allFix, Eigen::MatrixXd hCone_allFix,
        Eigen::VectorXd F_G,
        const double kContactForce, const double kFrictionE, const double kFrictionH,
        const double kCharacteristicLength,
        Eigen::MatrixXd G, const Eigen::VectorXd &b_G,
        const Eigen::MatrixXi &e_modes, const Eigen::MatrixXi &h_modes,
        const Eigen::VectorXi &e_modes_goal, const Eigen::VectorXi &h_modes_goal,
        HFVC *action);

    // Geometrical parameters:
    //  Jac_e, Jac_h
    //  eCone_allFix_r, hCone_allFix_r: each row denotes a generator
    //  F_G: gravity vector
    // Magnitude parameters:
    //  kContactForce
    //  kCharacteristicLength
    // List of all velocity-consistent contact modes
    //  e_cs_modes, h_cs_modes: number of cs modes x number of contacts
    //  e_ss_modes, h_ss_modes: a std vector, each element is a (number of modes x number of total sliding directions) matrix describing the stick-sliding modes for a particular cs modes.
    // Optional arguments (use an empty matrix to denote an unused argument):
    //  A, b: additional force constraint
    //  G, b_G: goal velocity description
    //  e_mode_goal, h_mode_goal: a particular goal mode
    //  return: <Geometrical stability margin, control stability margin>
    std::pair<double, double> wrenchStamping(Eigen::MatrixXd Jac_e, Eigen::MatrixXd Jac_h,
        Eigen::MatrixXd eCone_allFix, Eigen::MatrixXd hCone_allFix,
        Eigen::VectorXd F_G,
        const double kContactForce, const double kFrictionE, const double kFrictionH,
        const double kCharacteristicLength, const int kNumSlidingPlanes,
        const Eigen::MatrixXi &e_cs_modes,
        const std::vector<Eigen::MatrixXi> &e_ss_modes,
        const Eigen::MatrixXi &h_cs_modes,
        const std::vector<Eigen::MatrixXi> &h_ss_modes,
        Eigen::MatrixXd G, const Eigen::VectorXd &b_G,
        const Eigen::MatrixXi &e_cs_modes_goal,
        const std::vector<Eigen::MatrixXi> &e_ss_modes_goal,
        const Eigen::MatrixXi &h_cs_modes_goal,
        const std::vector<Eigen::MatrixXi> &h_ss_modes_goal,
        HFVC *action);

    void updateConstants(
        double kFrictionE, double kFrictionH, double kNumSlidingPlanes,
        double kContactForce, double kObjWeight, double kCharacteristicLength);

    /**
     * @brief      { function_description }
     *
     * @param[in]  kNumEContacts  The k number e contacts
     * @param[in]  kNumHContacts  The k number h contacts
     * @param[in]  CP_H_e         3 x n
     * @param[in]  CN_H_e         3 x n
     * @param[in]  CP_H_h         3 x n
     * @param[in]  CN_H_h         3 x n
     * @param[in]  CP_H_G         The cp h g
     * @param[in]  v_HG           The v hg
     */
    void updateContactGeometry(int kNumEContacts, int kNumHContacts,
        const Eigen::MatrixXd &CP_H_e, const Eigen::MatrixXd &CN_H_e,
        const Eigen::MatrixXd &CP_H_h, const Eigen::MatrixXd &CN_H_h,
        const Eigen::VectorXd &CP_H_G, const Eigen::VectorXd &v_HG);

    void computeContactModes();

    std::pair<double, double> wrenchStampingWrapper(const Eigen::MatrixXd &G, const Eigen::VectorXd &b_G,
      const Eigen::MatrixXi &e_cs_modes_goal,
      const std::vector<Eigen::MatrixXi> &e_ss_modes_goal,
      const Eigen::MatrixXi &h_cs_modes_goal,
      const std::vector<Eigen::MatrixXi> &h_ss_modes_goal, HFVC &action);

    std::pair<double, double> computeStabilityMargin(
      double kFrictionE, double kFrictionH, double kNumSlidingPlanes,
      double kContactForce, double kObjWeight, double kCharacteristicLength,
      int kNumEContacts, int kNumHContacts,
      const Eigen::MatrixXd &CP_H_e, const Eigen::MatrixXd &CN_H_e,
      const Eigen::MatrixXd &CP_H_h, const Eigen::MatrixXd &CN_H_h,
      const Eigen::VectorXd &CP_H_G, const Eigen::VectorXd &v_HG,
      const Eigen::MatrixXd &G, const Eigen::VectorXd &b_G,
      const Eigen::MatrixXi &e_cs_modes_goal,
      const std::vector<Eigen::MatrixXi> &e_ss_modes_goal,
      const Eigen::MatrixXi &h_cs_modes_goal,
      const std::vector<Eigen::MatrixXi> &h_ss_modes_goal);

    /**
     * Calculates the control from a motion plan. Must call updateConstants() before
     * calling this function.
     *
     * All contact normals points inside the object.
     *
     * @param[in]  obj_traj     7 x N, the object pose traj
     * @param[in]  finger_traj  6n x N, the finger contact traj, p1n1p2n2...
     * @param[in]  CP_W_e_traj  The environmental contact point location trajectory.
     *                          Each element is a 3xn matrix.
     * @param[in]  CN_W_e_traj  The environmental contact normal trajectory. Each
     *                          element is a 3xn matrix.
     * @param[in]  p_OG         Object COM location in the object frame.
     * @param[in]  e_ss_modes   Environmental sticking/sliding modes. Each element
     *                          is MatrixXi(1, kNumContactsE), 0: sticking, 1:
     *                          sliding.
     *
     * @return     True if success
     */
    bool computeControlFromMotionPlan(const Eigen::MatrixXd &obj_traj,
        const Eigen::MatrixXd &finger_traj, const std::vector<Eigen::MatrixXd> &CP_W_e_traj,
        const std::vector<Eigen::MatrixXd> &CN_W_e_traj, const Eigen::Vector3d &p_OG,
        const std::vector<Eigen::MatrixXi> &e_cs_modes,
        std::vector<HFVC> &action_traj);
private:
    bool modeCleaning(const Eigen::MatrixXi &cs_modes, const std::vector<Eigen::MatrixXi> &ss_modes, int kNumSlidingPlanes,
        Eigen::MatrixXi *sss_modes, std::vector<Eigen::MatrixXi> *s_modes);

    // sss mode:  -1: sticking   0: sliding   1: separation
    // return: each row is a generator
    Eigen::MatrixXd getConeOfTheMode(const Eigen::MatrixXd &cone_allFix,
        const Eigen::VectorXi &sss_mode, int kNumSlidingPlanes);

    Eigen::MatrixXd getConeOfTheMode_2d(const Eigen::MatrixXd &cone_allFix,
        const Eigen::VectorXi &modes);

    // Nu v >= 0
    void getConstraintOfTheMode_2d(
        const Eigen::MatrixXd &J_e_AF, const Eigen::MatrixXd &J_h_AF,
        const Eigen::VectorXi &mode_e, const Eigen::VectorXi &mode_h,
        Eigen::MatrixXd *N, Eigen::MatrixXd *Nu);
    void getConstraintOfTheMode(
        const Eigen::MatrixXd &J_e_AF, const Eigen::MatrixXd &J_h_AF,
        const Eigen::VectorXi &sss_mode_e, const Eigen::VectorXi &sss_mode_h,
        Eigen::MatrixXd *N, Eigen::MatrixXd *Nu);

    /**
     * Gets the wrench generators for sliding contact points.
     *
     * @param[in]  vel_samples_on_contact  The velocity samples in contact frame
     * @param[in]  Jt                      2xkDim rows of Jt for this contact
     * @param[in]  Jn                      1xkDim row of Jn for this contact
     * @param[in]  friction                The friction
     * @param      g_sampled               The sliding generators from one contact
     *
     * @return     true if no error.
     */
    bool getSlidingGeneratorsFromOneContact(const Eigen::MatrixXd &vel_samples_on_contact,
        const Eigen::MatrixXd &Jt, const Eigen::MatrixXd &Jn, double friction,
        std::vector<Eigen::MatrixXd> *g_sampled);

    /**
     * @brief      Finds identifier in modes.
     *
     * @param[in]  target_mode  The target mode
     * @param[in]  modes        The modes, each row is a mode
     *
     * @return     Return the row id of modes that is the same as target_mode.
     *             Return -1 if not found.
     */
    int findIdInModes(const Eigen::VectorXi &target_mode, const Eigen::MatrixXi &modes);

    /**
     * Force control for wrench stamping. Find a wrench that satisfies:
     *      * It is inside of the polyhedron defined by (pp_goal_A, pp_goal_b)
     *      * Its magnitude is smaller or equal to @p kContactForce
     * and maximize the distance towards every polyhedron in (pps_A, pps_b).
     */
    double forceControl(double kContactForce, int n_af,
        const Eigen::MatrixXd &pp_goal_A, const Eigen::VectorXd &pp_goal_b,
        const std::vector<Eigen::MatrixXd> &pps_A,
        const std::vector<Eigen::VectorXd> &pps_b,
        Eigen::VectorXd *wrench_best);

    /**
     * Algorithm Parameters
     */
    int print_level_;
    std::vector<int> num_seeds_;
    std::vector<int> hitAndRun_num_points_;
    std::vector<int> hitAndRun_num_discard_;
    std::vector<int> hitAndRun_num_runup_;
    std::vector<int> ransac_num_;
    std::vector<int> opt_ins_num_iter_;
    double opt_min_dist_improvement_;

    /**
     * Problem description parameters
     */
    double kFrictionE_;
    double kFrictionH_;
    double kContactForce_;
    double kObjWeight_;
    double kCharacteristicLength_;
    int kNumSlidingPlanes_;

    /**
     * Geometrical information parameters
     */
    int kNumEContacts_;
    int kNumHContacts_;

    Eigen::MatrixXd N_e_, T_e_, N_h_, T_h_;
    Eigen::MatrixXd eCone_allFix_, hCone_allFix_;
    Eigen::VectorXd F_G_;

    /**
     * Contact modes
     */
    Eigen::MatrixXi e_cs_modes_, h_cs_modes_;
    std::vector<Eigen::MatrixXi> e_ss_modes_, h_ss_modes_;
};

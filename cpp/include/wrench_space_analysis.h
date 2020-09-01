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



#define MIN_DIST_IMPROVEMENT 0.005
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
double wrenchSpaceAnalysis_2d(Eigen::MatrixXd Jac_e, Eigen::MatrixXd Jac_h,
    Eigen::MatrixXd eCone_allFix, Eigen::MatrixXd hCone_allFix,
    Eigen::VectorXd F_G,
    const double kContactForce, const double kFrictionE, const double kFrictionH,
    const double kCharacteristicLength,
    Eigen::MatrixXd G, const Eigen::VectorXd &b_G,
    const Eigen::MatrixXi &e_modes, const Eigen::MatrixXi &h_modes,
    const Eigen::VectorXi &e_modes_goal, const Eigen::VectorXi &h_modes_goal,
    int print_level);

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
void wrenchSpaceAnalysis(Eigen::MatrixXd Jac_e, Eigen::MatrixXd Jac_h,
    Eigen::MatrixXd eCone_allFix, Eigen::MatrixXd hCone_allFix,
    Eigen::VectorXd F_G,
    const double kContactForce, const double kFrictionE, const double kFrictionH,
    const double kCharacteristicLength, const int kNumSlidingPlanes,
    const Eigen::MatrixXi &e_cs_modes, const std::vector<Eigen::MatrixXi> &e_ss_modes,
    const Eigen::MatrixXi &h_cs_modes, const std::vector<Eigen::MatrixXi> &h_ss_modes,
    Eigen::MatrixXd G, const Eigen::VectorXd &b_G,
    const Eigen::MatrixXi &e_cs_modes_goal, const std::vector<Eigen::MatrixXi> &e_ss_modes_goal,
    const Eigen::MatrixXi &h_cs_modes_goal, const std::vector<Eigen::MatrixXi> &h_ss_modes_goal,
    int print_level);


bool modeCleaning(const Eigen::MatrixXi &cs_modes, const std::vector<Eigen::MatrixXi> &ss_modes, int kNumSlidingPlanes,
    Eigen::MatrixXi *sss_modes, std::vector<Eigen::MatrixXi> *s_modes);


Eigen::MatrixXd getConeOfTheMode_2d(const Eigen::MatrixXd &cone_allFix,
    const Eigen::VectorXi &modes);

// sss mode:  -1: sticking   0: sliding   1: separation
// return: each row is a generator
// Eigen::MatrixXd getConeOfTheMode(const Eigen::MatrixXd &cone_allFix,
//     const Eigen::VectorXi &sss_mode, const Eigen::VectorXi &s_mode, int kNumSlidingPlanes);
Eigen::MatrixXd getConeOfTheMode(const Eigen::MatrixXd &cone_allFix,
    const Eigen::VectorXi &sss_mode, int kNumSlidingPlanes);

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
 * Gets the sliding generators from one contact.
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


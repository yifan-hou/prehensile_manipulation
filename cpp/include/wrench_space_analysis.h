#include <vector>

#include <Eigen/Dense>
#include <Eigen/Geometry>

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
void wrenchSpaceAnalysis_2d(Eigen::MatrixXd Jac_e, Eigen::MatrixXd Jac_h,
    Eigen::MatrixXd eCone_allFix, Eigen::MatrixXd hCone_allFix,
    const Eigen::VectorXd &F_G,
    const double kContactForce, const double kFrictionE, const double kFrictionH,
    const double kCharacteristicLength, const int kNumSlidingPlanes,
    Eigen::MatrixXd G, const Eigen::VectorXd &b_G,
    const Eigen::MatrixXi &e_modes, const Eigen::MatrixXi &h_modes,
    const Eigen::MatrixXi &e_modes_goal, const Eigen::MatrixXi &h_modes_goal);

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
    const Eigen::VectorXd &F_G,
    const double kContactForce, const double kFrictionE, const double kFrictionH,
    const double kCharacteristicLength, const int kNumSlidingPlanes,
    const Eigen::MatrixXi &e_cs_modes, const std::vector<Eigen::MatrixXi> &e_ss_modes,
    const Eigen::MatrixXi &h_cs_modes, const std::vector<Eigen::MatrixXi> &h_ss_modes,
    Eigen::MatrixXd G, const Eigen::VectorXd &b_G,
    const Eigen::MatrixXi &e_cs_modes_goal, const std::vector<Eigen::MatrixXi> &e_ss_modes_goal,
    const Eigen::MatrixXi &h_cs_modes_goal, const std::vector<Eigen::MatrixXi> &h_ss_modes_goal);


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
 *
 * @return     The sliding generators from one contact.
 */
Eigen::MatrixXd getSlidingGeneratorsFromOneContact(const Eigen::MatrixXd &vel_samples_on_contact,
    const Eigen::MatrixXd &Jt, const Eigen::MatrixXd &Jn, double friction);

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


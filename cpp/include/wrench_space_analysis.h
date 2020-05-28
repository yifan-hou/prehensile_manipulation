#include <vector>

#include <Eigen/Dense>
#include <Eigen/Geometry>


bool modeCleaning(const Eigen::MatrixXi &cs_modes, const std::vector<Eigen::MatrixXi> &ss_modes, int kNumSlidingPlanes,
    Eigen::MatrixXi *sss_modes, std::vector<Eigen::MatrixXi> *s_modes);

// return: each row is a generator
Eigen::MatrixXd getConeOfTheMode(const Eigen::MatrixXd &cone_allFix,
    const Eigen::VectorXi &sss_mode, const Eigen::VectorXi &s_mode, int kNumSlidingPlanes);

void wrenchSpaceAnalysis(Eigen::MatrixXd Jac_e, Eigen::MatrixXd Jac_h,
    Eigen::MatrixXd eCone_allFix, Eigen::MatrixXd hCone_allFix,
    const Eigen::VectorXd &F_G,
    const double kContactForce, const double kCharacteristicLength, const int kNumSlidingPlanes,
    const Eigen::MatrixXi &e_cs_modes, const std::vector<Eigen::MatrixXi> &e_ss_modes,
    const Eigen::MatrixXi &h_cs_modes, const std::vector<Eigen::MatrixXi> &h_ss_modes,
    Eigen::MatrixXd G, const Eigen::VectorXd &b_G,
    const Eigen::VectorXi &e_mode_goal, const Eigen::VectorXi &h_mode_goal);
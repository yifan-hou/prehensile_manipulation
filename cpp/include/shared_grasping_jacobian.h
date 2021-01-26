#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>

typedef Eigen::Matrix<double, 6, 1> Vector6d;
using namespace Eigen;
/**
 * @brief      { function_description }
 *
 * @param[in]  kFrictionE    The environment friction coef
 * @param[in]  kFrictionH    The hand friction coef
 * @param[in]  CP_H_e        The list of environmental contact locations in the
 *                           hand frame.
 * @param[in]  CN_H_e        The list of environmental contact normals in the
 *                           hand frame.
 * @param[in]  CP_H_h        The list of hand contact locations in the hand
 *                           frame.
 * @param[in]  CN_H_h        The list of hand contact normals in the hand frame.
 * @param      N_e           Contact wrenches for environmental contact normals
 * @param      T_e           Contact wrenches for environmental contact
 *                           tangentials
 * @param      N_h           Contact wrenches for hand contact normals
 * @param      T_h           Contact wrenches for hand contact tangentials
 * @param      eCone_allFix  Contact wrenches for environmental friction cone
 *                           edges
 * @param      hCone_allFix  Contact wrenches for hand friction cone edges
 *
 * @return     { description_of_the_return_value }
 */
bool sharedGraspingGeometryProcessing2d(
    double kFrictionE, double kFrictionH,
    const MatrixXd &CP_H_e, const MatrixXd &CN_H_e,
    const MatrixXd &CP_H_h, const MatrixXd &CN_H_h,
    MatrixXd &N_e, MatrixXd &T_e, MatrixXd &N_h, MatrixXd &T_h,
    MatrixXd &eCone_allFix, MatrixXd &hCone_allFix);

bool sharedGraspingGeometryProcessing3d(
    double kFrictionE, double kFrictionH, double kNumSlidingPlanes,
    const MatrixXd &CP_H_e, const MatrixXd &CN_H_e,
    const MatrixXd &CP_H_h, const MatrixXd &CN_H_h,
    MatrixXd &N_e, MatrixXd &T_e, MatrixXd &N_h, MatrixXd &T_h,
    MatrixXd &eCone_allFix, MatrixXd &hCone_allFix);

Vector3d getGravityVector2d(double kObjWeight, const Vector2d &CP_H_G,
    const Vector2d &v_HG);

Vector6d getGravityVector3d(double kObjWeight, const Vector3d &CP_H_G,
    const Vector3d &v_HG);

/**
 * @brief      Computes the jacobian for 2d.
 *
 * @param[in]  kFriction  The friction coef
 * @param[in]  CP         2 x n
 * @param[in]  CN         2 x n
 * @param      N          n x 3, each row is a wrench
 * @param      T          n x 3, each row is a wrench
 * @param      Cone       2n x 3, each row is a wrench
 */
void getJacobian2d(double kFriction,
    const MatrixXd &CP, const MatrixXd &CN,
    MatrixXd &N, MatrixXd &T, MatrixXd &Cone);

/**
 * @brief      Computes the jacobian for 3d problems.
 *
 * @param[in]  kFriction          The friction coef
 * @param[in]  kNumSlidingPlanes  The k number sliding planes
 * @param[in]  CP                 3 x n
 * @param[in]  CN                 3 x n
 * @param      N                  n x 6, each row is a wrench
 * @param      T                  n x 6, each row is a wrench
 * @param      Cone               2*kNumSlidingPlanes*n x 3, each row is a
 *                                wrench
 */
void getJacobian3d(double kFriction, int kNumSlidingPlanes,
    const MatrixXd &CP, const MatrixXd &CN,
    MatrixXd &N, MatrixXd &T, MatrixXd &Cone);

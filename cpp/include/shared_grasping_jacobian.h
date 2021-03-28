/*
 *  Newton's Law:
 *      J_e'*tau_e - J_h'*tau_h = f_HO
 * Define J = [J_e, 0; -J_h, J_h]: 1 normal, kNumSlidingPlanes tangential; used by contact mode enumeration
 * where: J_e = [N_e; T_e]
 *        J_h = [N_h; T_h]   N: normal, T: tangential (XY for 3D, left for 2D)
 * Each contact contributes 1 normal, 2 (1 for planar problem) tangential constraints.
 * eCone, hCone: each row is a wrench space generator created by an edge of a friction cone.
 *   3D: Each contact contributes 2d + 1 edges; the last one is a copy of the first one
 *   2D: Each contact contributes 2 edges (left, right)
 * If planar, kNumSlidingPlanes must = 1, the computation puts everything on XY plane
 *
 * Left/right convention in 2d:
 *    Object motion w.r.t. the environment
 *
 *    cone edge 1,      cone edge 2
 *         -----\-------/------
 *         |     \  N  /      |
 *         |      \ ^ /       |
 * T_e <---|       \|/        | ---> right sliding
 *     =============|==================
 *
 */
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
 * @param      T                  2n x 6, each row is a wrench
 * @param      Cone               2*kNumSlidingPlanes*n x 3, each row is a
 *                                wrench
 */
void getJacobian3d(double kFriction, int kNumSlidingPlanes,
    const MatrixXd &CP, const MatrixXd &CN,
    MatrixXd &N, MatrixXd &T, MatrixXd &Cone);

/**
 * @brief      A simpler interface. Not considering friction
 *
 * @param[in]  CP    { parameter_description }
 * @param[in]  CN    { parameter_description }
 * @param      N     { parameter_description }
 * @param      T     { parameter_description }
 */
void getJacobian3d(const MatrixXd &CP, const MatrixXd &CN, MatrixXd &N,
    MatrixXd &T);
#include "shared_grasping_jacobian.h"
#include <iostream>

using Eigen::Vector3d, Eigen::Matrix3d;
using Eigen::VectorXd, Eigen::MatrixXd;

#define PI 3.1415926

bool sharedGraspingGeometryProcessing2d(
    double kFrictionE, double kFrictionH,
    const MatrixXd &CP_H_e, const MatrixXd &CN_H_e,
    const MatrixXd &CP_H_h, const MatrixXd &CN_H_h,
    MatrixXd &N_e, MatrixXd &T_e, MatrixXd &N_h, MatrixXd &T_h,
    MatrixXd &eCone_allFix, MatrixXd &hCone_allFix) {
  std::cout << "[sharedGraspingGeometryProcessing2d] Calling." << std::endl;
  getJacobian2d(kFrictionE, CP_H_e, CN_H_e, N_e, T_e, eCone_allFix);
  getJacobian2d(kFrictionH, CP_H_h, CN_H_h, N_h, T_h, hCone_allFix);
  std::cout << "[sharedGraspingGeometryProcessing2d] Done." << std::endl;
}

bool sharedGraspingGeometryProcessing3d(
    double kFrictionE, double kFrictionH, double kNumSlidingPlanes,
    const MatrixXd &CP_H_e, const MatrixXd &CN_H_e,
    const MatrixXd &CP_H_h, const MatrixXd &CN_H_h,
    MatrixXd &N_e, MatrixXd &T_e, MatrixXd &N_h, MatrixXd &T_h,
    MatrixXd &eCone_allFix, MatrixXd &hCone_allFix) {

  getJacobian3d(kFrictionE, kNumSlidingPlanes, CP_H_e, CN_H_e, N_e, T_e, eCone_allFix);
  getJacobian3d(kFrictionH, kNumSlidingPlanes, CP_H_h, CN_H_h, N_h, T_h, hCone_allFix);
}


Vector3d getGravityVector2d(double kObjWeight, const Vector2d &CP_H_G,
    const Vector2d &v_HG) {
  Vector3d F_G;
  // H_G = [v_HG; cross(CP_H_G, v_HG)];
  // F_G = kObjWeight*(H_G');
  F_G << v_HG, CP_H_G(0)*v_HG(1) - CP_H_G(1)*v_HG(0);
  F_G = kObjWeight * F_G;
  return F_G;
}

Vector6d getGravityVector3d(double kObjWeight, const Vector3d &CP_H_G,
    const Vector3d &v_HG) {
  Vector6d F_G;
  // H_G = [v_HG; cross(CP_H_G, v_HG)];
  // F_G = kObjWeight*(H_G');
  F_G << v_HG, CP_H_G.cross(v_HG);
  F_G = kObjWeight * F_G;
  return F_G;
}

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
    MatrixXd &N, MatrixXd &T, MatrixXd &Cone) {
  std::cout << "[getJacobian2d] calling." << std::endl;
  int kNumContacts = CP.cols();
  assert(kNumContacts != 0);
  assert(CP.rows() == 2);
  assert(CN.rows() == 2);
  assert(CN.cols() == kNumContacts);

  int kWrenchDim = 3, kEdgesPerContact = 2;
  N = MatrixXd::Zero(kNumContacts, kWrenchDim);
  T = MatrixXd::Zero(kNumContacts, kWrenchDim);

  Cone = MatrixXd::Zero(kEdgesPerContact*kNumContacts, kWrenchDim);

  Vector3d vr = Vector3d::UnitZ();

  double mu_norm = sqrt(1.0 + kFriction*kFriction);
  for (int i = 0; i < kNumContacts; ++i) {
    // contact normal
    VectorXd CN_i = CN.middleCols(i, 1).normalized();
    VectorXd CP_i = CP.middleCols(i, 1);
    // normal contact screw
    Vector3d contact_screw_n_i;
    contact_screw_n_i.head(2) = CN_i;
    contact_screw_n_i(2) = CP_i(0)*CN_i(1) - CP_i(1)*CN_i(0);
    N.middleRows(i, 1) = contact_screw_n_i.transpose();
    // contact tangential and friction cones
    Vector3d CN_i_extended;
    CN_i_extended << CN_i, 0;
    Vector3d CX = vr.cross(CN_i_extended).normalized();
    MatrixXd CT_W(2,2);
    CT_W << CX.head(2).transpose(), -CX.head(2).transpose(); // left, right
    MatrixXd CXY = CX.head(2).transpose();
    // friction cone edges
    MatrixXd CCone = (kFriction*CT_W + MatrixXd::Ones(2, 1)*CN_i.transpose())/mu_norm;
    // contact tangential directions in world frame
    CT_W = CT_W.topRows(1);
    // tangential contact screw (X,Y direction)
    MatrixXd contact_screw_t_i(1,3);
    contact_screw_t_i.leftCols(2) = CXY;
    contact_screw_t_i(2) = CP_i(0)*CXY(1) - CP_i(1)*CXY(0);
    T.middleRows(i, 1) = contact_screw_t_i;
    // cone contact screws (2 directions)
    MatrixXd Conei(2, 3);
    Conei.leftCols(2) = CCone;
    Conei(0, 2) = CP_i(0)*CCone(0,1) - CP_i(1)*CCone(0,0);
    Conei(1, 2) = CP_i(0)*CCone(1,1) - CP_i(1)*CCone(1,0);
    Cone.middleRows(kEdgesPerContact*i, kEdgesPerContact) = Conei;
  }
  std::cout << "[getJacobian2d] Done." << std::endl;
}

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
    MatrixXd &N, MatrixXd &T, MatrixXd &Cone) {
  int kNumContacts = CP.cols();
  assert(kNumContacts != 0);
  assert(CP.rows() == 3);
  assert(CN.rows() == 3);
  assert(CN.cols() == kNumContacts);

  int kWrenchDim = 6;
  int kEdgesPerContact = 2*kNumSlidingPlanes;
  Vector3d z = Vector3d::UnitZ();

  MatrixXd CT;
  // contact wrench for normals and tangentials
  N = MatrixXd::Zero(kNumContacts, kWrenchDim); // normals
  T = MatrixXd::Zero(kNumContacts*2, kWrenchDim); // tangentials

  // sliding planes, discretize sliding directions
  MatrixXd SP = MatrixXd::Zero(kNumSlidingPlanes, 3);
  for (int i = 0; i < kNumSlidingPlanes; ++i) {
    SP(i, 0) = cos(PI/kNumSlidingPlanes*(i+1));
    SP(i, 1) = sin(PI/kNumSlidingPlanes*(i+1));
  }
  MatrixXd SP_ALL = MatrixXd::Zero(kEdgesPerContact, 3);
  SP_ALL << SP, -SP;
  // contact tangential directions in contact frame
  CT = MatrixXd::Zero(kEdgesPerContact, 3);
  for (int i = 0; i < kEdgesPerContact; ++i) {
    Vector3d SP_i = SP_ALL.middleRows(i, 1).transpose();
    CT.middleRows(i, 1) = -z.cross(SP_i).transpose();
  }

  Cone = MatrixXd::Zero(kEdgesPerContact*kNumContacts, kWrenchDim);
  Vector3d vr = Vector3d::Random().normalized();
  std::cout << "vr: \n" << vr << std::endl;
  double mu_norm = sqrt(1.0 + kFriction*kFriction);
  for (int i = 0; i < kNumContacts; ++i) {
    // contact normal
    Vector3d CN_i = CN.middleCols(i,1).normalized();
    Vector3d CP_i = CP.middleCols(i,1);
    // normal contact screw
    Vector6d contact_screw_n_i;
    contact_screw_n_i.head(3) = CN_i;
    contact_screw_n_i.tail(3) = CP_i.cross(CN_i);
    N.middleRows(i, 1) = contact_screw_n_i.transpose();
    // contact tangential and friction cones
    Vector3d CX = vr.cross(CN_i).normalized();
    Vector3d CY = CN_i.cross(CX).normalized();
    MatrixXd CXY(2,3);
    CXY << CX.transpose(), CY.transpose();
    MatrixXd CT_W = CT.leftCols(2) * CXY; // tangential directions
    std::cout << "CXY: \n" << CXY << std::endl;
    std::cout << "CT_W: \n" << CT_W << std::endl;
    // friction cone edges
    MatrixXd CCone = (kFriction*CT_W + MatrixXd::Ones(kEdgesPerContact, 1)*CN_i.transpose())/mu_norm;
    // contact tangential directions in world frame
    CT_W = CT_W.topRows(kNumSlidingPlanes);
    // tangential contact screw (X,Y direction)
    MatrixXd contact_screw_t_i(2,6);
    contact_screw_t_i.leftCols(3) = CXY;
    contact_screw_t_i.block<1,3>(0,3) = CP_i.cross(CX).transpose();
    contact_screw_t_i.block<1,3>(1,3) = CP_i.cross(CY).transpose();
    T.middleRows(2*i, 2) = contact_screw_t_i;
    // cone contact screws (kEdgesPerContact directions)
    MatrixXd Conei(kEdgesPerContact, 6);
    Conei.leftCols(3) = CCone;
    for (int j = 0; j < kEdgesPerContact; ++j) {
      Vector3d CCone_i = CCone.middleRows(j,1).transpose();
      Conei.block<1,3>(j,3) = CP_i.cross(CCone_i).transpose();
    }
    Cone.middleRows(i*kEdgesPerContact, kEdgesPerContact) = Conei;
  }
}


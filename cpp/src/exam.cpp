#include "cnpy.h"
#include <iostream>
#include <Eigen/Dense>
#include "wrench_space_analysis.h"

using namespace Eigen;

int main(void) {
  Eigen::MatrixXd Jac_e, Jac_h, eCone_allFix, hCone_allFix;
  Eigen::VectorXd F_G;
  double kContactForce, kFrictionE, kFrictionH, kCharacteristicLength;
  int kNumSlidingPlanes;
  Eigen::MatrixXi e_cs_modes;
  std::vector<Eigen::MatrixXi> e_ss_modes;
  Eigen::MatrixXi h_cs_modes;
  std::vector<Eigen::MatrixXi> h_ss_modes;
  Eigen::MatrixXd G;
  Eigen::VectorXd b_G;
  Eigen::MatrixXi e_cs_modes_goal;
  std::vector<Eigen::MatrixXi> e_ss_modes_goal;
  Eigen::MatrixXi h_cs_modes_goal;
  std::vector<Eigen::MatrixXi> h_ss_modes_goal;
  int print_level;

  // read from python
  std::string filename = std::string("data/J_e.npy");
  std::cout << "reading: " << filename << std::endl;
  cnpy::NpyArray arr = cnpy::npy_load(filename);
  double* data_ptr_d = arr.data<double>();
  int rows = arr.shape[0];
  int cols = arr.shape[1];
  assert(cols == 6);
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  Jac_e = Eigen::MatrixXd(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      Jac_e(i,j) = data_ptr_d[j];
    }
    data_ptr_d += cols;
  }
  std::cout << "J_e:\n" << Jac_e << std::endl;

  filename = std::string("data/J_h.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_d = arr.data<double>();
  rows = arr.shape[0];
  cols = arr.shape[1];
  assert(cols == 6);
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  Jac_h = Eigen::MatrixXd(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      Jac_h(i,j) = data_ptr_d[j];
    }
    data_ptr_d += cols;
  }
  std::cout << "J_h:\n" << Jac_h << std::endl;

  filename = std::string("data/eCone_allFix.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_d = arr.data<double>();
  rows = arr.shape[0];
  cols = arr.shape[1];
  assert(cols == 6);
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  eCone_allFix = Eigen::MatrixXd(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      eCone_allFix(i,j) = data_ptr_d[j];
    }
    data_ptr_d += cols;
  }
  std::cout << "eCone_allFix:\n" << eCone_allFix << std::endl;

  filename = std::string("data/hCone_allFix.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_d = arr.data<double>();
  rows = arr.shape[0];
  cols = arr.shape[1];
  assert(cols == 6);
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  hCone_allFix = Eigen::MatrixXd(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      hCone_allFix(i,j) = data_ptr_d[j];
    }
    data_ptr_d += cols;
  }
  std::cout << "hCone_allFix:\n" << hCone_allFix << std::endl;

  filename = std::string("data/F_G.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_d = arr.data<double>();
  rows = arr.shape[0];
  cols = arr.shape[1];
  assert(rows == 6);
  assert(cols == 1);
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  F_G = Eigen::VectorXd(rows);
  for (int i = 0; i < rows; ++i) {
    F_G(i) = data_ptr_d[i];
  }
  std::cout << "F_G:\n" << F_G << std::endl;

  filename = std::string("data/kContactForce.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_d = arr.data<double>();
  kContactForce = data_ptr_d[0];
  std::cout << "kContactForce:\n" << kContactForce << std::endl;

  filename = std::string("data/kFrictionE.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_d = arr.data<double>();
  kFrictionE = data_ptr_d[0];
  std::cout << "kFrictionE:\n" << kFrictionE << std::endl;

  filename = std::string("data/kFrictionH.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_d = arr.data<double>();
  kFrictionH = data_ptr_d[0];
  std::cout << "kFrictionH:\n" << kFrictionH << std::endl;

  filename = std::string("data/kCharacteristicLength.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_d = arr.data<double>();
  kCharacteristicLength = data_ptr_d[0];
  std::cout << "kCharacteristicLength:\n" << kCharacteristicLength << std::endl;

  filename = std::string("data/kNumSlidingPlanes.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  int *data_ptr_i = arr.data<int>();
  kNumSlidingPlanes = data_ptr_i[0];
  std::cout << "kNumSlidingPlanes:\n" << kNumSlidingPlanes << std::endl;

  filename = std::string("data/e_cs_modes.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_i = arr.data<int>();
  rows = arr.shape[0];
  cols = arr.shape[1];
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  e_cs_modes = Eigen::MatrixXi(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      e_cs_modes(i,j) = data_ptr_i[j];
    }
    data_ptr_i += cols;
  }
  std::cout << "e_cs_modes:\n" << e_cs_modes << std::endl;

  filename = std::string("data/h_cs_modes.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_i = arr.data<int>();
  rows = arr.shape[0];
  cols = arr.shape[1];
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  h_cs_modes = Eigen::MatrixXi(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      h_cs_modes(i,j) = data_ptr_i[j];
    }
    data_ptr_i += cols;
  }
  std::cout << "h_cs_modes:\n" << h_cs_modes << std::endl;

  filename = std::string("data/e_cs_modes_goal.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_i = arr.data<int>();
  rows = arr.shape[0];
  cols = arr.shape[1];
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  e_cs_modes_goal = Eigen::MatrixXi(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      e_cs_modes_goal(i,j) = data_ptr_i[j];
    }
    data_ptr_i += cols;
  }
  std::cout << "e_cs_modes_goal:\n" << e_cs_modes_goal << std::endl;

  filename = std::string("data/h_cs_modes_goal.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_i = arr.data<int>();
  rows = arr.shape[0];
  cols = arr.shape[1];
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  h_cs_modes_goal = Eigen::MatrixXi(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      h_cs_modes_goal(i,j) = data_ptr_i[j];
    }
    data_ptr_i += cols;
  }
  std::cout << "h_cs_modes_goal:\n" << h_cs_modes_goal << std::endl;

  filename = std::string("data/G.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_d = arr.data<double>();
  rows = arr.shape[0];
  cols = arr.shape[1];
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  G = Eigen::MatrixXd(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      G(i,j) = data_ptr_d[j];
    }
    data_ptr_d += cols;
  }
  std::cout << "G:\n" << G << std::endl;

  filename = std::string("data/b_G.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_d = arr.data<double>();
  rows = arr.shape[0];
  cols = arr.shape[1];
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  b_G = Eigen::VectorXd(rows);
  for (int i = 0; i < rows; ++i) {
    b_G(i) = data_ptr_d[i];
  }
  std::cout << "b_G:\n" << b_G << std::endl;


  filename = std::string("data/print_level.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_i = arr.data<int>();
  print_level = data_ptr_i[0];
  std::cout << "print_level:\n" << print_level << std::endl;


  filename = std::string("data/e_ss_modes_rows.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_i = arr.data<int>();
  cols = arr.shape[1];
  std::cout << "arr.shape[0]: " << arr.shape[0] << ", arr.shape[1]: " << arr.shape[1] << std::endl;
  Eigen::VectorXi e_ss_modes_rows(cols);
  for (int i = 0; i < cols; ++i) {
    e_ss_modes_rows(i) = data_ptr_i[i];
  }
  std::cout << "e_ss_modes_rows:\n" << e_ss_modes_rows << std::endl;
  filename = std::string("data/e_ss_modes_data.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_i = arr.data<int>();
  rows = arr.shape[0];
  cols = arr.shape[1];
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  Eigen::MatrixXi e_ss_modes_data(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      e_ss_modes_data(i,j) = data_ptr_i[j];
    }
    data_ptr_i += cols;
  }
  std::cout << "e_ss_modes_data:\n" << e_ss_modes_data << std::endl;
  int count = 0;
  for (int c = 0; c < e_ss_modes_rows.rows(); ++c) {
    Eigen::MatrixXi e_ss_modes_c = e_ss_modes_data.middleRows(count, e_ss_modes_rows(c));
    count += e_ss_modes_rows(c);
    e_ss_modes.push_back(e_ss_modes_c);
  }


  filename = std::string("data/h_ss_modes_rows.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_i = arr.data<int>();
  rows = arr.shape[0];
  cols = arr.shape[1];
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  Eigen::VectorXi h_ss_modes_rows(cols);
  for (int i = 0; i < cols; ++i) {
    h_ss_modes_rows(i) = data_ptr_i[i];
  }
  std::cout << "h_ss_modes_rows:\n" << h_ss_modes_rows << std::endl;
  filename = std::string("data/h_ss_modes_data.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_i = arr.data<int>();
  rows = arr.shape[0];
  cols = arr.shape[1];
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  Eigen::MatrixXi h_ss_modes_data(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      h_ss_modes_data(i,j) = data_ptr_i[j];
    }
    data_ptr_i += cols;
  }
  std::cout << "h_ss_modes_data:\n" << h_ss_modes_data << std::endl;
  count = 0;
  for (int c = 0; c < h_ss_modes_rows.rows(); ++c) {
    Eigen::MatrixXi h_ss_modes_c = h_ss_modes_data.middleRows(count, h_ss_modes_rows(c));
    count += h_ss_modes_rows(c);
    h_ss_modes.push_back(h_ss_modes_c);
  }



  filename = std::string("data/e_ss_modes_goal_rows.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_i = arr.data<int>();
  rows = arr.shape[0];
  cols = arr.shape[1];
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  Eigen::VectorXi e_ss_modes_goal_rows(cols);
  for (int i = 0; i < cols; ++i) {
    e_ss_modes_goal_rows(i) = data_ptr_i[i];
  }
  std::cout << "e_ss_modes_goal_rows:\n" << e_ss_modes_goal_rows << std::endl;
  filename = std::string("data/e_ss_modes_goal_data.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_i = arr.data<int>();
  rows = arr.shape[0];
  cols = arr.shape[1];
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  Eigen::MatrixXi e_ss_modes_goal_data(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      e_ss_modes_goal_data(i,j) = data_ptr_i[j];
    }
    data_ptr_i += cols;
  }
  std::cout << "e_ss_modes_goal_data:\n" << e_ss_modes_goal_data << std::endl;
  count = 0;
  for (int c = 0; c < e_ss_modes_goal_rows.rows(); ++c) {
    Eigen::MatrixXi e_ss_modes_goal_c = e_ss_modes_goal_data.middleRows(count, e_ss_modes_goal_rows(c));
    count += e_ss_modes_goal_rows(c);
    e_ss_modes_goal.push_back(e_ss_modes_goal_c);
  }


  filename = std::string("data/h_ss_modes_goal_rows.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_i = arr.data<int>();
  rows = arr.shape[0];
  cols = arr.shape[1];
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  Eigen::VectorXi h_ss_modes_goal_rows(cols);
  for (int i = 0; i < cols; ++i) {
    h_ss_modes_goal_rows(i) = data_ptr_i[i];
  }
  std::cout << "h_ss_modes_goal_rows:\n" << h_ss_modes_goal_rows << std::endl;
  filename = std::string("data/h_ss_modes_goal_data.npy");
  std::cout << "\nreading: " << filename << std::endl;
  arr = cnpy::npy_load(filename);
  data_ptr_i = arr.data<int>();
  rows = arr.shape[0];
  cols = arr.shape[1];
  std::cout << "arr.shape[0]: " << rows << ", arr.shape[1]: " << cols << std::endl;
  Eigen::MatrixXi h_ss_modes_goal_data(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      h_ss_modes_goal_data(i,j) = data_ptr_i[j];
    }
    data_ptr_i += cols;
  }
  std::cout << "h_ss_modes_goal_data:\n" << h_ss_modes_goal_data << std::endl;
  count = 0;
  for (int c = 0; c < h_ss_modes_goal_rows.rows(); ++c) {
    Eigen::MatrixXi h_ss_modes_goal_c = h_ss_modes_goal_data.middleRows(count, h_ss_modes_goal_rows(c));
    count += h_ss_modes_goal_rows(c);
    h_ss_modes_goal.push_back(h_ss_modes_goal_c);
  }

  WrenchSpaceAnalysis wsa;
  wsa.wrenchStamping(Jac_e, Jac_h, eCone_allFix, hCone_allFix,
      F_G, kContactForce, kFrictionE, kFrictionH, kCharacteristicLength, kNumSlidingPlanes,
      e_cs_modes, e_ss_modes, h_cs_modes, h_ss_modes, G, b_G,
      e_cs_modes_goal, e_ss_modes_goal, h_cs_modes_goal, h_ss_modes_goal, print_level);
}

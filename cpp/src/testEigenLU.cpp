#include <iostream>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/LU>

using namespace Eigen;

int main() {

// MatrixXd A(3,5);
// A << 4, 2, 1, 4, 5,
//      4, 2, 1, 4, 5,
//      4, 2, 1, 4, 5;
MatrixXd A(MatrixXd::Random(17,10));
A.setRandom();
// A.middleCols(12, 1) = A.middleCols(11, 1);
std::cout << "A:\n" << A << std::endl;

// FullPivLU
Eigen::FullPivLU<MatrixXd> lu_full(A);
std::cout << "FullPivLU" << std::endl;
std::cout << "  col space:\n" << lu_full.image(A) << std::endl;
std::cout << "  kernel:\n" << lu_full.kernel() << std::endl;
std::cout << "  kernel'*kernel:\n" << lu_full.kernel().transpose() * lu_full.kernel() << std::endl;

// QR
HouseholderQR<MatrixXd> qr(A);
MatrixXd Q = qr.householderQ();
std::cout << "Q:\n" << Q << "\n";
std::cout << "Q^T * A:\n" << Q.transpose() * A << "\n";
std::cout << "qr.matrixQR():\n" << qr.matrixQR() << "\n";

// see how Identity works
MatrixXd thinQ(MatrixXd::Identity(5,3));
std::cout << "thinQ:" << thinQ << "\n";


return 0;
}
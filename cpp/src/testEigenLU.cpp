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
MatrixXd A(MatrixXd::Random(4,3));
A.setRandom();
// A.middleCols(12, 1) = A.middleCols(11, 1);
std::cout << "A:\n" << A << std::endl;

VectorXd b(4, 1);
b.setRandom();
std::cout << "b:\n" << b << std::endl;
// FullPivLU
Eigen::FullPivLU<MatrixXd> lu_full(A);
std::cout << "FullPivLU" << std::endl;
std::cout << "  col space:\n" << lu_full.image(A) << std::endl;
std::cout << "  kernel:\n" << lu_full.kernel() << std::endl;
std::cout << "  kernel'*kernel:\n" << lu_full.kernel().transpose() * lu_full.kernel() << std::endl;
VectorXd sol = lu_full.solve(b);
std::cout << "  solve:\n" << sol << std::endl;
std::cout << "  A*sol - b:\n" << A*sol - b << std::endl;


// QR
HouseholderQR<MatrixXd> qr(A);
MatrixXd Q = qr.householderQ();
std::cout << "Q:\n" << Q << "\n";
std::cout << "Q^T * A:\n" << Q.transpose() * A << "\n";
std::cout << "qr.matrixQR():\n" << qr.matrixQR() << "\n";

// see how Identity works
MatrixXd identity(MatrixXd::Identity(5,3));
std::cout << "identity:" << identity << "\n";


return 0;
}
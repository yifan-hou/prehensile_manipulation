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
// MatrixXd A(MatrixXd::Random(4,3));
MatrixXd A(MatrixXd::Random(17,12));
A.setRandom();
A.middleCols(11, 1) = A.middleCols(10, 1);
std::cout << "A:\n" << A << std::endl;

VectorXd b(17, 1);
b.setRandom();
std::cout << "b:\n" << b << std::endl;
// FullPivLU
std::cout << "\n\nFullPivLU test";
Eigen::FullPivLU<MatrixXd> lu_full(A);
std::cout << "  col space:\n" << lu_full.image(A) << std::endl;
std::cout << "  kernel:\n" << lu_full.kernel() << std::endl;
std::cout << "  kernel'*kernel:\n" << lu_full.kernel().transpose() * lu_full.kernel() << std::endl;
VectorXd sol = lu_full.solve(b);
std::cout << "  solve:\n" << sol << std::endl;
std::cout << "  A*sol - b:\n" << A*sol - b << std::endl;


// QR
std::cout << "\n\nQR test";
HouseholderQR<MatrixXd> qr(A);
MatrixXd Q = qr.householderQ();
std::cout << "Q:\n" << Q << "\n";
std::cout << "Q^T * A:\n" << Q.transpose() * A << "\n";
std::cout << "qr.matrixQR():\n" << qr.matrixQR() << "\n";

// Orthogonal decomposition
std::cout << "\n\nCompleteOrthogonalDecomposition test";
Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(A);
Eigen::MatrixXd image = cod.matrixQ();
std::cout << "  image:\n" << image << std::endl;
std::cout << "  image'*image:\n" << image.transpose() * image << std::endl;
std::cout << "  image'*A:\n" << image.transpose() * A << std::endl;
std::cout << "  ||image'*image - I||:\n" << (image.transpose() * image - MatrixXd::Identity(image.rows(), image.rows())).norm() << std::endl;
std::cout << "  image'*A, norm of rows:\n";
Eigen::MatrixXd Aproj = image.transpose() * A;
for (int i = 0; i < Aproj.rows(); ++i) {
    std::cout << i << ":\t" << Aproj.middleRows(i, 1).norm() << std::endl;
}
std::cout << std::endl;

// test threshold in rank()
Eigen::MatrixXd B(2,3);
B << 1, 0, 0,
     1, 1e-15, 0;
lu_full.compute(B);
std::cout << "B:\n" << B << std::endl;
std::cout << "lu rank: " << lu_full.rank() << std::endl;
lu_full.setThreshold(1e-8);
std::cout << "lu rank 2: " << lu_full.rank() << std::endl;
Eigen::ColPivHouseholderQR<MatrixXd> qr_col;
qr_col.compute(B);
std::cout << "qr rank: " << qr_col.rank() << std::endl;
qr_col.setThreshold(1e-8);
qr_col.compute(B);
std::cout << "qr rank 2: " << qr_col.rank() << std::endl;
cod.compute(B);
std::cout << "cod rank: " << cod.rank() << std::endl;

return 0;
}
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <iostream>
#include <Eigen/Dense>
namespace py = pybind11;

using namespace Eigen;

int mtimes(MatrixXd m1, MatrixXd m2) {
    std::cout << "m1: " << m1.rows() << " by " << m1.cols() << std::endl;
    std::cout << "m2: " << m2.rows() << " by " << m2.cols() << std::endl;
    return 0;
}

void manipulate(Eigen::Ref<Eigen::MatrixXd> data) {
    data = data*2;
}

PYBIND11_MODULE(example, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("mtimes", &mtimes, "A function which adds two numbers");
    m.def("manipulate", &manipulate, "A function which adds two numbers");
    m.def("diagonal", [](const Eigen::Ref<const Eigen::MatrixXd> &x) { return x.diagonal(); });
}

#include <Eigen/Dense>

namespace Poly {

/**
 * Vertex enumeration. Given a polyhedron (A, b):
 *      {x: A*x <= b},
 * Find its generator representation R = [m M]:
 *      {y: y = M_ray'*x1 + M_vertex'*x2, x1 >= 0, 0 <= x2 <= 1}
 * where M_ray and M_vertex are rows of M whose corresponding m equals 0 and 1.
 *
 *   Examples:
 *       A = [0 1 0; 1 0 0]; b = [0; 1];
 *       R = VertexEnumeration(A, b, &R);
 *
 * @param[in]  A     Each row represents an inequality.
 * @param[in]  b     Each row represents an inequality.
 * @param      R     Pointer to output R, each row of R represents a generator.
 *
 * @return     True if no error occurs.
 */
bool vertexEnumeration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, Eigen::MatrixXd *R);

/**
 * Facet enumeration. Given a polyhedron R = [m M]:
 *      {y: y = M_ray'*x1 + M_vertex'*x2, x1 >= 0, 0 <= x2 <= 1}
 * where M_ray and M_vertex are rows of M whose corresponding m equals 0 and 1;
 * Find its inequality representation (A, b):
 *      {x: A*x <= b},
 *
 * @param[in]  R     Each row represents a generator. First element of each row represents the type of generator (0: ray 1: vertex)
 * @param      A     Pointer of output A, each row represents an inequality.
 * @param      b     Pointer of output b, each row represents an inequality.
 *
 * @return     True if no error occurs.
 */
bool facetEnumeration(const Eigen::MatrixXd &R, Eigen::MatrixXd *A, Eigen::VectorXd *b);

/**
 * Compute the intersection of two polyhedra. Both input and output variables
 * are generator representations:  R = [m M]:
 *      {y: y = M_ray'*x1 + M_vertex'*x2, x1 >= 0, 0 <= x2 <= 1}
 * where M_ray and M_vertex are rows of M whose corresponding m equals 0 and 1;
 *
 * @param[in]  R1    Input polyhedron 1
 * @param[in]  R2    Input polyhedron 2
 * @param      R     The intersection of R1 and R2
 *
 * @return     True if no error occurs.
 */
bool intersection(const Eigen::MatrixXd &R1, const Eigen::MatrixXd &R2, Eigen::MatrixXd *R);

/**
 * Compute the intersection of two homogeneous cones. Both input and output
 * variables are generator representations: C
 *      {y: y = C'*x, x >= 0}
 *
 * @return     True if no error occurs.
 */
bool coneIntersection(const Eigen::MatrixXd &C1, const Eigen::MatrixXd &C2, Eigen::MatrixXd *C);

bool minkowskiSum(const Eigen::MatrixXd)

}

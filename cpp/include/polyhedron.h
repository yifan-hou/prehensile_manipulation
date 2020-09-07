#include <vector>
#include <Eigen/Dense>
#include <ppl.hh> // use PPL instead of CDD

#define KCONVHULL_ROUND 3
#define PPL_MULTIPLIER_BEFORE_ROUNDING 10000 // floating point will get round to integer after multiplied with this number

namespace PPL = Parma_Polyhedra_Library;

namespace Poly {

double angBTVec(const Eigen::VectorXd &x, const Eigen::VectorXd &b);

/**
 * Distance from a point p to the line {kn|k\in R}
 *
 * @param[in]  p     a point in N dimensional space
 * @param[in]  n     a ray in N dimensional space
 *
 * @return     Distance from p to the line of n. Value is negative if p.dot(n) is negative
 */
double distP2Line(const Eigen::VectorXd &p, const Eigen::VectorXd &n);

/**
 * Project a point to a line
 *
 * @param[in]  p     a point
 * @param[in]  n     Direction vector of the line, don't have to be unit
 *
 * @return     The projection.
 */
Eigen::VectorXd projectP2Line(const Eigen::VectorXd &p, const Eigen::VectorXd &n);

/**
* Compute the projection of a point p onto a hyperplane defined by
*      ax = b
*
* @param[in]  p     a point in n-D space
* @param[in]  a     a co-vector in n-D space
* @param[in]  b     a scalar
*
* @return     A point on the hyperplane
*/
Eigen::VectorXd projectP2Hyperplane(const Eigen::VectorXd &p, const Eigen::VectorXd &a, double b);

/**
 * Compute the angular distance between a ray and a cone. Both the ray and the
 * cone start from the origin (are homogeneous). The vector must not be in the
 * relative interior of the cone.
 *
 * @param[in]  p     a non-zero point on the ray
 * @param[in]  A     Hyperplane representation of the cone: Ax <= 0
 * @param[in]  R     Vertex representation of the cone, each row of R is a
 *                   generator
 *
 * @return     The angle
 */
double distRay2ConeFromOutside(const Eigen::VectorXd &p, const Eigen::MatrixXd &A, const Eigen::MatrixXd &R);

/**
 * Compute the distance between a point and a polyhedron using QP.
 *
 * @param[in]  p          The query point
 * @param[in]  A          Polyhedron description, A x <= b
 * @param[in]  b          Polyhedron description, A x <= b
 * @param[in]  x0         Initial guess for QP
 * @param      x_closest  The point in the polyhedron that is the closest to p
 *
 * @return     { description_of_the_return_value }
 */
double distP2Polyhedron(const Eigen::VectorXd &p, const Eigen::MatrixXd &A,
    const Eigen::VectorXd &b, const Eigen::VectorXd &x0, Eigen::VectorXd *x_closest = nullptr);

/**
 * Finds a point in P1 that is furthest away from P2. If P2 is empty, returns the
 * center of P1.
 *
 * @param[in]  A2    A 2
 * @param[in]  b2    The b 2
 * @param[in]  A1    A 1
 * @param[in]  b1    The b 1
 * @param[in]  xl    { parameter_description }
 * @param[in]  xu    { parameter_description }
 * @param      x     { parameter_description }
 *
 * @return     The away from polyhedrons.
 */
double getAwayFromPolyhedrons(
    const std::vector<Eigen::MatrixXd> &A2, const std::vector<Eigen::VectorXd> &b2,
    const Eigen::MatrixXd &A1, const Eigen::VectorXd &b1,
    const Eigen::VectorXd &xl, const Eigen::VectorXd &xu,
    Eigen::VectorXd *x);

/**
 * Hit and Run sampling in a polytope. It converges to a uniform sampling
 *
 * @param[in]  A        Polytope description, Ax <= b.
 * @param[in]  b        Polytope description, Ax <= b.
 * @param[in]  x0       A point inside of the polytope
 * @param[in]  No       Number of samples to output.
 * @param[in]  discard  The number of first discard samples are discarded
 *
 * @return     { Return the samples }
 */
Eigen::MatrixXd hitAndRunSampleInPolytope(const Eigen::MatrixXd &A,
        const Eigen::VectorXd &b, const Eigen::VectorXd &x0, int N, int discard = 10, int runup = 10, double max_radius = -1);

/**
 * Sample points in Polyhedron 1, outside of Polyhedron 2. This function samples
 * N random direction from a given initial point x0, then find the largest
 * feasible segment along each of those directions, return their center points.
 * There is one P1 and multiple P2.
 *
 * @param[in]  A1          Describes P1: A1 x <= b1. Rows must be normalized
 * @param[in]  b1          Describes P1: A1 x <= b1
 * @param[in]  A2          Describes P2: A2[i] x <= b2[i]. Rows must be normalized
 * @param[in]  b2          Describes P2: A2[i] x <= b2[i]
 * @param[in]  x0          The initial point. Must be within P1
 * @param[in]  N           Number of directions to try
 * @param[in]  max_radius  The maximum radius
 *
 * @return     At most N points
 */
std::vector<Eigen::VectorXd> sampleInP1OutOfP2(const Eigen::MatrixXd &A1,
    const Eigen::VectorXd &b1, const std::vector<Eigen::MatrixXd> &A2,
    const std::vector<Eigen::VectorXd> &b2, const Eigen::VectorXd &x0,
    int N, double max_radius = -1);

// 0: ray, 1: point, 2: line
bool constructPPLGeneratorsFromV(const Eigen::MatrixXd &R_input,
    PPL::Generator_System *gs);

bool constructPPLConstraintsFromH(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
    const Eigen::MatrixXd &Ae, const Eigen::VectorXd &be, PPL::Constraint_System *cs);

// 0: ray, 1: point, 2: line
bool constructPPLPolyFromV(const Eigen::MatrixXd &R_input, PPL::C_Polyhedron *ph);

bool constructPPLPolyFromH(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
    const Eigen::MatrixXd &Ae, const Eigen::VectorXd &be, PPL::C_Polyhedron *ph);

// 0: ray, 1: point, (Todo): implement 2: line
bool getVertexFromPPL(const PPL::C_Polyhedron &ph, Eigen::MatrixXd *R);

bool getFacetFromPPL(const PPL::C_Polyhedron &ph, Eigen::MatrixXd *A, Eigen::VectorXd *b);
/**
 * Vertex enumeration. Given a polyhedron (A, b, Ae, be):
 *      {x: A*x <= b, Ae*x = be},
 * Find its generator representation R = [m M]:
 *      {y: y = M_ray'*x1 + M_vertex'*x2, x1 >= 0, 0 <= x2 <= 1}
 * where M_ray and M_vertex are rows of M whose corresponding m equals 0 and 1.
 *
 *   Examples:
 *       A = [0 1 0; 1 0 0]; b = [0; 1];
 *       R = VertexEnumeration(A, b, Eigen::MatrixXd::Zero(0,0), Eigen::VectorXd::Zero(0), &R);
 *
 * @param[in]  A     Each row represents an inequality. It's good for numerical stability to normalize these rows
 * @param[in]  b     Each row represents an inequality.
 * @param[in]  Ae    Each row represents an equality. It's good for numerical stability to normalize these rows
 * @param[in]  be    Each row represents an equality.
 * @param      R     Pointer to output R, each row of R represents a generator.
 *
 * @return     True if no error occurs.
 */



bool vertexEnumeration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
        const Eigen::MatrixXd &Ae, const Eigen::VectorXd &be, Eigen::MatrixXd *R);
bool vertexEnumeration(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,  Eigen::MatrixXd *R);

/**
 * Facet enumeration. Given a polyhedron R = [m M]:
 *      {y: y = M_ray'*x1 + M_vertex'*x2, x1 >= 0, 0 <= x2 <= 1}
 * where M_ray and M_vertex are rows of M whose corresponding m equals 0 and 1;
 * Find its inequality representation (A, b):
 *      {x: A*x <= b},
 *
 * @param[in]  R     Each row represents a generator. First element: (0: ray 1: vertex)
 * @param      A     Pointer of output A, each row represents an inequality.
 * @param      b     Pointer of output b, each row represents an inequality.
 *
 * @return     True if no error occurs.
 */
bool facetEnumeration(const Eigen::MatrixXd &R, Eigen::MatrixXd *A, Eigen::VectorXd *b);

/**
 * Facet enumeration for a homogeneous cone. Given a polyhedron M:
 *      {y: y = M'*x, x >= 0}
 * Find its inequality representation (A):
 *      {x: A*x <= 0},
 *
 * @param[in]  M     Each row represents a generator
 * @param      A     Pointer of output A, each row represents an inequality.
 *
 * @return     True if no error occurs.
 */
bool coneFacetEnumeration(const Eigen::MatrixXd &M, Eigen::MatrixXd *A);

/**
 * Facet enumeration for a polytope M:
 *      {y: y = M'*x, 0 <= x <= 1}
 * Find its inequality representation (A, b):
 *      {x: A*x <= b},
 *
 * @param[in]  M     Each row represents a vertex
 * @param      A     Pointer of output A, each row represents an inequality.
 * @param      b     Pointer of output b, each row represents an inequality.
 *
 * @return     True if no error occurs.
 */
bool polytopeFacetEnumeration(const Eigen::MatrixXd &M, Eigen::MatrixXd *A, Eigen::VectorXd *b);

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
/**
 * Compute the intersection of two polytopes. Both input and output
 * variables are generator representations: P
 *      {y: y = P'*x, x \in (0, 1)}
 *
 * @return     True if no error occurs.
 */
bool polytopeIntersection(const Eigen::MatrixXd &P1, const Eigen::MatrixXd &P2, Eigen::MatrixXd *P);


/**
 * Offset a polytope by a vector.
 *
 * @param      polytope  Each row denotes a vertex of the polytope.
 * @param[in]  offset    The offset vector.
 *
 * @return     True if no error occurs.
 */
bool offsetPolytope(Eigen::MatrixXd *polytope, const Eigen::VectorXd &offset);

/**
 * Compute the convex hull of points.
 *
 * @param[in]  points  Each row denotes a point
 *
 * @return     The convex hull. Each row denotes a point
 */
Eigen::MatrixXd convhull(const Eigen::MatrixXd &points);
/**
 * Compute the convex hull of points
 *
 * @param[in]  vectors  Input, a vector of concatenated points
 * @param[in]  dim      The dimension of the points
 * @param[in]  num      The number of points
 * @param      results  Output, a vector of concatenated points
 *
 * @return     True if no error
 */
bool convhull(const std::vector<double> &vectors, int dim, int num, std::vector<double> *results);
/**
 * Compute the center of the largest inscribed sphere to the polytope.
 * Ax <= b may not be bounded. xu, xl boundary and A1, b1 only constraint the
 * sphere center, not the sphere itself.
 *
 * @param[in]  A     Ax <= b
 * @param[in]  b     Ax <= b
 * @param[in]  xl    xl <= x <= xu
 * @param[in]  xu    xl <= x <= xu
 * @param      xc    The found center
 *
 * @return     Radius of the found inscribed sphere; -1 if no solution
 */
double inscribedSphere(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
    const Eigen::VectorXd &xl, const Eigen::VectorXd &xu,
    const Eigen::MatrixXd &A1, const Eigen::VectorXd &b1, Eigen::VectorXd *xc);

/**
 * Compute the Minkowski sum of a bunch of vectors.
 *
 * @param[in]  vectors  Each row represents a vector (the end point. The
 *                      starting point is assumed to be the origin)
 * @param      results  The polytope generated by the Minkowski sum
 *
 * @return     True if no error
 */
bool minkowskiSumOfVectors(const Eigen::MatrixXd &vectors, Eigen::MatrixXd *results);
bool minkowskiSum(const Eigen::MatrixXd &poly1, const Eigen::MatrixXd &poly2, Eigen::MatrixXd *results);

/**
 * Check if a system of linear inequalities
 *   Ax <= b
 * has a solution.
 * TODO: write a full LP
 *
 * @return     True if a solution exists
 */
bool lpfeasibility(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
    Eigen::VectorXd *xs);

/**
 * Linear programming.
 *   min C'x
 *   s.t. Ax <= b
 *        Ae x == be
 *        xl <= x <= xu
 * Elements of xl and xu can be inf or NaN to indicate no constraint. Only
 * when isfinite() return true would the bound be considered.
 *
 * @param[in]  C             cost to be minimized
 * @param[in]  A             Inequality constraints, can be empty
 * @param[in]  b
 * @param[in]  Ae            Equality constraints, can be empty
 * @param[in]  be
 * @param[in]  xl            Lower bound of variable, can be empty.
 * @param[in]  xu            Upper bound of variable, can be empty.
 * @param      xs            Stores the solution
 * @param      optimal_cost  The optimal cost
 *
 * @return     True if the problem is feasible
 */
bool lp(const Eigen::VectorXd &C, const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
    const Eigen::MatrixXd &Ae, const Eigen::VectorXd &be,
    const Eigen::VectorXd &xl, const Eigen::VectorXd &xu, Eigen::VectorXd *xs,
    double *optimal_cost);
}
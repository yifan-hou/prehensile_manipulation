#include <iostream>
#include <stdio.h>            /* C input/output                       */
#include <stdlib.h>           /* C standard library                   */
#include <glpk.h>             /* GNU GLPK linear/mixed integer solver */
#include <Eigen/Dense>


int main(void)
{

  Eigen::MatrixXd A(4, 2);
  Eigen::VectorXd b(4);
  A << 1, 0,
       0, 1,
       -1, 0,
       0, -1;
  b << -1, 2, 0, 0;

  /* declare variables */
  glp_prob *lp;
  glp_smcp parm;
  glp_init_smcp(&parm);
  parm.presolve = GLP_OFF;
  parm.msg_lev = GLP_MSG_ERR; // error and warning only

  int *ia, *ja;
  double *ar;
  int rows = A.rows();
  int cols = A.cols();
  int size = rows * cols;
  ia = new int[size + 1000];
  ja = new int[size + 1000];
  ar = new double[size + 1000];
  double z, x1, x2;

  /* create problem */
  lp = glp_create_prob();
  glp_set_prob_name(lp, "short");
  glp_set_obj_dir(lp, GLP_MAX);

  /* fill problem */
  glp_add_rows(lp, rows);
  for (int r = 1; r <= rows; ++r) {
    glp_set_row_bnds(lp, r, GLP_UP, 0.0, b(r-1));
  }

  // set variable bounds and cost function
  glp_add_cols(lp, cols);
  for (int c = 1; c <= cols; ++c) {
    glp_set_col_bnds(lp, c, GLP_FR, 0.0, 0.0); // no boundary
    glp_set_obj_coef(lp, c, 0.0); // no cost function, so coef = 0
  }

  // fill in coefficient matrix
  int id = 0;
  for (int r = 1; r <= rows; ++r) {
    for (int c = 1; c <= cols; ++c) {
      id = (r-1)*cols + c;
      ia[id] = r, ja[id] = c, ar[id] = A(r-1, c-1);
    }
  }
  glp_load_matrix(lp, id, ia, ja, ar);
  /* solve problem */
  glp_simplex(lp, &parm);
  int result = glp_get_status(lp);
  std::cout << "Result code: " << result << std::endl;
  std::cout << "GLP_OPT: " << GLP_OPT << std::endl;
  std::cout << "GLP_FEAS: " << GLP_FEAS << std::endl;
  std::cout << "GLP_INFEAS: " << GLP_INFEAS << std::endl;
  std::cout << "GLP_NOFEAS: " << GLP_NOFEAS << std::endl;
  std::cout << "GLP_UNBND: " << GLP_UNBND << std::endl;
  std::cout << "GLP_UNDEF: " << GLP_UNDEF << std::endl;

  if ((result == GLP_OPT) || (result == GLP_FEAS)) {
    // feasible
    /* recover and display results */
    z = glp_get_obj_val(lp);
    Eigen::VectorXd xs(rows);
    for (int d = 0; d < cols; ++d) {
      xs(d) = glp_get_col_prim(lp, d + 1);
    }
    std::cout << "z: " << z << std::endl;
    std::cout << "solution: " << xs.transpose() << std::endl;
  } else {
    std::cout << "Infeasible." << std::endl;
  }

  /* housekeeping */
  glp_delete_prob(lp);
  glp_free_env();
  delete [] ia;
  delete [] ja;
  delete [] ar;

  return 0;
}
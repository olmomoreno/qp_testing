/*
CQP is an implementation of the Mehrotra's Predictor-Corrector
interior point method for the convex quadratic programming.
CQP solves problems of the form:

  min (x^t)Qx+(q^t)x
  s.t. Ax = b,
       Cx >= d.

  Dimensions of the problem data:
  Q is a (n x n)-matrix
  A is a (me x n)-matrix
  C is a (mi x n)-matrix
  q,b,d are vectors of appropriate dimensions.
  http://www.network-theory.co.uk/download/gslextras/Bundle/CQP-1.2/
*/

#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl_cqp.h>

using namespace std;

int s_GSL(){

    gsl_cqp_data *p_problem= new gsl_cqp_data;
    gsl_cqp_data problem;

    problem.Q = gsl_matrix_alloc (6,6);
    problem.q = gsl_vector_alloc (6);
    problem.A = gsl_matrix_alloc (1,6);
    problem.b = gsl_vector_calloc(1);
    problem.C = gsl_matrix_calloc(1,6);
    problem.d = gsl_vector_calloc(1);
    /*M*/
    gsl_matrix_set(problem.Q,0,0,0.0627+1e-9); gsl_matrix_set(problem.Q,0,1,-0.0183); gsl_matrix_set(problem.Q,0,2,-0.0603); gsl_matrix_set(problem.Q,0,3,0.0); gsl_matrix_set(problem.Q,0,4,0.0); gsl_matrix_set(problem.Q,0,5,0.0);
    gsl_matrix_set(problem.Q,1,0,-0.0183); gsl_matrix_set(problem.Q,1,1,0.068+1e-9); gsl_matrix_set(problem.Q,1,2,0.0176); gsl_matrix_set(problem.Q,1,3,-0.0183); gsl_matrix_set(problem.Q,1,4,-0.0603); gsl_matrix_set(problem.Q,1,5,0.0);
    gsl_matrix_set(problem.Q,2,0,-0.0603); gsl_matrix_set(problem.Q,2,1,0.0176); gsl_matrix_set(problem.Q,2,2,0.1207+1e-9); gsl_matrix_set(problem.Q,2,3,0.0); gsl_matrix_set(problem.Q,2,4,-0.0183); gsl_matrix_set(problem.Q,2,5,-0.0603);
    gsl_matrix_set(problem.Q,3,0,0.0); gsl_matrix_set(problem.Q,3,1,-0.0183); gsl_matrix_set(problem.Q,3,2,0.0); gsl_matrix_set(problem.Q,3,3,0.053+1e-9); gsl_matrix_set(problem.Q,3,4,0.0176); gsl_matrix_set(problem.Q,3,5,0.0);
    gsl_matrix_set(problem.Q,4,0,0.0); gsl_matrix_set(problem.Q,4,1,-0.0603); gsl_matrix_set(problem.Q,4,2,-0.0183); gsl_matrix_set(problem.Q,4,3,0.0176); gsl_matrix_set(problem.Q,4,4,0.0633+1e-9); gsl_matrix_set(problem.Q,4,5,0.0176);
    gsl_matrix_set(problem.Q,5,0,0.0); gsl_matrix_set(problem.Q,5,1,0.0); gsl_matrix_set(problem.Q,5,2,-0.0603); gsl_matrix_set(problem.Q,5,3,0.0); gsl_matrix_set(problem.Q,5,4,0.0176); gsl_matrix_set(problem.Q,5,5,0.0580+1e-9);
    /*q*/
    for (int i=0;i<=5;i++){
        gsl_vector_set(problem.q,i,0.0);
    }
    /*A*/
    for (int i =0;i<=5;i++){
        gsl_matrix_set(problem.A,0,i,0.0);
    }
    /*b*/
    gsl_vector_set(problem.b,0,0.0);
    /*CI*/
    gsl_matrix_set(problem.C,0,0,0.0627); gsl_matrix_set(problem.C,0,1,-0.0366); gsl_matrix_set(problem.C,0,2,0.1206);
    gsl_matrix_set(problem.C,0,3,0.0053); gsl_matrix_set(problem.C,0,4,0.0352); gsl_matrix_set(problem.C,0,5,0.0580);
    /*ci0*/
    gsl_vector_set(problem.d,0,0.6436);

    const size_t max_iter = 1000;
    size_t iter=1;
    const gsl_cqpminimizer_type * T;
    gsl_cqpminimizer *s;

    T = gsl_cqpminimizer_mg_pdip;
    //                                n         me          mi
    s = gsl_cqpminimizer_alloc(T, (size_t) 6,(size_t)  1,(size_t)   1);
    int status;
    status = gsl_cqpminimizer_set(s,&problem);
    status = gsl_cqpminimizer_iterate(s);

    printf("********************  Test  ********************\n\n");
    printf("== Itn ======= f ======== ||gap|| ==== ||residual||\n\n");

    do{
        status = gsl_cqpminimizer_iterate(s);
        status = gsl_cqpminimizer_test_convergence(s, 1e-15, 1e-15);
        printf("%4d   %14.8f  %13.6e  %13.6e\n", iter, gsl_cqpminimizer_f(s), gsl_cqpminimizer_gap(s), gsl_cqpminimizer_residuals_norm(s));
        if(status == GSL_SUCCESS){
            size_t j;
            printf("\nMinimum is found at\n");
            for(j=0; j<gsl_cqpminimizer_x(s)->size; j++)
                 printf("%9.6f ",gsl_vector_get(gsl_cqpminimizer_x(s), j));
             printf("\n\n");

             printf("\nLagrange-multipliers for Ax=b\n");
             for(j=0; j<gsl_cqpminimizer_lm_eq(s)->size; j++)
                 printf("%9.6f ",gsl_vector_get(gsl_cqpminimizer_lm_eq(s), j));
             printf("\n\n");

             printf("\nLagrange-multipliers for Cx>=d\n");
             for(j=0; j<gsl_cqpminimizer_lm_ineq(s)->size; j++)
                 printf("%9.6f ",gsl_vector_get(gsl_cqpminimizer_lm_ineq(s), j));
             printf("\n\n");
         }
         else{
             iter++;
         }

    }
    while(status == GSL_CONTINUE && iter<=max_iter);
    gsl_cqpminimizer_free(s);

    printf("End of it\n");
    return 0;
}

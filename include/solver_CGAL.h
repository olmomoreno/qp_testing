//http://doc.cgal.org/latest/QP_solver/index.html

#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/Gmpzf.h>
#include <CGAL/MP_Float.h>
using namespace std;
//typedef CGAL::Gmpzf ET;
typedef CGAL::MP_Float ET;
// program and solution types
typedef CGAL::Quadratic_program<double> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

void s_CGAL(){
        //QP Solver
        //Setting nonnegative QP with Ax >= b
        Program qp (CGAL::SMALLER, false, 0, false   , 0);

        //M
        //Setting D: Symmetric positive-semidefinite n×n matrix (the quadratic objective function)
        qp.set_d(0,0,0.0627 + 1.0e-9); //qp.set_d(0,1,-0.0183); qp.set_d(0,2,-0.0603); qp.set_d(0,3,0.0); qp.set_d(0,4,0.0); qp.set_d(0,5,0.0);
        qp.set_d(1,0,-0.0183); qp.set_d(1,1,0.068 + 1.0e-9); //qp.set_d(1,2,0.0176); qp.set_d(1,3,-0.0183); qp.set_d(1,4,-0.0603); qp.set_d(1,5,0.0);
        qp.set_d(2,0,-0.0603); qp.set_d(2,1,0.0176); qp.set_d(2,2,0.1207 + 1.0e-9); //qp.set_d(2,3,0.0); qp.set_d(2,4,-0.0183); qp.set_d(2,5,-0.0603);
        qp.set_d(3,0,0.0); qp.set_d(3,1,-0.0183); qp.set_d(3,2,0.0); qp.set_d(3,3,0.053 + 1.0e-9); //qp.set_d(3,4,0.0176); qp.set_d(3,5,0.0);
        qp.set_d(4,0,0.0); qp.set_d(4,1,-0.0603); qp.set_d(4,2,-0.0183); qp.set_d(4,3,0.0176); qp.set_d(4,4,0.0633 + 1.0e-9); //qp.set_d(4,5,0.0176);
        qp.set_d(5,0,0.0); qp.set_d(5,1,0.0); qp.set_d(5,2,-0.0603); qp.set_d(5,3,0.0); qp.set_d(5,4,0.0176); qp.set_d(5,5,0.0580 + 1.0e-9);

        //CI
        //Setting a: Constrains matrix
        qp.set_a(0,0,-0.0627); qp.set_a(0,1,0.0366); qp.set_a(0,2,-0.1206);
        qp.set_a(0,3,-0.0053); qp.set_a(0,4,-0.0352); qp.set_a(0,5,-0.0580);

        //ci0;
        //Setting b: Constrain vector
        qp.set_b(0,-0.6436);


        Solution s = CGAL::solve_quadratic_program(qp, ET());
    //    assert (s.solves_quadratic_program(qp));
        std::cout << s;

        double x[6];
        if(s.is_optimal()){
            int i=0;
            printf("N: %d\n",s.number_of_basic_variables());
            Solution::Variable_value_iterator it = s.variable_values_begin();
            Solution::Variable_value_iterator end = s.variable_values_end();
            for (; it != end; ++it) {
                x[i]=CGAL::to_double(*it);
                i++;
            }
        }
}

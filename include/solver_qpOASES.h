#include <iostream>
#include <qpOASES.hpp>
using namespace std;

    USING_NAMESPACE_QPOASES

    void s_qpOASES(){
    /* Setup data of first QP. */
    real_t H[6*6] = { 0.0627,   -0.0183,   -0.0603,         0,         0,         0,
                      -0.0183,    0.0680,    0.0176,   -0.0183,   -0.0603,         0,
                      -0.0603,    0.0176,    0.1207,         0,   -0.0183,   -0.0603,
                            0,   -0.0183,         0,    0.0053,    0.0176,         0,
                            0,   -0.0603,   -0.0183,    0.0176,    0.0633,    0.0176,
                            0,         0,   -0.0603,         0,    0.0176,    0.0580 };

    for (int i=0;i<36;i=i+7){
        H[i]+=1e-9;
    }

    real_t A[6*1] = { -0.0627,    0.0366,    0.1206,   -0.0053,   -0.0352,   -0.0580};
    real_t g[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    real_t lb[6] = { -200, -200.0, -200.0, -200.0, -200.0, -200.0};
    real_t ub[6] = { 200.0, 200.0, 200.0, 200.0, 200.0, 200.0};
    real_t lbA[1] = {0};
    real_t ubA[1] = {-0.6436};

/* Setting up QProblem object. */
//QProblem example( 3,0);
    QProblem example( 6,1, HST_SEMIDEF);
    Options options;
    example.setOptions( options );

    /* Solve first QP. */
    int_t nWSR = 100;
    example.init( H,g,A,lb,ub,NULL,ubA, nWSR );

    /* Get and print solution of first QP. */
    real_t xOpt[6];
    example.getPrimalSolution( xOpt );
    printf( "\nxOpt = [ %e, %e, %e, %e, %e, %e ];  objVal = %e\n\n",
            xOpt[0],xOpt[1],xOpt[2],xOpt[3],xOpt[4],xOpt[5],example.getObjVal() );

    example.printOptions();
    example.printProperties();
}

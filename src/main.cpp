#include <optimizer_all.h>

int main(){

    printf("| - - - - - Solver: MATLAB\n");
    printf("x[0]: 1.7515, x[1]: -0.7796, x[2]: -3.2520, x[3]: 0.0783, x[4]: 0.7157, x[5]: 1.5075, fv: 1.6431\n\n");

    printf("|---------- | Solver: qpOASES | ----------|\n");
    s_qpOASES();

    printf("|---------- | Solver: CGAL | ----------|\n");
    s_CGAL();

    printf("|---------- | Solver: GSL | ----------|\n");
    s_GSL();
    return 0;
}

///----------------------| MATLAB set-up |-----------------------------------------///
///Matlab Set-up
///f[0] = 1.3159; f[1] = -0.5082; f[2] = 1.0619;
///v[0] = -0.2504; v[1] = 0.0731; v[2] = 0.2408;
///desired_level = 1; tank_level = 0.2455;
/// T=eye(3);
/// Q=[v(1), v(2), v(3),0 ,0 ,0; 0, v(1), 0, v(2), v(3), 0;0, 0, v(1), 0, v(2), v(3)];
/// Qt = Q';
/// qdot = [v(1), v(2), v(3)]';
/// tau = [f(1), f(2), f(3)]';
/// M = Qt*T*Q+eye(6)*1e-9
/// CI = -(qdot'*Q)
/// ci0 = -[1-tl+tau'*qdot]
/// qp_options=optimoptions('quadprog','Algorithm','interior-point-convex','Display','iter-detailed');
/// [x, fv,exitflag]=quadprog(M,[],CI, ci0,[],[],[],[],[],qp_options)

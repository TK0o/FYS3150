#include "../header/function.hpp"

int main()
{
    int n;
    int N = 6;
    double a, d, h;

    set_initial_conditions(N, n, a, d, h, "N"); 
    arma::mat A = create_tri_mat(a, d, a, N);

    arma::vec arm_eig_val;
    arma::mat arm_eig_vec;
    solve_armadillo(A, arm_eig_vec, arm_eig_val);
    
    arma::vec ana_eig_val;
    arma::mat ana_eig_vec;    
    solve_analytic(ana_eig_vec, ana_eig_val, a, d, N, "0");

    compare_methods(arm_eig_vec, arm_eig_val, 
                    ana_eig_vec, ana_eig_val,
                    "armadillo", "analytical");
    return 0;
}

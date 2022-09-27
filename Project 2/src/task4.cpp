#include "../header/function.hpp"

int main()
{
    int n;
    int N = 6; 
    double a, d, h;
    int iterations;
    int maxiter = 5000;
    bool converged = false;
    double eps = std::pow(10, -8);
    
    set_initial_conditions(N, n, a, d, h, "N");

    arma::mat A = create_tri_mat(a, d, a, N);
    
    arma::vec eig_val;
    arma::mat eig_vec;
    jacobi_eigensolver(A, eps, eig_val, eig_vec, maxiter,
                       iterations, converged, "1", "1", "1");
    return 0;
}

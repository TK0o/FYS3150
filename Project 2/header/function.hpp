// Include guard
#ifndef __function_hpp__
#define __function_hpp__

#include <armadillo>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

/*-----------------------solver methods------------------------*/
void solve_armadillo(const arma::mat& A, arma::mat& eigenvec, 
                     arma::vec& eigenval);

void solve_analytic(arma::mat& eigenvec, arma::vec& eigenval, 
                    double a, double d, int N, std::string autosort);

double max_offdiag_symmetric(const arma::mat& A, int& k, int& l);

void jacobi_rotation(arma::mat& A, arma::mat& R, int k, int l);

void jacobi_eigensolver(const arma::mat A, double eps, arma::vec& eigenvalues,
        arma::mat& eigenvectors, const int maxiter, int& iterations,
        bool& converged, std::string compare_method_, std::string generic_info,
        std::string extra_info);


/*-----------------------specific task-------------------------*/
void run_scaling_data(int N_end, int maxiter, double eps);

void lowest3(arma::mat analytic_vec, arma::vec analytic_val,
             arma::mat jacobian_vec, arma::mat jacobian_val,
             arma::mat& data_holder);

void problem_6(int n_pow_to, int maxiter, double eps);


/*-------------------------utility-----------------------------*/
arma::mat create_tri_mat(double a_1, double d, double a_2, int N);

void set_initial_conditions(int& N, int& n, double& a, double& d, 
                            double& h, std::string based_on_n_or_N);

void compare_methods(const arma::mat& M1, const arma::vec& V1,
                     const arma::mat& M2, const arma::vec& V2,
                     std::string method1, std::string method2);

int print_statement(std::string print_result);

void test_max_offdiag_symmetric();
arma::mat test_matrix();

#endif
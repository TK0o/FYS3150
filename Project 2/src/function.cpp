#include "../header/function.hpp" 

void solve_armadillo(const arma::mat& A, arma::mat& eigenvec, 
                     arma::vec& eigenval)
{
    arma::vec eig_val;
    arma::eig_sym(eig_val, eigenvec, A);
    eigenval = arma::sort(eig_val);
}

void solve_analytic(arma::mat& eigenvec, arma::vec& eigenval, 
                    double a, double d, int N, std::string autosort)   
{
    double s = 0;
    arma::vec eig_val(N);
    arma::mat eig_vec(N, N);
    arma::vec temp_vec(N, arma::fill::zeros);

    for (double i = 1; i <= N; i++)
    {
        eig_val(i - 1) = d + 2 * a * std::cos(i * M_PI / (N + 1));        
        for (double j = 1; j <= N; j++)
        {
            temp_vec(j - 1) = std::sin(j * i * M_PI / (N + 1));  
            s += std::pow(temp_vec(j - 1), 2);
        }
        s = std::sqrt(s);
        for (int j = 1; j <= N; j++)
        {
            eig_vec(j - 1, i - 1) =  temp_vec(j - 1) / s;
        }
        s = 0;
    }
    eigenvec = eig_vec;
    if (print_statement(autosort))
    {
        eigenval = arma::sort(eig_val);
    }
    else
    {
        eigenval = eig_val;
    }
}

double max_offdiag_symmetric(const arma::mat& A, int& k, int& l)
{
    double max_val = 0;
    double temp_val;

    for (int i = 0; i < A.n_cols; i++)
    {
        for (int j = i + 1; j < A.n_cols; j++)
        {
            temp_val = std::abs(A(i, j));
            if (temp_val > max_val)
            {
                max_val = temp_val;
                k = i;
                l = j;
            }
        }  
    }
    return max_val;
}

void jacobi_rotation(arma::mat& A, arma::mat& R, int k, int l)
{
    int n_size = A.n_cols;
    double a_kk = A(k, k);
    double a_ll = A(l, l);
    double a_kl = A(k, l);
    double tau, c, s, t;

    if (a_kl != 0)
    {
        tau = (a_ll - a_kk) / (2 * a_kl);
        if (tau > 0)
        {
            t = 1 / (std::sqrt(1 + std::pow(tau, 2)) + tau);
        }
        else
        {
            t = -1 / (std::sqrt(1 + std::pow(tau, 2)) - tau);
        }
        c = 1 / std::sqrt(1 + std::pow(t, 2));
        s = c * t;
    }
    else
    {
        c = 1;
        s = 0;
        t = 0;
    }
    double c2 = std::pow(c, 2);
    double s2 = std::pow(s, 2);

    A(k, k) = a_kk * c2 - 2 * a_kl * c * s + a_ll * s2;
    A(l, l) = a_ll * c2 + 2 * a_kl * c * s + a_kk * s2;
    A(k, l) = 0;
    A(l, k) = 0;

    for (int i = 0; i < n_size; i++)
    {
        double r_ik = R(i, k);
        double r_il = R(i, l);
        R(i, k) = r_ik * c - r_il * s;
        R(i, l) = r_il * c + r_ik * s;

        if (i != k && i != l)
        {
            double a_ik = A(i, k);
            double a_il = A(i, l);
            
            A(i, k) = a_ik * c - a_il * s;
            A(k, i) = A(i, k);
            A(i, l) = a_il * c + a_ik * s;
            A(l, i) = A(i, l);
        }
    }
}

void jacobi_eigensolver(const arma::mat A, double eps, arma::vec& eigenvalues,
        arma::mat& eigenvectors, const int maxiter, int& iterations,
        bool& converged, std::string compare_method_, std::string generic_info,
        std::string extra_info)
{
    int k, l;
    double max_val;
    int N = A.n_cols;

    arma::mat A_copy = A;
    arma::mat R(N, N, arma::fill::eye);
    
    while (converged == false && iterations < maxiter)
    {
        max_val = max_offdiag_symmetric(A_copy, k, l);
        if (max_val < eps)
        {   
            converged = true;
            break;
        }
        jacobi_rotation(A_copy, R, k, l);
        iterations++;
    }

    arma::vec eigenval = A_copy.diag(0);
    arma::mat eigenvec = R.cols(arma::span::all);
    
    eigenvalues = eigenval;
    eigenvectors = eigenvec;

    if (print_statement(generic_info))
    {
        std::cout << "<----------general info---------->" << std::endl;
        std::cout << "eps = " << eps << std::endl;
        std::cout << "maxiter = " << maxiter << std::endl;
        std::cout << "max_val = " << std::setprecision(5) << max_val << '\n';
        std::cout << "amount of iterations = " << iterations << '\n';
        std::cout << "converged if the following value is (1 --> true) = "
                  << converged << std::endl;
        std::cout << "<------ end of general info------>\n" << std::endl;
    }

    if (print_statement(extra_info))
    {
        std::cout << "<----------extra info---------->" << std::endl;
        std::cout << "A_copy:\n" << A_copy << std::endl;
        std::cout << "eigenvectors:\n" << eigenvec << std::endl;
        std::cout << "eigenvalues:\n" << eigenval 
                  << "<------ end of extra info------>\n" << std::endl;
    }
    
    if (print_statement(compare_method_))
    {
        std::cout << "<-----comparison between jacobian and analytical----->"
                  << std::endl;
        arma::vec ana_eig_val;
        arma::mat ana_eig_vec; 
        solve_analytic(ana_eig_vec, ana_eig_val, 
                       A(1, 0), A(0, 0), N, "1");
        
        compare_methods(eigenvec, eigenval, 
                        ana_eig_vec, ana_eig_val,
                        "jacobian", "analytical");   
    }
}



void run_scaling_data(int N_end, int maxiter, double eps)
{
    int n;
    double a, d, h;
    arma::mat data(N_end, 2);
    
    std::ofstream MyFile("data/task5.txt");
    for (int N = 1; N <= N_end; N++)
    {
        int iterations = 0;
        bool converged = false;

        set_initial_conditions(N, n, a, d, h, "N");
        arma::mat A = create_tri_mat(a, d, a, N);
             
        arma::vec eig_val;
        arma::mat eig_vec;
        jacobi_eigensolver(A, eps, eig_val, eig_vec, maxiter,
                        iterations, converged, "0", "0", "0");

        if (iterations >= maxiter)
        {
            std::cout << "amount of iterations exceeds maxiter, at N = " 
                      << N << " \nincerase maxiter" << std::endl;
            std::cout << "stopping the loop" << std::endl;
            break;
        }
        data(N - 1, 1) = iterations;
        data(N - 1, 0) = N;  
    }
    MyFile << data << std::endl;
    MyFile.close();
}

void lowest3(arma::mat analytic_vec, arma::vec analytic_val,
             arma::mat jacobian_vec, arma::mat jacobian_val,
             arma::mat& data_holder)
{
    int N = analytic_val.n_rows;
    
    arma::vec x_hat = arma::linspace(0, 1, N + 2);
    arma::mat data(N + 2, 7, arma::fill::zeros);
    data.col(0) = x_hat;
    
    arma::uvec analytic_index = arma::sort_index(analytic_val);
    arma::uvec jacobian_index = arma::sort_index(jacobian_val);

    for (int i = 1; i < 4; i++)
    {
        int I = i - 1;
        data(arma::span(1, N), i) = analytic_vec.col(analytic_index(I));
        data(arma::span(1, N), i + 3) = jacobian_vec.col(jacobian_index(I));
    }
    data_holder = data;
}

void problem_6(int n_pow_to, int maxiter, double eps)
{
    int N;
    double a, d, h;
    
    for (int p = 1 ; p <= n_pow_to; p++)
    {
        int iterations;
        bool converged = false;
        int n = std::pow(10, p);
        
        std::string filename = "vectors_" + std::to_string(n) + ".txt";
        std::ofstream MyFile("data/" + filename);
        
        set_initial_conditions(N, n, a, d, h, "n");
        arma::mat A = create_tri_mat(a, d, a, N);
        arma::mat data_holder;
        
        arma::vec jac_eig_val;
        arma::mat jac_eig_vec;
        jacobi_eigensolver(A, eps, jac_eig_val, jac_eig_vec, maxiter,
                           iterations, converged, "0", "0", "0");

        arma::vec ana_eig_val;
        arma::mat ana_eig_vec;  
        solve_analytic(ana_eig_vec, ana_eig_val, a, d, N, "0");
        
        lowest3(ana_eig_vec, ana_eig_val, 
                jac_eig_vec, jac_eig_val, data_holder);
        
        MyFile << data_holder << std::endl;
        MyFile.close();
    }   
}



arma::mat create_tri_mat(double a_1, double d, double a_2, int N)
{
    arma::mat A = arma::mat(N, N, arma::fill::zeros);
    for (int i = 0; i < N; i++)
    {
        A(i, i) = d;
        if (i != 0)
        {
            A(i - 1, i) = a_1;
        }
        if (i != N - 1)
        {
            A(i + 1, i) = a_2;
        }
    }
    return A;
}

void set_initial_conditions(int& N, int& n, double& a, double& d, 
                            double& h, std::string based_on_n_or_N)
{
    if (based_on_n_or_N == "N")
    {   
        n = N + 1;
        h = 1 / double (n);
    }
    else if (based_on_n_or_N == "n")
    {
        N = n - 1;
        h = 1 / double (n);
    }
    else
    {
        std::cout << "the string based_on_n_or_N must be either " <<
                     "lower case 'n', or upper case 'N'" << std::endl;
    }
    d = 2 / std::pow(h, 2);
    a = -1 / std::pow(h, 2);
}

void compare_methods(const arma::mat& M1, const arma::vec& V1,
                     const arma::mat& M2, const arma::vec& V2,
                     std::string method1, std::string method2)
{
    std::string shortened_name_1;
    std::string shortened_name_2;
    for (int i = 0; i < 3; i++)
    {
        shortened_name_1 += method1[i];
        shortened_name_2 += method2[i];
    }

    std::cout << "\n";
    std::cout << "Shortened name for the two methods: " << std::endl;
    std::cout << method1 << " --> " << shortened_name_1 << std::endl;              
    std::cout << method2 << " --> " << shortened_name_2 << std::endl;
    std::cout << "\n";

    std::cout << shortened_name_1 << "_eigenval:" << std::setw(6) 
              << shortened_name_2 << "_eigenval:" << std::endl;

    for (int i = 0; i < V1.n_rows; i++)
    {
        std::cout << std::setprecision(4) << std:: scientific 
                  << std::setw(12) << V1(i) 
                  << std::setw(16) << V2(i) << "\n";     
    }
    std::cout << "\n";
    std::cout << shortened_name_1 << "_eigenvec:\n" << M1 << std::endl;
    std::cout << shortened_name_2 << "_eigenvec:\n" << M2 << std::endl;
}

int print_statement(std::string print_result)
{
    for (int i = 0; i < print_result.length(); i++)
    {
        print_result[i] = std::tolower(print_result[i]);
    }
    return 1 * (print_result == "1" || print_result == "true");
}

void test_max_offdiag_symmetric()
{
    int k, l;
    int N = 4;
    arma::mat A = test_matrix();
    double max_val = max_offdiag_symmetric(A, k, l);
    std::cout << "value of max_val = " << max_val << "\n"
         << "index of max_val in the matrix (k, l) = A(" 
         << k << ", " << l << ")" << std::endl;
}

arma::mat test_matrix()
{
    arma::mat A(4, 4, arma::fill::eye);
    arma::vec d2 = {0.5, -0.7, -0.7, 0.5};
    for (int i = 0; i < 4; i++)
    {
        A(3 - i, i) = d2(i);
    }
    return A;
}


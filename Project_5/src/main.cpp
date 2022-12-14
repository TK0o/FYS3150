#include "../header/readfile.hpp"

class Quantum_Box_Class
{
    public:
        arma::mat value_mat;
        arma::cx_mat A;
        arma::cx_mat B;
        arma::cx_mat b;
        arma::cx_mat V;
        arma::cx_mat U;
        arma::cx_mat res;
        arma::cx_mat temp_U;

        arma::cx_double _r_;

        std::vector<std::string> File_Names;

        double h;
        double x_c, y_c;
        double sigma_x, sigma_y;
        double p_x, p_y;
        double dt, t_stop;
        double v_0, wall_potential;
        double wall_thickness_x, wall_pos_centre_x;
        double slit_dist_y, slit_size_y, slit_centred_y;
        
        double nt;
        
        int n_slits;
        
        int M, N, n;
        int n_simulations_to_do;
        int current_simulation_id = 0;

        Quantum_Box_Class(std::string input_file)
        {
            readfile(value_mat, File_Names, input_file);
            check_default_values(value_mat);
            n_simulations_to_do = value_mat.size() / 17;
        }

        void create_wall()
        {             
            arma::vec wall_pos(n_slits);
            arma::vec y_val(n_slits + 2, arma::fill::zeros);
            y_val(n_slits + 1) = 1;

            for (int i = 0; i < n_slits; i++)
            {
                double slit_pos = (i - slit_centred_y * (n_slits - 1));
                double slit_span = (slit_size_y + slit_dist_y);
                wall_pos(i) = slit_pos * slit_span + slit_centred_y;
            }

            for (int i = 0; i < n_slits; i++)
            {
                y_val(i + 1) = wall_pos(i);
            }

            double s_size = slit_size_y / 2;
            double v = 1;
            double x, y;
            double x_min = 0;
            double x_max = 1;
            double y_min = 0;
            double y_max = 1;

            for (int i = 0; i < n_slits + 1; i++)
            {
                if (i != 0)
                {
                    y_min = y_val(i) + s_size;
                }
                else
                {
                    y_min = y_val(i);
                }

                if (i != n_slits)
                {
                    y_max = y_val(i + 1) - s_size;
                }
                else
                {
                    y_max = y_val(i + 1);
                }

                x_min = wall_pos_centre_x - wall_thickness_x / 2;
                x_max = wall_pos_centre_x + wall_thickness_x / 2;
       
                for (int j = 0; j < n; j++)
                {
                    for (int l = 0; l < n; l++)
                    {
                        x = double (l * h);
                        y = double (j * h);
 
                        int check = (x_min <= x) * (x <= x_max) *
                                    (y_min <= y) * (y <= y_max);
                        
                        if (check == 1)
                        {   
                            V(j, l) = arma::cx_double (wall_potential, 0);
                        }
                    }
                }
            }
        }

        void set_up_A_mat()
        {
            int n_ = n - 1;
            int n_r = N - n;
            arma::cx_mat A_(N, N, arma::fill::zeros);

            int s = 0;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    A_(s, s + (0 < ((s + 1) % n))) = -_r_;
                    A_(s, s - (0 < (s % n))) = -_r_;
                    
                    A_(s, s + n * (s < n_r)) = -_r_;
                    A_(s, s - n * (n_ < s)) = -_r_;
                    
                    A_(s, s) = arma::cx_double(1) + arma::cx_double(0, 4) 
                    * _r_ + arma::cx_double(0, dt * V(j, i).real() / 2);
                
                    s++;
                }
            }
            A = A_;
        }

        void set_up_B_mat()
        {
            int n_ = n - 1;
            int n_r = N - n;
            arma::cx_mat B_(N, N, arma::fill::zeros);
            int s = 0;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    B_(s, s + (0 < ((s + 1) % n))) = _r_;
                    B_(s, s - (0 < (s % n))) = _r_;

                    B_(s, s + n * (s < n_r)) = _r_;
                    B_(s, s - n * (n_ < s)) = _r_;

                    B_(s, s) = arma::cx_double(1) - arma::cx_double(0, 4) 
                    * _r_ - arma::cx_double(0, dt * V(j, i).real() / 2);      
                }                        
            }
            B = B_;
        }

        void initial_state()
        {
            arma::cx_mat U_(n, n, arma::fill::zeros);
            double X, Y;
            double real, imag;

            for (int j = 1; j < M - 1; j++)
            {
                for (int i = 1; i < M - 1; i++)
                {
                    X = double(i * h - x_c);
                    Y = double(j * h - y_c);

                    real = - (X * X) / (2 * std::pow(sigma_x, 2))
                           - (Y * Y) / (2 * std::pow(sigma_y, 2));
                    imag = (p_x * X) + (p_y * Y);        
                    U_(j - 1, i - 1) = std::exp(arma::cx_double(real, imag));
                }
            }
            U = U_ / std::sqrt(arma::accu(arma::conj(U_) % U_));
            for (int i = 0; i < n; i++)
            {
                res(0, arma::span(i * n, (i + 1) * n - 1)) = trans(U.col(i));
            }        
        }

        void mat_prod(int run_idx)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    res(run_idx + 1, j) += B(j, i) * res(run_idx, i);
                } 
            }
        }

        void solve_v2(arma::cx_mat A_mat, int run_idx)
        {
            int I;
            arma::cx_double diagonal_val, su_diag_val;
            for (int i = 0; i < N - 1; i++)
            {
                I = i + 1;

                res(run_idx + 1, i) /= A_mat(i, i);

                for (int j = i; j < N; j++)
                {
                    A_mat(i, j) = A_mat(i, j) / A_mat(i, i);
                }
                    
                
                for (int I_ = i + 1; I_ < N; I_++)
                {
                    if (A_mat(I_, i) != arma::cx_double(0, 0))
                    {
                        for (int j_ = i; j_ < N; j_++)
                        {
                            A_mat(I_, j_) = A_mat(I_, j_) - 
                                            A_mat(I_, i) * A_mat(i, j_);
                        }
                        res(run_idx + 1, I_) = res(run_idx + 1, I_) - 
                                A_mat(I_, i) * res(run_idx + 1, i);
                    }
                }
            }
            res(run_idx + 1, N - 1) /= A_mat(N - 1, N - 1);
            A_mat(N - 1, N - 1) /= A_mat(N - 1, N - 1);
            
            for (int i = 1; i < N; i++)
            {
                int N_ = N - i;    
                for (int I = (N_ - n) * (0 < N_ - n); I < N_; I++)
                {
                    res(run_idx + 1, I) = res(run_idx + 1, I) - A_mat(I, N_) 
                                        * res(run_idx + 1, N_);
                    A_mat(I, N_) = A_mat(I, N_) - A_mat(I, N_);
                }
            }
        }

        void simulate(std::string specific_filename)
        {
            std::ofstream MyFile("../data/" + specific_filename);
            for (int i = 0; i < nt - 1; i++)
            {
                this->mat_prod(i);
                this->solve_v2(A, i);
                // break;
                // this->solve(A, b, i);
                // temp_U.row(0) = b.row(0);
                // this->update_result(i);    
            }
            MyFile << res << std::endl;
            MyFile.close(); 
        }
        
        void execute_from_file()
        {
            for (int i = 0; i < n_simulations_to_do; i++)
            {
                this->setup_individual_simulations();
                this->simulate(File_Names[i]);
                current_simulation_id++;
            }
        }

        void setup_individual_simulations()
        {
            this->set_simulation_value_single_run(current_simulation_id);
            M = int(1 / h + 1);
            n = (M - 2);
            N = std::pow(n, 2);
            nt = t_stop / dt + 1;
            _r_ = arma::cx_double(0, dt / (std::pow(h, 2) * 2));
            arma::cx_mat Result(int(nt), N);
            res = Result.fill(0);
            arma::cx_mat single_run(1, N);
            b = single_run;
            arma::cx_mat potential(n, n, arma::fill::zeros); 
            potential.fill(v_0);
            V = potential;

            this->set_up_A_mat();
            this->set_up_B_mat();
            this->initial_state();
        }

        void set_simulation_value_single_run(int run_nr)
        {
            h = value_mat(run_nr, 0);
            x_c = value_mat(run_nr, 1);
            y_c = value_mat(run_nr, 2);
            sigma_x = value_mat(run_nr, 3);
            sigma_y = value_mat(run_nr, 4);
            p_x = value_mat(run_nr, 5);
            p_y = value_mat(run_nr, 6);
            dt = value_mat(run_nr, 7);
            t_stop = value_mat(run_nr, 8);
            v_0 = value_mat(run_nr, 9);
            wall_potential = value_mat(run_nr, 10);
            wall_thickness_x = value_mat(run_nr, 11);
            wall_pos_centre_x = value_mat(run_nr, 12);
            slit_dist_y = value_mat(run_nr, 13);
            slit_size_y = value_mat(run_nr, 14);
            slit_centred_y = value_mat(run_nr, 15);
            n_slits = value_mat(run_nr, 16);
            std::string filename = File_Names[run_nr];
        }        
};

int main ()
{ 
    Quantum_Box_Class qbc("input.txt");
    qbc.setup_individual_simulations();
    qbc.simulate();
    
    // qbc.execute_from_file();
    
    return 0;
}


/*
to be implemented if there is time for it

void dynamic_vec()
{
    std::vector<std::vector<int>> vec;
    for (int i = 0; i < 10; i++)
    {
        std::vector<int> temp_vec;
        for (int j = i; j < 10; j++)
        {
            temp_vec.push_back(j);
        }
        vec.push_back(temp_vec);
    }
}
*/

/* to be fixed
        void update_result(int run_idx)
        {
            for (int i = 0; i < n; i++)
            {
                res.row(run_idx) = temp_U;
            }
        }
        
        void solve(arma::cx_mat A_mat, arma::cx_mat B_mat, int run_idx)
        {
            int I;
            arma::cx_double diagonal_val, su_diag_val;
            for (int i = 0; i < N - 1; i++)
            {
                I = i + 1;

                B_mat(0, i) /= A_mat(i, i);

                for (int j = i; j < N; j++)
                {
                    A_mat(i, j) = A_mat(i, j) / A_mat(i, i);
                }
                    
                
                for (int I_ = i + 1; I_ < N; I_++)
                {
                    if (A_mat(I_, i) != arma::cx_double(0, 0))
                    {
                        // su_diag_val = A_mat(I_, i);
                        for (int j_ = i; j_ < N; j_++)
                        {
                            A_mat(I_, j_) = A_mat(I_, j_) - 
                                            A_mat(I_, i) * A_mat(i, j_);
                        }
                        B_mat(0, I_) = B_mat(0, I_) - 
                                       A_mat(I_, i) * B_mat(0, i);
                    }
                }
            }
            B_mat(0, N - 1) /= A_mat(N - 1, N - 1);
            A_mat(N - 1, N - 1) /= A_mat(N - 1, N - 1);
            
            for (int i = 1; i < N; i++)
            {
                int N_ = N - i;    
                for (int I = (N_ - n) * (0 < N_ - n); I < N_; I++)
                {
                    B_mat(0, I) = B_mat(0, I) - A_mat(I, N_) * B_mat(0, N_);
                    A_mat(I, N_) = A_mat(I, N_) - A_mat(I, N_);
                }
            }
           
            // temp_U.row(0) = B_mat.row(0);
            // temp_U = B_mat;
            res.row(run_idx) = B_mat;

        }
        */
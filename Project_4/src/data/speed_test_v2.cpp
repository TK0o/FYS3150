#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <chrono>
#include <omp.h>

class Speed_Class
{
    public:
        int L, Seed, N_flips;
        arma::vec Indexes;
        arma::mat result;
        arma::mat T;

    public:
        Speed_Class (int lattice_size, int seed)
        {
            L = lattice_size;
            Seed = seed;
            N_flips = std::pow(L, 2);
            
            arma::vec idx(L + 2, arma::fill::zeros);
            for (int i = 0; i < L; i++)
            {
                idx(i + 1) = i;
            }
            idx(0) = L - 1;
            idx(L) = 0;
            Indexes = idx;
        }

        arma::mat make_lattice()
        {   
            arma::arma_rng::set_seed(Seed);
            arma::mat lattice(L, L, arma::fill::randu);
            
            for (int i = 0; i < L; i++)
            {
                for (int j = 0; j < L; j++)
                {
                    lattice(i, j) = 1 - 2 * (lattice(i, j) < 0.5);
                }
            }
            return lattice;
        }

        double energy(arma::mat lattice)
        {
            double E = 0;
            int I, J, N = L - 1;
            for (int i = 0; i < L; i++)
            {
                I = i + L * ((N - i) / N) - 1;
                for (int j = 0; j < L; j++)
                {
                    J = (1 - j / N) * (j + 1);
                    E -= lattice(i, j) * (lattice(I, j) + lattice(i, J));                       
                }  
            }
            return E;
        }

        double get_mag(arma::mat lattice)
        {
            double S = arma::sum(arma::sum(lattice));
            return S;
        }

        arma::vec make_delta_E_vector(double T)
        {
            arma::vec dE_vec(5);
            for (int i = 0; i < 5; i++)
            {
                dE_vec(i) = std::exp(-(i * 4 - 8) / T);
            }
            return dE_vec;
        }

        void method_normal(arma::mat& lattice, int idx)
        {
            double dE; 
            int N = L - 1;
            int i, j;
            int I_up, I_down, J_left, J_right;
            
            std::mt19937 generator;
            std::uniform_real_distribution<> uniform(0, 1);
            std::uniform_int_distribution<int> uniform_dist(0, N);
            generator.seed(Seed + idx);

            for (int attempt = 0; attempt < N_flips; attempt++)
            {
                i = uniform_dist(generator);
                j = uniform_dist(generator);
                
                if (i != 0) {
                    I_up = i - 1;}
                else {I_up = N;}

                if (i != N) {
                    I_down = i + 1;}
                else {I_down = 0;}
                
                if (j != 0) {
                    J_left = j - 1;}
                else {J_left = N;}

                if (j != N) {
                    J_right = j + 1;}
                else {
                    J_right = 0;}    

                dE = 2 * lattice(i, j) * (lattice(i, J_left) + 
                                          lattice(I_up, j) + 
                                          lattice(I_down, j) + 
                                          lattice(i, J_right));
            }
        }

        void method_modulo(arma::mat& lattice, int idx)
        {
            double dE; 
            int N = L - 1;
            int i, j;
            int I_up, I_down, J_left, J_right;
            
            std::mt19937 generator;
            std::uniform_real_distribution<> uniform(0, 1);
            std::uniform_int_distribution<int> uniform_dist(0, N);
            generator.seed(Seed + idx);

            for (int attempt = 0; attempt < N_flips; attempt++)
            {
                i = uniform_dist(generator);
                j = uniform_dist(generator);

                I_up = (i + 1) % L;
                I_down = (i + L - 1) % L;

                J_left = (j + 1) % L;
                J_right = (j + L - 1) % L;

                dE = 2 * lattice(i, j) * (lattice(i, J_left) + 
                                          lattice(I_up, j) + 
                                          lattice(I_down, j) + 
                                          lattice(i, J_right));
            }
        }

        void method_vector(arma::mat& lattice, int idx)
        {
            double dE; 
            int N = L - 1;
            int i, j;
            int I_up, I_down, J_left, J_right;
            
            std::mt19937 generator;
            std::uniform_real_distribution<> uniform(0, 1);
            std::uniform_int_distribution<int> uniform_dist(0, N);
            generator.seed(Seed + idx);

            for (int attempt = 0; attempt < N_flips; attempt++)
            {
                i = uniform_dist(generator);
                j = uniform_dist(generator);

                I_up = Indexes(i);
                I_down = Indexes(i + 1);

                J_left = Indexes(j);
                J_right = Indexes(j + 1);

                dE = 2 * lattice(i, j) * (lattice(i, J_left) + 
                                          lattice(I_up, j) + 
                                          lattice(I_down, j) + 
                                          lattice(i, J_right));
            }
        }

        void method_int(arma::mat& lattice, int idx)
        {
            double dE; 
            int N = L - 1;
            int i, j;
            int I_up, I_down, J_left, J_right;
            
            std::mt19937 generator;
            std::uniform_real_distribution<> uniform(0, 1);
            std::uniform_int_distribution<int> uniform_dist(0, N);
            generator.seed(Seed + idx);

            for (int attempt = 0; attempt < N_flips; attempt++)
            {
                i = uniform_dist(generator);
                j = uniform_dist(generator);

                I_up = i + L * ((N - i) / N) - 1;
                I_down = (1 - i / N) * (i + 1);

                J_left = j + L * ((N - j) / N) - 1;
                J_right = (1 - j / N) * (j + 1);

                dE = 2 * lattice(i, j) * (lattice(i, J_left) + 
                                          lattice(I_up, j) + 
                                          lattice(I_down, j) + 
                                          lattice(i, J_right));
            }
        }

        void test_methods(int n)
        {   
            arma::mat time(n, 4);
            for (int k = 1; k <= n; k++)
            {
                int n_runs = std::pow(10, k); 
                arma::mat lattice = this->make_lattice();
                
                #pragma omp parallel
                #pragma omp single
                #pragma omp taskloop num_tasks(4)
                for (int j = 0; j < 4; j++)
                {   
                    if (j == 0)
                    {
                        auto t1 = std::chrono::high_resolution_clock::now();
                        for (int i = 0; i < n_runs; i++)
                        {
                            this->method_normal(lattice, i);
                        }
                        auto t2 = std::chrono::high_resolution_clock::now();
                        double duration_seconds = 
                             std::chrono::duration<double>(t2 - t1).count();
        
                        time(k - 1, j) = duration_seconds;
                    }

                    if (j == 1)
                    {
                        auto t1 = std::chrono::high_resolution_clock::now();
                        for (int i = 0; i < n_runs; i++)
                        {
                            this->method_modulo(lattice, i);
                        }
                        auto t2 = std::chrono::high_resolution_clock::now();
                        double duration_seconds = 
                             std::chrono::duration<double>(t2 - t1).count();
        
                        time(k - 1, j) = duration_seconds;
                    }
                    
                    if (j == 2)
                    {
                        auto t1 = std::chrono::high_resolution_clock::now();
                        for (int i = 0; i < n_runs; i++)
                        {
                            this->method_vector(lattice, i);
                        }
                        auto t2 = std::chrono::high_resolution_clock::now();
                        double duration_seconds = 
                             std::chrono::duration<double>(t2 - t1).count();
        
                        time(k - 1, j) = duration_seconds;
                    }
                    
                    if (j == 3)
                    {   
                        auto t1 = std::chrono::high_resolution_clock::now();
                        for (int i = 0; i < n_runs; i++)
                        {
                            this->method_int(lattice, i);
                        }
                        auto t2 = std::chrono::high_resolution_clock::now();
                        double duration_seconds = 
                             std::chrono::duration<double>(t2 - t1).count();
        
                        time(k - 1, j) = duration_seconds;             
                    }
                }
            }
            std::ofstream MyFile("data/test_methods.txt");
            MyFile << time << std::endl;
            MyFile.close();
        }

        void flip(arma::mat& lattice, arma::vec dE_vec, int idx)         
        {   
            bool case_1, case_2;
            
            double dE; 
            
            int N = L - 1;
            int i, j, dE_idx;
            int I_up, I_down, J_left, J_right;
            
            std::mt19937 generator;
            std::uniform_real_distribution<> uniform(0, 1);
            std::uniform_int_distribution<int> uniform_dist(0, N);
            generator.seed(Seed + idx);
            
            for (int attempt = 0; attempt < N_flips; attempt++)
            {
                i = uniform_dist(generator);
                j = uniform_dist(generator);

                I_up = i + L * ((N - i) / N) - 1;
                I_down = (1 - i / N) * (i + 1);

                J_left = j + L * ((N - j) / N) - 1;
                J_right = (1 - j / N) * (j + 1);

                dE = 2 * lattice(i, j) * (lattice(i, J_left) + 
                                          lattice(I_up, j) + 
                                          lattice(I_down, j) + 
                                          lattice(i, J_right));

                dE_idx = (dE + 8) / 4;

                case_1 = dE <= 0;
                case_2 = uniform(generator) < dE_vec(dE_idx);

                lattice(i, j) *= 1 - 2 * (case_1 + case_2 * (1 - case_1));
            }    
        }

        void temp_vec(arma::mat &T, int T_step, double T_0, double T_1)                 
        {
            double dT = 0;
            arma::mat T_vec(T_step + 1, 1, arma::fill::zeros);
            bool is_zero = (T_step != 0);
            
            dT = (T_1 - T_0) / (T_step * is_zero  + (1 - is_zero));

            for (double i = 0; i <= T_step; i++)
            {
                T_vec(i) = T_0 + i * dT;
            }
            T = T_vec;
        }

        void simulate(arma::mat& result, arma::mat T_vec, int Cycle)       
        {
            int n_T = T_vec.size();   

            double n1 = 1.0 / (Cycle * std::pow(L, 2));
            double n2 = 1.0 / (std::pow(Cycle, 2) * std::pow(L, 2));
            
            arma::mat res(n_T, 5, arma::fill::zeros);
            res.col(0) = T_vec;
            
            #pragma omp parallel
            #pragma omp single
            #pragma omp taskloop num_tasks(4)
            for (int I = 0; I < n_T; I++) 
            {   
                double T = T_vec(I);

                arma::mat lattice = this->make_lattice();
                arma::vec dE_vec = this->make_delta_E_vector(T);;

                double E, M;
                double E1 = 0;
                double E2 = 0;
                double M1 = 0;
                double M2 = 0;
                
                for (int i = 0; i < Cycle; i++)
                {
                    this->flip(lattice, dE_vec, i);
                    M = this->get_mag(lattice);
                    E = this->energy(lattice);
                    
                    M1 = M1 + M;
                    M2 = M2 + M * M;

                    E1 = E1 + E;
                    E2 = E2 + E * E; 
                }
                res(I, 1) = n1 * E1;
                res(I, 2) = n1 * M1;
                res(I, 3) = (n1 * E2 - n2 * std::pow(E1, 2)) / T;
                res(I, 4) = (n1 * M2 - n2 * std::pow(M1, 2)) / std::pow(T, 2);
            }
            result = res; 
        }
        
        void run_temp_setup(int T_step, double T0, double T1)
        {
            this->temp_vec(T, T_step, T0, T1);
        }
        
        void timer_simulate(int n)
        {
            arma::mat time(n, 2);
            for (int k = 1; k <= n; k++)
            {
                int Cycle = std::pow(10, k);
                
                auto t1 = std::chrono::high_resolution_clock::now();
                this->simulate(result, T, Cycle);
                auto t2 = std::chrono::high_resolution_clock::now();
                
                double duration_seconds = 
                        std::chrono::duration<double>(t2 - t1).count();
                time(k - 1, 0) = duration_seconds / Cycle;
                time(k - 1, 1) = Cycle;
            }

            #ifdef _OPENMP
            {
                std::ofstream MyFile("data/time_p_1.txt");
                MyFile << time << std::endl;
                MyFile.close();
            }
            #else
            {
                std::ofstream MyFile("data/time_p_0.txt");
                MyFile << time << std::endl;
                MyFile.close();
            }
            #endif
        }
};

int main()
{
    Speed_Class speed(20, 10);
    speed.test_methods(6);
    // speed.run_temp_setup(10, 1, 2);
    // speed.timer_simulate(5);
    return 0;
}



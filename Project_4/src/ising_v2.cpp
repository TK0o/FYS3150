#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <cmath>
#include <string>


void problem_5(int n);
void problem_8();

int main()
{
    problem_5(16);
    // problem_8();
    
    return 0;
}

class Ising_Class
{
    public:
        int L, Seed, Cycle, N_flips;
        arma::mat result;
        arma::mat T;

    public:
        Ising_Class (int lattice_size, int seed, int cycle)
        {
            L = lattice_size;
            Seed = seed;
            Cycle = cycle;
            N_flips = std::pow(L, 2);
        }

        arma::mat make_lattice(int ordered)
        {   
            arma::arma_rng::set_seed(Seed);
            arma::mat lattice(L, L, arma::fill::randu);
            
            for (int i = 0; i < L; i++)
            {
                for (int j = 0; j < L; j++)
                {
                    lattice(i, j) = 1 - 2 * (lattice(i, j) < 0.5) * 
                                            (1 - ordered);
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
            if (T_step != 0)
            {
            dT = (T_1 - T_0) / T_step;
            }

            for (double i = 0; i <= T_step; i++)
            {
                T_vec(i) = T_0 + i * dT;
            }
            T = T_vec;
        }

        void simulate_normal(arma::mat& result, arma::mat T_vec, int ordered)       
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

                arma::mat lattice = this->make_lattice(ordered);
                arma::vec dE_vec = this->make_delta_E_vector(T);

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

        void run(int ordered)
        {
            this->simulate_normal(result, T, ordered);
        }

        arma::mat get_result()
        {
            return result;
        }

        void write_to_file(std::string filename)
        {
            std::ofstream MyFile(filename);
            MyFile << result << std::endl;
            MyFile.close();
        }
};

void problem_5(int n)
{
    arma::mat res(1 + 5 * 4, n);
    arma::mat temp_res(5, 2);
    arma::mat cycle(1, n);
    for (int k = 1; k <= n; k++)
    {
        int cycle = std::pow(2, k);
        for (int i = 0; i < 2; i++)
        {
            Ising_Class ising (20, 100, cycle);
            ising.run_temp_setup(1, 1., 2.4);
            ising.run(i);
            temp_res = trans(ising.get_result());
            res(arma::span(1 + i * 10, 5 + i * 10), k - 1) = temp_res.col(0);
            res(arma::span(6 + i * 10, 10 + i * 10), k - 1) = temp_res.col(1); 
        }
        res(0, k - 1) = cycle;
    }
    std::ofstream MyFile("data/problem_5.txt");
    MyFile << res << std::endl;
    MyFile.close();
}

void problem_8()
{
    arma::mat res;
    for (int i = 0; i < 4; i++)
    {
        int L = 40 + i * 20;
        Ising_Class ising (L, 100, 10000);
        ising.run_temp_setup(50, 2.1, 2.4);
        ising.run(0);
        res = ising.get_result();

        std::string filename = "data/problem8_v2_N" + std::to_string(L) + ".txt";                  
        std::ofstream MyFile(filename);
        MyFile << res << std::endl;
        MyFile.close();
    }
}

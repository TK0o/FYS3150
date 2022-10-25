#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

class Particle;
class PenningTrap;

std::vector<Particle> single_particle_setup();
std::vector<Particle> double_particle_setup();

void problem_8();
void problem_8_single_particle(double B_0, double V_0, double d);
void problem_8_double_particle(double B_0, double V_0, double d);
void problem_8_double_particle_phase(double B_0, double V_0, double d);
void problem_8_single_relative_error(double B_0, double V_0, double d, 
                                     std::string method);

void problem_9_time_dependent(int n_particles, double axis_step_size,
                              double axis_step_0, double axis_step_1);

void problem_9_fine_search(int n_particles, double axis_step_size,
                           double axis_step_0, double axis_step_1, 
                           std::string filename);


int main()
{   
    double B_0 = 9.65e1;
    double V_0 = 2.41e6;
    double d = 5e2;
    // uncomment to run these
    // problem_8();

    // problem_9_time_dependent(100, 0.02, 0.2, 2.5);
    
    // first attempt
    // problem_9_fine_search(100, 0.001, 0.9, 1.1, 
    //                     "problem_9/p9_fine_search_interaction_.txt");
    
    problem_9_fine_search(100, 0.0002, 0.97, 1., 
                         "problem_9/p9_fine_search_ 2_interaction_.txt");
    return 0;
}


class Particle
{
    public:
        double Charge;
        double Mass;
        arma::vec pos;
        arma::vec vel;
    
    public:
        Particle (arma::vec v={0, 0, 0}, arma::vec r={0, 0, 0}, 
                  double q = 1, double m = 40.077)
        {
            Charge = q;
            Mass = m;
            vel = v;
            pos = r;      
        }
        void check_val()
        {
            std::cout << Charge << "\t" << Mass << "\n" 
                      << pos << "\n" << vel << '\n';
        }
    friend class PenningTrap;
};

class PenningTrap
{
    public:
        double B;
        double V;
        double d;
        double V0;
        
        double h_k; 
        double A_p, A_m;
        double w_0, w_z2, freq;
        
        int n, N, n_combi;
        
        std::vector<Particle> arr;
        arma::mat analytic_arr;
        arma::mat numeric_arr;
        arma::mat numeric_vel;
        arma::vec in_trap;
        arma::mat combi;
        arma::mat f;
    
    public:
        PenningTrap (double B_0=0, double V_0=0, double d_=0)
        {
            B = B_0;
            V = V_0;
            V0 = V_0;
            d = d_;
        }

        void add_particle(std::vector<Particle> A)
        {   
            arr = A;
            N = arr.size();
        }
        
        void generate_particles(int N_particles, int seed=10)
        {
            arma::arma_rng::set_seed(seed);
            std::vector<Particle> randomly_generated;
            randomly_generated.resize(N_particles);

            for (int i = 0; i < N_particles; i++)
            {
                arma::vec v = arma::vec(3).randn() * 0.1 * d;
                arma::vec r = arma::vec(3).randn() * 0.1 * d;
                randomly_generated[i] = Particle(v, r);
            }
            arr = randomly_generated;
            N = N_particles;
        }

        void auto_double_particle()
        {
            std::vector<Particle> Arr;
            Arr.resize(2);

            arma::vec v1 = {0, 25, 0};
            arma::vec p1 = {20, 0, 20};
            Arr[0] = Particle(v1, p1, 1, 40.077);

            arma::vec v2 = {0, 40, 5};
            arma::vec p2 = {25, 25, 0};
            Arr[1] = Particle(v2, p2, 1, 40.077);
            arr = Arr;
            N = 2;
        }

        void setup(int t, int n_step)
        {   
            arma::mat num_arr(3 * N, n_step);
            
            numeric_arr = num_arr;
            numeric_vel = num_arr;

            arma::mat F(3, int(N));
            f = F.zeros();
            n = n_step;
            
            h_k = double (t) / double (n_step);
            
            arma::vec in_the_trap(N);
            in_trap = in_the_trap.fill(1);
            this->make_combination();
        }

        void numerical(std::string method="0", int interaction=0, 
                       int external_field=0, int time_dependent_V=0,
                       double A_f=0, double w_v=0)
        {
            double d2 = std::pow(d, 2);
            if (N < 2)
            {
                interaction = 0;
            }

            if (method == "0" || method == "euler" || 
                method == "forward euler")    
            {
                for (int I = 0; I < n; I++)
                {
                    if (interaction == 1 || external_field == 1)
                    {
                        this->force_interaction(external_field, d2,
                                                interaction);
                    }

                    if (time_dependent_V == 1)
                    {
                        V = V0 * (1 + A_f * std::cos(w_v * I * h_k));
                    }

                    for (int i = 0; i < N; i++)
                    {
                        this->euler(i, I);
                    }
                    f = f.zeros();
                }
                numeric_arr = trans(numeric_arr);
                numeric_vel = trans(numeric_vel);
            }
            
            if (method == "1" || method == "rk4")
            {
                for (int I = 0; I < n; I++)
                {
                    if (interaction == 1 || external_field == 1)
                    {
                        this->force_interaction(external_field, d2,
                                                interaction);
                    }
                    
                    if (time_dependent_V == 1)
                    {
                        V = V0 * (1 + A_f * std::cos(w_v * I * h_k));
                    }

                    for (int i = 0; i < N; i++)
                    {
                        this->rk4(i, I);
                    }
                    arma::mat temp_f(3, int(N));
                    f = temp_f.zeros();
                }
                numeric_arr = trans(numeric_arr);
                numeric_vel = trans(numeric_vel);
            }
        }


        void force_interaction(int external_field, double d2, 
                               int interaction=0)
        {  
            double R_norm;
            int ind_i, ind_j;
            arma::vec R, temp_force, Pos_i;
            for (int J = 0; J < n_combi; J++)
            {
                ind_i = combi(J, 0);
                ind_j = combi(J, 1);
                
                Pos_i = arr[ind_i].pos;

                if (external_field == 1 && arma::cdot(Pos_i, Pos_i) > d2)
                {
                    in_trap(ind_i) = 0;
                }
                else if (R_norm != 0 && interaction == 1)
                {
                    R = Pos_i - arr[ind_j].pos;
                    R_norm = std::pow(arma::dot(R, R), 1.5);
                    temp_force = arr[ind_i].Charge * R / R_norm;

                    f.col(ind_i) += temp_force;
                    f.col(ind_j) -= temp_force;
                }

            }
        }

        arma::vec force(arma::vec vel, arma::vec pos, 
                        int i, arma::mat f={0, 0, 0})
        {   
            double fx, fy, fz, w0, wz;
            double m = arr[i].Mass;
            int q = arr[i].Charge;

            w0 = (q * B * in_trap(i)) / m;
            wz = (2. * q * V * in_trap(i)) / (m * std::pow(d, 2));

            fx = vel(1) * w0 + 0.5 * wz * pos(0) + f(0);
            fy = -vel(0) * w0 + 0.5 * wz * pos(1) + f(1);
            fz = - 0.5 * wz * pos(2) + f(2);
            arma::vec new_force;
            return new_force = {fx, fy, fz};
        }

        void euler(int i, int j)
        {   
            arma::vec Vel, Pos;
            Vel = arr[i].vel;
            Pos = arr[i].pos;

            arma::vec acceleration = force(Vel, Pos, i, f.col(i));
            arr[i].vel = Vel + acceleration * h_k;
            arr[i].pos = Pos + Vel * h_k;
            
            numeric_vel(arma::span(i * 3, i * 3 + 2), j) = arr[i].vel;
            numeric_arr(arma::span(i * 3, i * 3 + 2), j) = arr[i].pos;
        }

        void rk4(int i, int j)
        {
            arma::vec k_v1, k_v2, k_v3, k_v4;
            arma::vec k_r1, k_r2, k_r3, k_r4;
            arma::vec Vel, Pos;
            Vel = arr[i].vel;
            Pos = arr[i].pos;

            k_v1 = Vel;  
            k_r1 = force(k_v1, Pos, i, f.col(i));
            
            k_v2 = Vel + 0.5 * h_k * k_r1;
            k_r2 = force(k_v2, Pos + 0.5 * h_k * k_v2, i, f.col(i));
            
            k_v3 = Vel + 0.5 * h_k * k_r2;
            k_r3 = force(k_v3, Pos + 0.5 * h_k * k_v3, i, f.col(i));
            
            k_v4 = Vel + h_k * k_r3;
            k_r4 = force(k_v4, Pos + h_k * k_v4, i, f.col(i));

            arr[i].vel = Vel + h_k / 6. * (k_r1 + 2. * k_r2 +  
                                           2. * k_r3 + k_r4);
            arr[i].pos = Pos + h_k / 6. * (k_v1 + 2. * k_v2 +  
                                           2. * k_v3 + k_v4);

            numeric_vel(arma::span(i * 3, i * 3 + 2), j) = arr[i].vel;
            numeric_arr(arma::span(i * 3, i * 3 + 2), j) = arr[i].pos;
        }

        void make_combination()
        {         
            n_combi = (N * (N - 1)) / 2;
            arma::mat combination(n_combi, 2);
            
            int I = 0;
            for (int i = 0; i < N - 1; i++)
            {
                for (int j = i + 1; j < N; j++)
                {   
                    combination(I, 0) = i;
                    combination(I, 1) = j;
                    I++;
                }
            }
            combi = combination;   
        }

        void analytical(double x_0, double v_0, double z_0, 
                        double t, double i_steps,
                        double t_0=0, double q=1, double m=40.077)
        {
            double w_p, w_m, w_z, dt;
            std::complex<double> solution;
            arma::mat analytical_arr(int(i_steps), 3);

            w_0 = (q * B) / m;
            w_z2 = (2 * q * V) / (m * std::pow(d, 2));
            double W = std::sqrt(std::pow(w_0, 2) - 2 * w_z2);
            
            w_p = (w_0 + W) / 2;
            w_m = (w_0 - W) / 2;
            A_p = (v_0 + w_m * x_0) / (w_m - w_p);
            A_m = -(v_0 + w_p * x_0) / (w_m - w_p);

            w_z = std::sqrt(w_z2);
            dt = (t - t_0) / (i_steps);

            for (int i = t_0 + 1; i <= i_steps; i++)
            {
                std::complex<double> exponent1(0, -(w_p * dt * i));
                std::complex<double> exponent2(0, -(w_m * dt * i));
                solution = A_p * std::exp(exponent1) + 
                           A_m * std::exp(exponent2);
                
                analytical_arr((i - 1), 0) = std::real(solution);
                analytical_arr((i - 1), 1) = std::imag(solution);
                analytical_arr((i - 1), 2) = z_0 *std::cos(w_z * dt * i);
            }
            analytic_arr = analytical_arr;
        }
        
        void check_vals()
        {
            std::cout << "B = " << B << '\n'
                      << "V = " << V << '\n'
                      << "d = " << d << '\n';
        }

        void write_to_file(arma::mat &A, std::string filename)
        {
            std::ofstream MyFile("data/" + filename);
            MyFile << A << std::endl;
            MyFile.close();
        }

        void test()
        {
            for (int i = 0; i < N; i++)
            {
                arr[i].check_val();
            }
        }
};

std::vector<Particle> single_particle_setup()
{
    std::vector<Particle> arr;
    arr.resize(1);

    arma::vec v = {0, 25, 0};
    arma::vec p = {20, 0, 20};
    arr[0] = Particle(v, p, 1, 40.077);
    return arr;
}

std::vector<Particle> double_particle_setup()
{
    std::vector<Particle> arr;
    arr.resize(2);

    arma::vec v1 = {0, 25, 0};
    arma::vec p1 = {20, 0, 20};
    arr[0] = Particle(v1, p1, 1, 40.077);

    arma::vec v2 = {0, 40, 5};
    arma::vec p2 = {25, 25, 0};
    arr[1] = Particle(v2, p2, 1, 40.077);

    return arr;
}

void problem_8()
{
    double B_0 = 9.65e1;
    double V_0 = 2.41e6;
    double d = 5e2;
    
    problem_8_single_particle(B_0, V_0, d);
    problem_8_double_particle(B_0, V_0, d);
    problem_8_double_particle_phase(B_0, V_0, d);
    problem_8_single_relative_error(B_0, V_0, d, "0");
    problem_8_single_relative_error(B_0, V_0, d, "1");
}

void problem_8_single_particle(double B_0, double V_0, double d)
{  
    std::vector<Particle> arr = single_particle_setup();
    PenningTrap p_trap(B_0, V_0, d);
    p_trap.add_particle(arr);
    p_trap.setup(50, 4000);
    p_trap.numerical("1");
    p_trap.write_to_file(p_trap.numeric_arr, 
                         "problem_8/rk4_p8_1_particle.txt");
}

void problem_8_double_particle(double B_0, double V_0, double d)
{
    std::vector<Particle> arr = double_particle_setup();
    PenningTrap p_trap(B_0, V_0, d);
    
    for (int i = 0; i < 2; i++)
    {
        p_trap.add_particle(arr);
        p_trap.setup(50, 4000);
        p_trap.numerical("1", i);
        
        std::string filename;
        filename = "problem_8/rk4_p8_2_particle_interact_" 
                    + std::to_string(i) + "_.txt";

        p_trap.write_to_file(p_trap.numeric_arr, filename);
    }
}

void problem_8_double_particle_phase(double B_0, double V_0, double d)
{   
    int I, k;
    int J = 0;
    std::vector<Particle> arr = double_particle_setup();
    PenningTrap p_trap(B_0, V_0, d);
    for (int i = 0; i < 2; i++)
    {   
        p_trap.add_particle(arr);
        p_trap.auto_double_particle();
        p_trap.setup(50, 4000);
        p_trap.numerical("1", i);
        
        arma::mat phase(4000, 8);
        for (int j = 0; j < 4; j++)
        {   
            I = 2 * (j % 2);
            k = j * 2;
            phase(arma::span(0, 3999), k) = p_trap.numeric_arr.col(J + I);
            phase(arma::span(0, 3999), 1 + k) = p_trap.numeric_vel.col(J + I);
            
            if (j % 2 == 1)
            {
                J = 3;
            }
        }
        J = 0;
        std::string filename;
        filename = "problem_8/rk4_p8_phase_interact_" 
                    + std::to_string(i) + "_.txt";
        p_trap.write_to_file(phase, filename);
    }
}

void problem_8_single_relative_error(double B_0, double V_0, double d, 
                                     std::string method)
{
    std::vector<Particle> arr = single_particle_setup();
    PenningTrap p_trap(B_0, V_0, d);
    std::string method_string;
    for (int i = 0; i < 4; i++)
    {
        int i_step = 4000 * std::pow(2, i);
        p_trap.add_particle(arr);
        p_trap.setup(50, i_step);

        if (method == "0")
        {
            method_string = "euler";
            p_trap.numerical("0");
            
        }
        else if(method == "1")
        {
            method_string = "rk4";
            p_trap.numerical("1");
        }
        arma::mat m(i_step, 6);
        p_trap.analytical(20, 25, 20, 50, i_step);

        for (int j = 0; j < 3; j++)
        {   
            m(arma::span(0, i_step - 1), j) = p_trap.numeric_arr.col(j);
            m(arma::span(0, i_step - 1), j + 3) = p_trap.analytic_arr.col(j);
        }
        
        std::string filename;
        
        filename = "problem_8/" + method_string + "_p8_relative_error" 
                        + std::to_string(i_step) + "_.txt";
        p_trap.write_to_file(m, filename);
    }
}

void problem_9_time_dependent(int n_particles, double axis_step_size,
                              double axis_step_0, double axis_step_1)
{
    double n0 = axis_step_0 / axis_step_size;
    double n1 = axis_step_1 / axis_step_size;
    double n_tot = n1 - n0 + 1;

    double B_0 = 9.65e1;
    double V_0 = 2.41e6;
    double d = 5e2;

    arma::mat counter(n_tot, 4);
    double Angular_f, w_v;
    PenningTrap p_trap(B_0, V_0, d);

    std::vector<Particle> Arr;
    p_trap.generate_particles(100);
    Arr = p_trap.arr;

    for (int j = 0; j < 3; j++)
    {   
        Angular_f = 0.1 + j * 0.3;
        for (int i = 0; i < n_tot; i++)
        {   
            w_v = (i + n0) * axis_step_size;

            p_trap.add_particle(Arr);
            p_trap.setup(500, 40000);
            p_trap.numerical("1", 0, 1, 1, Angular_f, w_v);
            
            counter(i, j + 1) = arma::sum(p_trap.in_trap);
            counter(i, 0) = w_v;
            std::cout << i << '\t' << j << '\t' << arma::sum(p_trap.in_trap) 
                      << std::endl;
        }
        p_trap.write_to_file(counter, 
        "problem_9/p9_time_dependent.txt");
    }
}


void problem_9_fine_search(int n_particles, double axis_step_size,
                           double axis_step_0, double axis_step_1, 
                           std::string filename)
{
    double n0 = axis_step_0 / axis_step_size;
    double n1 = axis_step_1 / axis_step_size;
    double n_tot = n1 - n0 + 1;

    double B_0 = 9.65e1;
    double V_0 = 2.41e6;
    double d = 5e2;

    arma::mat counter(n_tot, 4);
    double Angular_f, w_v;
    PenningTrap p_trap(B_0, V_0, d);

    std::vector<Particle> Arr;
    p_trap.generate_particles(100);
    Arr = p_trap.arr;

    Angular_f = 0.1;
    for (int j = 0; j < 2; j++)
    {   
        for (int i = 0; i < n_tot; i++)
        {               
            w_v = (i + n0) * axis_step_size;

            p_trap.add_particle(Arr);
            p_trap.setup(500, 40000);
            p_trap.numerical("1", j, 1, 1, Angular_f, w_v);
            
            counter(i, j + 1) = arma::sum(p_trap.in_trap);
            counter(i, 0) = w_v;
            std::cout << i << '\t' << j << '\t' << arma::sum(p_trap.in_trap) 
                      << std::endl;
        }
        p_trap.write_to_file(counter, filename);
    }
}

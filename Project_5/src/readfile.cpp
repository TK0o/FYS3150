#include "../header/readfile.hpp"

#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <cstring>


void remove_spaces(std::string& parameter_name)
{
    for (int i = 0; i < parameter_name.length(); i++)
    {
        if (parameter_name[i] == ' ')
        {
            parameter_name.erase(parameter_name.begin() + i);
            i--;
        }
    }
}
void allocate_values(std::string parameter_name, double value, 
                     std::vector<std::string> parameters_vec,
                     arma::mat& temp)
{
    remove_spaces(parameter_name);
    int i = 0;
    for (auto j : parameters_vec)
    {
        if (parameter_name == j)
        {
            temp(0, i) = value;
        }
        i++;
    }
}

void make_parameter_vec(std::vector<std::string>& parameters)
{
    parameters = {
    "h",
    "x_c",
    "y_c",
    "sigma_x",
    "sigma_y",
    "p_x",
    "p_y",
    "delta_t",
    "t_stop",
    "v_0",
    "wall_potential",
    "wall_thickness_x",
    "wall_pos_centre_x",
    "slit_dist_y",
    "slit_size_y",
    "slit_centred_y",
    "n_slits"};
}

void readfile(arma::mat& parameter_values_mat, 
              std::vector<std::string>& file_name_vec,
              std::string input_file)
{    
    std::vector<std::string> name_vec;
    make_parameter_vec(name_vec);

    std::string line; 
    std::ifstream MyFile("../" + input_file);
    std::size_t double_slash;
    std::size_t l_comment_s;
    std::size_t l_comment_e;
    std::size_t eq;
    std::size_t file_name;
    std::size_t n_simulations;
    int multi_line_comment = 0;
    
    arma::mat simulation_values (1, 17, arma::fill::zeros);
    arma::mat temp_values (1, 17, arma::fill::zeros);

    int sim_idx = 0;

    while (std::getline(MyFile, line))
    { 
        l_comment_s = line.find("/*");
        l_comment_e = line.find("*/");

        int slash_asterisk = (l_comment_s != std::string::npos);
        int asterisk_slash = (l_comment_e != std::string::npos);

        if (slash_asterisk == 1 && asterisk_slash == 1)       
        {
            line.erase(l_comment_s, l_comment_e - l_comment_s + 2);
        }
        else if (multi_line_comment == 0 && slash_asterisk == 1)
        {
            multi_line_comment = 1;
            line.erase(l_comment_s);
        }
        else if (multi_line_comment == 1 && asterisk_slash == 1)
        {
            multi_line_comment = 0;
            line.erase(0, l_comment_e + 2);
        }
        else if (multi_line_comment == 1)
        {
            continue;
        }

        double_slash = line.find("//");
        if (double_slash != std::string::npos)
        {
            line.erase(double_slash);
        }        

        for (int i = 0; i < line.length(); i++)
        {
            line[i] = std::tolower(line[i]);
        }
        
        eq = line.find("=");
        file_name = line.find("file_name");
        if (file_name != std::string::npos)
        {
            std::string parameter = line.substr(0, file_name);            
            std::string name = line.substr(eq + 1, line.length());
            remove_spaces(name);

            file_name_vec[sim_idx - 1] = name;
            continue;
        } 

        if (eq != std::string::npos)
        {
            std::string parameter = line.substr(0, eq);            
            std::string num = line.substr(eq + 1, line.length());
            std::istringstream str_num(num);
            double value;
            str_num >> value;
            allocate_values(parameter, value, name_vec, temp_values);
        }
         
        simulation_values(sim_idx - (sim_idx > 0), arma::span(0, 16)) = 
                          temp_values.row(0);

        n_simulations = line.find("simulate");
        if (sim_idx == 0 && n_simulations != std::string::npos)
        {   
            sim_idx++;
            file_name_vec.push_back("simulation_" + std::to_string(sim_idx) +
                                    ".txt");
            continue;
        }
        else if (n_simulations != std::string::npos)
        {
            arma::mat sim_values (sim_idx + 1, 17, arma::fill::zeros);
            sim_values(sim_idx - 1, arma::span(0, 16)) = temp_values.row(0);
            
            sim_values(arma::span(0, sim_idx - 1), arma::span(0, 16)) = 
                        simulation_values;
            temp_values.fill(0);
            simulation_values = sim_values;
            sim_idx++;

            file_name_vec.push_back("simulation_" + std::to_string(sim_idx) +
                                    ".txt");
        }
    }
    simulation_values(sim_idx - (sim_idx > 0), arma::span(0, 16)) = 
                          temp_values.row(0);
    MyFile.close();
    parameter_values_mat = simulation_values;
}

void check_default_values(arma::mat& A)
{   
    /*
    index for a given parameter
    "h"                 |     0
    "x_c"               |     1
    "y_c"               |     2
    "sigma_x"           |     3
    "sigma_y"           |     4
    "p_x"               |     5
    "p_y"               |     6
    "delta_t"           |     7
    "t_stop"            |     8
    "v_0"               |     9
    "wall_potential"    |    10
    "wall_thickness_x"  |    11
    "wall_pos_centre_x" |    12
    "slit_dist_y"       |    13
    "slit_size_y"       |    14
    "slit_centred_y"    |    15
    "n_slits"           |    16
    */

    // index for a given parameter
    arma::vec Default_Values =
    {
        0.005,
        0.25,
        0.5,
        0.005,
        0.005,
        200,
        0,
        2.5e-5,
        0.008,
        1e10,
        1e15,
        0.02,
        0.5,
        0.05,
        0.05,
        0.5,
        2.0
    };
    
    int n_simulations = A.size() / 17;
    for (int i = 0; i < n_simulations; i++)
    {
        for (int j = 0; j < 17; j++)
        {   
            if (A(i, j) == 0)
            {
                A(i, j) = Default_Values(j);
            }
        }
    }
}
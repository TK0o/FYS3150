// Include guard
#ifndef __readfile_hpp__
#define __readfile_hpp__

#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <cstring>

void remove_spaces(std::string& parameter_name);

void allocate_values(std::string parameter_name, double value, 
                     std::vector<std::string> parameters_vec,
                     arma::mat& temp);

void make_parameter_vec(std::vector<std::string>& parameters);

void readfile(arma::mat& parameter_values_mat, 
              std::vector<std::string>& file_name_vec,
              std::string input_file);

void check_default_values(arma::mat& A);
#endif

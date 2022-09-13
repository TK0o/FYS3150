#include <iostream>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <math.h>
#include <fstream>
using namespace std;

vector<float> get_x_values(int n);
vector<float> get_f_values(vector<float>& x_val);

vector<float> analytical (vector<float>& x_val);

void special_algorithm(vector<float>& f_val, int n);

void absolute_error(string file_name, int n, int n_up_to);
void relative_error(string file_name, int n, int n_up_to);

void print_m(vector<vector<float>>& vec);
void print_v(vector<float>& vec);

int main() {
  int n = 10;

  string file_name = "absolute_error.txt";
  absolute_error (file_name, n, 6);

  string file_name2 = "relative_error.txt";
  relative_error (file_name2, n, 7);
  return 0;
}

void absolute_error(string file_name, int n, int n_up_to){
  ofstream MyFile("data/" + file_name);
  float delta_val;
  int N;

  for (int i = 1; i < n_up_to + 1; i++){
    N = pow (n, i);

    MyFile << scientific << setprecision(1) << "n = " << float(N) << endl;
    MyFile << setw(6) << "x" << setw(8) << "|" << setw(10) << "log10" << endl;

    vector<float> x_val = get_x_values(N);
    vector<float> f_val = get_f_values(x_val);
    vector<float> exact = analytical(x_val);
    special_algorithm(f_val, N);

    for (int j = 0; j < exact.size(); j++){
      delta_val = abs(exact[j] - f_val[j + 1]);

      if (delta_val == 0) {
        delta_val = NAN;
      }
      else{
        delta_val = log10(delta_val);
      }
      MyFile << setw(7) << setprecision (5)
          << scientific << x_val[j] << setw(17) << delta_val << endl;
    }
    MyFile << "-----------------------------\n" << endl;
  }
  MyFile.close();
}

void relative_error(string file_name, int n, int n_up_to){
  ofstream MyFile("data/" + file_name);
  float delta_val;
  vector<float> relative_val(n_up_to, 0);
  int N;

  MyFile << setw(8) << "n^x" << setw(6) << "|"
         << setw(10) << "log10" << endl;

  for (int i = 1; i < n_up_to + 1; i++){
    N = pow (n, i);

    vector<float> x_val = get_x_values(N);
    vector<float> f_val = get_f_values(x_val);
    vector<float> exact = analytical(x_val);
    special_algorithm(f_val, N);

    for (int j = 0; j < exact.size(); j++){
      delta_val = abs((exact[j] - f_val[j + 1]) / exact[j]);

      if (delta_val == 0) {
        delta_val = NAN;
      }
      else{
        if (relative_val[i] < delta_val){
          relative_val[i] = delta_val;
        }
        delta_val = log10(delta_val);
      }
    }
    relative_val[i] = log10(relative_val[i]);
    MyFile << defaultfloat << setw(8) << i << setw(20) << scientific
    << setprecision(5) << relative_val[i] << endl;
  }
  MyFile.close();
}

void special_algorithm(vector<float>& f_val, int n){
  int N = n - 1;
  float h = 1 / float(n);
  float const_h = 100 * pow (h, 2);

  for (float i = 0; i < N; i++){
    f_val[i + 1] = (f_val[i + 1] + f_val[i]) * ((i + 1) / (i + 2));
  }

  for (float i = N - 1; i > 0; i--){
    f_val[i] = f_val[i] + f_val[i + 1] * (i / (i + 1));
  }

  for (int i = 1; i < f_val.size() - 1; i++){
    f_val[i] = f_val[i] * const_h;
  }
}

vector<float> get_x_values(int n){
  float h = 1 / float (n);
  vector<float> vec_for_x(n, 0);
  for (int i = 0; i < vec_for_x.size(); i++){
    vec_for_x[i] = (i + 1) * h;
  }
  return vec_for_x;
}

vector<float> get_f_values(vector<float>& x_val){
  int N = x_val.size();
  vector<float> vec_for_f(N + 1, 0);

  for (int i = 0; i < N - 1; i++){
    vec_for_f[i + 1] = exp (-10 * x_val[i]);
  }
  return vec_for_f;
}

vector<float> analytical (vector<float>& x_val){
  vector<float> vec_analytical(x_val.size() - 1, 0);

  for (int i = 0; i < x_val.size() - 1; i++) {
      vec_analytical[i] = 1 - (1 - exp(-10)) * x_val[i] - exp(-10 * x_val[i]);
  }
  return vec_analytical;
}

void print_v(vector<float>& vec){
  for(int i = 0; i < vec.size(); i++){
    cout << scientific << setprecision (5) << setw(12) << vec[i];
  }
  cout << endl;
}

void print_m(vector<vector<float>>& vec){
  for(int i = 0; i < vec.size(); i++){
    for(int j = 0; j < vec[0].size(); j++){
      cout << defaultfloat <<
      setprecision (5) << setw(15) << vec[i][j] << setw(10);
    }
  	cout << endl;
  }
}

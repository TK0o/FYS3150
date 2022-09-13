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

void write_to_file(vector<float>& x_val, vector<vector<float>>& u_val,
  vector<float>& u_analytic, string file_name, int n);

void print_m(vector<vector<float>>& vec);
void print_v(vector<float>& vec);

int main() {
  int n = 10;

  string file_name = "task9.txt";

  vector<float> x_val = get_x_values(n);
  vector<float> f_val = get_f_values(x_val);
  vector<float> exact = analytical(x_val);

  special_algorithm(f_val, n);
  return 0;
}

// implementation of the special algorithm
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

void write_to_file(vector<float>& x_val, vector<vector<float>>& u_val,
  vector<float>& u_analytic, string file_name, int n){
  ofstream MyFile("data/" + file_name);

  MyFile << "n + 1 = " << n + 1 << endl;
  MyFile << setw(6) << "x" << setw(9) << "|" << setw(15)
         << "u(x) analytic" << setw(3) << "|" << setw(15)
         << "u(x) numeric" << endl;

  MyFile << 0 << setw(17) << 0 << setw(19) << 0 << endl;

  for (int i = 0; i < u_analytic.size(); i++){
    MyFile << setw(7) << setprecision (5)
        << scientific << x_val[i] << setw(17) <<
        u_analytic[i] << setw(19) << u_val[i][n - 1] << endl;
  }
  MyFile << 1 << setw(17) << 0 << setw(19) << 0 << endl;
  MyFile.close();
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

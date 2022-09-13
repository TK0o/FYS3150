#include <chrono>
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
vector<vector<float>> make_matrix_v2 (int n, int a, int b, int c);

void special_algorithm(vector<float>& f_val, int n);
void general_algorithm(vector<vector<float>>& A, int n);

void take_time(string file_name, int n, int repetitions, int n_up_to);

void print_m(vector<vector<float>>& vec);
void print_v(vector<float>& vec);

int main() {
  int n = 10;
  int repetitions = 10;
  string file_name = "task10.txt";

  take_time(file_name, n, repetitions, 4);
  return 0;
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

vector<vector<float>> make_matrix_v2 (int n, int a, int b, int c){
  vector<float> row(n, 0);
  vector<vector<float>> matrix_A(n - 1, row);
  int N = matrix_A.size();
  float h = 1 / float(n);

  for (int i = 0; i < N; i++){
    matrix_A[i][i] = b;
    matrix_A[i][n - 1] = exp(-10 * (i + 1) * h);

    if (i != N - 1){
      matrix_A[i][i + 1] = c;
    }
    if (i != 0) {
      matrix_A[i][i - 1] = a;
    }
  }
  return matrix_A;
}

void general_algorithm(vector<vector<float>>& A, int n){
  float diagonal_val, su_diagonal_val;
  int index, I;
  int N = n - 1;
  int iterations = 2 * (N - 1);

  for (int i = 0; i < iterations; i++){
    if (i < N - 1){
      I = i;
      index = I + 1;
    }
    else {
      I = 2 * (N - 1) - i;
      index = I - 1;
    }

    diagonal_val = A[I][I];
    su_diagonal_val = A[index][I];

    for (int j = 0; j < n; j++){
      A[I][j] = A[I][j] / diagonal_val;
      A[index][j] = A[index][j] - su_diagonal_val * A[I][j];
    }
  }

  float h = 1 / float(n);
  float const_h = 100 * pow (h, 2);
  for (int i = 0; i < N; i++){
    A[i][N] = A[i][N] * const_h;
  }
}

vector<float> analytical (vector<float>& x_val){
  vector<float> vec_analytical(x_val.size() - 1, 0);

  for (int i = 0; i < x_val.size() - 1; i++) {
      vec_analytical[i] = 1 - (1 - exp(-10)) * x_val[i] - exp(-10 * x_val[i]);
  }
  return vec_analytical;
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

void take_time(string file_name, int n, int repetitions, int n_up_to){
  float duration_seconds_general;
  float duration_seconds_special;
  int N;

  vector<float> row(2, 0);
  vector<vector<float>> Time_vector (n_up_to, row);

  for (int i = 1; i < n_up_to + 1; i++){
    N = pow (n, i);
    for (int j = 0; j < repetitions; j++){
      vector<float> x_val = get_x_values(N);
      vector<float> f_val = get_f_values(x_val);
      vector<vector<float>> A = make_matrix_v2(N, -1, 2, -1);

      auto t1 = std::chrono::high_resolution_clock::now();
      general_algorithm(A, N);
      auto t2 = std::chrono::high_resolution_clock::now();

      duration_seconds_general =
            std::chrono::duration<float>(t2 - t1).count();

      auto t3 = std::chrono::high_resolution_clock::now();
      special_algorithm(f_val, N);
      auto t4 = std::chrono::high_resolution_clock::now();

      duration_seconds_special =
            std::chrono::duration<float>(t4 - t3).count();

      A.clear();
      x_val.clear();
      f_val.clear();

      Time_vector[i - 1][0] += duration_seconds_general;
      Time_vector[i - 1][1] += duration_seconds_special;
    }
    Time_vector[i - 1][0] /= float(repetitions);
    Time_vector[i - 1][1] /= float(repetitions);
  }

  ofstream MyFile("data/" + file_name);
  MyFile << " base n = " << n << setw(2) << "|" << setw(19)
         << "N repetitions = " << repetitions << setw(4)
         << "|" << setw(23) << "Time given in seconds" << endl;

  MyFile << setw(8) << "n^x" << setw(6) << "|" << setw(22)
         << "average time general" << setw(3) << "|"
         << setw(22) << "average time special" << endl;

  for (int i = 0; i < n_up_to; i++){
    MyFile << defaultfloat << setw(7) << i + 1 << scientific
           << setprecision(5) << setw(24) << Time_vector[i][0]
           << setw(25) << Time_vector[i][1] << endl;
  }
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

#include <iostream>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <math.h>
#include <fstream>
using namespace std;

vector<float> get_x_values(int n);
vector<float> analytical (vector<float>& x_val);
vector<float> get_f_values(vector<float>& x_val);

vector<vector<float>> make_matrix_v2 (int n, int a, int b, int c);

void general_algorithm(vector<vector<float>>& A, int n);
void write_to_file(string file_name, int n, int n_up_to);

void print_m(vector<vector<float>>& vec);
void print_v(vector<float>& vec);


int main() {
  int n = 10;
  string file_name = "task7_b.txt";

  // Segmentation fault (core dumped) if int n_up_to is more than 4
  write_to_file(file_name, n, 4);
  return 0;
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

void write_to_file(string file_name, int n, int n_up_to){
  ofstream MyFile("data/" + file_name);
  int N = pow(n, n_up_to);

  for (int i = 1; i < n_up_to + 1; i++){
    N = pow (n, i);

    vector<float> x_val = get_x_values(N);
    vector<vector<float>> A = make_matrix_v2(N, -1, 2, -1);
    vector<float> exact = analytical(x_val);
    general_algorithm(A, N);

    MyFile << "n + 1 = " << N + 1 << endl;
    MyFile << setw(6) << "x" << setw(9) << "|"
           << setw(15) << "u(x) analytic"
           << setw(3) << "|" << setw(15)
           << "u(x) numeric" << endl;

    MyFile << 0 << setw(17) << 0 << setw(19) << 0 << endl;

    for (int j = 0; j < exact.size(); j++){
      MyFile << setw(7) << setprecision (5)
             << scientific << x_val[j] << setw(17)
             << exact[j] << setw(19) << A[j][N - 1] << endl;
    }

    MyFile << 1 << setw(17) << 0 << setw(19) << 0 << endl;
    MyFile << "-----------------------------------------------\n" << endl;
    A.clear();
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

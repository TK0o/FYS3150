#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <iomanip>
using namespace std;

vector<float> get_x_values(int n);
vector<float> analytical (vector<float>& x_val);
void write_to_file(vector<float> x_val, vector<float> u_val, string file_name);

int main() {
  string file_name = "task2.txt";
  int n = 100;

  vector<float> x_val = get_x_values(n);
  vector<float> exact = analytical(x_val);

  write_to_file(x_val, exact, file_name);
  return 0;
}


vector<float> get_x_values(int n){
  float h = 1 / float (n);
  vector<float> vec_for_x(n, 0);
  for (int i = 0; i < vec_for_x.size(); i++){
    vec_for_x[i] = (i + 1) * h;
  }
  return vec_for_x;
}

vector<float> analytical (vector<float>& x_val){
  vector<float> vec_analytical(x_val.size() - 1, 0);

  for (int i = 0; i < x_val.size() - 1; i++) {
      vec_analytical[i] = 1 - (1 - exp(-10)) * x_val[i]
                                 - exp(-10 * x_val[i]);
  }
  return vec_analytical;
}

void write_to_file(vector<float> x_val, vector<float> u_val, string file_name){
  ofstream MyFile("data/" + file_name);

  MyFile << setw(6) << "x" << setw(9) << "|" << setw(10)
         << "u(x)" << endl;
  MyFile << 0 << setw(17) << 0 << endl;

  for (int i = 0; i < u_val.size(); i++){
    MyFile << setw(7) << setprecision (5)
        << scientific << x_val[i] << setw(17) << u_val[i] << endl;
  }
  MyFile << 1 << setw(17) << 0 << endl;
  MyFile.close();
}

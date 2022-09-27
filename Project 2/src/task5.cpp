#include "../header/function.hpp"

int main()
{
    double eps = std::pow(10, -8);
    run_scaling_data(100, 20000, eps);
    return 0;
}

#include <iostream>
#include "matlib.h"
//========================================================================
// This example demonstrates solving the following simultaneous equations
// 2 * x1 + 8 * x2 + 5 * x3 = 5,
// 1 * x1 + 1 * x2 + 1 * x3 = -2,
// 1 * x1 + 2 * x2 - 1 * x3 = 2.
//========================================================================
int main() {
Matrix<3, 3, double> A;
Vector<3, double> x;
Vector<3, double> b;
// Write in Ax = b form
A = 2, 8, 5,
1, 1, 1,
1, 2, -1;
b = 5, -2, 2;
// solve for x
x = inverse(A) * b;
std::cout << "solution: " << transpose(x) << std::endl;
return 0;
}

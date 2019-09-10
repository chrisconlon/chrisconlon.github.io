#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

void swapmat(mat &X, mat &Y)
{
  mat tmp = ones<mat>(4,5);
  X=tmp;
}
int main(int argc, char** argv)
  {

  mat A ;
  mat B = randu<mat>(4,5);

  swapmat(A,B);
  cout << A << endl;
  cout << B << endl;

  return 0;
}



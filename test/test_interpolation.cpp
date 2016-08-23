/* 
  Testing the class Tinterpolation in interpolation.h
      
  Author: Martin Horvat, May 2013,  August 2016
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

#include "../src/matrix.h"
#include "../src/interpolation.h"

typedef double myreal;

int main(){
  
  const char *fileIN = "dataIN.dat";
  
  std::ofstream f(fileIN);
  
  f.precision(std::numeric_limits<myreal>::digits10 + 1);
  f << std::scientific;
  
  int i, 
      order = 6,
      N = 1000;
  
  myreal t, 
         T = 4*M_PI, 
         dt = T/N,
         **X = matrix<myreal>(N + 1, 3);

  for (i = 0; i <= N; ++i) {
    t = dt*i;
    X[i][0] = t;
    X[i][1] = std::sin(t);
    X[i][2] = std::cos(t);
    f << X[i][0] << ' ' << X[i][1] << ' ' << X[i][2] << '\n';
  }
  f.close();
  
  Tinterpolation<myreal> V(order, 2, N + 1, X, Tinterpolation<myreal>::General);
  const char *fileOUT = "dataOUT.dat";

  //Tinterpolation<myreal> V(order, 2, N + 1, X, Equid);
  //const char *fileOUT = "dataOUT2.dat";

  //Tinterpolation <myreal> V(order, 2, fileIN, General);
  //const char *fileOUT = "dataOUT3.dat";

  //Tinterpolation <myreal> V(order, 2, fileIN, Equid_General);
  //const char *fileOUT = "dataOUT4.dat";
  
  //std::cerr << "N=" << V.get_size() << '\n';
   
  free_matrix(X);

  int M = 1234;
  
  dt = T/M;
  
  f.open(fileOUT);
    
  // for OUT1,2,3,4
  for (t = dt; t < T; t+= dt) {
    f << t << ' ' 
      << V.get_value(t, 0, 0) << ' '                    // values X[0]
      << V.get_value(t, 1, 0) << ' '                    // values X[1] 
      << V.get_value(t, 0, 0) - std::sin(t) << ' '      // values X[0] - X_exact[0]
      << V.get_value(t, 1, 0) - std::cos(t) << '\n';    // values X[1] - X_exact[1]
  }   
  
  f.close();
  
  Tinterpolation <myreal> U(order, 2, fileIN, Tinterpolation<myreal>::Equid_General);
  const char *fileERR = "dataERR.dat";
  
  f.open(fileERR);
    
  // for dataERR
  for (t = dt; t < T; t += dt) {
    f << t << ' ' 
      << U.get_value(t, 0, 0) << ' '                    // values X[0]
      << U.get_value(t, 0, 0) - std::sin(t) << ' '      // values X[0] - X_exact[0]
      << U.get_value(t, 0, 1) - std::cos(t) << ' '      // values X'[0] - X_exact'[0]
      << U.get_value(t, 0, 2) + std::sin(t) << '\n';    // values X''[0] - X_exact''[0]
  }
  
  f.close();
  
  /* Checking with k = 6
    ./test_interpolation

    # Error in values of sin(t):
    gnuplot> plot 'dataERR.dat' u 1:3 w l
    #Finding: |abs. error| < 5e-16

    # Error in first derivative of sin(t):
    gnuplot> plot 'dataERR.dat' u 1:4 w l
    #Finding: outside boundaries |abs. error| < 2e-13

    # Error in second derivative of sin(t)
    gnuplot> plot 'dataERR.dat' u 1:5 w l
    #Finding: |abs. error| < 5e-11
  */
    
  return 0;
}

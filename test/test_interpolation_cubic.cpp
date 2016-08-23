/* 
  Testing the class Tinterpolation in interpolation.h
     
  Author: Martin Horvat, April 2013, August 2016
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

#include "../src/matrix.h"
#include "../src/interpolation_cubic.h"

typedef double myreal;

int main(){
  
  //
  // Prepare file in which we store latter interpolated vector function
  //
  
  const char *fileIN = "dataIN_cubic.dat";
  
  std::ofstream f(fileIN);
  
  f.precision(std::numeric_limits<myreal>::digits10 + 1);
  f << std::scientific;
  
  
  //
  // Define vector function
  //
  int i, N = 1000;
  
  myreal t, 
         T = 4*M_PI, 
         dt = T/N,
         **X = matrix<myreal>(N + 1, 3);

  for (i = 0; i <= N; ++i) {
    X[i][0] = (t = dt*i);
    X[i][1] = std::sin(t);
    X[i][2] = std::cos(t);
    f << X[i][0] << ' ' << X[i][1] << ' ' << X[i][2] << '\n';
  }
  f.close();
  
  //Tinterpolation_cubic <myreal> V(2, N + 1, X, General);
  //const char *fileOUT = "dataOUT1_cubic.dat";

  //Tinterpolation_cubic <myreal> V(2, N + 1, X, Equid);
  //const char *fileOUT = "dataOUT2_cubic.dat";

  //Tinterpolation_cubic <myreal> V(2, fileIN, General);
  //const char *fileOUT = "dataOUT3_cubic.dat";

  Tinterpolation_cubic <myreal> 
    V(2, fileIN, Tinterpolation_cubic <myreal>::Equid_General);
       
  const char *fileOUT = "dataOUT_cubic.dat";
  
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
  
  Tinterpolation_cubic <myreal> 
    U(2, fileIN, Tinterpolation_cubic <myreal>::Equid_General);
    
  const char *fileERR = "dataERR_cubic.dat";
  
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
  
  /* Checking
    ./test_interpolation_cubic

    # Error in values of sin(t):
    gnuplot> plot 'dataERR.dat' u 1:3 w l
    #Finding: |abs. error| < 6e-11

    # Error in first derivative of sin(t):
    gnuplot> plot 'dataERR.dat' u 1:4 w l
    #Finding: outside boundaries |abs. error| < 2e-8

    # Error in second derivative of sin(t)
    gnuplot> plot 'dataERR.dat' u 1:5 w l
    #Finding: |abs. error| < 1.5e-5
  */
    
  return 0;
}

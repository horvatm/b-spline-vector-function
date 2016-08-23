/* 
  Testing solver for the system of linear equations 
  
    w.x = b
  
  involving banded matrix w

  Author: Martin Horvat, April 2013, August 2016
*/


#include <cmath>
#include <iostream>
#include <cstdlib>

#include "../src/matrix.h"
#include "../src/solver_banded.h"

int main(){
  
  
  int n = 5,
      kl = 2,
      ku = 1,
      m = kl + ku + 1,
      i, j, iflag;
  
         
  //      
  // Define the right side
  //
  double *b = new double[n];
  for (i = 0; i < n; ++i) b[i] = i+1;
 
  //
  // Define the matrix
  //
  double **w = matrix<double>(m, n), a; 
  
  for (i = 0; i < m; ++i)
    for (j = 0; j < n; ++j) 
      w[i][j] = 0;
  
  
  for (i = 0; i < n; ++i) 
    for (j = 0; j < n; ++j) {
      if ( i == j ) 
        a = 1;
      else if (std::abs(i-j) == 1) 
        a = 0.5;
      else if (i-j == 2)   
        a = 0.25;
      else a = 0;
        
      if (a != 0) w[ku + i - j][j] = a;
    }
  
  //
  // Print matrix
  //
  for (i = 0; i < m; ++i){
    for (j = 0; j < n; ++j) std::cout <<  w[i][j] << '\t';
    std::cout << '\n';
  }    
  
  //
  // Solving banded matrix
  //
    
  banfac (w, m, n, kl, ku, iflag);
  
  std::cout << "iflag:" << iflag << '\n';
  
  banslv (w, m, n, kl, ku, b);
  
  //
  // Print solutions
  //
  for (i = 0; i < n; ++i) 
    std::cout << "b[" <<  i << "]= "<< b[i] << '\n';
  
  delete [] b;
  
  free_matrix(w);
    
  return EXIT_FAILURE;
}

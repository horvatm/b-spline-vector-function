/* 
  Testing bspline routines
  
  Compile:
    g++ test_bsplines.cpp -o  test_bsplines -O2 -Wall
     
  Author: Martin Horvat, May 2013
*/


#include <iostream>
#include <cmath>
#include <vector>

//#define LAPACK
//#define SPLINE_EVAL_ARRAYS

#include "../src/bsplines.h"

int main(){
  
  // **************************************************************************
  // **************************************************************************
  
  #if 0 // testing locate_knot
  
  std::vector<double> L;
  L.push_back(1);
  L.push_back(1);
  L.push_back(1);
  L.push_back(2);
  L.push_back(3);
  L.push_back(4);
  L.push_back(4);
  L.push_back(4);
  L.push_back(5);
  L.push_back(6);
  L.push_back(7);
  L.push_back(8);
  L.push_back(8);
  
  int i;
  
  std::cout.precision(16);
  std::cout << std::scientific;
  
  for (double x = 0; x <= 10; x += 0.2) {
    
    i = locate_knot(x, L);
    
    std::cout 
      << x << ' ' 
      << i << ' ' 
      << (i == -1 ? -1 : L[i]) << ' '
      << (i == -1 ? -1 : x-L[i]) << ' '
      << (i == -1 ? -1 : L[i] <= x && x < L[i+1]) << '\n';
  }
  
  #endif
  
  // **************************************************************************
  // **************************************************************************
  
  
  #if 0 // testing spline_eval, k (= order ) = 4, results: a.dat
  
  int k = 4;
  
  std::vector<double> L;
  L.push_back(1);
  L.push_back(1);
  L.push_back(1);
  L.push_back(1);
  L.push_back(2);
  L.push_back(3);
  L.push_back(4);
  L.push_back(5);
  L.push_back(6);
  L.push_back(7);
  L.push_back(8);
  L.push_back(8);
  L.push_back(8);
  L.push_back(8);
  
  std::vector<double> W ( L.size() - k);
  
  for (unsigned i = 0; i < W.size(); ++i) W[i] = 1;
  
  int l;
  
  double s;
  
  std::cout.precision(16);
  std::cout << std::scientific;
  
  
  double x, xmin, xmax;
  
  xmin = xmax = L[1];
  
  for (unsigned i = 1; i < L.size(); ++i) {
    xmin = std::min(xmin, L[i]);
    xmax = std::max(xmax, L[i]);
  }
  
  for (x = xmin; x <= xmax; x += 0.1) {
  
    l = locate_knot(x, L);
    
    if (l < 0) continue;
    
    s = spline_eval(l, x, k, W, L);
    
    std::cout << x << ' ' << s << '\n';
  }
  
  #endif
  
  // **************************************************************************
  // **************************************************************************
 
  #if 0 // testing spline_eval, k (= order ) = 4, results: b.dat
  
  int k = 4;
  
  std::vector<double> L;
  L.push_back(1);
  L.push_back(1);
  L.push_back(1);
  L.push_back(1);
  L.push_back(2);
  L.push_back(3);
  L.push_back(4);
  L.push_back(5);
  L.push_back(6);
  L.push_back(7);
  L.push_back(8);
  L.push_back(8);
  L.push_back(8);
  L.push_back(8);
  
  std::vector<double> W ( L.size() - k);
  
  int l;
  
  unsigned i, j;  
  
  double s, x, xmin, xmax;
  
  xmin = xmax = L[1];
  
  for (i = 1; i < L.size(); ++i) {
    xmin = std::min(xmin, L[i]);
    xmax = std::max(xmax, L[i]);
  }
  
  std::cout.precision(16);
  std::cout << std::scientific;
   
  for (i = 0; i < W.size(); ++i) {
    
    for (j = 0; j <  W.size(); ++j) W[j] = 0;
    W[i] = 1;
  
    for (x = xmin; x <= xmax; x += 0.01) {
    
      l = locate_knot(x, L);
      
      if (l < 0) {
        std::cerr << "Out of range\n";
        continue;
      }  
      
      s = spline_eval(l, x, k, W, L);
      
      std::cout << x << ' ' << s << '\n';
    }
    
    std::cout << '\n';
  }
  
  /* GNUPLOT: Plotting separate B-splines
   
   f(x1,n1,x2,n2,x3,n3,x4,n4,n)=(n1==n?x1 : (n2==n? x2 : (n3==n? x3 : (n4==n?x4 :1/0))))
   plot for [n=0:9] 'c.dat' u 1:(f($2,$3,$4,$5,$6,$7,$8,$9,n)) w l
  
  */ 
  #endif

  // **************************************************************************
  // **************************************************************************

  #if 0 // testing bspline_table, k (= order ) = 4, results: c.dat
  
  unsigned k = 4;
  
  std::vector<double> L;
  L.push_back(1);
  L.push_back(1);
  L.push_back(1);
  L.push_back(1);
  L.push_back(2);
  L.push_back(3);
  L.push_back(4);
  L.push_back(5);
  L.push_back(6);
  L.push_back(7);
  L.push_back(8);
  L.push_back(8);
  L.push_back(8);
  L.push_back(8);
  
  std::vector<double> B;
  
  int l;
  
  unsigned i, len = (k*(k+1))/2;  
  
  B.resize(len);
  
  double x, xmin, xmax;
  
  xmin = xmax = L[0];
  
  for (i = 1; i < L.size(); ++i) {
    xmin = std::min(xmin, L[i]);
    xmax = std::max(xmax, L[i]);
  }
  
  std::cout.precision(16);
  std::cout << std::scientific;
  
  for (x = xmin; x <= xmax; x += 0.01) {
  
    l = locate_knot(x, L);
    
    if (l < 0) {
      std::cerr << "Out of range\n";
      continue;
    }  

    bspline_table(l, x, k, L, B);
    
    std::cout << x ;
    
    // plotting point of N_{l-k+1,k}, .., N_{l,k}
    for (i = 1; i <= k; ++i) 
      std::cout << ' ' << B[len-i] << ' ' << l + 1 - i << ' ';
    
    std::cout << '\n';
    
  }
  
  #endif  
  
  
  // **************************************************************************
  // **************************************************************************
  
  #if 0 // testing spline_eval, k (= order ) = 4, results: d.dat
  
  unsigned i, k = 4;
  
  // L
  std::vector<double> L;
  L.push_back(1);
  L.push_back(1);
  L.push_back(1);
  L.push_back(1);
  L.push_back(2);
  L.push_back(3);
  L.push_back(4);
  L.push_back(5);
  L.push_back(6);
  L.push_back(7);
  L.push_back(8);
  L.push_back(8);
  L.push_back(8);
  L.push_back(8);
  unsigned nL = L.size();
  
  // W
  unsigned nW = nL - k;
  std::vector<double> W (nW);
  for (i = 0; i < nW; ++i) W[i] = i + 1;
  
  // B
  unsigned nB = (k*(k+1))/2;
  std::vector<double> B(nB);
  
  // x in [xmin, xmax]
  double x, xmin, xmax;
  xmin = xmax = L[0];
  for (i = 1; i < nL; ++i) {
    xmin = std::min(xmin, L[i]);
    xmax = std::max(xmax, L[i]);
  }
  
  // setting derivative orders
  std::vector <unsigned> P;
  P.push_back(0);
  P.push_back(1);
  P.push_back(2);
  unsigned nP = P.size();

  // calc. function s(x), s'(x), s''(x)
  std::vector <double> S(nP);
  
  std::cout.precision(16);
  std::cout << std::scientific;
  
  int l;
  
  for (x = xmin; x <= xmax; x += 0.01) {
  
    l = locate_knot(x, L);
    
    if (l < 0) {
      std::cerr << "Out of range\n";
      continue;
    }  

    bspline_table(l, x, k, L, B);
    
    spline_eval(l, k, B, W, L, P, S);
    
    // plotting s(x), s'(x), s''(x)
    std::cout << x ;
    for (i = 0; i < nP; ++i) std::cout << ' ' << S[i];
    std::cout << '\n';
    
  }
  
  #endif  
  
  
  // **************************************************************************
  // **************************************************************************
  
  #if 1 // testing spline_interpolate, k (= order ) = 6
  

  unsigned i, j, k = 6;
  
  // D
  unsigned dim = 1,
           nD = 20;
           
  double **D  = matrix <double>(nD, dim + 1);
  
  for (i = 0; i < nD; ++i) {
    D[i][0] = i;
    D[i][1] = std::sin(M_PI*i/15.);
  }     
  
  // L
  unsigned nL = nD + k;
  double *L = new double [nL];

  // W
  unsigned nW = nD;
  double **W = matrix<double> (dim, nW);
  
  // calculate the weights
  //spline_interpolate_lapack(k, nD, dim, D, W, L);   // results in e.dat
  spline_interpolate(k, nD, dim, D, W, L);          // results in f.dat
  
  // x in [xmin, xmax]
  double x, xmin, xmax;
  
  xmin = xmax = D[0][0];
  
  for (i = 1; i < nD; ++i) {
    xmin = std::min(xmin, D[i][0]);
    xmax = std::max(xmax, D[i][0]);
  }
 
  // B
  unsigned nB = (k*(k+1))/2;
  double *B = new double [nB];

  // setting derivative orders
  #if  defined(SPLINE_EVAL_ARRAYS)
  unsigned nP = 3,
      *P = new unsigned[nP];
  #else
  int nP = 3,
      *P = new int[nP];
  #endif
      
  P[0] = 0;
  P[1] = 1;
  P[2] = 2;
 
  // calc. function s(x), s'(x), s''(x)
  double *S = new double [nP];

  std::cout.precision(16);
  std::cout << std::scientific;
  
  int l;
  
  for (x = xmin; x <= xmax; x += 0.01) {
  
    l = locate_knot(x, nL, L);
    
    if (l < 0) {
      std::cerr << "Out of range\n";
      continue;
    }  

    bspline_table(l, x, k, L, B);
    
    // plotting s(x), s'(x), s''(x)
    std::cout << x ;
    
    for (i = 0; i < dim; ++i) {
      
      spline_eval(l, k, B, W[i], L, nP, P, S);
      
      for (j = 0; j < nP; ++j) std::cout << ' ' << S[j];
    } 
     
    std::cout << '\n';
  }
  
  delete [] L;
  free_matrix (D);
  free_matrix (W);
  delete [] B;
  delete [] P;
  delete [] S;
  
  #endif  
  
  
  
  return 0;
}

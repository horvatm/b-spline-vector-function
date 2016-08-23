#if !defined(__bsplines_h)
#define __bsplines_h

/*
  Library for spline interpolation of univariate vector functions using B-splines.
 
  Refs:
  * P. Dierckx, Curve and Surface Fitting with Splines, Clarendon Press, 1995
  * C. De Boor, On Calculating with B-splines, J. Approx. Theory 6, 50-62 (1972)

  Author: Martin Horvat, May 2013, August 2016
*/ 


#include <cmath>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cassert>

// allocation and freeing C-style matrix
#include "matrix.h"

// solving linear banded systems by DeBoor's algorithm or 
// by LAPACK using dgbsv_
#if !defined(LAPACK)
#include "solver_banded.h"
#endif

// enabling assert()
//#define NDEBUG 

// enabling bound check in std::vector
//#define _GLIBCXX_DEBUG


/*
  Locate knot index of the interval in which lies x. Performing a binary search.

  Input:
    x - argument
    L - knots = {lambda_i : i = 1, ... , n_L = n_W + k}
  
  Returns:
    l - 1  :  l in [lambda_l, lambda_{l+1})
    -1  :  otherwise
*/ 

template <class T> 
int locate_knot(const T &x, std::vector<T> &L){

  if (x < L.front()) return -1;
  
  if (x >= L.back()) return -1;
  
  static int lo = 0, hi = 1;
  
  if (!(L[lo] <= x && x < L[hi])) {
    
    int k;
    lo = 0;
    hi = L.size() - 1;

    while (hi - lo > 1){
      k = (hi + lo) >> 1;
      if (L[k] > x) hi = k; else lo = k;
    }
    
    assert(L[lo] <= x && x < L[hi]);
  }
  
  return lo;
}

template <class T> 
int locate_knot(const T &x, int len, T *L){

  if (x < L[0]) return -1;
  
  if (x >= L[len-1]) return -1;
  
  static int lo = 0, hi = 1;
  
  if (!(L[lo] <= x && x < L[hi])) {
    
    int k;
    lo = 0;
    hi = len - 1;

    while (hi - lo > 1){
      k = (hi + lo) >> 1;
      if (L[k] > x) hi = k; else lo = k;
    }
    
    assert(L[lo] <= x && x < L[hi]);
  }
  
  return lo;
}

/*
  Evaluate spline given in B-spline form.
    
    s(x) = sum_i w_i N_{i,k,L}(x)
  
  where w_i are weight to normalized B-splines N_{i,k,L} using de Casteljau 
  (de Boor) algorithm:
  
    s(x) = w_l^{[k-1]}(x)
  
  with recursion
  
    w_i^{[j]}= [ w_i      : j = 0
               [ u_i w_{i}^{[j-1]} + (1 - u_i) w_{i-1}^{[j-1]} : j > 0
  
    u_i = (x - lambda_i)/(lambda_{i+k-j} -lambda_{i}) 
  
  that generates triangular scheme (Dierckx, p. 14, Table 1.4; deBoor1972, p. 56) 
  
    (i,j) := w_i^{[j]}
    (i,0) := w_i
  
    (l-k+1,0)
    (l-k+2,0)   (l-k+2,1)
    (l-k+3,0)   (l-k+3,1)   .
                              .
        .           .           .
        .           .             
        .           .             (l-1,k-2)
    (l,0)       (l,1)             (l,k-2)     (l,k-1)
  
  and
    
    s(x) = (l,k-1)
    
  Input:
    l - index of the knot interval, x in [lambda_l, lambda_{l+1})
    x - argument
    k - order (= degree + 1)
    W - weights = {w_i : i = 1, ... , n_W}
    L - knots = {lambda_i : i = 1, ... , n_L = n_W + k}
    
  Return:
    s(x) = sum_i w_i N_{i,k,L}
  
  Following algorithm given in DeBoor1972, p. 54, 56 .
  
*/
 
template <class T> 
T spline_eval(
  const int &l,
  const T &x, 
  const int &k, 
  std::vector<T> &W, 
  std::vector<T> &L
){
  
  int r, s, 
      q,        // temporary variable
      d = k -1; // degree = order  - 1
  
  T *A0 = new T[k], 
    *A1 = new T[k], 
    u;
  
  q = l - d;
    
  for (r = 0; r < k; ++r) A0[r] = W[q + r];
    
  for (s = 0; s < d; ++s) {
    
    for (r = s + 1; r <= d; ++r) {
      
      q = l + r;
      
      u = (x - L[q - d])/(L[q - s] - L[q - d]);
      
      A1[r] = u*A0[r] + (1 - u)*A0[r - 1];
    }
    
    for (r = s + 1; r <= d; ++r) A0[r] = A1[r];
  }
  
  u = A0[k-1];
  
  delete [] A0;
  delete [] A1;
  
  return u;
}


template <class T> 
T spline_eval(
  const int &l,
  const T &x, 
  const int &k, 
  T *W, 
  T *L
){
  
  int r, s, 
      q,        // temporary variable
      d = k -1; // degree = order  - 1
  
  T *A0 = new T[k], 
    *A1 = new T[k], 
    u;
  
  q = l - d;
    
  for (r = 0; r < k; ++r) A0[r] = W[q + r];
    
  for (s = 0; s < d; ++s) {
    
    for (r = s + 1; r <= d; ++r) {
      
      q = l + r;
      
      u = (x - L[q - d])/(L[q - s] - L[q - d]);
      
      A1[r] = u*A0[r] + (1 - u)*A0[r - 1];
    }
    
    for (r = s + 1; r <= d; ++r) A0[r] = A1[r];
  }
  
  u = A0[k-1];
  
  delete [] A0;
  delete [] A1;
  
  return u;
}


/*
  Generate triangular scheme of non-zero normalized B-spline using the Cox - deBoor recursion:
  
    N_{i,j+1}(x) =  (x - lambda_{i}    )/(lambda_{i+j} - lambda_{i}    ) N_{i,j}(x)
                  + (lambda_{i+j+1} - x)/(lambda_{i+j+1} - lambda_{i+1}) N_{i+1,j}(x)
  
  with
  
    N_{i,1} (x) = [  1 : x in [lambda_i, lambda_{i+1})
                  [  0 : otherwise
  and
    
    N_{i,j}(x) = 0  if x not in  [lambda_i, lambda_{i+k}]
  
  The recursion generates the triangular scheme, which for l, such that x in [lambda_l, lambda_{l+1}) reads
  
    (i,j) := N_{i,j}    
  
    (l,1)   (l-1,2)   ...   (l-k+2,k-1)   (l-k+1,k)
            (l,2)     ...   (l-k+3,k-1)   (l-k+2,k)
                               
                  .             .             .
                    .           .             .
                      .         .             .
                        .    
                           
                              (l,k-1)     (l-1,k)
                                          (l,k)
  
  with  (l,1) = 1 (Dierckx, p. 13, Table 1.3). Other B-splines are zero! 
  Non-zero are only
    
    (i,j) : j = 1, ..., k;  i = l, ..., k - j +1
   
    
  Input:
    l - index of the knot interval, x in [lambda_l, lambda_{l+1})
    x - argument
    k - order = degree + 1
    L - knots = {lambda_i : i = 1, ... , n_L = n_W + k}
  
  Output:
    B[ - linearized table of B-spline values
      = { (l,1), 
          (l-1,2), (l,2), 
          (l-2,3), (l-1,3), (l,3), 
          ....,
          (l-k+1,k), (l-k+1,k), ... ,(l-1,k), (l,k) }
      
      = {B_i : i = 1, ..., n_B = k(k+1)/2 }
      
    B-splines of order i are available at
      
      B[j + i(i-1)/2] j = 0, ... , i -1
  
  Following algorithm given in DeBoor1972, p.58 .
*/

template <class T> 
void bspline_table(
  const int &l, 
  const T &x, 
  const unsigned &k,
  std::vector<T> &L,
  std::vector<T> &B
){
  
  int 
    s, r, 
    d = k - 1;  // degree = k - 1
    
  unsigned  o = 1,      // offset in B array  
            len = (k*(k+1))/2;
             
  T *N0 = new T [k],
    *N1 = new T [k],
    M, a, b;
  
  
  assert (B.size() >= len);
  
  B[0] = N0[0] = 1;
    
  for (s = 1; s <= d; ++s) {
    
    N1[0] = 0;
    
    for (r = 1; r <= s; ++r) {
      
      a = L[l + r];
      
      b = L[l + r - s];
      
      M = N0[r - 1]/(a - b);
      
      N1[r - 1] += (a - x)*M;
      
      N1[r] = (x - b)*M;
    } 
    
    for (r = 0; r <= s; ++r) B[o++] = N0[r] = N1[r];
  } 
  
  assert (o == len);
  
  delete [] N0;
  delete [] N1;
}

//#define TWO_ARRAY_CASE
#if defined(TWO_ARRAY_CASE)

template <class T> 
void bspline_table(
  const int &l, 
  const T &x, 
  const unsigned &k,
  T *L,
  T *B
){
  
  int 
    s, r, 
    d = k - 1;  // degree = k - 1
    
  unsigned  o = 1;      // offset in B array  
             
  T *N0 = new T [k],
    *N1 = new T [k],
    M, a, b;
  
  B[0] = N0[0] = 1;
    
  for (s = 1; s <= d; ++s) {
    
    N1[0] = 0;
    
    for (r = 1; r <= s; ++r) {
      
      a = L[l + r];
      
      b = L[l + r - s];
      
      M = N0[r - 1]/(a - b);
      
      N1[r - 1] += (a - x)*M;
      
      N1[r] = (x - b)*M;
    } 
    
    for (r = 0; r <= s; ++r) B[o++] = N0[r] = N1[r];
  } 
  
  assert (o == (k*(k+1))/2);
  
  delete [] N0;
  delete [] N1;
}

#elif defined(ONE_ARRAY)

template <class T> 
void bspline_table(
  const int &l, 
  const T &x, 
  const unsigned &k,
  T *L,
  T *B
){
  
  int s, r, 
      d = k - 1;                // degree = k - 1
    
  unsigned o = 1, o_ = 0;       // offset in B array  
             
  T M, a, b;
  
  B[0] = 1;
    
  for (s = 1; s <= d; ++s) {
    
    B[o + 0] = 0;
    
    for (r = 1; r <= s; ++r) {
      
      a = L[l + r];
      
      b = L[l + r - s];
      
      M = B[o_ + r - 1]/(a - b);
      
      B[o + r - 1] += (a - x)*M;
      
      B[o + r] = (x - b)*M;
    } 
    
    o_ += s;
    o  += s + 1;
  }
}
#else

template <class T> 
void bspline_table(
  const int &l, 
  const T &x, 
  const unsigned &k,
  T *L,
  T *B
){
  
  int s, r, 
      d = k - 1;                // degree = k - 1
    
  T M, a, b, *P = B + 1;
  
  *B = 1;     
  
  L += l;     // shifting point of knots
    
  for (s = 1; s <= d; ++s) {
    
    *P = 0;
    
    for (r = 1; r <= s; ++r) {
      
      a = L[r];
      
      b = L[r - s];
      
      M = *(B++)/(a - b);
      
      *(P++) += (a - x)*M;
      
      *P = (x - b)*M;
    }
    ++P;
  }
}


#endif


/*
  Calculate non-zero normalized B-splines for given x in [lambda_l, lambda_{l+1}).
  
    N_{l-k+1,k}, N_{l-k+1,k}, ... , N_{l-1,k}, N_{l,k}
     
  Input:
    l - index of the knot interval, x in [lambda_l, lambda_{l+1})
    x - argument
    k - order = degree + 1
    L - knots = {lambda_i : i = 1, ... , n_L = n_W + k}
  
  Output:
    N[k] = { N_{l-k+1,k}, N_{l-k+1,k}, ... , N_{l-1,k}, N_{l,k} }
        
  Following algorithm given in DeBoor1972, p.58.
*/ 

template <class T> 
void bspline_vals(
  const int &l, 
  const T &x, 
  const int &k,
  std::vector<T> &L,
  std::vector<T> &N
){
  
  int s, r, 
      q,          // temporary variable
      d = k - 1;  // degree = k - 1
      
  T *N1 = new T [k], M, a, b;
      
  assert (N.size() >= k);
  
  N[0] = 1;
  
  for (s = 1; s <= d; ++s) {
    
    N1[0] = 0;
    
    for (r = 1; r <= s; ++r) {
      
      a = L[l + r];
      b = L[l + r - s];
      
      M = N[r - 1]/(a -b);
      
      N1[r - 1] += (a - x)*M;
      
      N1[r] = (x - b)*M;
    } 
    
    for (r = 0; r <= s; ++r) N[r] = N1[r];
  } 

  delete [] N1;
}

template <class T> 
void bspline_vals(
  const int &l, 
  const T &x, 
  const int &k,
  T * L,
  T * N
){
  
  int s, r, 
      d = k - 1;  // degree = k - 1
  
      
  T *N1 = new T [k], 
    M, 
    a, b;
  
  
  L += l;     // shifting points of knots
  
  N[0] = 1;
  
  for (s = 1; s <= d; ++s) {
    
    N1[0] = 0;
    
    for (r = 1; r <= s; ++r) {
      
      a = L[r];
      b = L[r - s];
      
      M = N[r - 1]/(a - b);
      
      N1[r - 1] += (a - x)*M;
      
      N1[r] = (x - b)*M;
    } 
    
    for (r = 0; r <= s; ++r) N[r] = N1[r];
  } 

  delete [] N1;
}



/*
  Evaluate spline derivatives given in B-spline form
    
    s^{(p_i)}    p_i  in N_0 : i = 1, ..., n_P
    
  where w_i are weight to normalized B-splines N_{i,k,L} using formula
  
    s^{(nu)}(x) = [ prod_{i=1}^nu (k-i)]
                  sum_{i=l-k+1+nu} w_i^{(nu)} N_{i,k-nu, L}(x)
  
  where x in [lambda_l, lambda_{l+1}) and w_j^{(nu)} are obtained via 
  recursion

                [  w_j  :  i=0
    w_j^{(i)} = [
                [ (w_j^{(i-1)} - w_{j-1}^{(i-1)})/(lambda_{j+k-i} - lambda_j) : i > 0
  
  that generates a trapezoidal scheme (Dierckx, p. 15, Table 1.5):
  
    (j, i) := w_j^{(i)}
    (i, 0) := w_i
    
    (l-k+1,0)
    (l-k+2,0)   (l-k+2,1)
    (l-k+3,0)   (l-k+3,1)     .
                                .
      .             .             (l-k+1+nu, nu)
      .             .               .
      .             .               .
        
    (l,0)       (l,1)             (l,nu)
    
  Input:
    l - index of the knot interval, x in [lambda_l, lambda_{l+1})
    k - order (= degree + 1)
    B - table of B-splines evaluated x (from bspline_table () )
    W - weights = {w_i : i = 1, ... , n_W}
    L - knots = {lambda_i : i = 1, ... , n_L = n_W + k}
    P - set of orders of derivatives = {p_i : i = 1, ... , n_P}
        p_i < p_{i+1}
  Output:
  
    S = { s^{(p_i)}(x) : i = 1, ... , n_S = n_P}
  
  Following algorithm given in DeBoor1972, p.53.
*/

template <class T> 
void spline_eval(
  const int &l,
  const unsigned &k, 
  std::vector<T> &B,
  std::vector<T> &W, 
  std::vector<T> &L,
  std::vector<unsigned> &P,
  std::vector<T> &S
){
  
  assert(S.size() >= P.size());
  
  //
  // get max derivative order, Pmax
  //
  
  unsigned i, 
           Pmax = P[0],
           nP = P.size();     // number of values
      
  for (i = 1; i < nP; ++i) Pmax = std::max(Pmax, P[i]);
  
  //
  // determine the size vector A and reserve space 
  //
  unsigned nA  = 0;
  
  for (i = 0; i <= Pmax; ++i) nA += k - i; 
    
  T *A = new T [nA];
  
  
  // 
  // calculate coefficients A == w_j^{(i)}
  //
  unsigned d = k - 1,         // degree = k - 1  
           q, p, r, s, z;
  
  q = l - d;           
  for (i = 0; i < k; ++i) A[i] = W[q + i];
  
  s = 0;
  q = k;
  for (p = 1; p <= Pmax; ++p){
    
    r = d - p;
    z = l + r;
    
    for (i = l; i <= z; ++i, ++q, ++s)
      A[q]  = (A[s + 1] - A[s])/(L[i + 1] - L[i - r]);
    
    ++s;
  }  
  
  
  //
  // calculate values of derivatives
  //
  unsigned oA = 0,            // offset in A
           oB = (k*(k-1))/2;  // offset in B
           
  T sum,
    f = 1;                    // factor
  
  for (i = p = 0; i < nP; ++i) {

    // calculate factor and offsets
    for (r = P[i]; p < r; ++p) {
      q = d - p;
      
      f *= q;      
      oA += q + 1;  
      oB -= q;
    }
    
    assert(p == P[i]);
    
    // SUM = sum_{q = l - k + p}^l A_q^{(p)} N_{q, k - p}(x)
    sum = 0;
    for (q = 0, r = d - p; q <= r; ++q) sum += A[oA + q]*B[oB + q];
    
    // result: factor*SUM 
    S[i] = f*sum;
  }
  
  delete [] A;
}


#if defined(SPLINE_EVAL_ARRAYS)

template <class T> 
void spline_eval(
  const int &l,
  const unsigned &k, 
  T *B,
  T *W, 
  T *L,
  const unsigned &nP,
  unsigned *P,
  T *S
){
  
  //
  // get max derivative order, Pmax
  //
  
  unsigned i, 
           Pmax = P[0];
      
  for (i = 1; i < nP; ++i) Pmax = std::max(Pmax, P[i]);
  
  //
  // determine the size vector A and reserve space 
  //
  //unsigned nA  = 0;
  //for (i = 0; i <= Pmax; ++i) nA += k - i; 
  unsigned nA = ((Pmax + 1)*(2*k - Pmax))/2;
  T *A = new T [nA];
  
  
  // 
  // calculate coefficients A == w_j^{(i)}
  //
  unsigned d = k - 1,         // degree = k - 1  
           q, p, r, s, z;
           
  q = l - d;           
  for (i = 0; i < k; ++i) A[i] = W[q + i];
    
  s = 0;
  q = k;
  for (p = 1; p <= Pmax; ++p){
    
    r = d - p;
    z = l + r;
    
    for (i = l; i <= z; ++i, ++q, ++s)
      A[q]  = (A[s + 1] - A[s])/(L[i + 1] - L[i - r]);
    
    ++s;
  }  
  
  
  //
  // calculate values of derivatives
  //
  unsigned oA = 0,            // offset in A
           oB = (k*(k-1))/2;  // offset in B
  
  T sum,
    f = 1;                    // factor
    
  for (i = p = 0; i < nP; ++i) {

    // calculate factor and offsets
    for (r = P[i]; p < r; ++p) {
      q = d - p;
      
      f *= q;      
      oA += q + 1;  
      oB -= q;
    }
    
    assert(p == P[i]);
    
    // SUM = sum_{q = l - k + p}^l A_q^{(p)} N_{q, k - p}(x)
    sum = 0;
    for (q = 0, r = d - p; q <= r; ++q) sum += A[oA + q]*B[oB + q];
    
    // result: factor*SUM 
    S[i] = f*sum;
  }
  
  delete [] A;
}

#else   // using pointers to address A, B

template <class T> 
void spline_eval(
  const int &l,
  const int &k, 
  T *B,
  T *W, 
  T *L,
  const int &nP,
  int *P,
  T *S
){
  
  //
  // get max derivative order, Pmax
  //
  
  int i, 
      Pmax = P[0];
                
  for (i = 1; i < nP; ++i) Pmax = std::max(Pmax, P[i]);
  
  //
  // determine the size vector A and reserve space 
  //
  int nA = ((Pmax + 1)*(2*k - Pmax))/2;
  T *A = new T [nA];
  
  // 
  // calculate coefficients A == w_j^{(i)}
  //
  int d = k - 1,         // degree = k - 1  
      q, p, r, s;
           
  q = l - d;           
  for (i = 0; i < k; ++i) A[i] = W[q + i];
  
  L += l;         // shifting points of knots
    
  s = 0;
  q = k;
  for (p = 1; p <= Pmax; ++p){
    
    r = d - p;
    
    for (i = 0; i <= r; ++i, ++q, ++s)
      A[q] = (A[s + 1] - A[s])/(L[i + 1] - L[i - r]);
    
    ++s;
  }  
  
  
  //
  // calculate values of derivatives
  //
  
  T sum,
    f = 1,                        // factor
    *pA = A,                      // pointer to A
    *pB = B + (k*(k-1))/2;        // pointer to B  
      
  for (i = p = 0; i < nP; ++i) {

    // calculate factor and offsets
    for (r = P[i]; p < r; ++p) {
      q = d - p;
      
      f *= q;
      
      pA += q + 1;  
      pB -= q;
    }
    
    assert(p == P[i]);
    
    // SUM = sum_{q = l - k + p}^l A_q^{(p)} N_{q, k - p}(x)
    sum = 0;
    for (q = 0, r = d - p; q <= r; ++q) sum += pA[q]*pB[q];
    
    // result: factor*SUM 
    S[i] = f*sum;
  }
  
  delete [] A;
}

#endif

/*
  Calculate a not-a-knot interpolation based on given data
    
    D = X x Y = { (x_i, y_i) in T x T^{dim} :  i = 1, ..., n_D}
  
  for EVEN ORDER k = 2m !!!!
  
  Input:
    k - order = degree + 1
        
    D - data = {(x_i, y_{i,j} : i = 1, ..., n_D, j = 1, ... , dim}
      = Y[n_D][dim + 1]
  
  Output:
    W - weights = {w_{i,j} : i = 1,..,dim, i = 1, ... , n_W = n_D}
      = W[dim][n_W]
      
    L - knots = {lambda_i: i = 1, ... , n_L = n_W + k}
      = L[n_L]
  Note:  
    L = {x_1, ..., x_1, x_{m+1}, ..., x_{n_D-m}, x_{n_D}, ..., x_{n_D}}
         -------------  ----------------------   ---------------------
            k times             n_D - k              k times

  We formulate spline
  
    S(x) = sum_{i=0}^{nD-1} w_i N_i(x)    w_i in T^{dim}

  The interpolation condition
    
    S(x_i) = y_i    i = 0, .., nD-1
  
  results in system of equations
    
    sum_{j=0}^{nD-1} w_j N_j(x_i)  = y_i    y_i in T^{dim}
  
  that can be written in matrix form as
  
      B W^T = Y     Y=[y_{i,j}]_{i,j}   B = [N_j(x_i)]_{i,j}
  
  Collocation matrix B is stored in banded format 
  
    http://www.netlib.org/lapack/lug/node124.html
  
  and the system is solved using  
  
    http://www.netlib.org/lapack/double/dgbsv.f
  
  Input:
    k -- order
    nD -- nr. of data points
    dim -- dimension of data points
    D[nD][dim+1] -- data
    
  Output:
    W[dim][nD] -- weights
    L[nL = nD + k] -- knots 
*/

#if defined(LAPACK)
extern "C" int dgbsv_(
  int *n, int *kl, int *ku, int *nrhs, 
  double *ab, int *ldab, int *ipiv, double *b, 
	int *ldb, int *info);
  
template <class T>
bool spline_interpolate_lapack(
  int k, 
  int nD,
  int dim,
  T ** D,
  T ** W,
  T * L
) {

  // interpolation supported for even orders 
  assert(k % 2 == 0);
  
  // create knots L for not-a-know boundary condition
  int 
    i,
    m = k/2,
    d = k - 1;    

  T a = D[0][0],
    b = D[nD-1][0];
    
  for (i = 0; i < k; ++i){
    L[i] = a;
    L[nD + i] = b;
  }
  
  for (i = k; i < nD; ++i) L[i] = D[i - m][0];
   //
  // defining collocation matrix B as required by LAPACK routine dgbsv
  //  
  int l, r, s, j,
      kl = k - 1,           // nr. of subdiagonals in B
      ku = k - 1,           // nr. of superdiagonals in B
      ldab = 2*kl + ku + 1; // leading dimension of A
  
  T **A = matrix <T> (nD, ldab); // banded format of B in column-first notation
  
  for (i = 0; i < nD; ++i)
    for (j = 0; j < ldab; ++j) 
      A[i][j] = 0;
  
  //
  // calculating the collocation matrix
  //
  T *N = new T [k];
  
  for (i = 0; i < nD; ++i) {
  
    // calc. l such that x_i in [L[l], L[l+1]) 
    // except if(x = x_{nD-1}) then l = nD - 1
    //
    // L = {x_0, ..., x_0, x_{m}, ..., x_{n_D-m-1}, x_{n_D-1}, ..., x_{n_D-1}}
    //       0 , ..., k-1,   k,   ...,    nD -1,       nD,     ...,  nD + k -1 
    //       -------------  ----------------------   ------------------------
    //          k times             n_D - k                   k times
    
    if (i < m)          
      l = d;
    else if (i < nD - m)  
      l = i + m;
    else 
      l = nD -1;
    
    // calc { N_{l-k+1,k}, N_{l-k+1,k}, ... , N_{l-1,k}, N_{l,k} }
    bspline_vals(l, D[i][0], k, L, N);
    
    s = l - d; 
    
    assert(s >= 0);
    
    // store collocation matrix B in column-first banded format 
    //    B[i][s + r] = N[r]   r = 0, ..., k - 1
    // from full to banded format (in C)
    //    A [kl + ku + ii - jj][jj] = B[ii][jj]  ii = i;  jj = s + r
    // and transposing to comply with Fortran
    //    "A = B^T"
    j = ku + kl + i;
    for (r = 0 ; r < k; ++r) A[s + r][j - s - r] = N[r];
  }
  
  delete [] N;
    
  // right side of equation B W = Y in column first notation
  for (i = 0; i < dim; ++i) 
    for (j = 0; j < nD; ++j) 
      W[i][j] = D[j][i+1];

  // calculate weights W by solving B W^T = Y
  int *ipiv = new int [nD], info;
 
  dgbsv_(&nD, &kl, &ku, &dim, A[0], &ldab, ipiv, W[0], &nD, &info);
  
  free_matrix(A);  
  
  delete [] ipiv;
  
  if (info < 0){
    std::cerr 
    << "spline_interpolate_lapack: the " << -info 
    << "-th argument had an illegal value!\n";
  } else if (info > 0){
    std::cerr 
    << "spline_interpolate_lapack: the factor U is exactly singular.\n";
  }
  return (info==0);
}
#endif // if defined(LAPACK)

//
// For solving the system  B W^T = Y we use deBoor's method in solver_banded.h
//

template <class T>
bool spline_interpolate(
  int k, 
  int nD,
  int dim,
  T ** D,
  T ** W,
  T * L
) {

  // interpolation supported for even orders 
  assert(k % 2 == 0);
  
  // create knots L for not-a-know boundary condition
  int 
    i,
    m = k/2,
    d = k - 1;    

  T a = D[0][0],
    b = D[nD-1][0];
    
  for (i = 0; i < k; ++i){
    L[i] = a;
    L[nD + i] = b;
  }
  
  for (i = k; i < nD; ++i) L[i] = D[i - m][0];
 
  // defining collocation matrix B as required by banfac and banslv
 
  
  int l, r, s, j,
      kl = k - 1,               // nr. of subdiagonals in B
      ku = k - 1,               // nr. of superdiagonals in B
      nB = kl + ku + 1;         // leading dimension of A
  
  T **A = matrix <T> (nB, nD);  // banded format of B 
  
  for (i = 0; i < nB; ++i)
    for (j = 0; j < nD; ++j) 
      A[i][j] = 0;
  
  //
  // calculating the collocation matrix
  //
  T *N = new T [k];
  
  for (i = 0; i < nD; ++i) {
  
    // calc. l such that x_i in [L[l], L[l+1]) 
    // except if(x = x_{nD-1}) then l = nD - 1
    //
    // L = {x_0, ..., x_0, x_{m}, ..., x_{n_D-m-1}, x_{n_D-1}, ..., x_{n_D-1}}
    //       0 , ..., k-1,   k,   ...,    nD -1,       nD,     ...,  nD + k -1 
    //       -------------  ----------------------   ------------------------
    //          k times             n_D - k                   k times
    
    if (i < m)          
      l = d;
    else if (i < nD - m)  
      l = i + m;
    else 
      l = nD -1;
    
    // calc { N_{l-k+1,k}, N_{l-k+1,k}, ... , N_{l-1,k}, N_{l,k} }
    bspline_vals(l, D[i][0], k, L, N);
    
    s = l - d; 
    
    assert(s >= 0);
    
    // store collocation matrix B in banded format 
    //    B[i][s + r] = N[r]   r = 0, ..., k - 1
    // from full to banded format (in C)
    //    A [ku + ii - jj][jj] = B[ii][jj]  ii = i;  jj = s + r
  
    j = ku + i;
    for (r = 0 ; r < k; ++r) A[j - s - r][s + r] = N[r];
  }
  
  delete [] N;
  
  int info;
  
  //
  // making the LU decomposition
  //
  
  banfac (A, nB, nD, kl, ku, info);
  
  if (info != 1){
    std::cerr << "spline_interpolate: LU decomposition did not succeed\n";
    return  false;
  }  
    
  // using LU decomposition to solve B W_i = Y_i
  for (i = 0; i < dim; ++i) {
    
    for (j = 0; j < nD; ++j) W[i][j] = D[j][i+1];
    
    banslv (A, nB, nD, kl, ku, W[i]);
  }
  
  free_matrix(A);
  
  return true;
}


#endif // #if !defined(__bsplines_h)

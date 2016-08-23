#if !defined(__interpolation_h)
#define __interpolation_h

/*
  A class for calculation of univariate spline of vector functions

     y(t) in T^dim      (e.g. T = Real )
  
  and its derivatives.

  The data is given in the general form
    
    (t_i, y(t))           i = 0, ..., N-1
  
  with t_{i+1} - t_i = dt_i > 0 or on a lattice 
  
    y(t_i)                i = 0, ..., N-1
    t_i = tmin + i*dt
  
  
  Analysis is based on univariate B-splines of degree n.
  
  
Ref:
  * http://en.wikipedia.org/wiki/B-spline
    http://mathworld.wolfram.com/B-Spline.html
  
  * Carl De Boor - A Practical Guide to Splines (Springer, 2001)
    http://pages.cs.wisc.edu/~deboor/
    http://pages.cs.wisc.edu/~deboor/pgs/
    
  * Paul Dierckx - Curve and Surface Fitting With Splines (Clarendon, 1993)
    http://www.netlib.org/dierckx/
  
  Author: Martin Horvat, May 2013, August 2016
*/ 

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <cassert>
#include <utility>
#include <limits>

// B-spline interpolation
#include "bsplines.h"

// allocation and freeing of C-style matrices
#include "matrix.h"

// enabling assert()
//#define NDEBUG

template <class T>
class Tinterpolation {
  
  public:
  
  // supported definitions of timestamps
  enum Ttime {
      Equid,          // Equidistant: t_i = t_0 + h*i
      Equid_General,  // Equidistant general: t_{i+1} - t_i = h.
      General         // General
    };

  protected:
  
  bool init_status;
  
  int k,    // order (=degree+1) == even of polynomial used in interpolation
      m,    // = k / 2
      d,    // dimension of points
      nD,   // number of data points
      
      lo,   // interval variables of locate_knot
      hi,   // interval variables of locate_knot
      
      nB,   // number of entries in B
      nL;   // number of knots L = nD + k
      
  Ttime time_type;
  
  
  T tmin,   // minimal time/dependent variable
    tmax,   // maximal time/dependent variable
    dt,     // in case of equidistant data: t_{1} - t_{0}
    
    x_old,  // last used argument
    
    *ofs,   // offset to all dimension of data ofs[d+1]
    
    *B,     // B-spline table
    *L,     // B-spline knots  L[nL]
    **W;    // B-spline weights W[d][nD]
    
  int  *P;     // orders of derivatives
  
  int locate_knot(const T &t){
    
    if (t < tmin) return -1;

    if (t > tmax) return -1;
    
    if (t == tmax) {
      lo = nD - 1;  // considering comment by Dierckx
      hi = nD;
      return lo;
    }
    
    if (!(L[lo] <= t && t < L[hi])) {    
    
      int i;
      
      switch (time_type){
          
        case General:
      
          lo = 0;
          hi = nL - 1;

          while (hi - lo > 1){
            i = (hi + lo) >> 1;
            if (L[i] > t) hi = i; else lo = i;
          }
      
          assert(L[lo] <= t && t < L[hi]);
                
          break;
        
        default:

          i = int ((t - tmin)/dt);
          
          if (i < m)
            lo = k - 1;
          else if (i < nD - m)  
            lo = i + m;
          else 
            lo = nD - 1;
            
          hi = lo + 1;
      }
    }
    
    return lo;
  }

  void set_offset(int d, T *offset) {
    
    int i; 
    
    ofs = new T [d + 1];
    
    if (offset) {
      for (i = 0; i <= d; ++i) ofs[i] = offset[i];
    } else {
      for (i = 0; i <= d; ++i) ofs[i] = 0;
    }
  }
  
  void prepare_spline(int k, int d, int nD,  T **D){
    
    // k is even, m is half of order
    assert (k % 2 == 0);
    
    m = k/2;
    
    // L -- knots
    nL = nD + k;
    L = new T [nL];
    
    // W -- weights 
    W = matrix<T> (d, nD);
    
    if (!spline_interpolate(k, nD, d, D, W, L)){
      std::cerr 
        << "TInterpolation::prepare_spline: Spline interpolate failed!\n";
      exit(EXIT_FAILURE);
    }
    
    // B -- B-spline tables
    nB = (k*(k+1))/2;
    B = new T [nB]; 
    
    // P -- order of derivatives
    P = new int [d];
    
    // old arguments to unreachable
    x_old = std::numeric_limits<T>::max();
  }
  
  public:
  
  Tinterpolation() : init_status(false) {}
  
  /* Interpolation of different types of data stored in a file.  
   
     time_type = General:
    
      In file <filename> we have columns
      
        t_i X[i][0] ...X[i][d-1]        i = 0, ..., N-1
      
      where t_{i+1} > t_i.
      
     time_type = Equid_General:
    
     In file <filename> we have columns
      
        t_i X[i][0] ...X[i][d-1]        i = 0, ..., N-1
      
      where t_{i+1} - t_i = h. h is calculated from the first two rows.
    
     time_type = Equad:
 
      In file <filename> we have in the first column
        
        t0  h
        
      and continue next rows with columns
      
        X[i][0] ...X[i][d-1]                        i = 0, ... , N-1
      
      where X(t_i) = X[i] and t_i = t_0 + i*dt.
      
    Input:
      k -- order(=degree + 1)== even of polynomials
      d -- dimension of vector function
      filename -- filename of the file with data
      time_type --  type of time/independent variable
      offset[d+1]  -- offset of all variables       
  */
  
  Tinterpolation(
    const int &k,             
    const int &d,             
    const char *filename,    
    const Tinterpolation<T>::Ttime &time_type,  
    T* offset = 0             
  ) : k(k), d(d), lo(0), hi(1), time_type(time_type) {
    
    // setting ofs 
    set_offset(d, offset);   
    
     // read data
    std::ifstream f(filename);
    
    if (!f.is_open()) {
      std::cerr << "Error opening file " << filename << '\n';
      exit(EXIT_FAILURE);
    }
    
    int i, j, d1 = d + 1;

    T t;
    
    std::vector<T> buf;
    
    switch (time_type) {
      
      case Equid:
        f >> tmin >> dt;
        i = j = 0;
        while (f.good())
          if (f >> t) {
            if (++i == d1) {
              buf.push_back(tmin + dt*j);
              i = 0;
              ++j;
            }
            buf.push_back(t);
          }
          
      break;
      
      default:        
        while (f.good()) if (f >> t) buf.push_back(t);
    }
    
    // determine length of data
    nD = buf.size();
    assert (nD % d1 ==0);
    nD /= d1;

    // format data matrix D
    T **D = new T* [nD];
    D[0] = &buf[0];
    for (i = 1; i < nD; ++i) D[i] = D[i-1] + d1;

    // apply offset to data
    if (offset)
      for (i = 0; i < nD; ++i)
        for (j = 0; j <= d; ++j)
          D[i][j] -= ofs[j];
    
    // time axis attributes
    tmin = D[0][0];
    tmax = D[nD-1][0];
    dt = D[1][0] - tmin;
     
    // prepare spline data
    prepare_spline(k, d, nD, D);  // setting nL, L, W, nB, B
   
    delete [] D;
      
    init_status = true;
  }
  
  /* 
    Interpolation of equidistant data ( time_type = Equid)
      
      X(t_i) = (X_0(t_i) , ... , X_{d-1}(t_i)) 
             = X[i] in T^d                        i = 0, ... , nD-1
       
    with t_i = t_0 + i*dt.
     
    Input:
      k -- order (=degree+1) == even of polynomials 
      d -- dimension of vector function
      nD -- number of data points
      X[nD][d] -- data
      t0 -- initial time
      dt -- time step
      offset[d+1] -- offset of all variables
  */
  
  Tinterpolation(
    const int &k,
    const int &d,
    const int &nD,
    T **X,
    const T &t0,
    const T &dt,
    T* offset = 0
  ) : k(k), d(d), nD(nD), lo(0), hi(1), 
      time_type (Tinterpolation<T>::Equid), tmin(t0), dt(dt) {
    
    // setting ofs 
    set_offset(d, offset);   
    
    // read data  
    int i, j;
    
    T **D = matrix<T> (nD, d + 1);
    
    for (i = 0; i < nD; ++i) {
      D[i][0] = tmin + i*dt;
      for (j = 0; j < d; ++j) D[i][j+1] = X[i][j];
    }
    
    // apply offset to data
    if (offset) 
      for (i = 0; i < nD; ++i)
        for (j = 0; j <= d; ++j)
          D[i][j] -= ofs[j];
    
    // time axis attributes (could be changes by offsets)
    if (offset) tmin = D[0][0];
    tmax = D[nD - 1][0];
    
    // prepare spline data
    prepare_spline(k, d, nD, D);  // setting nL, L, W, nB, B
   
    free_matrix(D);
    
    init_status = true;   
  }
  
  /* 
    Interpolation of data stored in general format 
    
       Y[N][d+1]
    
    where    
    
       Y[i] = {t_i, X(t_i)} in T^(d+1)          i = 0, ... , N-1
    
    with t_{i} < t_{i+1}. 
    
    time_type = Equid_General
    
      Equidistant time dt = t_{i+1} - t_i
    
    time_type = General
    
      General time t_{i+1} > t_{i}
      
    Input:
      n -- order (=degree+1) == even of polynomials 
      d -- dimension of vector function
      nD -- number of points
      Y[nD][d+1] -- data
      time_type -- type of time/independent variable
      offset[d+1]  -- offset of all variables
  */ 
  
  Tinterpolation(
    const int &k,
    const int &d,
    const int &nD,
    T **Y,
    const Tinterpolation<T>::Ttime & time_type,
    T* offset = 0
  ) : k(k), d(d), nD(nD), lo(0), hi(1), time_type(time_type) {
    
    // setting ofs 
    set_offset(d, offset);
    
    T **D;
    
    // copy data and apply offset to data
    if (offset) {
      
      int i, j;
      
      D = matrix<T>(nD, d + 1);
      
      for (i = 0; i < nD; ++i)
        for (j = 0; j <= d; ++j)
          D[i][j] = Y[i][j] - ofs[j];
      
    } else D = Y;
    
    // time axis attributes
    tmin = D[0][0];
    tmax = D[nD - 1][0];
    dt = D[1][0] - tmin;
        
    // prepare spline data
    prepare_spline(k, d, nD, D); // setting nL, L, W, nB, B
    
    if (offset) free_matrix(D);  
     
    init_status = true;
  }
  
  ~Tinterpolation(){
    if (init_status){
      delete [] ofs;
      delete [] L;
      delete [] B;
      delete [] P;
      free_matrix(W);
    }
  }
  
  /*
    Return order
  */ 
  int get_order() const {return k;}
  
  /* 
    Returns dimension
  */
  int get_dim() const {return d;}
  
  /* 
    Returns the number of points 
  */
  int get_size() const {return nD;}
  
  /* 
    Returns minimal time
  */
  T get_min_time() const { return tmin + ofs[0]; }
  
  /* 
    Returns maximal time
  */
  T get_max_time() const { return tmax + ofs[0]; }
  
  /* 
    Calculate interpolated component of a vector
    
    Input:
      t - time
      i - index of the component
      j - order of the derivative
      
    Return: 
      X_i^{(j)}(t)    i in [0, d-1],  j in [0,n-1]
  */
  
  T get_value(const T & t, const int & i, const int &j){
    
    if (init_status){
    
      T x  = t - ofs[0], s;

      int l = locate_knot(x);
      
      if (l < 0){
        std::cerr 
          << "Tinterpolation::get_value: Index out of bounds\n"
          << "value=" << t << '\n';
        return 0;
      }
      
      if (j == 0){  // calculating only values
        
        s = spline_eval(l, x, k, W[i], L);
        
      } else {
        
        // values of B-splines
        bspline_table(l, x, k, L, B);
                
        // orders of derivatives
        P[0] = j;
                
        // evaluate derivatives      
        spline_eval(l, k, B, W[i], L, 1, P, &s);
        
      }
      
      return s;
      
    } else {
      std::cerr << "Tinterpolation::get_value: Not initialized\n";
      exit(EXIT_FAILURE);
    }
    return 0;
  }
  
  /* 
    Calculate interpolated vector functions
    
    Input:
      t - time
      j - order of the derivative
      
    Output: 
      v[d]:  v[i] = y_i^{(j)}(t)    i in [0, d-1],  j in [0, k-1]
  */
  
  void  get_value(const T & t, const int &j, T *v){
    
    if (init_status){
      
      T x  = t - ofs[0];
      
      int l = locate_knot(x);
      
      if (l < 0){
        std::cerr 
          << "Tinterpolation::get_value: Index out of bounds\n"
          << "value=" << t << '\n';
        return;
      }
           
      if (j == 0){  // calculating only values
        
        for (int i = 0; i < d; ++i)
          v[i] = ofs[i+1] + spline_eval(l, x, k, W[i], L);
        
      } else {
        
        // values of B-splines
        if (x_old != x) bspline_table(l, x, k, L, B);
        x_old = x;
                
        // orders of derivatives
        P[0] = j;
                
        // evaluate derivatives      
        for (int i = 0; i < d; ++i)
          spline_eval(l, k, B, W[i], L, 1, P, v+i);
        
      }
      
    } else {
      std::cerr << "Tinterpolation::get_value: Not initialized\n";
      exit(EXIT_FAILURE);
    }
  }
 
   /* 
    Calculate interpolated vector functions
    
    Input:
      t - time
      jmax - order of the derivative
      
    Output: 
      v - leading row matrix
      v[d][jmax+1]:  v[i][j] = y_i^{(j)}(t)    i in [0, d-1],  j in [0,jmax]
  */
  
  void  get_value(const T & t, const int &jmax, T **v){
    
    if (init_status){
      
      T x  = t - ofs[0];
      
      int l = locate_knot(x);
      
      if (l < 0){
        std::cerr 
          << "Tinterpolation::get_value: Index out of bounds\n"
          << "value=" << t << '\n';
        return;
      }
           
      if (jmax == 0){  // calculating only values
        
        for (int i = 0; i < d; ++i) 
          v[i][0] = ofs[i+1] + spline_eval(l, x, k, W[i], L);
        
      } else {
        
        // values of B-splines
        if (x_old != x) bspline_table(l, x, k, L, B);
        x_old = x;
        
        int i, j;
        
        // orders of derivatives
        for (j = 0; j <= jmax; ++j) P[j] = j;
        
        // evaluate derivatives
        T *u = v[0];     
        
        j = jmax + 1;
              
        for (i = 0; i < d; ++i) {
          
          spline_eval(l, k, B, W[i], L, j, P, u);
          
          u += j;
        }
        
        for (i = 0; i < d; ++i) v[i][0] += ofs[i+1];
      }
      
    } else {
      std::cerr << "Tinterpolation::get_value: Not initialized\n";
      exit(EXIT_FAILURE);
    }
  }
 
   /* 
    Calculate interpolated vector functions for a set of indices
    
    Input:
      t - time
      set = {i_1,..., i_len}      - index of the component
      j - order of the derivative
      
    Output: 
      v[n]:  v[i] = y_{i_k}^{(j)}(t)    k in [0, len-1],  j in [0,n-1]
  */
  
  void  get_value(const T & t, const int &j, std::vector<int> & set, T *v){
    
    if (init_status){
      
      T x  = t - ofs[0];
      
      int l = locate_knot(x);
      
      if (l < 0){
        std::cerr 
          << "Tinterpolation::get_value: index out of bounds\n"
          << "value=" << t << '\n';
        return;
      }
             
      if (j == 0){  // calculating only values
        
        for (int i = 0, len = set.size(); i < len; ++i)
          v[i] = ofs[set[i]+1] + spline_eval(l, x, k, W[set[i]], L);
        
      } else {
        
        // values of B-splines
        if (x_old != x) bspline_table(l, x, k, L, B);
        x_old = x;
                
        // orders of derivatives
        P[0] = j;
                
        // evaluate derivatives      
        for (int i = 0, len = set.size(); i < len; ++i)
          spline_eval(l, k, B, W[set[i]], L, 1, P, v+i);
      }
      
    } else {
      std::cerr << "Tinterpolation::get_value: Not initialized\n";
      exit(EXIT_FAILURE);
    }  
  }
  
  void  get_value(const T & t, const int &j, const int &len, int* set, T *v){
    
    if (init_status){
      
      T x  = t - ofs[0];
      
      int l = locate_knot(x);
      
      if (l < 0){
        std::cerr 
          << "Tinterpolation::get_value: index out of bounds\n"
          << "value=" << t << '\n';
        return;
      }
             
      if (j == 0){  // calculating only values
        
        for (int i = 0; i < len; ++i)
          v[i] = ofs[set[i]+1] + spline_eval(l, x, k, W[set[i]], L);
        
      } else {
        
        // values of B-splines
        if (x_old != x) bspline_table(l, x, k, L, B);
        x_old = x;
                
        // orders of derivatives
        P[0] = j;
                
        // evaluate derivatives      
        for (int i = 0; i < len; ++i)
          spline_eval(l, k, B, W[set[i]], L, 1, P, v+i);
      }
      
    } else {
      std::cerr << "Tinterpolation::get_value: Not initialized\n";
      exit(EXIT_FAILURE);
    }  
  }
  
  
  /* 
    Calculate interpolated vector functions for a set of indices
    
    Input:
      t - time
      set[len][2] - index and derivative order of the component 
                 = {(i_1, j_1), ...,  (i_q, j_q) }  
                  
                  i_k <= i_{k+1}
                  
                  i_k = i_{k+1} : j_{k} < j_{k+1}
                  
                  i in [0, d-1],  j in [0, k-1]
                
    Output: 
      v[k]:  v[k] = y_{i_k}^{(j_k)}(t)    k in [0, len-1]
  */

  void  get_value(const T & t, const int &jmax, const int &len, int **set, T *v){
    
    if (init_status){
      
      T x  = t - ofs[0];
      
      int l = locate_knot(x);
      
      if (l < 0){
        std::cerr 
          << "Tinterpolation::get_value: index out of bounds\n"
          << "value=" << t << '\n';
        return;
      }
      
      int i;
      
      if (jmax == 0){  // calculating only values
        
        for (i = 0; i < len; ++i) 
          v[i] = ofs[set[i][0]+1] + spline_eval(l, x, k, W[set[i][0]], L);
        
      } else {
        
        // values of B-splines
        if (x_old != x) bspline_table(l, x, k, L, B);
        x_old = x;
                
        // orders of derivatives
        i = 0;
        
        int j, q, nP;
        
        while (i < len){
          
          // collecting orders for q-th components
          nP = 1;
          q = set[i][0];
          P[0] = set[i][1];
          
          j = i + 1;
          while (j < len && q == set[j][0]) P[nP++] = set[j++][1];
                   
          // evaluate q-th component 
          spline_eval(l, k, B, W[q], L, nP, P, v+i);
          
          if (P[0] == 0) v[i] += ofs[q+1];
          
          i += nP;
        }
      }
      
    } else {
      std::cerr << "Tinterpolation::get_value: Not initialized\n";
      exit(EXIT_FAILURE);
    }  
  }
};


#endif // if !defined(__interpolation_h)

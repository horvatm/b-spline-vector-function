#if !defined(__interpolation_cubic_h)
#define __interpolation_cubic_h

/*
  A class for calculation of natural splines  y''(a) = y''(b) = 0 of vector 
  functions

     y(t) in T^d        (e.g. T = Real )
  
  and its derivatives. We follow
   
    J. Stoer, Numerische mathematik, 1 (Springer-Verlag, 2005)
  
  and
  
    Numerical recipes in C

  The data is given in the form
    
    (t_i, y(t))           i = 0, ..., N-1
  
  with t_{i+1} - t_i = dt_i > 0 or 
  
    y(t_i)                i = 0, ..., N-1
    t_i = tmin + i*dt
    
  The computation is explained in Stoer p. 108-114 and NR p.113.
  
  The error according to Stoer p.114 - 118 is bounded as
  
   |f(t)- S^{(k)} (t) < c_k L Delta^{(4-k)}  k = 0,1,2
    
  where 
  
    Delta = max_i |dt_i|
    
    L = max_{t in [a,b]} |f^{(4)}(t)| 
    
    c_0 = 5/384, c_1 = 1/24, c_2 = 3/8
    
  Author: Martin Horvat, April 2013, August 2016
*/ 


#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <cassert>

// allocation and freeing of C-style matrices
#include "matrix.h"

// enabling assert()
//#define NDEBUG

template <class T>
class Tinterpolation_cubic {
  
  public:
  
  // supported definitions of timestamps
  enum Ttime {
      Equid,          // Equidistant: t_i = t_0 + h*i
      Equid_General,  // Equidistant general: t_{i+1} - t_i = h.
      General         // General
    };

  protected:
  
  bool init_status;
  
  int d, lo, hi, N;
  
  Ttime time_type;
  
  T tmin,   // minimal time
    tmax,   // maximal time 
    dt,     // in case of equidistant data: t_{1} - t_{0}
    a, 
    b,
    
    *x,     // t_0 , t_1  , ... , t_{d-1}     x_i in T 
    *ofs,   // t0, y0_0   , ... , y0_{d-1}    t_0, y0_i in T 

    **y,    // y_0 , y_1   , ... , y_{d-1}    y_i in T^N
    **y2;   // y_0'', y_1'', ... , y_{d-1}''  y_i'' in T^N


  void prepare_spline(int d, int N, T **D){
    
    int	i, j;
  
    T	p, sig, *u = new T [N-1];
    
    x = new T [N];
    y =  matrix <T> (d, N);
    y2 = matrix <T> (d, N);
        
    for (i = 0; i < N; ++i) {
      x[i] = D[i][0];
      for(j = 0; j < d; ++j) y[j][i] = D[i][j+1];
    }
         
    switch (time_type){
      
      case General:
             
        for (i = 0; i < d; ++i) {
          
          y2[i][0] = u[0] = 0;
        
          for(j = 1; j < N-1; ++j){
            //sig = (x[i] - x[i-1])/(x[i+1] - x[i-1]);
            sig = (x[j] - x[j-1])/(x[j+1] - x[j-1]);
            
            //p = sig*y2[i-1] + 2.0;
            p = sig*y2[i][j-1] + 2;
            
            //y2[i] = (sig - 1.0)/p;
            y2[i][j] = (sig - 1)/p;
            
            //u[i] = (y[i+1] - y[i])/(x[i+1] - x[i]) 
            //       - (y[i] - y[i-1])/(x[i] - x[i-1]);
            u[j] = (y[i][j+1]-y[i][j])/(x[j+1]-x[j])
                   -(y[i][j]-y[i][j-1])/(x[j]-x[j-1]);
            
            //u[i] = (6.0*u[i]/(x[i+1] - x[i-1]) - sig*u[i-1])/p;
            u[j] = (6*u[j]/(x[j+1] - x[j-1]) - sig*u[j-1])/p;
          }
          
          y2[i][N-1] = 0; 
          
          //for(k = n-2; k >= 0; k--)y2[k] = y2[k]*y2[k+1] + u[k];
          for (j = N-2; j >= 0; --j) y2[i][j] = y2[i][j]*y2[i][j+1] + u[j];
        }
      break;
      
      default:
      
        for (i = 0; i < d; ++i) {
          
          y2[i][0] = u[0] = 0;
        
          for(j = 1; j < N-1; ++j){
            //sig = (x[i] - x[i-1])/(x[i+1] - x[i-1]) = 1/2
            
            //p = sig*y2[i-1] + 2.0; *= 2
            p = y2[i][j-1] + 4;
            
            //y2[i] = (sig - 1.0)/p;
            y2[i][j] = -1/p;
            
            //u[i] = (y[i+1] - y[i])/(x[i+1] - x[i]) 
            //        - (y[i] - y[i-1])/(x[i] - x[i-1]);
            u[j] = (y[i][j+1] - 2*y[i][j] + y[i][j-1])/dt;
            
            //u[i] = (6.0*u[i]/(x[i+1] - x[i-1]) - sig*u[i-1])/p;
            u[j] = (6*u[j]/dt - u[j-1])/p;
          }
          
          y2[i][N-1] = 0; 
          for (j = N-2; j >= 0; --j) y2[i][j] = y2[i][j]*y2[i][j+1] + u[j];
        }
    }
      
    delete [] u;    
  }

  // for time_type = General
  bool calc_interval(const T & t){
    
    if (t < tmin) return false;

    if (t > tmax) return false;

    if (t == tmin) {
      lo = 0;
      hi = 1;
      a = 1;
      b = 0;
      return true;
    } 
   
    if (t == tmax) {
      lo = N - 2;
      hi = N - 1;
      a = 0;
      b = 1;
      return true;
    } 
    
    int i;
        
    switch (time_type){
    
      case General:
      
        if (!(x[lo] <= t && t < x[hi])) {

          lo = 0;
          hi = N - 1;
          
          while (hi - lo > 1){
            i = (hi + lo) >> 1;
            if (x[i] > t) hi = i; else  lo = i;
          }
          
          dt = x[hi] - x[lo];
        }
          
        a = (x[hi] - t)/dt;
        b = (t - x[lo])/dt;
      
      break;
    
      default:
  
        b = (t - tmin)/dt;
        
        lo = int(b);
        hi = lo + 1;
        
        a = hi - b;
        b -= lo;
    }
        
    return true;
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
    
  public:
  
  Tinterpolation_cubic(): init_status(false) {}
  
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
      d -- dimension of vector function
      filename -- filename of the file with data
      time_type --  type of time/independent variable
      offset[d+1]  -- offset of all variables       
    
  */
  
  Tinterpolation_cubic(
    const int &d, 
    const char *filename, 
    const Tinterpolation_cubic<T>::Ttime & time_type, 
    T* offset = 0
  
  ) : d(d), lo(0), hi(1), time_type(time_type) {
     
     // read data
    std::ifstream f(filename);
    
    if (!f.is_open()) {
      std::cerr << "Error opening file " << filename << '\n';
      exit(EXIT_FAILURE);
    }
    
    set_offset(d, offset);    // setting ofs
    
    int i, j, d1 = d + 1;
     
    T t;
      
    std::vector <T> buf;

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
    N = buf.size();
    assert (N % d1 == 0);
    N /= d1;

    // format data matrix D
    T **D = new T* [N];
    D[0] = &buf[0];
    for (i = 1; i < N; ++i) D[i] = D[i-1] + d1;

    // apply offset to data
    for (i = 0; i < N; ++i)
      for (j = 0; j <= d; ++j)
        D[i][j] -= ofs[j];
    
    // time axis attributes
    tmin = D[0][0];
    tmax = D[N-1][0];
    dt = D[1][0] - tmin;
    
    // prepare spline data, setting y and y2   
    prepare_spline(d, N, D);  
    
    delete [] D;
    
    init_status = true;
  }
  
    
  /* 
    Interpolation of equidistant data
      
      X(t_i) = (X_0(t_i) , ... , X_{d-1}(t_i)) 
             = X[i] in T^d                        i = 0, ... , N-1
       
    with t_i = t_0 + i*dt.
    
    Input:
      d -- dimension of vector function
      nD -- number of data points
      X[nD][d] -- data
      t0 -- initial time
      dt -- time step
      offset[d+1] -- offset of all variables
  */
   
  Tinterpolation_cubic(
    const int &d, 
    const int &N, 
    T **X, 
    T t0, 
    T dt, 
    T* offset = 0
  ) : d(d), lo(0), hi(1), N(N), 
      time_type (Tinterpolation_cubic<T>::Equid), 
      tmin(t0), dt(dt) {
    
    // read data  
    int i, j;
    
    // settings ofs
    set_offset(d, offset);
    
    T **D = matrix<T> (N, d+1);
    
    for (i = 0; i < N; ++i){
      D[i][0] = tmin + i*dt;
      for (j = 0; j < d; ++j) D[i][j+1] = X[i][j];
    }
    
    // apply offset to data
    for (i = 0; i < N; ++i)
      for (j = 0; j <= d; ++j)
        D[i][j] -= ofs[j];
    
    // time axis attributes
    // tmin = D[0][0];
    tmax = D[N-1][0];
    //dt = D[1][0] - tmin;
    
    // prepare spline data, setting y and y2   
    prepare_spline(d, N, D);  
    
    free_matrix(D);
    
    init_status = true;   
  }
  
  /* 
    Interpolation of data stored in general format 
    
       Y[N][d+1]
    
    where    
    
       Y[i] = {t_i, X(t_i)} in T^(d+1)          i = 0, ... , N-1
    
    with t_{i} < t_{i+1}. 
    
    time_type = Equid_General:    Equidistant time dt = t_{i+1} - t_i
    
    time_type: General:   General time t_{i+1} > t_{i}, 
  */ 
  
  Tinterpolation_cubic(
    const int &d, 
    const int &N, 
    T **Y, 
    const Tinterpolation_cubic<T>::Ttime & time_type, 
    T* offset = 0
   ): d(d), lo(0), hi(1), N(N), time_type(time_type) {
    
    set_offset(d, offset); //  using d, set ofs    
    
    int i, j;
    
    T **D;

    // copy data and apply offset to data
    if (offset) {
      
      int i, j;
      
      D = matrix<T>(N, d + 1);
      
      for (i = 0; i < N; ++i)
        for (j = 0; j <= d; ++j)
          D[i][j] = Y[i][j] - ofs[j];
      
    } else D = Y;
    
    // time axis attributes
    tmin = D[0][0];
    tmax = D[N-1][0];
    dt = D[1][0] - tmin;
         
    // prepare spline data, setting y and y2
    prepare_spline(d, N, Y); 
    
    if (offset) free_matrix(D);  
    
    init_status = true;
  }
  
  ~Tinterpolation_cubic(){
    
    if (init_status) {
      
      delete [] x; 
      delete [] ofs;

      free_matrix(y);
      free_matrix(y2);
    }
  }
  
  
  /* 
    Returns dimension
  */
  int get_dim() const {return d;}
  
  /* 
    Returns the number of points 
  */
  int get_size() const {return N;}
  
  /* 
    Returns minimal time
  */
  T get_min_time() const { return tmin+ofs[0]; }
  
  /* 
    Returns maximal time
  */
  T get_max_time() const { return tmax+ofs[0]; }
  
  /* 
    Calculate interpolated component of a vector
    
    Input:
      t - time
      i - index of the component
      j - order of the derivative
      
    Return: 
      X_i^{(j)}(t)    i in [0, d-1],  j in [0,2]
  */
  
  T get_value(const T & t, const int & i, const int &j){
    
    if (init_status){
    
      if (!calc_interval(t-ofs[0])){   // sets lo, hi, a, b
        std::cerr 
          << "Tinterpolation_cubic::get_value: The index is out of bounds\n";
        return 0;
      }
      
      switch (j) { 
        case 0:
          // calculate value of a spline: y_i(t)
          // a*y[lo]+b*y[hi]+((a*a*a-a)*y2[lo]+(b*b*b-b)*y2[hi])*(h*h)/6
          return ofs[i+1] + a*y[i][lo] + b*y[i][hi] + 
                 (a*(a*a-1)*y2[i][lo] + b*(b*b-1)*y2[i][hi])*(dt*dt)/6;
        case 1:
          // calculate first derivative of a spline: y_i'(t)
          // -y2[lo]*a^2*dt/2+y2[hi]*b^2*dt/2 +(y[hi]-y[lo])/dt+(y2[lo]-y2[hi])*dt/6
          return (y[i][hi] - y[i][lo])/dt +
                 ((1-3*a*a)*y2[i][lo] + (3*b*b-1)*y2[i][hi])*dt/6;
        case 2:
          // calculate second derivative of a spline : y_i''(t)
          // y2[lo]*a + y2[hi]*b
          return y2[i][lo]*a + y2[i][hi]*b;
      }
    } else {
      std::cerr << "Tinterpolation not initialized\n";
      exit(EXIT_FAILURE);
    }
    return 0;
  }
  
  /* 
    Calculate interpolated vector functions
    
    Input:
      t - time
      i - index of the component
      j - order of the derivative
      
    Output: 
      v[d]:  v[i] = y_i^{(j)}(t)    i in [0, d-1],  j in [0,2]
  */
  
  void  get_value(const T & t, const int &j, T *v){
    
    if (init_status){
       
      if (!calc_interval(t-ofs[0])){  // sets lo, hi, a, b
        std::cerr 
          << "Tinterpolation_cubic::get_value: The index is out of bounds\n";
        return;
      }
      
      int i;
      
      T A, B, F;
      
      switch (j) { 
        case 0:
          // calculate value of a spline: y_i(t)
          // a*y[lo]+b*y[hi]+((a*a*a-a)*y2[lo]+(b*b*b-b)*y2[hi])*(h*h)/6
          F = dt*dt/6;
          A = a*(a*a-1)*F, 
          B = b*(b*b-1)*F;
            
          for (i = 0; i < d; ++i) 
            v[i] = ofs[i+1] + a*y[i][lo] + b*y[i][hi] + A*y2[i][lo] + B*y2[i][hi];
            
        break;
        
        case 1:
          // calculate first derivative of a spline: y_i'(t)
          // -y2[lo]*a^2*dt/2+y2[hi]*b^2*dt/2 +(y[hi]-y[lo])/dt+(y2[lo]-y2[hi])*dt/6
          F = dt/6;
          A = (3*a*a-1)*F, 
          B = (3*b*b-1)*F;
              
          for (i = 0; i < d; ++i) 
            v[i] = (y[i][hi] - y[i][lo])/dt - A*y2[i][lo] + B*y2[i][hi];
        break;
        
        case 2:
          // calculate second derivative of a spline : y_i''(t)
          // y2[lo]*a + y2[hi]*b
          for (i = 0; i < d; ++i) v[i] = y2[i][lo]*a + y2[i][hi]*b;
        break;
            
        default:
          for (i = 0; i < d; ++i) v[i] = 0;
        break;
      }
    } else {
      std::cerr << "Tinterpolation not initialized\n";
      exit(EXIT_FAILURE);
    }  
  }
 
 
   /* 
    Calculate interpolated vector functions for a set of indices
    
    Input:
      t - time
      set = {i_1, i_n}      - index of the component
      j - order of the derivative
      
    Output: 
      v[n]:  v[i] = y_{i_k}^{(j)}(t)    k in [0, n-1],  j in [0,2]
  */
  
  void  get_value(const T & t, const int &j, std::vector<int> & set, T *v){
    
    if (init_status){
       
      if (!calc_interval(t-ofs[0])){  // sets lo, hi, a, b
        std::cerr 
          << "Tinterpolation_cubic::get_value: The index is out of bounds\n";
        return;
      }
                  
      int i, ii, len = set.size();
      
      T A, B, F;
      
      switch (j) { 
        case 0:
          // calculate value of a spline: y_i(t)
          // a*y[lo]+b*y[hi]+((a*a*a-a)*y2[lo]+(b*b*b-b)*y2[hi])*(h*h)/6
          F = dt*dt/6;
          A = a*(a*a-1)*F;
          B = b*(b*b-1)*F;
          
          for (i = 0; i < len; ++i) {
            ii = set[i];
            v[i] = ofs[ii+1] + 
                  a*y[ii][lo] + b*y[ii][hi] + A*y2[ii][lo] + B*y2[ii][hi];
          }
        break;
        
        case 1:
          // calculate first derivative of a spline: y_i'(t)
          // -y2[lo]*a^2*dt/2+y2[hi]*b^2*dt/2 +(y[hi]-y[lo])/dt+(y2[lo]-y2[hi])*dt/6
          F = dt/6;
          A = (3*a*a-1)*F; 
          B = (3*b*b-1)*F;
          
          for (i = 0; i < len; ++i) {
            ii = set[i];
            v[i] = (y[ii][hi] - y[ii][lo])/dt - A*y2[ii][lo] + B*y2[ii][hi];
          }    
        break;
        
        case 2:
          // calculate second derivative of a spline : y_i''(t)
          // y2[lo]*a + y2[hi]*b 
          for (i = 0; i < len; ++i) {
            ii = set[i]; 
            v[i] = y2[ii][lo]*a + y2[ii][hi]*b;
          }
          
        break;
            
        default:
          for (i = 0; i < len; ++i) v[i] = 0;
        break;
      }
    } else {
      std::cerr << "Tinterpolation not initialized\n";
      exit(EXIT_FAILURE);
    }  
  }
 
};


#endif // if !defined(__interpolation_h)

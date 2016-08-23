#if !defined(__matrix_h)
#define __matrix_h

/*
  Library defining C-style matrix
  
  Ref: 
  * Numerical recipes in C

  Author: Martin Horvat, August 2016
*/

/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
template <class T> T** matrix(int nrow, int ncol) {

  T **m = new T* [nrow];
  
  m[0] = (T *) new char [nrow*ncol*sizeof(T)];
        
  for( int i = 1; i < nrow; i++) m[i] = m[i-1]+ ncol;

  return m;       
}


/* free a double matrix allocated by dmatrix() */
template <class T> void free_matrix(T **&m) {
  
  delete [] (char*) m[0];
  
  delete [] m;
  
  m = 0;
}

/* tensor 3. order                              */
template <class T> T*** tensor3 (const int &n1, const int &n2, const int &n3){
  
  int i, j;
  
  T ***m = new T** [n1], *p = new T[n1*n2*n3];
  
  for (i = 0; i < n1; ++i) {
  
    m[i] = new T* [n2];
  
    for (j = 0; j < n2; ++j){
      m[i][j] = p;
      p += n3;
    }
  }

  return m;
}

template <class T> void free_tensor3(T ***&m, const int &n1, const int &n2=0, const int &n3=0){
  
  delete [] m[0][0];
  
  for (int i = 0; i < n1; ++i) delete [] m[i];
  
  delete [] m;
  
  m = 0;
}

#endif //#if !defined(__matrix_h)

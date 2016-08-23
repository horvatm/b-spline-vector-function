#if !defined(__solver_banded_h)
#define __solver_banded_h

/*
  Library for solving linear systems with banded matrices via Gauss elimination 
  without pivoting. Code is based on de Boor's code (from the book) found in 
  
  http://people.sc.fsu.edu/~jburkardt/f_src/pppack/pppack.html
  
 
  Refs:
  * Carl DeBoor, A Practical Guide to Splines, Springer, 2001

  Author: Martin Horvat, May 2013
*/ 


#include <iostream>
#include <cmath>

/*****************************************************************************

   BANFAC factors a banded matrix without pivoting.

  Discussion:

    BANFAC returns in W the LU-factorization, without pivoting, of 
    the banded matrix A of order NROW with (NBANDL+1+NBANDU) bands 
    or diagonals in the work array W.
 
    Gauss elimination without pivoting is used.  The routine is 
    intended for use with matrices A which do not require row 
    interchanges during factorization, especially for the totally 
    positive matrices which occur in spline calculations.

    The matrix storage mode used is the same one used by LINPACK 
    and LAPACK, and results in efficient innermost loops.
 
    Explicitly, A has 
 
      NBANDL bands below the diagonal
      1     main diagonal
      NBANDU bands above the diagonal

    and thus, with MIDDLE=NBANDU+1,
    A(I+J,J) is in W(I+MIDDLE,J) for I=-NBANDU,...,NBANDL, J=1,...,NROW.

    For example, the interesting entries of a banded matrix
    matrix of order 9, with NBANDL=1, NBANDU=2:

      11 12 13  0  0  0  0  0  0
      21 22 23 24  0  0  0  0  0
       0 32 33 34 35  0  0  0  0
       0  0 43 44 45 46  0  0  0
       0  0  0 54 55 56 57  0  0
       0  0  0  0 65 66 67 68  0
       0  0  0  0  0 76 77 78 79
       0  0  0  0  0  0 87 88 89
       0  0  0  0  0  0  0 98 99

    would appear in the first 1+1+2=4 rows of W as follows:

       0  0 13 24 35 46 57 68 79
       0 12 23 34 45 56 67 78 89
      11 22 33 44 55 66 77 88 99
      21 32 43 54 65 76 87 98  0
 
    All other entries of W not identified in this way with an
    entry of A are never referenced.
 
    This routine makes it possible to solve any particular linear system 
    A*X=B for X by the call

      call banslv ( w, nroww, nrow, nbandl, nbandu, b )

    with the solution X contained in B on return.
 
    If IFLAG=2, then one of NROW-1, NBANDL, NBANDU failed to be nonnegative, 
    or else one of the potential pivots was found to be zero 
    indicating that A does not have an LU-factorization.  This 
    implies that A is singular in case it is totally positive.

  Modified:

    14 February 2007

  Author:

    Carl DeBoor

  Reference:

    Carl DeBoor,
    A Practical Guide to Splines,
    Springer, 2001,
    ISBN: 0387953663,
    LC: QA1.A647.v27.

  Parameters:
 
    Input/output, real ( kind = 8 ) W(NROWW,NROW).
    On input, W contains the "interesting" part of a banded 
    matrix A, with the diagonals or bands of A stored in the
    rows of W, while columns of A correspond to columns of W. 
    On output, W contains the LU-factorization of A into a unit 
    lower triangular matrix L and an upper triangular matrix U 
    (both banded) and stored in customary fashion over the 
    corresponding entries of A.  

    Input, integer ( kind = 4 ) NROWW, the row dimension of the work array W.
    NROWW must be at least NBANDL+1 + NBANDU.
 
    Input, integer ( kind = 4 ) NROW, the number of rows in A.

    Input, integer ( kind = 4 ) NBANDL, the number of bands of A below 
    the main diagonal.
 
    Input, integer ( kind = 4 ) NBANDU, the number of bands of A above 
    the main diagonal.
 
    Output, integer ( kind = 4 ) IFLAG, error flag.
    1, success.
    2, failure, the matrix was not factored.

******************************************************************************/
template <class T>
void banfac (
  T **w, 
  const int & nroww, 
  const int & nrow, 
  const int & nbandl, 
  const int & nbandu, 
  int & iflag){
  
  T factor;
  
  int i, j, k, l, n, middle, nrow_, pivot;
  
  nrow_ = nrow - 1;
  
  iflag = 1;

  if ( nrow < 1 ) {
    iflag = 2;
    return;
  }

//
//  W(MIDDLE,*) contains the main diagonal of A.
//

  middle = nbandu + 1;
  
  if ( nrow_ == 0 ) {
    
    if ( w[nbandu][nrow_] == T(0) ) iflag = 2;
    
    return;
  }

//
//  A is upper triangular.  Check that the diagonal is nonzero.
//

  if ( nbandl <= 0 ) {

    
    for (i = 0; i < nrow_; ++i){
      
      if ( w[nbandu][i] == T(0) ) {
        
        iflag = 2;
        
        return;
      }
      
    }

    if ( w[nbandu][nrow_] == T(0) ) iflag = 2;
  
    return;
//
//  A is lower triangular.  Check that the diagonal is nonzero and
//  divide each column by its diagonal.
//

  } else if ( nbandu <= 0 ) {

    for (i = 0; i < nrow_; ++i) {
      
      pivot = w[nbandu][i];

      if ( pivot == T(0) ) {
        
        iflag = 2;
        
        return;
      }
      
      k = std::min(nbandl, nrow_ - i);
    
      for (j = 0; j < k; ++j) w[middle + j][i] /= pivot;
    }

    return;
  }
  
//
//  A is not just a triangular matrix.  
//  Construct the LU factorization.
//
  for (i = 0; i < nrow_; ++i){
  
//
//  W(MIDDLE,I) is the pivot for the I-th step.
//
    if ( w[nbandu][i] == T(0) ) {

      iflag = 2;

      std::cerr 
        << "BANFAC - Fatal error!\n"
        << "Zero pivot encountered in column " <<  i << '\n';

      return;
    }
//
//  Divide each entry in column I below the diagonal by PIVOT.
//
    k = std::min( nbandl, nrow_ - i);
    
    for (j = 0; j < k; ++j) w[middle + j][i] /= w[nbandu][i];
    
//
//  Subtract A(I,I+K)*(I-th column) from (I+K)-th column (below row I).
//
    l = std::min( nbandu, nrow_ - i);
    
    for (k = 1; k <= l; ++k) {
    
      factor = w[nbandu - k][i + k];
      
      n = std::min( nbandl, nrow_ - i);
      
      for (j = 0; j < n; ++j)
        w[middle - k + j][i + k] -= w[middle+j][i]*factor;
    }
  }
//
//  Check the last diagonal entry.
//
  if ( w[nbandu][nrow_] == T(0) ) iflag = 2;
  
}

/*****************************************************************************80

 BANSLV solves a banded linear system A * X = B factored by BANFAC.

  Modified:

    14 February 2007

  Author:

    Carl DeBoor

  Reference:

    Carl DeBoor,
    A Practical Guide to Splines,
    Springer, 2001,
    ISBN: 0387953663,
    LC: QA1.A647.v27.

  Parameters:

    Input, real ( kind = 8 ) W(NROWW,NROW).  W contains the banded matrix,
    after it has been factored by BANFAC.

    Input, integer ( kind = 4 ) NROWW, the row dimension of the work array W.
    NROWW must be at least NBANDL+1 + NBANDU.
 
    Input, integer ( kind = 4 ) NROW, the number of rows in A.

    Input, integer ( kind = 4 ) NBANDL, the number of bands of A below the 
    main diagonal.
 
    Input, integer ( kind = 4 ) NBANDU, the number of bands of A above the 
    main diagonal.
 
    Input/output, real ( kind = 8 ) B(NROW).
    On input, B contains the right hand side of the system to be solved.
    On output, B contains the solution, X.

********************************************************************************/ 

template <class T>
void banslv (
  T **w, 
  const int &nroww, 
  const int &nrow, 
  const int &nbandl, 
  const int &nbandu, 
  T *b) {

  int i, j, k, jmax, nrow_;
  
  nrow_ = nrow - 1;
  
  if ( nrow == 1 ){
    b[0] /= w[nbandu][0];
    return;
  }

//
//  Forward pass:
//
//  For I = 1, 2, ..., NROW-1, subtract RHS(I)*(I-th column of L) 
//  from the right hand side, below the I-th row.
//

  if ( 0 < nbandl ) {
    for (i = 0; i < nrow_; ++i){
      
      jmax = std::min ( nbandl, nrow_ - i);
      
      for (j = 1; j <= jmax; ++j)
        b[i+j] -= b[i] * w[nbandu + j][i];
    }
  }

//
//  Backward pass:
//
//  For I=NROW, NROW-1,...,1, divide RHS(I) by 
//  the I-th diagonal entry of U, then subtract 
//  RHS(I)*(I-th column of U) from right hand side, above the I-th row.
//

  for (i = nrow_; i > 0; --i){
     
    b[i] /= w[nbandu][i];
    
    k = std::min (nbandu, i);
    
    for (j = 1; j <= k; ++j) b[i-j] -= b[i] * w[nbandu - j][i];
  }

  b[0] /= w[nbandu][0];
}

#endif // #if !defined(__solver_banded_h)

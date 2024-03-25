      subroutine strsm ( side, uplo, transa, diag, m, n, alpha, a, lda,
     $                   b, ldb )
c
cc STRSM solves A*x=alpha*B with optional transposes.
c
*     .. scalar arguments ..
      character*1        side, uplo, transa, diag
      integer            m, n, lda, ldb
      real               alpha
*     .. array arguments ..
      real               a( lda, * ), b( ldb, * )
*     ..
*
*  purpose
*  =======
*
*  strsm  solves one of the matrix equations
*
*     op( a )*x = alpha*b,   or   x*op( a ) = alpha*b,
*
*  where alpha is a scalar, x and b are m by n matrices, a is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( a )  is one  of
*
*     op( a ) = a   or   op( a ) = a'.
*
*  the matrix x is overwritten on b.
*
*  parameters
*  ==========
*
*  side   - character*1.
*           on entry, side specifies whether op( a ) appears on the left
*           or right of x as follows:
*
*              side = 'l' or 'l'   op( a )*x = alpha*b.
*
*              side = 'r' or 'r'   x*op( a ) = alpha*b.
*
*           unchanged on exit.
*
*  uplo   - character*1.
*           on entry, uplo specifies whether the matrix a is an upper or
*           lower triangular matrix as follows:
*
*              uplo = 'u' or 'u'   a is an upper triangular matrix.
*
*              uplo = 'l' or 'l'   a is a lower triangular matrix.
*
*           unchanged on exit.
*
*  transa - character*1.
*           on entry, transa specifies the form of op( a ) to be used in
*           the matrix multiplication as follows:
*
*              transa = 'n' or 'n'   op( a ) = a.
*
*              transa = 't' or 't'   op( a ) = a'.
*
*              transa = 'c' or 'c'   op( a ) = a'.
*
*           unchanged on exit.
*
*  diag   - character*1.
*           on entry, diag specifies whether or not a is unit triangular
*           as follows:
*
*              diag = 'u' or 'u'   a is assumed to be unit triangular.
*
*              diag = 'n' or 'n'   a is not assumed to be unit
*                                  triangular.
*
*           unchanged on exit.
*
*  m      - integer.
*           on entry, m specifies the number of rows of b. m must be at
*           least zero.
*           unchanged on exit.
*
*  n      - integer.
*           on entry, n specifies the number of columns of b.  n must be
*           at least zero.
*           unchanged on exit.
*
*  alpha  - real            .
*           on entry,  alpha specifies the scalar  alpha. when  alpha is
*           zero then  a is not referenced and  b need not be set before
*           entry.
*           unchanged on exit.
*
*  a      - real             array of dimension ( lda, k ), where k is m
*           when  side = 'l' or 'l'  and is  n  when  side = 'r' or 'r'.
*           before entry  with  uplo = 'u' or 'u',  the  leading  k by k
*           upper triangular part of the array  a must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           a is not referenced.
*           before entry  with  uplo = 'l' or 'l',  the  leading  k by k
*           lower triangular part of the array  a must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           a is not referenced.
*           note that when  diag = 'u' or 'u',  the diagonal elements of
*           a  are not referenced either,  but are assumed to be  unity.
*           unchanged on exit.
*
*  lda    - integer.
*           on entry, lda specifies the first dimension of a as declared
*           in the calling (sub) program.  when  side = 'l' or 'l'  then
*           lda  must be at least  max( 1, m ),  when  side = 'r' or 'r'
*           then lda must be at least max( 1, n ).
*           unchanged on exit.
*
*  b      - real             array of dimension ( ldb, n ).
*           before entry,  the leading  m by n part of the array  b must
*           contain  the  right-hand  side  matrix  b,  and  on exit  is
*           overwritten by the solution matrix  x.
*
*  ldb    - integer.
*           on entry, ldb specifies the first dimension of b as declared
*           in  the  calling  (sub)  program.   ldb  must  be  at  least
*           max( 1, m ).
*           unchanged on exit.
*
*
*  level 3 blas routine.
*
*
*  -- written on 8-february-1989.
*     jack dongarra, argonne national laboratory.
*     iain duff, aere harwell.
*     jeremy du croz, numerical algorithms group ltd.
*     sven hammarling, numerical algorithms group ltd.
*
*
*     .. external functions ..
      logical            lsame
      external           lsame
*     .. external subroutines ..
      external           xerbla
*     .. intrinsic functions ..
      intrinsic          max
*     .. local scalars ..
      logical            lside, nounit, upper
      integer            i, info, j, k, nrowa
      real               temp
*     .. parameters ..
      real               one         , zero
      parameter        ( one = 1.0e+0, zero = 0.0e+0 )
*     ..
*     .. executable statements ..
*
*     test the input parameters.
*
      lside  = lsame( side  , 'l' )
      if( lside )then
         nrowa = m
      else
         nrowa = n
      end if
      nounit = lsame( diag  , 'n' )
      upper  = lsame( uplo  , 'u' )
*
      info   = 0
      if(      ( .not.lside                ).and.
     $         ( .not.lsame( side  , 'r' ) )      )then
         info = 1
      else if( ( .not.upper                ).and.
     $         ( .not.lsame( uplo  , 'l' ) )      )then
         info = 2
      else if( ( .not.lsame( transa, 'n' ) ).and.
     $         ( .not.lsame( transa, 't' ) ).and.
     $         ( .not.lsame( transa, 'c' ) )      )then
         info = 3
      else if( ( .not.lsame( diag  , 'u' ) ).and.
     $         ( .not.lsame( diag  , 'n' ) )      )then
         info = 4
      else if( m  .lt.0               )then
         info = 5
      else if( n  .lt.0               )then
         info = 6
      else if( lda.lt.max( 1, nrowa ) )then
         info = 9
      else if( ldb.lt.max( 1, m     ) )then
         info = 11
      end if
      if( info.ne.0 )then
         call xerbla( 'strsm ', info )
         return
      end if
*
*     quick return if possible.
*
      if( n.eq.0 )
     $   return
*
*     and when  alpha.eq.zero.
*
      if( alpha.eq.zero )then
         do 20, j = 1, n
            do 10, i = 1, m
               b( i, j ) = zero
   10       continue
   20    continue
         return
      end if
*
*     start the operations.
*
      if( lside )then
         if( lsame( transa, 'n' ) )then
*
*           form  b := alpha*inv( a )*b.
*
            if( upper )then
               do 60, j = 1, n
                  if( alpha.ne.one )then
                     do 30, i = 1, m
                        b( i, j ) = alpha*b( i, j )
   30                continue
                  end if
                  do 50, k = m, 1, -1
                     if( b( k, j ).ne.zero )then
                        if( nounit )
     $                     b( k, j ) = b( k, j )/a( k, k )
                        do 40, i = 1, k - 1
                           b( i, j ) = b( i, j ) - b( k, j )*a( i, k )
   40                   continue
                     end if
   50             continue
   60          continue
            else
               do 100, j = 1, n
                  if( alpha.ne.one )then
                     do 70, i = 1, m
                        b( i, j ) = alpha*b( i, j )
   70                continue
                  end if
                  do 90 k = 1, m
                     if( b( k, j ).ne.zero )then
                        if( nounit )
     $                     b( k, j ) = b( k, j )/a( k, k )
                        do 80, i = k + 1, m
                           b( i, j ) = b( i, j ) - b( k, j )*a( i, k )
   80                   continue
                     end if
   90             continue
  100          continue
            end if
         else
*
*           form  b := alpha*inv( a' )*b.
*
            if( upper )then
               do 130, j = 1, n
                  do 120, i = 1, m
                     temp = alpha*b( i, j )
                     do 110, k = 1, i - 1
                        temp = temp - a( k, i )*b( k, j )
  110                continue
                     if( nounit )
     $                  temp = temp/a( i, i )
                     b( i, j ) = temp
  120             continue
  130          continue
            else
               do 160, j = 1, n
                  do 150, i = m, 1, -1
                     temp = alpha*b( i, j )
                     do 140, k = i + 1, m
                        temp = temp - a( k, i )*b( k, j )
  140                continue
                     if( nounit )
     $                  temp = temp/a( i, i )
                     b( i, j ) = temp
  150             continue
  160          continue
            end if
         end if
      else
         if( lsame( transa, 'n' ) )then
*
*           form  b := alpha*b*inv( a ).
*
            if( upper )then
               do 210, j = 1, n
                  if( alpha.ne.one )then
                     do 170, i = 1, m
                        b( i, j ) = alpha*b( i, j )
  170                continue
                  end if
                  do 190, k = 1, j - 1
                     if( a( k, j ).ne.zero )then
                        do 180, i = 1, m
                           b( i, j ) = b( i, j ) - a( k, j )*b( i, k )
  180                   continue
                     end if
  190             continue
                  if( nounit )then
                     temp = one/a( j, j )
                     do 200, i = 1, m
                        b( i, j ) = temp*b( i, j )
  200                continue
                  end if
  210          continue
            else
               do 260, j = n, 1, -1
                  if( alpha.ne.one )then
                     do 220, i = 1, m
                        b( i, j ) = alpha*b( i, j )
  220                continue
                  end if
                  do 240, k = j + 1, n
                     if( a( k, j ).ne.zero )then
                        do 230, i = 1, m
                           b( i, j ) = b( i, j ) - a( k, j )*b( i, k )
  230                   continue
                     end if
  240             continue
                  if( nounit )then
                     temp = one/a( j, j )
                     do 250, i = 1, m
                       b( i, j ) = temp*b( i, j )
  250                continue
                  end if
  260          continue
            end if
         else
*
*           form  b := alpha*b*inv( a' ).
*
            if( upper )then
               do 310, k = n, 1, -1
                  if( nounit )then
                     temp = one/a( k, k )
                     do 270, i = 1, m
                        b( i, k ) = temp*b( i, k )
  270                continue
                  end if
                  do 290, j = 1, k - 1
                     if( a( j, k ).ne.zero )then
                        temp = a( j, k )
                        do 280, i = 1, m
                           b( i, j ) = b( i, j ) - temp*b( i, k )
  280                   continue
                     end if
  290             continue
                  if( alpha.ne.one )then
                     do 300, i = 1, m
                        b( i, k ) = alpha*b( i, k )
  300                continue
                  end if
  310          continue
            else
               do 360, k = 1, n
                  if( nounit )then
                     temp = one/a( k, k )
                     do 320, i = 1, m
                        b( i, k ) = temp*b( i, k )
  320                continue
                  end if
                  do 340, j = k + 1, n
                     if( a( j, k ).ne.zero )then
                        temp = a( j, k )
                        do 330, i = 1, m
                           b( i, j ) = b( i, j ) - temp*b( i, k )
  330                   continue
                     end if
  340             continue
                  if( alpha.ne.one )then
                     do 350, i = 1, m
                        b( i, k ) = alpha*b( i, k )
  350                continue
                  end if
  360          continue
            end if
         end if
      end if
*
      return
*
*     end of strsm .
*
      end

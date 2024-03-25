      subroutine sgemv ( trans, m, n, alpha, a, lda, x, incx,
     $                   beta, y, incy )
c
cc SGEMV computes y = alpha*A*x+beta*y with optional transpose.
c
*     .. scalar arguments ..
      real               alpha, beta
      integer            incx, incy, lda, m, n
      character*1        trans
*     .. array arguments ..
      real               a( lda, * ), x( * ), y( * )
*     ..
*
*  purpose
*  =======
*
*  sgemv  performs one of the matrix-vector operations
*
*     y := alpha*a*x + beta*y,   or   y := alpha*a'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and a is an
*  m by n matrix.
*
*  parameters
*  ==========
*
*  trans  - character*1.
*           on entry, trans specifies the operation to be performed as
*           follows:
*
*              trans = 'n' or 'n'   y := alpha*a*x + beta*y.
*
*              trans = 't' or 't'   y := alpha*a'*x + beta*y.
*
*              trans = 'c' or 'c'   y := alpha*a'*x + beta*y.
*
*           unchanged on exit.
*
*  m      - integer.
*           on entry, m specifies the number of rows of the matrix a.
*           m must be at least zero.
*           unchanged on exit.
*
*  n      - integer.
*           on entry, n specifies the number of columns of the matrix a.
*           n must be at least zero.
*           unchanged on exit.
*
*  alpha  - real            .
*           on entry, alpha specifies the scalar alpha.
*           unchanged on exit.
*
*  a      - real             array of dimension ( lda, n ).
*           before entry, the leading m by n part of the array a must
*           contain the matrix of coefficients.
*           unchanged on exit.
*
*  lda    - integer.
*           on entry, lda specifies the first dimension of a as declared
*           in the calling (sub) program. lda must be at least
*           max( 1, m ).
*           unchanged on exit.
*
*  x      - real             array of dimension at least
*           ( 1 + ( n - 1 )*abs( incx ) ) when trans = 'n' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( incx ) ) otherwise.
*           before entry, the incremented array x must contain the
*           vector x.
*           unchanged on exit.
*
*  incx   - integer.
*           on entry, incx specifies the increment for the elements of
*           x. incx must not be zero.
*           unchanged on exit.
*
*  beta   - real            .
*           on entry, beta specifies the scalar beta. when beta is
*           supplied as zero then y need not be set on input.
*           unchanged on exit.
*
*  y      - real             array of dimension at least
*           ( 1 + ( m - 1 )*abs( incy ) ) when trans = 'n' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( incy ) ) otherwise.
*           before entry with beta non-zero, the incremented array y
*           must contain the vector y. on exit, y is overwritten by the
*           updated vector y.
*
*  incy   - integer.
*           on entry, incy specifies the increment for the elements of
*           y. incy must not be zero.
*           unchanged on exit.
*
*
*  level 2 blas routine.
*
*  -- written on 22-october-1986.
*     jack dongarra, argonne national lab.
*     jeremy du croz, nag central office.
*     sven hammarling, nag central office.
*     richard hanson, sandia national labs.
*
*
*     .. parameters ..
      real               one         , zero
      parameter        ( one = 1.0e+0, zero = 0.0e+0 )
*     .. local scalars ..
      real               temp
      integer            i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
*     .. external functions ..
      logical            lsame
      external           lsame
*     .. external subroutines ..
      external           xerbla
*     .. intrinsic functions ..
      intrinsic          max
*     ..
*     .. executable statements ..
*
*     test the input parameters.
*
      info = 0
      if     ( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 1
      else if( m.lt.0 )then
         info = 2
      else if( n.lt.0 )then
         info = 3
      else if( lda.lt.max( 1, m ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      else if( incy.eq.0 )then
         info = 11
      end if
      if( info.ne.0 )then
         call xerbla( 'sgemv ', info )
         return
      end if
*
*     quick return if possible.
*
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
*
*     set  lenx  and  leny, the lengths of the vectors x and y, and set
*     up the start points in  x  and  y.
*
      if( lsame( trans, 'n' ) )then
         lenx = n
         leny = m
      else
         lenx = m
         leny = n
      end if
      if( incx.gt.0 )then
         kx = 1
      else
         kx = 1 - ( lenx - 1 )*incx
      end if
      if( incy.gt.0 )then
         ky = 1
      else
         ky = 1 - ( leny - 1 )*incy
      end if
*
*     start the operations. in this version the elements of a are
*     accessed sequentially with one pass through a.
*
*     first form  y := beta*y.
*
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, leny
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, leny
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, leny
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, leny
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero )
     $   return
      if( lsame( trans, 'n' ) )then
*
*        form  y := alpha*a*x + y.
*
         jx = kx
         if( incy.eq.1 )then
            do 60, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  do 50, i = 1, m
                     y( i ) = y( i ) + temp*a( i, j )
   50             continue
               end if
               jx = jx + incx
   60       continue
         else
            do 80, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  iy   = ky
                  do 70, i = 1, m
                     y( iy ) = y( iy ) + temp*a( i, j )
                     iy      = iy      + incy
   70             continue
               end if
               jx = jx + incx
   80       continue
         end if
      else
*
*        form  y := alpha*a'*x + y.
*
         jy = ky
         if( incx.eq.1 )then
            do 100, j = 1, n
               temp = zero
               do 90, i = 1, m
                  temp = temp + a( i, j )*x( i )
   90          continue
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
  100       continue
         else
            do 120, j = 1, n
               temp = zero
               ix   = kx
               do 110, i = 1, m
                  temp = temp + a( i, j )*x( ix )
                  ix   = ix   + incx
  110          continue
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
  120       continue
         end if
      end if
*
      return
*
*     end of sgemv .
*
      end

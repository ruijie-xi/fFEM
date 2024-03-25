 
        subroutine ums2su (nlu, npiv, n, lup, lui, lux, x, w)
c
cc UMS2SU is a utility which solves an upper triangular system.
c
        integer nlu, npiv, n, lup (nlu), lui (*)
        real
     $          lux (*), x (n), w (n)
 
c=== ums2su ============================================================
c
c  unsymmetric-pattern multifrontal package (umfpack). version 2.0s.
c  copyright (c) 1995, timothy a. davis, university of florida, usa.
c  joint work with iain s. duff, rutherford appleton laboratory, uk.
c  september 1995. work supported by the national science foundation
c  (dms-9223088 and dms-9504974) and the state of florida; and by cray
c  research inc. through the allocation of supercomputing resources.
 
c***********************************************************************
c* notice:  "the umfpack package may be used solely for educational,   *
c* research, and benchmarking purposes by non-profit organizations and *
c* the u.s. government.  commericial and other organizations may make  *
c* use of umfpack solely for benchmarking purposes only.  umfpack may  *
c* be modified by or on behalf of the user for such use but at no time *
c* shall umfpack or any such modified version of umfpack become the    *
c* property of the user.  umfpack is provided without warranty of any  *
c* kind, either expressed or implied.  neither the authors nor their   *
c* employers shall be liable for any direct or consequential loss or   *
c* damage whatsoever arising out of the use or misuse of umfpack by    *
c* the user.  umfpack must not be sold.  you may make copies of        *
c* umfpack, but this notice and the copyright notice must appear in    *
c* all copies.  any other use of umfpack requires written permission.  *
c* your use of umfpack is an implicit agreement to these conditions."  *
c*                                                                     *
c* the ma38 package in the harwell subroutine library (hsl) has        *
c* equivalent functionality (and identical calling interface) as       *
c* umfpack.  it is available for commercial use.   technical reports,  *
c* information on hsl, and matrices are available via the world wide   *
c* web at http://www.cis.rl.ac.uk/struct/arcd/num.html, or by          *
c* anonymous ftp at seamus.cc.rl.ac.uk/pub.  also contact john         *
c* harding, harwell subroutine library, b 552, aea technology,         *
c* harwell, didcot, oxon ox11 0ra, england.                            *
c* telephone (44) 1235 434573, fax (44) 1235 434340,                   *
c* email john.harding@aeat.co.uk, who will provide details of price    *
c* and conditions of use.                                              *
c***********************************************************************
 
c=======================================================================
c  not user-callable.
 
c=======================================================================
c  description:
c=======================================================================
c
c  solves ux = b, where u is the upper triangular factor of a matrix
c  (if btf not used) or a single diagonal block (if btf is used).
c  b is overwritten with the solution x.
 
c=======================================================================
c  input:
c=======================================================================
c
c       nlu:            number of lu arrowheads in the lu factors
c       npiv:           number of pivots found (normally n)
c       n:              order of matrix
c       lup (1..nlu):   pointer to lu arrowheads in lui
c       lui ( ... ):    integer values of lu arrowheads
c       lux ( ... ):    real values of lu arroheads
c       x (1..n):       the right-hand-side
 
c=======================================================================
c  workspace:
c=======================================================================
c
c       w (1..n)
 
c=======================================================================
c  output:
c=======================================================================
c
c       x (1..n):       the solution to ux=b
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutine:   ums2s2
c       subroutines called:     strsv, sgemv
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer j, k, s, luip, luxp, luk, ludegr, ludegc, lurp, uxp,
     $          lucp, col
        real
     $          one
        parameter (one = 1.0)
 
c  s:       an element, or lu arrowhead
c  k:       kth pivot
c  j:       jth column in u2 array in element s
c  luip:    s is in lui (luip...)
c  luxp:    real part of s is in lux (luxp...)
c  luk:     number of pivots in s
c  ludegc:  column degree of non-pivotal part of s
c  ludegr:  row degree of non-pivotal part of s
c  lucp:    pattern of column of s in lui (lucp...lucp+ludegc-1)
c  lurp:    pattern of row of s in lui (lurp...lurp+ludegr-1)
c  uxp:     the luk-by-ludegr u2 block of s is in lux (uxp...)
c  col:     column index
 
c=======================================================================
c  executable statments:
c=======================================================================
 
        k = npiv
        do 30 s = nlu, 1, -1
 
c          -------------------------------------------------------------
c          get s-th lu arrowhead (s = nlu..1, in reverse pivotal order)
c          -------------------------------------------------------------
 
           luip   = lup (s)
           luxp   = lui (luip)
           luk    = lui (luip+1)
           ludegr = lui (luip+2)
           ludegc = lui (luip+3)
           lucp   = (luip + 7)
           lurp   = lucp + ludegc
           uxp    = luxp + luk * (ludegc + luk)
 
           if (luk .eq. 1) then
 
c             ----------------------------------------------------------
c             only one pivot, stride-1 sparse dot product
c             ----------------------------------------------------------
 
cfpp$ nodepchk l
              do j = 1, ludegr
                 col = lui (lurp+j-1)
c                row: k, u (row,col): lux (uxp+j-1)
                 x (k) = x (k) - lux (uxp+j-1) * x (col)
              end do
c             divide by pivot, u (k,k): lux (luxp)
              x (k) = x (k) / lux (luxp)
              k = k - 1
 
           else
 
c             ----------------------------------------------------------
c             more than one pivot
c             ----------------------------------------------------------
 
              k = k - luk
              do 20 j = 1, ludegr
                 col = lui (lurp+j-1)
                 w (j) = x (col)
20            continue
              call sgemv ('n', luk, ludegr, -one,
     $           lux (uxp), luk, w, 1, one, x (k+1), 1)
              call strsv ('u', 'n', 'n', luk,
     $           lux (luxp), ludegc + luk, x (k+1), 1)
 
           endif
 
30      continue
        return
        end

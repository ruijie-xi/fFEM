 
        subroutine ums2sl (nlu, npiv, n, lup, lui, lux, x, w)
c
cc UMS2SL is a utility which solves a lower triangular system.
c
        integer nlu, npiv, n, lup (nlu), lui (*)
        real
     $          lux (*), x (n), w (n)
 
c=== ums2sl ============================================================
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
c  solves lx = b, where l is the lower triangular factor of a matrix
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
c       x (1..n):       the solution to lx=b
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutine:   ums2s2
c       subroutines called:     strsv, sgemv
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer i, k, s, luip, luxp, luk, ludegc, lucp, lxp, row
        real
     $          one
        parameter (one = 1.0)
 
c  s:       an element, or lu arrowhead
c  k:       kth pivot
c  i:       ith row in l2 array in element s
c  luip:    integer part of s is in lui (luip...)
c  luxp:    real part of s is in lux (luxp...)
c  luk:     number of pivots in s
c  ludegc:  column degree of non-pivotal part of s
c  lucp:    pattern of column of s in lui (lucp...lucp+ludegc-1)
c  lxp:     the ludegc-by-luk l2 block of s is in lux (lxp...)
c  row:     row index
 
c=======================================================================
c  executable statments:
c=======================================================================
 
        k = 0
        do 40 s = 1, nlu
 
c          -------------------------------------------------------------
c          get the s-th lu arrowhead (s = 1..nlu, in pivotal order)
c          -------------------------------------------------------------
 
           luip   = lup (s)
           luxp   = lui (luip)
           luk    = lui (luip+1)
           ludegc = lui (luip+3)
           lucp   = (luip + 7)
           lxp    = luxp + luk
 
           if (luk .eq. 1) then
 
c             ----------------------------------------------------------
c             only one pivot, stride-1 sparse saxpy
c             ----------------------------------------------------------
 
              k = k + 1
c             l (k,k) is one
cfpp$ nodepchk l
              do i = 1, ludegc
                 row = lui (lucp+i-1)
c                col: k, l (row,col): lux (lxp+i-1)
                 x (row) = x (row) - lux (lxp+i-1) * x (k)
              end do
 
           else
 
c             ----------------------------------------------------------
c             more than one pivot
c             ----------------------------------------------------------
 
              call strsv ('l', 'n', 'u', luk,
     $           lux (luxp), ludegc + luk, x (k+1), 1)
              do 20 i = 1, ludegc
                 row = lui (lucp+i-1)
                 w (i) = x (row)
20            continue
              call sgemv ('n', ludegc, luk, -one,
     $           lux (lxp), ludegc + luk, x (k+1), 1, one, w, 1)
              do 30 i = 1, ludegc
                 row = lui (lucp+i-1)
                 x (row) = w (i)
30            continue
              k = k + luk
           endif
40      continue
        return
        end

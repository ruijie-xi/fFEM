! quicksort_module.f90
! This module contains two recursive subroutines for sorting arrays of real(8) and integer
! values using the quicksort algorithm. The subroutines are overloaded by the interface
! quicksort, which allows the user to call the quicksort subroutine with an array of any
! type that has a quicksort subroutine defined for it. The quicksort algorithm is a
! divide-and-conquer algorithm that sorts an array by partitioning it into two subarrays
! and recursively sorting the subarrays. The algorithm is not stable, meaning that the
! order of equal elements may not be preserved. The quicksort algorithm has an average
! time complexity of O(n log n) and a worst-case time complexity of O(n^2), where n is the
! number of elements in the array. The quicksort algorithm is an in-place algorithm, meaning
! that it does not require additional memory beyond the array being sorted. The quicksort
! algorithm is widely used in practice due to its efficiency and simplicity.
! TODO: use select type to handle different data types
module quicksort_module
   implicit none

   private
   public quicksort
   
   interface quicksort
      module procedure quicksort_real8, quicksort_integer
   end interface quicksort
   
 contains
 
 ! Subroutine to sort an array of real(8) values using the quicksort algorithm
 ! a : input/output array of real(8) values to be sorted
 ! index: input/output array of integers to store the original index of the elements
 ! make sure that the index array is initialized to 1, 2, 3, ..., n before calling this subroutine
 recursive subroutine quicksort_real8(a, index)
   implicit none
   real(8) :: a(:)
   integer :: index(:)
   integer :: t_index
   real(8) x, t
   integer :: first = 1, last
   integer i, j
 
   last = size(a,1)
   x = a( (first+last) / 2 )
   i = first
   j = last
   
   do
      do while (a(i) < x)
         i=i+1
      end do
      do while (x < a(j))
         j=j-1
      end do
      if (i >= j) exit
      t = a(i);  a(i) = a(j);  a(j) = t
      t_index = index(i);  index(i) = index(j);  index(j) = t_index
      i=i+1
      j=j-1
   end do
   
   if (first < i - 1) call quicksort_real8(a(first : i - 1),index(first : i - 1))
   if (j + 1 < last)  call quicksort_real8(a(j + 1 : last),index(j + 1 : last))
 end subroutine quicksort_real8

 recursive subroutine quicksort_integer(a,index)
 implicit none
 integer :: a(:)
 integer :: index(:)
 integer :: t_index
 integer x, t
 integer :: first = 1, last
 integer i, j

 last = size(a,1)
 x = a( (first+last) / 2 )
 i = first
 j = last
 
 do
    do while (a(i) < x)
       i=i+1
    end do
    do while (x < a(j))
       j=j-1
    end do
    if (i >= j) exit
    t = a(i);  a(i) = a(j);  a(j) = t
    t_index = index(i);  index(i) = index(j);  index(j) = t_index
    i=i+1
    j=j-1
 end do
 
 if (first < i - 1) call quicksort_integer(a(first : i - 1),index(first : i - 1))
 if (j + 1 < last)  call quicksort_integer(a(j + 1 : last),index(j + 1 : last))
end subroutine quicksort_integer
   
 end module quicksort_module
 
 
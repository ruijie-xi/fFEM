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
module ffem_quicksort
   implicit none

   private
   public quicksort, unique, binarysearch
   
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

function unique(array) result(unique_array)
   integer, intent(in) :: array(:)
   integer, allocatable,dimension(:) :: unique_array

   integer :: i,count,left,right,N
   integer, allocatable, dimension(:) :: temp_array,index,temp_unique_array

   N = size(array)
   allocate(temp_array(N),index(N),temp_unique_array(N))
   temp_array = array
   index = [(i,i=1,N)]

   call quicksort_integer(temp_array, index)

   count = 1
   left = 1
   right = 1
   temp_unique_array(1) = temp_array(1)
   do while (right <= N)
      if (temp_array(right) /= temp_array(left)) then
         count = count + 1
         left = right
         temp_unique_array(count) = temp_array(right)
      end if
      right = right + 1
   end do

   allocate(unique_array(count))

   unique_array = temp_unique_array(1:count)

end function unique

function binarysearch(sorted_array, value) result(index)
   integer, intent(in) :: sorted_array(:)
   integer, intent(in) :: value
   integer :: index
   integer :: left, right, mid
   integer :: N

   N = size(sorted_array)
   left = 1
   right = N
   index = -1
   
   if (value < sorted_array(1) .or. value > sorted_array(N)) then 
      index = -1
      return
   end if

   do while (left <= right)
      mid = (left + right) / 2
      if (sorted_array(mid) == value) then
         index = mid
         return
      else if (sorted_array(mid) < value) then
         left = mid + 1
      else
         right = mid - 1
      end if
   end do

end function binarysearch
   
end module ffem_quicksort
 
 
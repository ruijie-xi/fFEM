program test_quicksort
  use quicksort_module
  implicit none
  real(8) :: a(1000),btmp,a_restore(1000)
  integer :: b(1000),b_restore(1000)
  integer :: i
  integer :: index(1000)

  index = [(i,i=1,1000)]
  call random_seed()
  do i = 1, size(a)
    call random_number(a(i))
  end do
  print *, 'Unsorted real8 array:'
  do i = 1, size(a)
    print *, a(i)
  end do
  call quicksort(a,index)
  print *, 'Sorted real8 array:'
  do i = 1, size(a)
    print *, a(i)
  end do
  print *, 'Index array:'
  do i = 1, size(index)
    print *, index(i)
  end do

  ! restore array
  do i = 1, size(a)
    a_restore(index(i)) = a(i)
  end do
  print *, "restored array"
  do i = 1, size(a)
    print *, a_restore(i)
  end do

  call random_seed()
  do i = 1, size(b)
    call random_number(btmp)
    b(i) = int(btmp*10000)
  end do
  print *, 'Unsorted integer array:'
  do i = 1, size(b)
    print *, b(i)
  end do
  index = [(i,i=1,1000)]
  call quicksort(b,index)
  print *, 'Sorted integer array:'
  do i = 1, size(b)
    print *, b(i)
  end do
  print *, 'Index array:'
  do i = 1, size(index)
    print *, index(i)
  end do

  ! restore array
  do i = 1, size(b)
    b_restore(index(i)) = b(i)
  end do
  print *, "restored array"
  do i = 1, size(b)
    print *, b_restore(i)
  end do



end program test_quicksort


module sigmoid_module

  use utils_module
  
  implicit none

  interface map_sigmoid
     module procedure map_sigmoid_rank1
     module procedure map_sigmoid_rank2
     module procedure map_sigmoid_rank3
     module procedure map_sigmoid_rank4
  end interface map_sigmoid

  interface unmap_sigmoid
     module procedure unmap_sigmoid_rank1
     module procedure unmap_sigmoid_rank2
     module procedure unmap_sigmoid_rank3
     module procedure unmap_sigmoid_rank4
  end interface unmap_sigmoid

contains

  function map_sigmoid_rank1(x, L) result (z)

    real*8 :: L, mu, lambda
    real*8, dimension(:) :: x
    real*8, dimension(:), allocatable :: z

    allocate(z(size(x)))
    
    mu = 0.d5 * L
    lambda = 0.d25 * L

    z = min(x, mu + lambda * log(x / (L - x)))

  end function map_sigmoid_rank1

  function map_sigmoid_rank2(x, L) result (z)

    real*8 :: L, mu, lambda
    real*8, dimension(:,:) :: x
    real*8, dimension(:,:), allocatable :: z
    integer*4 :: i

    allocate(z(size(x,1), size(x,2)))

    do i = 1,size(x,1)
       z(i,:) = map_sigmoid_rank1(x(i,:), L)
    end do

  end function map_sigmoid_rank2

  function map_sigmoid_rank3(x, L) result (z)

    real*8 :: L, mu, lambda
    real*8, dimension(:,:,:) :: x
    real*8, dimension(:,:,:), allocatable :: z
    integer*4 :: i, j

    allocate(z(size(x,1), size(x,2), size(x,3)))

    do j = 1,size(x,2)
       do i = 1,size(x,1)
          z(i,j,:) = map_sigmoid_rank1(x(i,j,:), L)
       end do
    end do

  end function map_sigmoid_rank3
  
  function map_sigmoid_rank4(x, L) result (z)

    real*8 :: L, mu, lambda
    real*8, dimension(:,:,:,:) :: x
    real*8, dimension(:,:,:,:), allocatable :: z
    integer*4 :: i, j, k

    allocate(z(size(x,1), size(x,2), size(x,3), size(x,4)))

    do k = 1,size(x,3)
       do j = 1,size(x,2)
          do i = 1,size(x,1)
             z(i,j,k,:) = map_sigmoid_rank1(x(i,j,k,:), L)
          end do
       end do
    end do

  end function map_sigmoid_rank4

  function unmap_sigmoid_rank1(z, L) result (x)

    real*8 :: L, mu, lambda
    real*8, dimension(:) :: z
    real*8, dimension(:), allocatable :: x

    allocate(x(size(z)))

    mu = 0.d5 * L
    lambda = 0.d25 * L

    x = dist_normal(z, mu, lambda)
    x = max(z, L * (x / (1 + x)))

  end function unmap_sigmoid_rank1

  function unmap_sigmoid_rank2(x, L) result (z)

    real*8 :: L, mu, lambda
    real*8, dimension(:,:) :: x
    real*8, dimension(:,:), allocatable :: z
    integer*4 :: i

    allocate(z(size(x,1), size(x,2)))

    do i = 1,size(x,1)
       z(i,:) = unmap_sigmoid_rank1(x(i,:), L)
    end do

  end function unmap_sigmoid_rank2

  function unmap_sigmoid_rank3(x, L) result (z)

    real*8 :: L, mu, lambda
    real*8, dimension(:,:,:) :: x
    real*8, dimension(:,:,:), allocatable :: z
    integer*4 :: i, j

    allocate(z(size(x,1), size(x,2), size(x,3)))

    do j = 1,size(x,2)
       do i = 1,size(x,1)
          z(i,j,:) = unmap_sigmoid_rank1(x(i,j,:), L)
       end do
    end do

  end function unmap_sigmoid_rank3
  
  function unmap_sigmoid_rank4(x, L) result (z)

    real*8 :: L, mu, lambda
    real*8, dimension(:,:,:,:) :: x
    real*8, dimension(:,:,:,:), allocatable :: z
    integer*4 :: i, j, k

    allocate(z(size(x,1), size(x,2), size(x,3), size(x,4)))

    do k = 1,size(x,3)
       do j = 1,size(x,2)
          do i = 1,size(x,1)
             z(i,j,k,:) = unmap_sigmoid_rank1(x(i,j,k,:), L)
          end do
       end do
    end do

  end function unmap_sigmoid_rank4

end module sigmoid_module

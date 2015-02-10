module output_module

  use constants_module
  use input_module
  
  implicit none

  type variables
     character(slen) :: name = ''
     integer*4 :: rank = -1
     real*8, pointer :: rank0 => null()
     real*8, dimension(:), pointer :: rank1 => null()
     real*8, dimension(:,:), pointer :: rank2 => null()
     real*8, dimension(:,:,:), pointer :: rank3 => null()
     integer*4 :: fid = -1
     integer*4 :: n = 0
     type(variables), pointer :: sum
  end type variables
  
  interface get_pointer
     module procedure get_pointer_rank0
     module procedure get_pointer_rank1
     module procedure get_pointer_rank2
     module procedure get_pointer_rank3
  end interface get_pointer
  
contains

  subroutine alloc_variable(var, name, dims)

    type(variables), dimension(:), allocatable, intent(inout) :: var
    type(variables), dimension(:), allocatable :: tmp
    character(*), intent(in) :: name
    integer*4, dimension(:), intent(in), optional :: dims
    integer*4 :: i

    if (allocated(var)) then
       allocate(tmp(size(var)))
       tmp = var
       deallocate(var)
       allocate(var(size(tmp) + 1))
       var(1:size(tmp)) = tmp
       deallocate(tmp)
    else
       allocate(var(1))
    end if

    i = size(var)

    var(i)%name = trim(name)

    if (present(dims)) then
       var(i)%rank = size(dims)

       select case (var(i)%rank)
       case (1)
          allocate(var(i)%rank1(dims(1)))
          var(i)%rank1 = 0.d0
       case (2)
          allocate(var(i)%rank2(dims(1), dims(2)))
          var(i)%rank2 = 0.d0
       case (3)
          allocate(var(i)%rank3(dims(1), dims(2), dims(3)))
          var(i)%rank3 = 0.d0
       end select
    else
       var(i)%rank = 0
       
       allocate(var(i)%rank0)
       var(i)%rank0 = 0.d0
    end if
    
  end subroutine alloc_variable

  subroutine alloc_variable_sum(var, name)

    type(variables), dimension(:), allocatable, intent(inout) :: var
    character(*), intent(in) :: name
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          select case (var(i)%rank)
          case (0)
             allocate(var(i)%sum%rank0)
             var(i)%sum%rank0 = 0.d0
          case (1)
             allocate(var(i)%sum%rank1(size(var(i)%rank1)))
             var(i)%sum%rank1 = 0.d0
          case (2)
             allocate(var(i)%sum%rank2( &
                  size(var(i)%rank2, 1), &
                  size(var(i)%rank2, 2)))
             var(i)%sum%rank2 = 0.d0
          case (3)
             allocate(var(i)%sum%rank3( &
                  size(var(i)%rank3, 1), &
                  size(var(i)%rank3, 2), &
                  size(var(i)%rank3, 3)))
             var(i)%sum%rank3 = 0.d0
          end select
       end if
    end do
    
  end subroutine alloc_variable_sum
  
  subroutine get_pointer_rank0(var, name, ptr)

    type(variables), dimension(:), allocatable, intent(inout) :: var
    character(*), intent(in) :: name
    real*8, pointer, intent(inout) :: ptr
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) .eq. trim(name)) then
          ptr => var(i)%rank0
          var(i)%rank = 0
          exit
       end if
    end do
    
  end subroutine get_pointer_rank0

  subroutine get_pointer_rank1(var, name, ptr)

    type(variables), dimension(:), allocatable, intent(inout) :: var
    character(*), intent(in) :: name
    real*8, dimension(:), pointer, intent(inout) :: ptr
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) .eq. trim(name)) then
          ptr => var(i)%rank1
          var(i)%rank = 1
          exit
       end if
    end do
    
  end subroutine get_pointer_rank1

  subroutine get_pointer_rank2(var, name, ptr)

    type(variables), dimension(:), allocatable, intent(inout) :: var
    character(*), intent(in) :: name
    real*8, dimension(:,:), pointer, intent(inout) :: ptr
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) .eq. trim(name)) then
          ptr => var(i)%rank2
          var(i)%rank = 2
          exit
       end if
    end do
    
  end subroutine get_pointer_rank2

  subroutine get_pointer_rank3(var, name, ptr)

    type(variables), dimension(:), allocatable, intent(inout) :: var
    character(*), intent(in) :: name
    real*8, dimension(:,:,:), pointer, intent(inout) :: ptr
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) .eq. trim(name)) then
          ptr => var(i)%rank3
          var(i)%rank = 3
          exit
       end if
    end do
    
  end subroutine get_pointer_rank3
  
  subroutine output_init(var, vars)

    type(variables), dimension(:), intent(inout) :: var
    character(*), dimension(:), intent(in) :: vars
    integer*4 :: i, j

    do i = 1,size(vars)
       do j = 1,size(var)
          if (trim(var(j)%name) == trim(vars(i))) then
             var(j)%fid = 100 + i
       
             open(unit=var(j)%fid, file=trim(vars(i))//".out", &
                  action="write", status="replace", form="unformatted")
          end if
       end do
    end do

  end subroutine output_init

  subroutine output_init_sum(var, var_sum)

    type(variables), dimension(:), intent(in) :: var
    type(variables), dimension(:), allocatable, intent(inout) :: var_sum
    integer*4 :: i

    allocate(var_sum(size(var)))
    
    do i = 1,size(var)
       var_sum(i)%name = var(i)%name
       var_sum(i)%rank = var(i)%rank
       if (var(i)%fid .gt. 0) then
          var_sum(i)%fid = var(i)%fid + 100
          select case (var(i)%rank)
          case (0)
             allocate(var_sum(i)%rank0)
             var_sum(i)%rank0 = 0.d0
          case (1)
             allocate(var_sum(i)%rank1( &
                  size(var(i)%rank1)))
             var_sum(i)%rank1 = 0.d0
          case (2)
             allocate(var_sum(i)%rank2( &
                  size(var(i)%rank2, 1), &
                  size(var(i)%rank2, 2)))
             var_sum(i)%rank2 = 0.d0
          case (3)
             allocate(var_sum(i)%rank3( &
                  size(var(i)%rank3, 1), &
                  size(var(i)%rank3, 2), &
                  size(var(i)%rank3, 3)))
             var_sum(i)%rank3 = 0.d0
          end select
          
          open(unit=var_sum(i)%fid, file=trim(var(i)%name)//".sum", &
               action="write", status="replace", form="unformatted")
       end if
    end do

  end subroutine output_init_sum

  subroutine output_update(var, var_sum)

    type(variables), dimension(:), intent(in) :: var
    type(variables), dimension(:), intent(inout) :: var_sum
    integer*4 :: i

    do i = 1,size(var)
       if (var(i)%fid > 0) then
          var(i)%n = var(i)%n + 1
          select case (var(i)%rank)
          case (0)
             var_sum(i)%rank0 = var_sum(i)%rank0 + var(i)%rank0
          case (1)
             var_sum(i)%rank1 = var_sum(i)%rank1 + var(i)%rank1
          case (2)
             var_sum(i)%rank2 = var_sum(i)%rank2 + var(i)%rank2
          case (3)
             var_sum(i)%rank3 = var_sum(i)%rank3 + var(i)%rank3
          end select
       end if
    end do
    
  end subroutine output_update

  subroutine output_write(var)

    type(variables), dimension(:), intent(in) :: var
    integer*4 :: i

    do i = 1,size(var)
       if (var(i)%fid > 0) then
          select case (var(i)%rank)
          case (0)
             write(var(i)%fid) var(i)%rank0
          case (1)
             write(var(i)%fid) var(i)%rank1
          case (2)
             write(var(i)%fid) var(i)%rank2
          case (3)
             write(var(i)%fid) var(i)%rank3
          end select
       end if
    end do
    
  end subroutine output_write

  subroutine output_clear(var)

    type(variables), dimension(:), intent(inout) :: var
    integer*4 :: i
    
    do i = 1,size(var)
       if (var(i)%fid > 0) then
          select case (var(i)%rank)
          case (0)
             var(i)%rank0 = 0.d0
          case (1)
             var(i)%rank1 = 0.d0
          case (2)
             var(i)%rank2 = 0.d0
          case (3)
             var(i)%rank3 = 0.d0
          end select
       end if
    end do

  end subroutine output_clear

  subroutine output_close(var)

    type(variables), dimension(:), intent(in) :: var
    integer*4 :: i
    
    do i = 1,size(var)
       if (var(i)%fid .gt. 0) then
          close(var(i)%fid)
       end if
    end do

  end subroutine output_close
  
  subroutine write_dimensions(par)

    integer*4, parameter :: fid=1
    type(parameters), intent(in) :: par
  
    open(unit=fid, file="dims.out", action="write", status="replace", form="unformatted")
    write(fid) par%nx
    write(fid) par%dx
    write(fid) par%nt
    write(fid) par%dt
    write(fid) par%nfractions
    write(fid) par%nlayers+2
    write(fid) par%ntout
    write(fid) par%tout
    close(fid)

  end subroutine write_dimensions

  function is_output(var, name) result (ret)

    type(variables), dimension(:) :: var
    character(*) :: name
    logical :: ret
    integer*4 :: i
    
    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          ret = (var(i)%fid > 0)
          exit
       end if
    end do

  end function is_output
  
end module output_module

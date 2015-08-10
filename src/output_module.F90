module output_module

  use constants_module
  use input_module
  use utils_module
  use bed_module
  
  implicit none

  type variables
     character(slen) :: name = ''
     integer*4 :: rank = -1
     integer*4 :: n = 0
     integer*4, dimension(:), pointer :: dims => null()
     type(variables_data), pointer :: val => null()
     type(variables_data), pointer :: sum => null()
     type(variables_data), pointer :: avg => null()
     type(variables_data), pointer :: var => null()
     type(variables_data), pointer :: min => null()
     type(variables_data), pointer :: max => null()
     logical :: isset = .false.
  end type variables

  type variables_data
     real*8, dimension(:), pointer :: rank1 => null()
     real*8, dimension(:,:), pointer :: rank2 => null()
     real*8, dimension(:,:,:), pointer :: rank3 => null()
     real*8, dimension(:,:,:,:), pointer :: rank4 => null()
     real*8 :: init
     integer*4 :: fid = -1
  end type variables_data

  interface get_pointer
     module procedure get_pointer_rank1
     module procedure get_pointer_rank2
     module procedure get_pointer_rank3
     module procedure get_pointer_rank4
  end interface get_pointer

contains

  subroutine alloc_variable(var, name, dims)

    type(variables), dimension(:), allocatable, intent(inout) :: var
    type(variables), dimension(:), allocatable :: tmp
    character(*), intent(in) :: name
    integer*4, dimension(:), intent(in) :: dims
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

    if (sum(dims) .eq. 0) then
       var(i)%rank = 0
    else
       var(i)%rank = size(dims)
       allocate(var(i)%dims(var(i)%rank)) ! MEMERR
       var(i)%dims = dims
    end if

    call alloc_variable_data(var(i)%val, dims, var(i)%rank) ! MEMERR
    call alloc_variable_data(var(i)%sum, dims, var(i)%rank)
    call alloc_variable_data(var(i)%var, dims, var(i)%rank)
    call alloc_variable_data(var(i)%min, dims, var(i)%rank,  1d10)
    call alloc_variable_data(var(i)%max, dims, var(i)%rank, -1d10)

    allocate(var(i)%avg) ! MEMERR

  end subroutine alloc_variable

  subroutine alloc_variable_data(var, dims, rank, val)

    type(variables_data), pointer, intent(inout) :: var
    integer*4, dimension(:), intent(in) :: dims
    integer*4, intent(in) :: rank
    real*8, optional, intent(in) :: val

    allocate(var) ! MEMERR

    if (present(val)) then
       var%init = val
    else
       var%init = 0.d0
    end if

    allocate(var%rank1(product(dims)))
    var%rank1 = var%init

    ! create remapping pointer
    select case (rank)
    case (2)
       var%rank2(1:dims(1), 1:dims(2)) => var%rank1
    case (3)
       var%rank3(1:dims(1), 1:dims(2), 1:dims(3)) => var%rank1
    case (4)
       var%rank4(1:dims(1), 1:dims(2), 1:dims(3), 1:dims(4)) => var%rank1
    end select

  end subroutine alloc_variable_data

  subroutine dealloc_variables(var)

    type(variables), dimension(:), intent(inout) :: var
    integer*4 :: i

    do i = 1,size(var)
       if (associated(var(i)%dims)) then
          deallocate(var(i)%dims)
       end if
       
       call dealloc_variable_data(var(i)%val)
       call dealloc_variable_data(var(i)%sum)
       call dealloc_variable_data(var(i)%avg)
       call dealloc_variable_data(var(i)%var)
       call dealloc_variable_data(var(i)%min)
       call dealloc_variable_data(var(i)%max)
       
    end do

  end subroutine dealloc_variables

  subroutine dealloc_variable_data(var)

    type(variables_data), pointer, intent(inout) :: var

    nullify(var%rank4)
    nullify(var%rank3)
    nullify(var%rank2)

    if (associated(var%rank1)) then
       deallocate(var%rank1)
    end if

    if (associated(var)) then
       deallocate(var)
    end if

  end subroutine dealloc_variable_data

  subroutine get_pointer_rank1(var, name, dims, ptr)

    type(variables), dimension(:), allocatable, intent(inout) :: var
    character(*), intent(in) :: name
    integer*4, dimension(:), intent(in) :: dims
    real*8, dimension(:), pointer, intent(inout) :: ptr
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          ptr(1:dims(1)) => var(i)%val%rank1
          exit
       end if
    end do
    
  end subroutine get_pointer_rank1

  subroutine get_pointer_rank2(var, name, dims, ptr)

    type(variables), dimension(:), allocatable, intent(inout) :: var
    character(*), intent(in) :: name
    integer*4, dimension(:), intent(in) :: dims
    real*8, dimension(:,:), pointer, intent(inout) :: ptr
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          ptr(1:dims(1), 1:dims(2)) => var(i)%val%rank1
          exit
       end if
    end do
    
  end subroutine get_pointer_rank2

  subroutine get_pointer_rank3(var, name, dims, ptr)

    type(variables), dimension(:), allocatable, intent(inout) :: var
    character(*), intent(in) :: name
    integer*4, dimension(:), intent(in) :: dims
    real*8, dimension(:,:,:), pointer, intent(inout) :: ptr
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          ptr(1:dims(1), 1:dims(2), 1:dims(3)) => var(i)%val%rank1
          exit
       end if
    end do
    
  end subroutine get_pointer_rank3
  
  subroutine get_pointer_rank4(var, name, dims, ptr)

    type(variables), dimension(:), allocatable, intent(inout) :: var
    character(*), intent(in) :: name
    integer*4, dimension(:), intent(in) :: dims
    real*8, dimension(:,:,:,:), pointer, intent(inout) :: ptr
    integer*4 :: i
    
    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          ptr(1:dims(1), 1:dims(2), 1:dims(3), 1:dims(4)) => var(i)%val%rank1
          exit
       end if
    end do

  end subroutine get_pointer_rank4

  subroutine output_init(par, var)

    type(parameters), intent(in) :: par
    type(variables), dimension(:), intent(inout) :: var
    integer*4 :: i, j

    do i = 1,size(par%outputvars)
       do j = 1,size(var)
          if (trim(var(j)%name) == trim(par%outputvars(i))) then
             call output_init_data(var(j)%val, par%output_dir, trim(par%outputvars(i))//".out", 100+i)
             if (isin(par%outputtypes, 'sum')) then
                call output_init_data(var(j)%sum, par%output_dir, trim(par%outputvars(i))//".sum.out", 200+i)
             end if
             if (isin(par%outputtypes, 'avg')) then
                call output_init_data(var(j)%avg, par%output_dir, trim(par%outputvars(i))//".avg.out", 300+i)
             end if
             if (isin(par%outputtypes, 'var')) then
                call output_init_data(var(j)%var, par%output_dir, trim(par%outputvars(i))//".var.out", 400+i)
             end if
             if (isin(par%outputtypes, 'min')) then
                call output_init_data(var(j)%min, par%output_dir, trim(par%outputvars(i))//".min.out", 500+i)
             end if
             if (isin(par%outputtypes, 'max')) then
                call output_init_data(var(j)%max, par%output_dir, trim(par%outputvars(i))//".max.out", 600+i)
             end if
          end if
       end do
    end do

  end subroutine output_init

  subroutine output_init_data(var, dir, name, fid)

    type(variables_data), intent(inout) :: var
    character(*), intent(in) :: name, dir
    integer*4, intent(in) :: fid

    var%fid = fid
    open(unit=var%fid, file=trim(dir) // trim(name), &
         action="write", status="replace", form="unformatted")

  end subroutine output_init_data

  subroutine output_update(var)

    type(variables), dimension(:), intent(inout) :: var
    integer*4 :: i

    do i = 1,size(var)
       if (var(i)%val%fid > 0) then
          var(i)%n = var(i)%n + 1
          var(i)%sum%rank1 = var(i)%sum%rank1 + var(i)%val%rank1
          var(i)%var%rank1 = var(i)%var%rank1 + var(i)%val%rank1**2
          var(i)%min%rank1 = min(var(i)%min%rank1, var(i)%val%rank1)
          var(i)%max%rank1 = max(var(i)%max%rank1, var(i)%val%rank1)
       end if
    end do
    
  end subroutine output_update

  subroutine output_write(var)

    type(variables), dimension(:), intent(in) :: var
    integer*4 :: i

    do i = 1,size(var)
       if (var(i)%val%fid > 0) &
            write(var(i)%val%fid) var(i)%val%rank1
       if (var(i)%sum%fid > 0) &
            write(var(i)%sum%fid) var(i)%sum%rank1
       if (var(i)%avg%fid > 0) &
            write(var(i)%avg%fid) var(i)%sum%rank1 / var(i)%n
       if (var(i)%var%fid > 0) &
            write(var(i)%var%fid) &
            (var(i)%var%rank1 - var(i)%sum%rank1**2 / var(i)%n) / (var(i)%n - 1)
       if (var(i)%min%fid > 0) &
            write(var(i)%min%fid) var(i)%min%rank1
       if (var(i)%max%fid > 0) &
            write(var(i)%max%fid) var(i)%max%rank1
    end do
    
  end subroutine output_write

  subroutine output_clear(var)

    type(variables), dimension(:), intent(inout) :: var
    integer*4 :: i

    do i = 1,size(var)
       var(i)%n = 0
       if (var(i)%sum%fid > 0 .or. var(i)%avg%fid > 0) &
            call output_clear_data(var(i)%sum, var(i)%rank)
       if (var(i)%var%fid > 0) &
            call output_clear_data(var(i)%var, var(i)%rank)
       if (var(i)%min%fid > 0) &
            call output_clear_data(var(i)%min, var(i)%rank)
       if (var(i)%max%fid > 0) &
            call output_clear_data(var(i)%max, var(i)%rank)
    end do

  end subroutine output_clear

  subroutine output_clear_data(var, rank)

    type(variables_data), intent(inout) :: var
    integer*4, intent(in) :: rank
    
    var%rank1 = var%init

  end subroutine output_clear_data

  subroutine output_close(var)

    type(variables), dimension(:), intent(in) :: var
    integer*4 :: i
    
    do i = 1,size(var)
       if (var(i)%val%fid .gt. 0) close(var(i)%val%fid)
       if (var(i)%sum%fid .gt. 0) close(var(i)%sum%fid)
       if (var(i)%var%fid .gt. 0) close(var(i)%var%fid)
       if (var(i)%min%fid .gt. 0) close(var(i)%min%fid)
       if (var(i)%max%fid .gt. 0) close(var(i)%max%fid)
    end do

  end subroutine output_close
  
  subroutine write_dimensions(par)

    integer*4, parameter :: fid=1
    type(parameters), intent(in) :: par
  
    open(unit=fid, file=trim(par%output_dir) // "dims.out", &
         action="write", status="replace", form="unformatted")
    write(fid) par%nx
    write(fid) par%ny
    write(fid) par%nt
    write(fid) par%dt
    write(fid) par%nfractions
    write(fid) par%nlayers
    write(fid) par%ntout
    write(fid) par%tout
    close(fid)

  end subroutine write_dimensions

  function get_rank(var, name) result (rank)

    type(variables), dimension(:) :: var
    character(*) :: name
    integer*4 :: rank
    integer*4 :: i

    rank = 0
    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          rank = var(i)%rank
          exit
       end if
    end do

  end function get_rank

  function get_shape(var, name) result (dims)

    type(variables), dimension(:) :: var
    character(*) :: name
    integer*4, dimension(:), allocatable :: dims
    integer*4 :: i, rank

    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          rank = var(i)%rank
          if (rank > 0) then
             allocate(dims(rank))
             dims = var(i)%dims
          end if
          exit
       end if
    end do

    if (.not. allocated(dims)) then
       allocate(dims(0))
       dims = 0.d0
    end if

  end function get_shape

  subroutine set_var_reset(var)

    type(variables), dimension(:), intent(inout) :: var
    integer*4 :: i

    do i = 1,size(var)
       var(i)%isset = .false.
    end do
    
  end subroutine set_var_reset
  
  function is_set(var, name) result(ret)

    type(variables), dimension(:), intent(inout) :: var
    character(*) :: name
    logical :: ret
    integer*4 :: i
  
    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          ret = var(i)%isset
          exit
       end if
    end do
    
  end function is_set

  function is_output(var, name) result (ret)

    type(variables), dimension(:) :: var
    character(*) :: name
    logical :: ret
    integer*4 :: i
    
    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          ret = (var(i)%val%fid > 0)
          exit
       end if
    end do

  end function is_output
  
end module output_module

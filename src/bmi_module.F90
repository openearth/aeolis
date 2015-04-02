module bmi_module

  use iso_c_binding
  use iso_c_utils
  use logging
  use output_module
  use init_module
  use moist_module
  use run_module

  implicit none

  type(parameters) :: par
  type(spaceparams) :: s
  type(spaceparams_linear) :: sl
  type(variables), dimension(:), allocatable :: var

  interface get_c_pointer
     module procedure get_c_pointer_rank0
     module procedure get_c_pointer_rank1
     module procedure get_c_pointer_rank2
     module procedure get_c_pointer_rank3
     module procedure get_c_pointer_rank4
  end interface get_c_pointer
    
contains

  integer(c_int) function initialize(c_configfile) result(ierr) bind(C, name="initialize")
    !DEC$ ATTRIBUTES DLLEXPORT::initialize

    implicit none

    ! variables
    character(kind=c_char), intent(in) :: c_configfile(slen)
    character(len=strlen(c_configfile)) :: configfile
    integer*4 :: i

    ! convert c string to fortran string
    configfile = char_array_to_string(c_configfile)
    
    ierr = 0
    write(msgbuf,*) 'Initializing with ', configfile
    call log(LEVEL_INFO, trim(msgbuf))

    ! read parameters file
    par = read_params(configfile)

    ! initialize structures
    call init(par, var, s, sl)
    par%t = 0

    write(*,*) 'Model run started...'

  end function initialize

  integer(c_int) function update(dt) result(ierr) bind(C, name="update")
    !DEC$ ATTRIBUTES DLLEXPORT::update

    real(c_double), value, intent(in) :: dt

    ierr = 0
    write(msgbuf,*) 'Updating with dt: ', dt
    call log(LEVEL_DEBUG, trim(msgbuf))

    if (dt > 0.d0) then
       par%dt = dt
    end if

    call step(par, s, sl, var)

  end function update

  integer(c_int) function finalize() result(ierr) bind(C, name="finalize")
    !DEC$ ATTRIBUTES DLLEXPORT::finalize

    ierr = 0
    call log(LEVEL_INFO, 'Finalize')

    call output_close(var)

    write(*,*) 'Done.'

  end function finalize

  subroutine get_var_type(c_var_name, c_type_name)  bind(C, name="get_var_type")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_var_type

    character(kind=c_char), intent(in) :: c_var_name(*)
    character(kind=c_char), intent(out) :: c_type_name(MAXSTRINGLEN)

    character(len=strlen(c_var_name)) :: var_name
    character(len=MAXSTRINGLEN) :: type_name

    var_name = char_array_to_string(c_var_name)
    type_name = 'double'

    c_type_name = string_to_char_array(trim(type_name))

  end subroutine get_var_type

  subroutine get_var_rank(c_var_name, rank) bind(C, name="get_var_rank")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_var_rank

    character(kind=c_char), intent(in) :: c_var_name(*)
    integer(c_int), intent(out) :: rank
    character(len=strlen(c_var_name)) :: var_name

    var_name = char_array_to_string(c_var_name)
    rank = get_rank(var, var_name)

  end subroutine get_var_rank

  subroutine get_var_shape(c_var_name, var_shape) bind(C, name="get_var_shape")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_var_shape

    character(kind=c_char), intent(in) :: c_var_name(*)
    integer(c_int), intent(inout) :: var_shape(MAXDIMS)
    character(len=strlen(c_var_name)) :: var_name
    integer :: rank

    var_name = char_array_to_string(c_var_name)
    var_shape = get_shape(var, var_name)

    rank = get_rank(var, var_name)
    call arrflip(var_shape, rank)
    
  end subroutine get_var_shape

  subroutine get_var(c_var_name, xptr) bind(C, name="get_var")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_var

    ! return a pointer to the variable
    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), intent(inout) :: xptr
    
    real(c_double), target, save :: x_0d_double_ptr
    real(c_double), allocatable, target, save :: x_1d_double_ptr(:)
    real(c_double), allocatable, target, save :: x_2d_double_ptr(:,:)
    real(c_double), allocatable, target, save :: x_3d_double_ptr(:,:,:)
    real(c_double), allocatable, target, save :: x_4d_double_ptr(:,:,:,:)

    character(len=strlen(c_var_name)) :: var_name
    integer*4, dimension(:), allocatable :: shp
    integer*4 :: rank

    var_name = char_array_to_string(c_var_name)

    rank = get_rank(var, var_name)
    allocate(shp(rank))
    shp = get_shape(var, var_name)

    select case(rank)
    case(0)
       call get_c_pointer(var, var_name, x_0d_double_ptr)
       xptr = c_loc(x_0d_double_ptr)
    case(1)
       if (allocated(x_1d_double_ptr)) then
          deallocate(x_1d_double_ptr)
       end if
       allocate(x_1d_double_ptr(shp(1)))
       call get_c_pointer(var, var_name, x_1d_double_ptr)
       xptr = c_loc(x_1d_double_ptr)
     case(2)
       if (allocated(x_2d_double_ptr)) then
          deallocate(x_2d_double_ptr)
       end if
       allocate(x_2d_double_ptr(shp(1), shp(2)))
       call get_c_pointer(var, var_name, x_2d_double_ptr)
       xptr = c_loc(x_2d_double_ptr)
     case(3)
       if (allocated(x_3d_double_ptr)) then
          deallocate(x_3d_double_ptr)
       end if
       allocate(x_3d_double_ptr(shp(1), shp(2), shp(3)))
       call get_c_pointer(var, var_name, x_3d_double_ptr)
       xptr = c_loc(x_3d_double_ptr)
     case(4)
       if (allocated(x_4d_double_ptr)) then
          deallocate(x_4d_double_ptr)
       end if
       allocate(x_4d_double_ptr(shp(1), shp(2), shp(3), shp(4)))
       call get_c_pointer(var, var_name, x_4d_double_ptr)
       xptr = c_loc(x_4d_double_ptr)
     end select

  end subroutine get_var

  subroutine set_var(c_var_name, xptr) bind(C, name="set_var")
    !DEC$ ATTRIBUTES DLLEXPORT :: set_var
    
    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), value, intent(in) :: xptr

    character(len=strlen(c_var_name)) :: var_name
    integer*4, dimension(:), allocatable :: shp
    integer*4 :: rank

    var_name = char_array_to_string(c_var_name)

    rank = get_rank(var, var_name)
    allocate(shp(rank))
    shp = get_shape(var, var_name)
    
    call set_c_pointer(var, var_name, xptr, rank, shp)

  end subroutine set_var
  
  subroutine get_current_time(time) bind(C, name="get_current_time")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_current_time

    real(c_double) :: time
    time = par%t

  end subroutine get_current_time

  subroutine get_start_time(time) bind(C, name="get_start_time")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_start_time

    real(c_double) :: time
    time = 0

  end subroutine get_start_time

  subroutine get_end_time(time) bind(C, name="get_end_time")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_end_time

    real(c_double) :: time
    time = par%tstop

  end subroutine get_end_time

  subroutine get_c_pointer_rank0(var, name, val)

    type(variables), dimension(:), intent(in) :: var
    character(*), intent(in) :: name
    real(c_double), intent(out) :: val
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          val = var(i)%val%rank0
          exit
       end if
    end do
    
  end subroutine get_c_pointer_rank0

  subroutine get_c_pointer_rank1(var, name, val)

    type(variables), dimension(:), intent(in) :: var
    character(*), intent(in) :: name
    real(c_double), intent(out) :: val(:)
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          val = var(i)%val%rank1
          exit
       end if
    end do
    
  end subroutine get_c_pointer_rank1

  subroutine get_c_pointer_rank2(var, name, val)

    type(variables), dimension(:), intent(in) :: var
    character(*), intent(in) :: name
    real(c_double), intent(out) :: val(:,:)
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          val = var(i)%val%rank2
          exit
       end if
    end do
    
  end subroutine get_c_pointer_rank2

  subroutine get_c_pointer_rank3(var, name, val)

    type(variables), dimension(:), intent(in) :: var
    character(*), intent(in) :: name
    real(c_double), intent(out) :: val(:,:,:)
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          val = var(i)%val%rank3
          exit
       end if
    end do
    
  end subroutine get_c_pointer_rank3

  subroutine get_c_pointer_rank4(var, name, val)

    type(variables), dimension(:), intent(in) :: var
    character(*), intent(in) :: name
    real(c_double), intent(out) :: val(:,:,:,:)
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          val = var(i)%val%rank4
          exit
       end if
    end do
    
  end subroutine get_c_pointer_rank4

  subroutine set_c_pointer(var, name, val, rank, shp)

    type(variables), dimension(:), intent(in) :: var
    character(*), intent(in) :: name
    type(c_ptr), intent(in) :: val
    integer*4, intent(in) :: rank
    integer*4, dimension(:), intent(in) :: shp
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          select case (rank)
          case(0)
             call c_f_pointer(val, var(i)%val%rank0, shp)
          case(1)
             call c_f_pointer(val, var(i)%val%rank1, shp)
          case(2)
             call c_f_pointer(val, var(i)%val%rank2, shp)
          case(3)
             call c_f_pointer(val, var(i)%val%rank3, shp)
          end select
          exit
       end if
    end do
    
  end subroutine set_c_pointer

end module bmi_module

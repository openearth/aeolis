module bmi_module

  use iso_c_binding
  use iso_c_utils
  use init_module
  use moist_module
  use run_module
  use utils_module
  use bed_module

  implicit none

  type(parameters) :: par
  type(spaceparams) :: s
  type(spaceparams_linear) :: sl
  type(variables), dimension(:), allocatable :: var
  logical, save :: clear_output = .false.

  interface get_c_pointer
     module procedure get_c_pointer_rank0
     module procedure get_c_pointer_rank1
     module procedure get_c_pointer_rank2
     module procedure get_c_pointer_rank3
     module procedure get_c_pointer_rank4
  end interface get_c_pointer

  interface set_c_pointer
     module procedure set_c_pointer_rank0
     module procedure set_c_pointer_rank1
  end interface set_c_pointer

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
    
    ! read parameters file
    par = read_params(configfile)

    ! initialize structures
    call init(par, var, s, sl)
    par%t = 0

    ! initialize output
    call output_init(par, var)

  end function initialize

  integer(c_int) function update(dt) result(ierr) bind(C, name="update")
    !DEC$ ATTRIBUTES DLLEXPORT::update

    real(c_double), value, intent(in) :: dt

    if (dt > 0.d0) then
       par%dt = dt
    end if

    if (clear_output) then
       call output_clear(var)
       clear_output = .false.
    end if

    call step(par, s, sl, var)
    call output_update(var)
    call set_var_reset(var)
    
  end function update

  integer(c_int) function finalize() result(ierr) bind(C, name="finalize")
    !DEC$ ATTRIBUTES DLLEXPORT::finalize

    call output_close(var)
    
    if (allocated(par%grain_size)) then
       deallocate(par%grain_size)
    end if
    if (allocated(par%grain_dist)) then
       deallocate(par%grain_dist)
    end if
    if (allocated(par%uw)) then
       deallocate(par%uw)
    end if
    if (allocated(par%zs)) then
       deallocate(par%zs)
    end if
    if (allocated(par%meteo)) then
       deallocate(par%meteo)
    end if
    if (allocated(par%moist)) then
       deallocate(par%moist)
    end if    

    deallocate(s%xu)
    deallocate(s%yu)
    deallocate(s%xv)
    deallocate(s%yv)
    deallocate(s%dsz)
    deallocate(s%dnz)
    deallocate(s%dsdnzi)
    deallocate(s%alfaz)
    deallocate(s%dsu)
    deallocate(s%dnu)
    deallocate(s%dsdnui)
    deallocate(s%alfau)
    deallocate(s%dsv)
    deallocate(s%dnv)
    deallocate(s%dsdnvi)
    deallocate(s%alfav)
    deallocate(s%dsc)
    deallocate(s%dnc)

    call dealloc_variables(var)
    deallocate(var)

    call dealloc_bedcomposition()

    nullify(s%uw)
    nullify(s%udir)
    nullify(s%zs)
    nullify(s%rho)
    nullify(s%dist)
    nullify(s%xz)
    nullify(s%yz)
    nullify(s%xu)
    nullify(s%yu)
    nullify(s%xv)
    nullify(s%yv)
    nullify(s%zb)
    nullify(s%uws)
    nullify(s%uwn)
    nullify(s%uth)
    nullify(s%moist)
    nullify(s%Cu)
    nullify(s%Ct)
    nullify(s%supply)
    nullify(s%thlyr)
    nullify(s%p)
    nullify(s%d10)
    nullify(s%d50)
    nullify(s%d90)
    nullify(s%mass)
    nullify(s%dsz)
    nullify(s%dnz)
    nullify(s%dsdnzi)
    nullify(s%alfaz)
    nullify(s%dsu)
    nullify(s%dnu)
    nullify(s%dsdnui)
    nullify(s%alfau)
    nullify(s%dsv)
    nullify(s%dnv)
    nullify(s%dsdnvi)
    nullify(s%alfav)
    nullify(s%dsc)
    nullify(s%dnc)

!    write(0,*) 'Done.'

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
    character(slen) :: name, type
    character(len=strlen(c_var_name)) :: var_name

    var_name = char_array_to_string(c_var_name)
    call split_var(var_name, name, type)
    
    rank = get_rank(var, name)

  end subroutine get_var_rank

  subroutine get_var_shape(c_var_name, var_shape) bind(C, name="get_var_shape")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_var_shape

    character(kind=c_char), intent(in) :: c_var_name(*)
    integer(c_int), intent(inout) :: var_shape(MAXDIMS)
    character(len=strlen(c_var_name)) :: var_name
    character(slen) :: name, type
    integer :: rank

    var_name = char_array_to_string(c_var_name)
    call split_var(var_name, name, type)

    var_shape = get_shape(var, name)
    rank = get_rank(var, name)
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

    character(10), dimension(:), allocatable :: var_name_parts
    character(slen) :: name, type

    var_name = char_array_to_string(c_var_name)
    call split_var(var_name, name, type)

    rank = get_rank(var, name)
    if (rank == 0) then
       call get_c_pointer(par, name, x_0d_double_ptr)
       xptr = c_loc(x_0d_double_ptr)
    else
       allocate(shp(rank))
       shp = get_shape(var, name)
       clear_output = .true.

       select case(rank)
       case(1)
          if (allocated(x_1d_double_ptr)) then
             deallocate(x_1d_double_ptr)
          end if
          allocate(x_1d_double_ptr(shp(1)))
          call get_c_pointer(var, name, type, x_1d_double_ptr)
          xptr = c_loc(x_1d_double_ptr)
       case(2)
          if (allocated(x_2d_double_ptr)) then
             deallocate(x_2d_double_ptr)
          end if
          allocate(x_2d_double_ptr(shp(1), shp(2)))
          call get_c_pointer(var, name, type, x_2d_double_ptr)
          xptr = c_loc(x_2d_double_ptr)
       case(3)
          if (allocated(x_3d_double_ptr)) then
             deallocate(x_3d_double_ptr)
          end if
          allocate(x_3d_double_ptr(shp(1), shp(2), shp(3)))
          call get_c_pointer(var, name, type, x_3d_double_ptr)
          xptr = c_loc(x_3d_double_ptr)
       case(4)
          if (allocated(x_4d_double_ptr)) then
             deallocate(x_4d_double_ptr)
          end if
          allocate(x_4d_double_ptr(shp(1), shp(2), shp(3), shp(4)))
          call get_c_pointer(var, name, type, x_4d_double_ptr)
          xptr = c_loc(x_4d_double_ptr)
       end select
    end if

  end subroutine get_var

  subroutine set_var(c_var_name, xptr) bind(C, name="set_var")
    !DEC$ ATTRIBUTES DLLEXPORT :: set_var
    
    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), value, intent(in) :: xptr

    character(len=strlen(c_var_name)) :: var_name
    integer*4, dimension(:), allocatable :: shp
    integer*4 :: rank

    real(c_double), pointer :: x_2d_double_ptr(:,:)
    real*8, dimension(par%nx+1,par%ny+1) :: zbx

    var_name = char_array_to_string(c_var_name)

    rank = get_rank(var, var_name)
    if (rank == 0) then
       call set_c_pointer(par, var_name, xptr)
    else
       allocate(shp(rank))
       shp = get_shape(var, var_name)
       if (var_name == 'zbx') then
          call c_f_pointer(xptr, x_2d_double_ptr, shp)
          zbx = x_2d_double_ptr
          call update_bedlevel(par, sl, zbx)
       else
          call set_c_pointer(var, var_name, xptr, rank, shp)
       endif
    endif

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

  subroutine get_c_pointer_rank0(par, name, val)

    type(parameters), intent(in) :: par
    character(*), intent(in) :: name
    real(c_double), intent(out) :: val

    if (trim(name) == 'dt') then
       val = par%dt
    elseif (trim(name) == 'accfac') then
       val = par%accfac
    end if
    
  end subroutine get_c_pointer_rank0

  subroutine get_c_pointer_rank1(var, name, type, val)

    type(variables), dimension(:), intent(in) :: var
    character(*), intent(in) :: name, type
    real(c_double), intent(out) :: val(:)
    integer*4 :: i
    
    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          select case (trim(type))
          case ('sum')
             val = var(i)%sum%rank1
          case ('avg')
             val = var(i)%sum%rank1 / var(i)%n
          case ('var')
             val = (var(i)%var%rank1 - var(i)%sum%rank1**2 / var(i)%n) / (var(i)%n - 1)
          case ('min')
             val = var(i)%min%rank1
          case ('max')
             val = var(i)%max%rank1
          case default
             val = var(i)%val%rank1
          end select
          exit
       end if
    end do
    
  end subroutine get_c_pointer_rank1

  subroutine get_c_pointer_rank2(var, name, type, val)

    type(variables), dimension(:), intent(in) :: var
    character(*), intent(in) :: name, type
    real(c_double), intent(out) :: val(:,:)
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          select case (trim(type))
          case ('sum')
             val = var(i)%sum%rank2
          case ('avg')
             val = var(i)%sum%rank2 / var(i)%n
          case ('var')
             val = (var(i)%var%rank2 - var(i)%sum%rank2**2 / var(i)%n) / (var(i)%n - 1)
          case ('min')
             val = var(i)%min%rank2
          case ('max')
             val = var(i)%max%rank2
          case default
             val = var(i)%val%rank2
          end select
          exit
       end if
    end do
    
  end subroutine get_c_pointer_rank2

  subroutine get_c_pointer_rank3(var, name, type, val)

    type(variables), dimension(:), intent(in) :: var
    character(*), intent(in) :: name, type
    real(c_double), intent(out) :: val(:,:,:)
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          select case (trim(type))
          case ('sum')
             val = var(i)%sum%rank3
          case ('avg')
             val = var(i)%sum%rank3 / var(i)%n
          case ('var')
             val = (var(i)%var%rank3 - var(i)%sum%rank3**2 / var(i)%n) / (var(i)%n - 1)
          case ('min')
             val = var(i)%min%rank3
          case ('max')
             val = var(i)%max%rank3
          case default
             val = var(i)%val%rank3
          end select
          exit
       end if
    end do
    
  end subroutine get_c_pointer_rank3

  subroutine get_c_pointer_rank4(var, name, type, val)

    type(variables), dimension(:), intent(in) :: var
    character(*), intent(in) :: name, type
    real(c_double), intent(out) :: val(:,:,:,:)
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          select case (trim(type))
          case ('sum')
             val = var(i)%sum%rank4
          case ('avg')
             val = var(i)%sum%rank4 / var(i)%n
          case ('var')
             val = (var(i)%var%rank4 - var(i)%sum%rank4**2 / var(i)%n) / (var(i)%n - 1)
          case ('min')
             val = var(i)%min%rank4
          case ('max')
             val = var(i)%max%rank4
          case default
             val = var(i)%val%rank4
          end select
          exit
       end if
    end do
    
  end subroutine get_c_pointer_rank4

  subroutine set_c_pointer_rank0(par, name, val)
    
    type(parameters), intent(inout) :: par
    character(*), intent(in) :: name
    type(c_ptr), intent(in) :: val
    
    real(c_double), pointer :: x_1d_double_ptr(:)

    call c_f_pointer(val, x_1d_double_ptr, (/ 1 /))
    
    if (trim(name) == 'dt') then
       par%dt = x_1d_double_ptr(1)
    elseif (trim(name) == 'accfac') then
       par%accfac = x_1d_double_ptr(1)
    end if

  end subroutine set_c_pointer_rank0
  
  subroutine set_c_pointer_rank1(var, name, val, rank, shp)

    type(variables), dimension(:), intent(inout) :: var
    character(*), intent(in) :: name
    type(c_ptr), intent(in) :: val
    integer*4, intent(in) :: rank
    integer*4, dimension(:), intent(in) :: shp
    integer*4 :: i

    real(c_double), pointer :: x_1d_double_ptr(:)
    real(c_double), pointer :: x_2d_double_ptr(:,:)
    real(c_double), pointer :: x_3d_double_ptr(:,:,:)
    
    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          select case (rank)
          case(0)
             exit
          case(1)
             call c_f_pointer(val, x_1d_double_ptr, shp)
             var(i)%val%rank1 = x_1d_double_ptr
          case(2)
             call c_f_pointer(val, x_2d_double_ptr, shp)
             var(i)%val%rank2 = x_2d_double_ptr
          case(3)
             call c_f_pointer(val, x_3d_double_ptr, shp)
             var(i)%val%rank3 = x_3d_double_ptr
          end select
          var(i)%isset = .true.
          exit
       end if
    end do
    
  end subroutine set_c_pointer_rank1

end module bmi_module

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
  end type variables

  type variables_data
     real*8, pointer :: rank0 => null()
     real*8, dimension(:), pointer :: rank1 => null()
     real*8, dimension(:,:), pointer :: rank2 => null()
     real*8, dimension(:,:,:), pointer :: rank3 => null()
     real*8, dimension(:,:,:,:), pointer :: rank4 => null()
     real*8 :: init
     integer*4 :: fid = -1
  end type variables_data

  interface get_pointer
     module procedure get_pointer_rank0
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

    select case (rank)
    case (0)
       allocate(var%rank0)
       var%rank0 = var%init
    case default
       allocate(var%rank1(product(dims)))
       var%rank1 = var%init
    end select

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
    if (associated(var%rank0)) then
       deallocate(var%rank0)
    end if

    if (associated(var)) then
       deallocate(var)
    end if

  end subroutine dealloc_variable_data

  subroutine get_pointer_rank0(var, name, dims, ptr)

    type(variables), dimension(:), allocatable, intent(inout) :: var
    character(*), intent(in) :: name
    integer*4, dimension(:), intent(in) :: dims
    real*8, pointer, intent(inout) :: ptr
    integer*4 :: i

    do i = 1,size(var)
       if (trim(var(i)%name) == trim(name)) then
          ptr => var(i)%val%rank0
          exit
       end if
    end do
    
  end subroutine get_pointer_rank0

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
          select case (var(i)%rank)
          case (0)
             var(i)%sum%rank0 = var(i)%sum%rank0 + var(i)%val%rank0
             var(i)%var%rank0 = var(i)%var%rank0 + var(i)%val%rank0**2
             var(i)%min%rank0 = min(var(i)%min%rank0, var(i)%val%rank0)
             var(i)%max%rank0 = max(var(i)%max%rank0, var(i)%val%rank0)
!          case (1)
          case default
             var(i)%sum%rank1 = var(i)%sum%rank1 + var(i)%val%rank1
             var(i)%var%rank1 = var(i)%var%rank1 + var(i)%val%rank1**2
             var(i)%min%rank1 = min(var(i)%min%rank1, var(i)%val%rank1)
             var(i)%max%rank1 = max(var(i)%max%rank1, var(i)%val%rank1)
!          case (2)
!             var(i)%sum%rank2 = var(i)%sum%rank2 + var(i)%val%rank2
!             var(i)%var%rank2 = var(i)%var%rank2 + var(i)%val%rank2**2
!             var(i)%min%rank2 = min(var(i)%min%rank2, var(i)%val%rank2)
!             var(i)%max%rank2 = max(var(i)%max%rank2, var(i)%val%rank2)
!          case (3)
!             var(i)%sum%rank3 = var(i)%sum%rank3 + var(i)%val%rank3
!             var(i)%var%rank3 = var(i)%var%rank3 + var(i)%val%rank3**2
!             var(i)%min%rank3 = min(var(i)%min%rank3, var(i)%val%rank3)
!             var(i)%max%rank3 = max(var(i)%max%rank3, var(i)%val%rank3)
!          case (4)
!             var(i)%sum%rank4 = var(i)%sum%rank4 + var(i)%val%rank4
!             var(i)%var%rank4 = var(i)%var%rank4 + var(i)%val%rank4**2
!             var(i)%min%rank4 = min(var(i)%min%rank4, var(i)%val%rank4)
!             var(i)%max%rank4 = max(var(i)%max%rank4, var(i)%val%rank4)
          end select
       end if
    end do
    
  end subroutine output_update

  subroutine output_write(var)

    type(variables), dimension(:), intent(in) :: var
    integer*4 :: i

    do i = 1,size(var)
       select case (var(i)%rank)
       case (0)
          if (var(i)%val%fid > 0) &
               write(var(i)%val%fid) var(i)%val%rank0
          if (var(i)%sum%fid > 0) &
               write(var(i)%sum%fid) var(i)%sum%rank0
          if (var(i)%avg%fid > 0) &
               write(var(i)%avg%fid) var(i)%sum%rank0 / var(i)%n
          if (var(i)%var%fid > 0) &
               write(var(i)%var%fid) &
               (var(i)%var%rank0 - var(i)%sum%rank0**2 / var(i)%n) / (var(i)%n - 1)
          if (var(i)%min%fid > 0) &
               write(var(i)%min%fid) var(i)%min%rank0
          if (var(i)%max%fid > 0) &
               write(var(i)%max%fid) var(i)%max%rank0
       case (1)
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
       case (2)
          if (var(i)%val%fid > 0) &
               write(var(i)%val%fid) var(i)%val%rank2
          if (var(i)%sum%fid > 0) &
               write(var(i)%sum%fid) var(i)%sum%rank2
          if (var(i)%avg%fid > 0) &
               write(var(i)%avg%fid) var(i)%sum%rank2 / var(i)%n
          if (var(i)%var%fid > 0) &
               write(var(i)%var%fid) &
               (var(i)%var%rank2 - var(i)%sum%rank2**2 / var(i)%n) / (var(i)%n - 1)
          if (var(i)%min%fid > 0) &
               write(var(i)%min%fid) var(i)%min%rank2
          if (var(i)%max%fid > 0) &
               write(var(i)%max%fid) var(i)%max%rank2
       case (3)
          if (var(i)%val%fid > 0) &
               write(var(i)%val%fid) var(i)%val%rank3
          if (var(i)%sum%fid > 0) &
               write(var(i)%sum%fid) var(i)%sum%rank3
          if (var(i)%avg%fid > 0) &
               write(var(i)%avg%fid) var(i)%sum%rank3 / var(i)%n
          if (var(i)%var%fid > 0) &
               write(var(i)%var%fid) &
               (var(i)%var%rank3 - var(i)%sum%rank3**2 / var(i)%n) / (var(i)%n - 1)
          if (var(i)%min%fid > 0) &
               write(var(i)%min%fid) var(i)%min%rank3
          if (var(i)%max%fid > 0) &
               write(var(i)%max%fid) var(i)%max%rank3
       case (4)
          if (var(i)%val%fid > 0) &
               write(var(i)%val%fid) var(i)%val%rank4
          if (var(i)%sum%fid > 0) &
               write(var(i)%sum%fid) var(i)%sum%rank4
          if (var(i)%avg%fid > 0) &
               write(var(i)%avg%fid) var(i)%sum%rank4 / var(i)%n
          if (var(i)%var%fid > 0) &
               write(var(i)%var%fid) &
               (var(i)%var%rank4 - var(i)%sum%rank4**2 / var(i)%n) / (var(i)%n - 1)
          if (var(i)%min%fid > 0) &
               write(var(i)%min%fid) var(i)%min%rank4
          if (var(i)%max%fid > 0) &
               write(var(i)%max%fid) var(i)%max%rank4
       end select
    end do
    
  end subroutine output_write

  subroutine output_clear(var)

    type(variables), dimension(:), intent(inout) :: var
    integer*4 :: i

    do i = 1,size(var)
       var(i)%n = 0
       call output_clear_data(var(i)%sum, var(i)%rank)
       call output_clear_data(var(i)%var, var(i)%rank)
       call output_clear_data(var(i)%min, var(i)%rank)
       call output_clear_data(var(i)%max, var(i)%rank)
    end do

  end subroutine output_clear

  subroutine output_clear_data(var, rank)

    type(variables_data), intent(inout) :: var
    integer*4, intent(in) :: rank
    
    if (var%fid > 0) then
       select case (rank)
       case (0)
          var%rank0 = var%init
!       case (1)
       case default
          var%rank1 = var%init
!       case (2)
!          var%rank2 = var%init
!       case (3)
!          var%rank3 = var%init
!       case (4)
!          var%rank4 = var%init
       end select
    end if

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

  subroutine write_output(par, sl, var)

    type(parameters), intent(inout) :: par
    type(spaceparams_linear), intent(inout) :: sl
    type(variables), dimension(:), intent(inout) :: var

    ! write output
    if (par%t .le. par%dt  .or. par%tout < par%dt .or. &
         mod(par%t, par%tout) < par%dt) then

       ! update derived variables
       if (is_output(var, 'mass')) sl%mass = get_layer_mass(par)
       if (is_output(var, 'd10')) sl%d10 = get_layer_percentile(par, 0.1d0)
       if (is_output(var, 'd50')) sl%d50 = get_layer_percentile(par, 0.5d0)
       if (is_output(var, 'd90')) sl%d90 = get_layer_percentile(par, 0.9d0)
       
       call output_write(var)
       call output_clear(var)

    end if
    
  end subroutine write_output
  
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
          allocate(dims(rank))
          dims = var(i)%dims
!          select case (rank)
!             case(0)
!                shp = shape(var(i)%val%rank0)
!             case(1)
!                shp = shape(var(i)%val%rank1)
!             case(2)
!                shp = shape(var(i)%val%rank2)
!             case(3)
!                shp = shape(var(i)%val%rank3)
!             case(4)
!                shp = shape(var(i)%val%rank4)
!             end select
          exit
       end if
    end do

  end function get_shape

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

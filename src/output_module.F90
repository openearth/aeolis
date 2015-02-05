module output_module

  use constants_module
  use input_module
  
  implicit none

  type output_files
     character(slen) :: var
     integer*4 :: fid
  end type output_files

  interface output_write
     module procedure output_write_rank0
     module procedure output_write_rank1
     module procedure output_write_rank2
     module procedure output_write_rank3
  end interface output_write

contains

  function output_init(vars) result (out)

    type(output_files), dimension(:), allocatable :: out
    character(*), dimension(:) :: vars
    integer*4 :: n, i

    n = size(vars)
    
    allocate(out(n))
    
    do i = 1, n
       out(i)%var = trim(vars(i))
       out(i)%fid = 100 + i
       
       open(unit=out(i)%fid, file=trim(vars(i))//".out", &
            action="write", status="replace", form="unformatted")
    end do

  end function output_init

  subroutine output_write_rank0(out, var, val)

    type(output_files), dimension(:), intent(in) :: out
    character(*) :: var
    real*8 :: val
    integer*4 :: fid

    fid = get_fid(out, var)

    if (fid > 0) then
       write(fid) val
    end if

  end subroutine output_write_rank0

  subroutine output_write_rank1(out, var, val)

    type(output_files), dimension(:), intent(in) :: out
    character(*) :: var
    real*8, dimension(:) :: val
    integer*4 :: fid

    fid = get_fid(out, var)

    if (fid > 0) then
       write(fid) val
    end if

  end subroutine output_write_rank1
  
  subroutine output_write_rank2(out, var, val)

    type(output_files), dimension(:), intent(in) :: out
    character(*) :: var
    real*8, dimension(:,:) :: val
    integer*4 :: fid

    fid = get_fid(out, var)

    if (fid > 0) then
       write(fid) val
    end if

  end subroutine output_write_rank2
  
  subroutine output_write_rank3(out, var, val)

    type(output_files), dimension(:), intent(in) :: out
    character(*) :: var
    real*8, dimension(:,:,:) :: val
    integer*4 :: fid

    fid = get_fid(out, var)

    if (fid > 0) then
       write(fid) val
    end if

  end subroutine output_write_rank3

  subroutine output_close(out)

    type(output_files), dimension(:), intent(inout) :: out
    integer*4 :: i
    
    do i = 1, size(out)
       close(out(i)%fid)
    end do

  end subroutine output_close

  function get_fid(out, var) result (fid)

    type(output_files), dimension(:) :: out
    character(*) :: var
    integer*4 :: fid, i

    fid = -1
    
    do i = 1, size(out)
       if (trim(out(i)%var) .eq. trim(var)) then
          fid = out(i)%fid
          exit
       end if
    end do

  end function get_fid
  
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
    close(fid)

  end subroutine write_dimensions

end module output_module

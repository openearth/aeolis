module init_module

  use input_module
  use output_module
  use bed_module
  use wind_module
  use moist_module
  
  implicit none

contains

  subroutine init(par, var, s, sl)

    type(parameters), intent(inout) :: par
    type(spaceparams), intent(out) :: s
    type(spaceparams_linear), intent(out) :: sl
    type(variables), dimension(:), allocatable, intent(out) :: var
    real*8, dimension(:,:), allocatable :: x_tmp, y_tmp, zb_tmp
    integer*4, parameter :: fid=20

    write(0,*) 'Initialization started...'
    
    ! bed
    write(0,*) 'Generating bed profile and composition...'
    call generate_bed(par, x_tmp, y_tmp, zb_tmp)

    call alloc_variable(var, 'x',  (/par%nx+1, par%ny+1/))
    call alloc_variable(var, 'y',  (/par%nx+1, par%ny+1/))
    call alloc_variable(var, 'zb', (/par%nx+1, par%ny+1/))

    call get_pointer(var, 'x',  (/par%nx+1, par%ny+1/), s%xz)
    call get_pointer(var, 'y',  (/par%nx+1, par%ny+1/), s%yz)
    call get_pointer(var, 'zb', (/par%nx+1, par%ny+1/), s%zb)
    
    call get_pointer(var, 'x',  (/par%nc/), sl%xz)
    call get_pointer(var, 'y',  (/par%nc/), sl%yz)
    call get_pointer(var, 'zb', (/par%nc/), sl%zb)

    s%xz = x_tmp
    s%yz = y_tmp
    s%zb = zb_tmp
    deallocate(x_tmp)
    deallocate(y_tmp)
    deallocate(zb_tmp)

    call generate_bedcomposition(par)
!    open(unit=fid, file=trim(par%output_dir) // "bed.in", &
!         action="write", status="replace", form="unformatted")
!    write(fid) s%xz
!    write(fid) s%yz
!    write(fid) s%zb
!    close(fid)

    ! wind
    write(0,*) 'Generating wind time series...'
    call generate_wind(par, par%uw)
!    open(unit=fid, file=trim(par%output_dir) // "wind.in", &
!         action="write", status="replace", form="unformatted")
!    write(fid) par%uw%t
!    write(fid) par%uw%u
!    write(fid) par%uw%dir
!    close(fid)

    ! time
    par%ntout = int(par%tstop / par%tout) + 1

    ! moist
    write(0,*) 'Generating moisture time series...'
    call generate_moist(par, par%moist)

    ! tide and meteo
    write(0,*) 'Generating tide and meteo time series...'
    call generate_tide(par, par%zs)
    call generate_meteo(par, par%meteo)

    ! fractions
    call alloc_variable(var, 'rho',    (/par%nfractions/))
    call alloc_variable(var, 'dist',   (/par%nfractions/))

    ! variables
    call alloc_variable(var, 'Cu',     (/par%nfractions, par%nx+1, par%ny+1/))
    call alloc_variable(var, 'Ct',     (/par%nfractions, par%nx+1, par%ny+1/))
    call alloc_variable(var, 'p',      (/par%nfractions, par%nx+1, par%ny+1/))
    call alloc_variable(var, 'uth',    (/par%nfractions, par%nx+1, par%ny+1/))
    call alloc_variable(var, 'mass',   (/par%nfractions, par%nlayers, par%nx+1, par%ny+1/))
    call alloc_variable(var, 'supply', (/par%nfractions, par%nx+1, par%ny+1/))
    call alloc_variable(var, 'moist',  (/par%nlayers, par%nx+1, par%ny+1/))
    call alloc_variable(var, 'thlyr',  (/par%nlayers, par%nx+1, par%ny+1/))

    ! extra output
    call alloc_variable(var, 'uw',     (/par%nx+1, par%ny+1/))
    call alloc_variable(var, 'uws',    (/par%nx+1, par%ny+1/))
    call alloc_variable(var, 'uwn',    (/par%nx+1, par%ny+1/))
    call alloc_variable(var, 'udir',   (/par%nx+1, par%ny+1/))
    call alloc_variable(var, 'zs',     (/par%nx+1, par%ny+1/))
!    call alloc_variable(var, 'd10',    (/par%nlayers, par%nx+1, par%ny+1/))
!    call alloc_variable(var, 'd50',    (/par%nlayers, par%nx+1, par%ny+1/))
!    call alloc_variable(var, 'd90',    (/par%nlayers, par%nx+1, par%ny+1/))

    ! create remapping pointers
    call get_pointer(var, 'rho',    (/par%nfractions/), s%rho)
    call get_pointer(var, 'dist',   (/par%nfractions/), s%dist)
    call get_pointer(var, 'Cu',     (/par%nfractions, par%nx+1, par%ny+1/), s%Cu)
    call get_pointer(var, 'Ct',     (/par%nfractions, par%nx+1, par%ny+1/), s%Ct)
    call get_pointer(var, 'p',      (/par%nfractions, par%nx+1, par%ny+1/), s%p)
    call get_pointer(var, 'uth',    (/par%nfractions, par%nx+1, par%ny+1/), s%uth)
    call get_pointer(var, 'mass',   (/par%nfractions, par%nlayers, par%nx+1, par%ny+1/), s%mass)
    call get_pointer(var, 'supply', (/par%nfractions, par%nx+1, par%ny+1/), s%supply)
    call get_pointer(var, 'moist',  (/par%nlayers, par%nx+1, par%ny+1/), s%moist)
    call get_pointer(var, 'thlyr',  (/par%nlayers, par%nx+1, par%ny+1/), s%thlyr)
    call get_pointer(var, 'uw',     (/par%nx+1, par%ny+1/), s%uw)
    call get_pointer(var, 'uws',    (/par%nx+1, par%ny+1/), s%uws)
    call get_pointer(var, 'uwn',    (/par%nx+1, par%ny+1/), s%uwn)
    call get_pointer(var, 'udir',   (/par%nx+1, par%ny+1/), s%udir)
    call get_pointer(var, 'zs',     (/par%nx+1, par%ny+1/), s%zs)
!    call get_pointer(var, 'd10',    (/par%nlayers, par%nx+1, par%ny+1/), s%d10)
!    call get_pointer(var, 'd50',    (/par%nlayers, par%nx+1, par%ny+1/), s%d50)
!    call get_pointer(var, 'd90',    (/par%nlayers, par%nx+1, par%ny+1/), s%d90)

    call get_pointer(var, 'rho',    (/par%nfractions/), sl%rho)
    call get_pointer(var, 'dist',   (/par%nfractions/), sl%dist)
    call get_pointer(var, 'Cu',     (/par%nfractions, par%nc/), sl%Cu)
    call get_pointer(var, 'Ct',     (/par%nfractions, par%nc/), sl%Ct)
    call get_pointer(var, 'p',      (/par%nfractions, par%nc/), sl%p)
    call get_pointer(var, 'uth',    (/par%nfractions, par%nc/), sl%uth)
    call get_pointer(var, 'mass',   (/par%nfractions, par%nlayers, par%nc/), sl%mass)
    call get_pointer(var, 'supply', (/par%nfractions, par%nc/), sl%supply)
    call get_pointer(var, 'moist',  (/par%nlayers, par%nc/), sl%moist)
    call get_pointer(var, 'thlyr',  (/par%nlayers, par%nc/), sl%thlyr)
    call get_pointer(var, 'uw',     (/par%nc/), sl%uw)
    call get_pointer(var, 'uws',    (/par%nc/), sl%uws)
    call get_pointer(var, 'uwn',    (/par%nc/), sl%uwn)
    call get_pointer(var, 'udir',   (/par%nc/), sl%udir)
    call get_pointer(var, 'zs',     (/par%nc/), sl%zs)
!    call get_pointer(var, 'd10',    (/par%nlayers, par%nc/), sl%d10)
!    call get_pointer(var, 'd50',    (/par%nlayers, par%nc/), sl%d50)
!    call get_pointer(var, 'd90',    (/par%nlayers, par%nc/), sl%d90)

    s%rho = par%rhom
    s%dist = par%grain_dist

    ! create spatial grid matrixes
    call gridprops(par, s)

    write(0,*) 'Starting simulation...'

  end subroutine init

  subroutine gridprops(par, s)

    ! THIS FUNCTION IS BASED ON THE XBEACH GRIDPROPS FUNCTION
    ! https://svn.oss.deltares.nl/repos/xbeach/trunk/src/xbeachlibrary/spaceparams.F90
    
    implicit none

    type(parameters), intent(in) :: par
    type(spaceparams), intent(inout) :: s

    integer                           :: i, j
    real*8,dimension(:,:),allocatable :: xc      ! x-coordinate c-points
    real*8,dimension(:,:),allocatable :: yc      ! y-coordinate c-points
    real*8                            :: dsdnu   ! surface of cell centered around u-point
    real*8                            :: dsdnv   ! surface of cell centered around v-point
    real*8                            :: dsdnz   ! surface of cell centered around z-point
    real*8                            :: x1,y1,x2,y2,x3,y3,x4,y4

    allocate(s%xu(par%nx+1, par%ny+1))
    allocate(s%yu(par%nx+1, par%ny+1))
    allocate(s%xv(par%nx+1, par%ny+1))
    allocate(s%yv(par%nx+1, par%ny+1))
    allocate(s%dsz(par%nx+1, par%ny+1))
    allocate(s%dnz(par%nx+1, par%ny+1))
    allocate(s%dsdnzi(par%nx+1, par%ny+1))
    allocate(s%alfaz(par%nx+1, par%ny+1))
    allocate(s%dsu(par%nx+1, par%ny+1))
    allocate(s%dnu(par%nx+1, par%ny+1))
    allocate(s%dsdnui(par%nx+1, par%ny+1))
    allocate(s%alfau(par%nx+1, par%ny+1))
    allocate(s%dsv(par%nx+1, par%ny+1))
    allocate(s%dnv(par%nx+1, par%ny+1))
    allocate(s%dsdnvi(par%nx+1, par%ny+1))
    allocate(s%alfav(par%nx+1, par%ny+1))
    allocate(s%dsc(par%nx+1, par%ny+1))
    allocate(s%dnc(par%nx+1, par%ny+1))

    allocate(xc(par%nx+1, par%ny+1))
    allocate(yc(par%nx+1, par%ny+1))

    s%xu = 0.d0
    s%yu = 0.d0
    s%xv = 0.d0
    s%yv = 0.d0
    s%dsz = 0.d0
    s%dnz = 0.d0
    s%dsdnzi = 0.d0
    s%alfaz = 0.d0
    s%dsu = 0.d0
    s%dnu = 0.d0
    s%dsdnui = 0.d0
    s%alfau = 0.d0
    s%dsv = 0.d0
    s%dnv = 0.d0
    s%dsdnvi = 0.d0
    s%alfav = 0.d0
    s%dsc = 0.d0
    s%dnc = 0.d0

    xc = 0.d0
    yc = 0.d0

    ! world coordinates of u-points
    do i=1,par%ny+1
       do j=1,par%nx
          s%xu(j,i)=.5d0*(s%xz(j,i)+s%xz(j+1,i))
          s%yu(j,i)=.5d0*(s%yz(j,i)+s%yz(j+1,i))
       enddo
       s%xu(par%nx+1,i)=1.5d0*s%xz(par%nx+1,i)-0.5d0*s%xz(par%nx,i)
       s%yu(par%nx+1,i)=1.5d0*s%yz(par%nx+1,i)-0.5d0*s%yz(par%nx,i)
    enddo

    ! world coordinates of v-points
    if (par%ny>0) then
       do j=1,par%nx+1
          do i=1,par%ny
             s%xv(j,i)=.5d0*(s%xz(j,i)+s%xz(j,i+1))
             s%yv(j,i)=.5d0*(s%yz(j,i)+s%yz(j,i+1))
          enddo
          s%xv(j,par%ny+1)=1.5d0*s%xz(j,par%ny+1)-0.5d0*s%xz(j,par%ny)
          s%yv(j,par%ny+1)=1.5d0*s%yz(j,par%ny+1)-0.5d0*s%yz(j,par%ny)
       enddo
    else
       s%xv=s%xz
       s%yv=s%yz
    endif

    ! world coordinates of corner points
    if (par%ny>0) then
       do i=1,par%ny
          do j=1,par%nx
             xc(j,i)=.25d0*(s%xz(j,i)+s%xz(j,i+1)+s%xz(j+1,i)+s%xz(j+1,i+1))
             yc(j,i)=.25d0*(s%yz(j,i)+s%yz(j,i+1)+s%yz(j+1,i)+s%yz(j+1,i+1))
          enddo
          xc(par%nx+1,i)=0.5d0*(s%xu(par%nx+1,i)+s%xu(par%nx+1,i+1))
          yc(par%nx+1,i)=0.5d0*(s%yu(par%nx+1,i)+s%yu(par%nx+1,i+1))
       enddo
       do j=1,par%nx
          xc(j,par%ny+1)=0.5d0*(s%xv(j,par%ny+1)+s%xv(j+1,par%ny+1))
          yc(j,par%ny+1)=0.5d0*(s%yv(j,par%ny+1)+s%yv(j+1,par%ny+1))
       enddo
       xc(par%nx+1,par%ny+1)=1.5d0*s%xu(par%nx+1,par%ny+1)-0.5*s%xu(par%nx+1,par%ny)
       yc(par%nx+1,par%ny+1)=1.5d0*s%yu(par%nx+1,par%ny+1)-0.5*s%yu(par%nx+1,par%ny)
    else
       xc=s%xu
       yc=s%yu
    endif

    ! s%dsu
    do i=1,par%ny+1
       do j=1,par%nx
          s%dsu(j,i)=((s%xz(j+1,i)-s%xz(j,i))**2+(s%yz(j+1,i)-s%yz(j,i))**2)**(0.5d0)
       enddo
    enddo
    s%dsu(par%nx+1,:)=s%dsu(par%nx,:)

    ! s%dsz
    do i=1,par%ny+1
       do j=2,par%nx+1
          s%dsz(j,i)=((s%xu(j,i)-s%xu(j-1,i))**2+(s%yu(j,i)-s%yu(j-1,i))**2)**(0.5d0)
       enddo
    enddo
    s%dsz(1,:)=s%dsz(2,:)

    ! s%dsv
    if (par%ny>0) then
       do i=1,par%ny+1
          do j=2,par%nx+1
             s%dsv(j,i)=((xc(j,i)-xc(j-1,i))**2+(yc(j,i)-yc(j-1,i))**2)**(0.5d0)
          enddo
       enddo
       s%dsv(par%nx+1,:)=s%dsv(par%nx,:)
    else
       s%dsv=s%dsz
    endif

    ! s%dsc
    if (par%ny>0) then
       do i=1,par%ny+1
          do j=1,par%nx
             s%dsc(j,i)=((s%xv(j+1,i)-s%xv(j,i))**2+(s%yv(j+1,i)-s%yv(j,i))**2)**(0.5d0)
          enddo
       enddo
       s%dsc(par%nx+1,:)=s%dsc(par%nx,:)
    else
       s%dsc=s%dsu
    endif

    ! s%dnu
    if (par%ny>0) then
       do i=2,par%ny+1
          do j=1,par%nx+1
             s%dnu(j,i)=((xc(j,i)-xc(j,i-1))**2+(yc(j,i)-yc(j,i-1))**2)**(0.5d0)
          enddo
       enddo
       s%dnu(:,1)=s%dnu(:,2)
    else
       s%dnu=100.d0
    endif

    ! s%dnz
    if (par%ny>0) then
       do i=2,par%ny+1
          do j=1,par%nx+1
             s%dnz(j,i)=((s%xv(j,i)-s%xv(j,i-1))**2+(s%yv(j,i)-s%yv(j,i-1))**2)**(0.5d0)
          enddo
       enddo
       s%dnz(:,1)=s%dnz(:,2)
    else
       s%dnz=100.d0
    endif

    ! s%dnv
    if (par%ny>0) then  
       do i=1,par%ny
          do j=1,par%nx+1
             s%dnv(j,i)=((s%xz(j,i+1)-s%xz(j,i))**2+(s%yz(j,i+1)-s%yz(j,i))**2)**(0.5d0)
          enddo
       enddo
       s%dnv(:,par%ny+1)=s%dnv(:,par%ny)
    else
       s%dnv=100.d0
    endif

    ! s%dnc
    if (par%ny>0) then  
       do i=1,par%ny
          do j=1,par%nx+1
             s%dnc(j,i)=((s%xu(j,i+1)-s%xu(j,i))**2+(s%yu(j,i+1)-s%yu(j,i))**2)**(0.5d0)
          enddo
       enddo
       s%dnc(:,par%ny+1)=s%dnc(:,par%ny)
    else
       s%dnc=100.d0
    endif


    if (par%ny>0) then 

       ! dsdnu
       do i=2,par%ny+1
          do j=1,par%nx
             x1=s%xv(j  ,i  ) - s%xv(j  ,i-1)
             x3=s%xv(j+1,i-1) - s%xv(j  ,i-1)
             x2=s%xv(j+1,i  ) - s%xv(j+1,i-1)
             x4=s%xv(j+1,i  ) - s%xv(j  ,i  )
             y1=s%yv(j  ,i  ) - s%yv(j  ,i-1)
             y3=s%yv(j+1,i-1) - s%yv(j  ,i-1)
             y2=s%yv(j+1,i  ) - s%yv(j+1,i-1)
             y4=s%yv(j+1,i  ) - s%yv(j  ,i  )
             dsdnu=0.5d0*(abs(x1*y3-x3*y1)+abs(x2*y4-x4*y2))
             s%dsdnui(j,i)=1.d0/dsdnu
          enddo
       enddo
       s%dsdnui(:,1)=s%dsdnui(:,2)
       s%dsdnui(par%nx+1,:)=s%dsdnui(par%nx,:)

       ! dsdnv
       do i=1,par%ny
          do j=2,par%nx+1
             x1=s%xu(j-1,i+1) - s%xu(j-1,i  )
             x3=s%xu(j  ,i  ) - s%xu(j-1,i  )
             x2=s%xu(j  ,i+1) - s%xu(j  ,i  )
             x4=s%xu(j  ,i+1) - s%xu(j-1,i+1)
             y1=s%yu(j-1,i+1) - s%yu(j-1,i  )
             y3=s%yu(j  ,i  ) - s%yu(j-1,i  )
             y2=s%yu(j  ,i+1) - s%yu(j  ,i  )
             y4=s%yu(j  ,i+1) - s%yu(j-1,i+1)
             dsdnv=0.5d0*(abs(x1*y3-x3*y1)+abs(x2*y4-x4*y2))
             s%dsdnvi(j,i)=1.d0/dsdnv
          enddo
       enddo
       s%dsdnvi(:,par%ny+1)=s%dsdnvi(:,par%ny)
       s%dsdnvi(1,:)=s%dsdnvi(2,:)

       ! dsdnz
       do i=2,par%ny+1
          do j=2,par%nx+1
             x1=xc(j-1,i  ) - xc(j-1,i-1)
             x3=xc(j  ,i-1) - xc(j-1,i-1)
             x2=xc(j  ,i  ) - xc(j  ,i-1)
             x4=xc(j  ,i  ) - xc(j-1,i  )
             y1=yc(j-1,i  ) - yc(j-1,i-1)
             y3=yc(j  ,i-1) - yc(j-1,i-1)
             y2=yc(j  ,i  ) - yc(j  ,i-1)
             y4=yc(j  ,i  ) - yc(j-1,i  )
             dsdnz=0.5d0*(abs(x1*y3-x3*y1)+abs(x2*y4-x4*y2))
             s%dsdnzi(j,i)=1.d0/dsdnz
          enddo
       enddo
       s%dsdnzi(:,1)=s%dsdnzi(:,2)
       s%dsdnzi(1,:)=s%dsdnzi(2,:)

    else

       s%dsdnui=1.d0/(s%dsu*s%dnu)
       s%dsdnvi=1.d0/(s%dsv*s%dnv)
       s%dsdnzi=1.d0/(s%dsz*s%dnz)

    endif

    ! s%alfaz, grid orientation in z-points
    do i=1,par%ny+1
       do j=2,par%nx
          s%alfaz(j,i)=atan2(s%yz(j+1,i)-s%yz(j-1,i),s%xz(j+1,i)-s%xz(j-1,i))
       enddo
       s%alfaz(1,i)=s%alfaz(2,i)
       s%alfaz(par%nx+1,i)=s%alfaz(par%nx,i)
    enddo

    ! s%alfau, grid orientation in u-points
    do i=1,par%ny+1
       do j=1,par%nx
          s%alfau(j,i)=atan2(s%yz(j+1,i)-s%yz(j,i),s%xz(j+1,i)-s%xz(j,i))
       enddo
       s%alfau(par%nx+1,i)=s%alfau(par%nx,i)
    enddo

    ! s%alfav, grid orientation in v-points
    if (par%ny>0) then
       do j=1,par%nx+1
          do i=1,par%ny
             s%alfav(j,i)=atan2(s%yz(j,i+1)-s%yz(j,i),s%xz(j,i+1)-s%xz(j,i))
          enddo
          s%alfav(j,par%ny+1)=s%alfav(j,par%ny)
       enddo
    else
       s%alfav=s%alfaz
    endif

    deallocate (xc)
    deallocate (yc)

  end subroutine gridprops

end module init_module

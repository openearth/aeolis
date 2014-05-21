module bed_module

  use constants_module
  use utils_module
  use input_module

  use precision
  use bedcomposition_module
  use message_module

  implicit none

contains

  subroutine generate_bed(par, x, z)

    type(parameters), intent(inout) :: par
    integer*4 :: fid, ierr, n, i
    real*8, dimension(:), allocatable, intent(out) :: x, z
    real*8, dimension(:), allocatable :: x_c, z_c

    fid = 99

    ! count lines
    n = -1
    ierr = 0
    open(fid, file=trim(par%bed_file))
    do while (ierr == 0)
       read(fid, *, iostat=ierr)
       n = n + 1
    end do
    rewind(fid)

    ! allocate arrays
    allocate(x_c(n))
    allocate(z_c(n))

    ! read data
    i = 1
    ierr = 0
    do while (ierr == 0)
       read(fid, *, iostat=ierr) x_c(i), z_c(i)
       i = i + 1
    end do
    close(fid)
        
    ! allocate arrays
    par%nx = nint(maxval(x_c) / par%dx)
    allocate(x(par%nx+1))
    allocate(z(par%nx+1))

    ! interpolate bed
    x = (/(i*par%dx, i=0, par%nx)/)
    do i=1, par%nx+1
       z(i) = linear_interp(x_c, z_c, x(i))
    end do
    
  end subroutine generate_bed

  subroutine generate_bedcomposition(par, x, z)

    type(parameters), intent(in) :: par
    real*8, dimension(:), intent(in) :: x, z

    integer :: iostat
    type(bedcomp_data), pointer :: morlyr
    type(message_stack), pointer :: messages
    logical, pointer :: exchlyr
    integer, pointer :: nmlb
    integer, pointer :: nmub
    integer, pointer :: flufflyr
    integer , pointer :: idiffusion
    integer , pointer :: iporosity
    integer , pointer :: iunderlyr
    integer , pointer :: maxwarn
    integer , pointer :: ndiff
    integer , pointer :: neulyr
    integer , pointer :: nfrac
    integer , pointer :: nlalyr
    integer , pointer :: nlyr
    integer , pointer :: updbaselyr
    real(fp) , pointer :: minmass
    real(fp) , pointer :: theulyr
    real(fp) , pointer :: thlalyr
    real(fp) , dimension(:) , pointer :: dpsed
    real(fp) , dimension(:) , pointer :: thtrlyr
    real(fp) , dimension(:) , pointer :: zdiff
    real(fp) , dimension(:,:) , pointer :: kdiff
    real(fp) , dimension(:,:) , pointer :: mfluff
    real(fp) , dimension(:,:) , pointer :: svfrac
    real(fp) , dimension(:,:) , pointer :: thlyr
    real(fp) , dimension(:,:,:) , pointer :: msed
    real(prec) , dimension(:,:) , pointer :: bodsed

    character(message_len) :: message
    integer :: i
    integer :: l
    integer :: nm
    integer :: nstep
    integer :: istat
    integer , dimension(:) , allocatable :: sedtyp
    real(fp) :: t
    real(fp) :: dt
    real(fp) :: g
    real(fp) :: morfac
    real(fp) :: rhow
    real(fp) :: tend
    real(fp) :: tstart
    real(fp) , dimension(:) , allocatable :: cdryb
    real(fp) , dimension(:) , allocatable :: chezy
    real(fp) , dimension(:) , allocatable :: dp
    real(fp) , dimension(:) , allocatable :: dz
    real(fp) , dimension(:) , allocatable :: h0
    real(fp) , dimension(:) , allocatable :: h1
    real(fp) , dimension(:) , allocatable :: logsedsig
    real(fp) , dimension(:) , allocatable :: rhosol
    real(fp) , dimension(:) , allocatable :: sedd50
    real(fp) , dimension(:) , allocatable :: sedd90
    real(fp) , dimension(:) , allocatable :: taub
    real(fp) , dimension(:) , allocatable :: umod
    real(fp) , dimension(:,:), allocatable :: mass
    real(fp) , dimension(:,:), allocatable :: massfluff
    real(fp) , dimension(:,:), allocatable :: r0
    real(fp) , dimension(:,:), allocatable :: r1
    real(fp) , dimension(:,:), allocatable :: rn
    real(fp) , dimension(:,:), allocatable :: sink
    real(fp) , dimension(:,:), allocatable :: sinkf
    real(fp) , dimension(:,:), allocatable :: sour
    real(fp) , dimension(:,:), allocatable :: sourf
    real(fp) , dimension(:,:), allocatable :: ws

    allocate (morlyr)
    allocate (messages)
    call initstack(messages)

    message = 'initializing bed composition module'
    call addmessage(messages, message)
    if (initmorlyr(morlyr) /= 0 ) call adderror(messages, message)

    iostat = 0
    message = 'initializing logical and scalar values'
    call addmessage(messages, message)
    istat = bedcomp_getpointer_integer(morlyr, 'flufflayer_model_type' , flufflyr)
    istat = bedcomp_getpointer_integer(morlyr, 'diffusion_model_type' , idiffusion)
    istat = bedcomp_getpointer_integer(morlyr, 'porosity_model_type' , iporosity)
    istat = bedcomp_getpointer_integer(morlyr, 'bed_layering_type' , iunderlyr)
    istat = bedcomp_getpointer_integer(morlyr, 'number_of_fractions' , nfrac)
    istat = bedcomp_getpointer_integer(morlyr, 'first_column_number' , nmlb)
    istat = bedcomp_getpointer_integer(morlyr, 'last_column_number' , nmub)
    istat = bedcomp_getpointer_integer(morlyr, 'MaxNumShortWarning' , maxwarn)
    istat = bedcomp_getpointer_integer(morlyr, 'number_of_diffusion_values' , ndiff)
    istat = bedcomp_getpointer_integer(morlyr, 'number_of_eulerian_layers' , neulyr)
    istat = bedcomp_getpointer_integer(morlyr, 'number_of_lagrangian_layers' , nlalyr)
    istat = bedcomp_getpointer_integer(morlyr, 'base_layer_updating_type' , updbaselyr)
    istat = bedcomp_getpointer_realfp (morlyr, 'MinMassShortWarning' , minmass)
    istat = bedcomp_getpointer_realfp (morlyr, 'thickness_of_eulerian_layers' , theulyr)
    istat = bedcomp_getpointer_realfp (morlyr, 'thickness_of_lagrangian_layers' , thlalyr)
    if (istat /= 0) call adderror(messages, message)

    nmlb    = 1                
    nmub    = size(x)                
    tstart  = 0.0              
    tend    = par%tstop           
    dt      = par%dt            
    morfac  = 1.0              
    nstep  = (tend-tstart)/dt; 

    nfrac       = 2
    iunderlyr   = 2            
    neulyr      = 3            
    nlalyr      = 0            
    theulyr     = 0.1_fp       
    thlalyr     = 0.2_fp       
    updbaselyr  = 1            
                               
    maxwarn     = 100          
    minmass     = 0.0_fp       
    idiffusion  = 0            
    ndiff       = 5            
    flufflyr    = 0            
    iporosity   = 0            
                               
    message = 'allocating bed composition module'
    call addmessage(messages, message)
    if (allocmorlyr(morlyr) /=0 ) call adderror(messages, message)

    message = 'initializing bed compositon arrays'     
    call addmessage(messages, message)
    istat = bedcomp_getpointer_integer(morlyr, 'number_of_layers' , nlyr)
    istat = bedcomp_getpointer_realfp (morlyr, 'layer_mass' , msed)
    istat = bedcomp_getpointer_realfp (morlyr, 'layer_thickness' , thlyr)
    istat = bedcomp_getpointer_realfp (morlyr, 'thickness_of_transport_layer' , thtrlyr)
    istat = bedcomp_getpointer_realfp (morlyr, 'solid_volume_fraction' , svfrac)
    if (istat/=0) call adderror(messages, message)

    allocate(cdryb(nfrac))
    allocate(logsedsig(nfrac))
    allocate(rhosol(nfrac))
    allocate(sedd50(nfrac))
    allocate(sedd90(nfrac))
    allocate(sedtyp(nfrac))
    
    allocate(chezy(nmlb:nmub))
    allocate(h0(nmlb:nmub))
    allocate(h1(nmlb:nmub))
    allocate(umod(nmlb:nmub))
    allocate(taub(nmlb:nmub))
    allocate(r0(nfrac,nmlb:nmub))
    allocate(r1(nfrac,nmlb:nmub))
    allocate(rn(nfrac,nmlb:nmub))
    allocate(ws(nfrac,nmlb:nmub))
    
    allocate(mass(nfrac,nmlb:nmub))
    allocate(massfluff(nfrac,nmlb:nmub))
    allocate(sink(nfrac,nmlb:nmub))
    allocate(sinkf(nfrac,nmlb:nmub))
    allocate(sour(nfrac,nmlb:nmub))
    allocate(sourf(nfrac,nmlb:nmub))
    allocate(dp(nmlb:nmub))
    allocate(dz(nmlb:nmub))

    sedtyp(1)   = SEDTYP_NONCOHESIVE_SUSPENDED 
    sedtyp(2)   = SEDTYP_COHESIVE              
    cdryb       = 1650.0_fp                    
    rhosol      = 2650.0_fp                    
    sedd50      = 0.0001_fp                    
    sedd90      = 0.0002_fp                    
    logsedsig   = log(1.34_fp)                 

    thtrlyr = 0.1_fp       
    thlyr   = 0.1_fp       
    svfrac  = 1.0_fp       
    msed = 0.0_fp          
    do l = 1, nfrac
       msed(l,:,:) = thlyr*cdryb(l)/nfrac 
    enddo

    call setbedfracprop(morlyr, sedtyp, sedd50, logsedsig, cdryb)

  end subroutine generate_bedcomposition

  subroutine compute_threshold_bedslope(par, x, z, u_th)

    type(parameters), intent(in) :: par
    real*8, dimension(:), intent(in) :: x, z
    real*8, dimension(:), intent(inout) :: u_th
    integer :: i
    real*8 :: phi, theta

    phi = par%phi / 180.d0 * pi

    do i=1,par%nx
       theta = -atan((z(i+1) - z(i)) / (x(i+1) - x(i)))
       u_th(i) = sqrt((tan(phi) - tan(theta) / tan(phi) + cos(theta)) / \
                 (tan(phi)+1)) * u_th(i)
    end do

  end subroutine compute_threshold_bedslope

end module bed_module

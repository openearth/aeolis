module bed_module

  use constants_module
  use utils_module
  use input_module
  use moist_module

  use precision
  use bedcomposition_module
  use message_module

  implicit none

  type(bedcomp_data), pointer :: morlyr
  type(message_stack), pointer :: messages
  character(message_len) :: message

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

    type(parameters), intent(inout) :: par
    real*8, dimension(:), intent(in) :: x, z

    integer :: iostat
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
    integer , pointer :: keuler
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

    integer :: i, j
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
    istat = bedcomp_getpointer_integer(morlyr, 'keuler' , keuler)
    istat = bedcomp_getpointer_realfp (morlyr, 'MinMassShortWarning' , minmass)
    istat = bedcomp_getpointer_realfp (morlyr, 'thickness_of_eulerian_layers' , theulyr)
    istat = bedcomp_getpointer_realfp (morlyr, 'thickness_of_lagrangian_layers' , thlalyr)
    if (istat /= 0) call adderror(messages, message)

    nmlb    = 1                
    nmub    = size(x)                
    tstart  = 0.0              
    tend    = par%tstop
    dt      = par%dt
    morfac  = 1.d0
    nstep   = (tend-tstart)/dt; 

    nfrac       = par%nfractions
    iunderlyr   = 2 ! graded sediments
    neulyr      = 1 !par%nlayers ! all layers are fixed
    nlalyr      = par%nlayers ! except one
    theulyr     = 10 !par%layer_thickness
    thlalyr     = par%layer_thickness
    updbaselyr  = 1 ! base layer is independent
                               
    maxwarn     = 100
    minmass     = 0.5_fp       
    idiffusion  = 0 ! no diffusion between layers
    ndiff       = 5
    flufflyr    = 0 ! no fluff layers            
    iporosity   = 0 ! porosity included in densities
                               
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
    cdryb       = par%rhom
    rhosol      = par%rhop
    sedd50      = 0.0001_fp                    
    sedd90      = 0.0002_fp                    
    logsedsig   = log(1.34_fp)                 
    svfrac  = 1.0_fp

    ! set layer thickness
    thtrlyr = thlalyr
    do j = 1,nlalyr+1
       thlyr(j,:) = thlalyr
    end do
    do j = nlalyr+2,nlyr
       thlyr(j,:) = theulyr
    end do

    ! set layer mass
    msed = 0.0_fp
    do j = 1,nlyr
       do i = 1,par%nfractions
          msed(i,j,:) = thlyr(j,:) * cdryb(i) * par%grain_dist(i) / &
               max(1e-10, sum(par%grain_dist))
       end do
    end do

    ! update parameters
    par%nlayers = nlyr

    call setbedfracprop(morlyr, sedtyp, sedd50, logsedsig, cdryb)

  end subroutine generate_bedcomposition

  function update_bed(par, z, mass, rho) result (z_new)

    type(parameters), intent(in) :: par
    real*8, dimension(:), intent(in) :: z
    real*8, dimension(:,:), intent(in) :: mass
    real*8, dimension(:), intent(in) :: rho
    real*8, dimension(:), allocatable :: z_new, dz

    allocate(z_new(size(z)))
    allocate(dz(size(z)))

    if ( updmorlyr(morlyr, mass, 0.d0 * mass, rho, par%dt, 1.d0, dz, messages) /= 0 ) then
       call adderror(messages, message)
    end if

    if (par%bedupdate) then
       z_new = z + dz
    else
       z_new = z
    end if
    
  end function update_bed

  function get_layer_mass() result (mass)

    integer :: istat
    real*8 , dimension(:,:,:), pointer :: mass

    istat = bedcomp_getpointer_realfp(morlyr, 'layer_mass', mass)
    if (istat/=0) call adderror(messages, message)

  end function get_layer_mass

  function get_layer_thickness() result (thlyr)

    integer :: istat
    real*8 , dimension(:,:), pointer :: thlyr

    istat = bedcomp_getpointer_realfp(morlyr, 'layer_thickness', thlyr)
    if (istat/=0) call adderror(messages, message)

  end function get_layer_thickness

  function get_layer_percentile(par, p) result (perc)

    type(parameters), intent(in) :: par
    integer :: istat, i, j, k
    real*8 :: sedtot
    real*8, intent(in) :: p
    real*8, dimension(:), allocatable :: frac
    real*8, dimension(:,:), allocatable :: perc
    real*8, dimension(:,:,:), allocatable :: msed_cs
    real*8, dimension(:,:,:), pointer :: msed

    allocate(frac(par%nfractions))
    allocate(perc(par%nlayers, par%nx+1))
    allocate(msed_cs(par%nfractions, par%nlayers, par%nx+1))

    istat = bedcomp_getpointer_realfp(morlyr, 'layer_mass', msed)
    if (istat/=0) call adderror(messages, message)

    do i = 1,par%nx+1
       do j = 1,par%nlayers
          sedtot = 0.0_fp
          do k = 1,par%nfractions
             sedtot = sedtot + msed(k,j,i)
             msed_cs(k,j,i) = sedtot
          enddo
          perc(j,i) = 10**linear_interp(msed_cs(:,j,i)/sedtot, log10(par%grain_size), p)
       enddo
    enddo

  end function get_layer_percentile

  function get_layer_massfraction(par) result (fractions)

    type(parameters), intent(in) :: par
    integer :: istat, i, j, k
    real*8 :: sedtot
    real*8, dimension(:,:,:), allocatable :: fractions
    real*8, dimension(:,:,:), pointer :: msed

    allocate(fractions(par%nfractions, par%nlayers, par%nx+1))

    istat = bedcomp_getpointer_realfp(morlyr, 'layer_mass', msed)
    if (istat/=0) call adderror(messages, message)

    do i = 1,par%nx+1
       do j = 1,par%nlayers
          sedtot = 0.0_fp
          do k = 1,par%nfractions
             sedtot = sedtot + msed(k,j,i)
          enddo
          do k = 1,par%nfractions
             fractions(k,j,i) = msed(k,j,i)/sedtot
          enddo
       enddo
    enddo

  end function get_layer_massfraction

  function get_layer_volumefraction(par) result (fractions)

    type(parameters), intent(in) :: par
    integer :: istat, i, j, k
    real*8, dimension(:,:,:), allocatable :: fractions
    real*8, dimension(:,:,:), pointer :: msed
    real*8, dimension(:,:), pointer :: svfrac
    real*8, dimension(:,:), pointer :: thlyr
    real*8, dimension(:), pointer :: dens

    allocate(fractions(par%nfractions, par%nlayers, par%nx+1))

    istat = bedcomp_getpointer_realfp(morlyr, 'layer_thickness', thlyr)
    istat = bedcomp_getpointer_realfp(morlyr, 'layer_mass', msed)
    istat = bedcomp_getpointer_realfp(morlyr, 'sediment_density', dens)
    istat = bedcomp_getpointer_realfp(morlyr, 'solid_volume_fraction', svfrac)
    if (istat/=0) call adderror(messages, message)

    do i = 1,par%nx+1
       do j = 1,par%nlayers
          do k = 1,par%nfractions+2
             fractions(k,j,i) = msed(k,j,i)/(dens(k)*svfrac(j,i)*thlyr(j,i))
          enddo
       enddo
    enddo

  end function get_layer_volumefraction

  subroutine compute_threshold_grainsize(par, u_th)

    type(parameters), intent(in) :: par
    real*8, dimension(:,:), intent(inout) :: u_th
    integer :: i

    if (.not. par%th_grainsize) return
    
    ! Bagnold

    do i=1,par%nfractions
       u_th(i,:) = par%A * sqrt(((par%rhop - par%rhoa) * &
                   par%g * par%grain_size(i)) / par%rhop)
    end do

  end subroutine compute_threshold_grainsize

  subroutine compute_threshold_bedslope(par, x, z, u_th)

    type(parameters), intent(in) :: par
    real*8, dimension(:), intent(in) :: x, z
    real*8, dimension(:,:), intent(inout) :: u_th
    integer :: i
    real*8 :: phi, theta

    if (.not. par%th_bedslope) return

    ! Dyer, 1986

    phi = par%phi / 180.d0 * pi

    do i=1,par%nx
       theta = -atan((z(i+1) - z(i)) / (x(i+1) - x(i)))
       u_th(:,i) = sqrt((tan(phi) - tan(theta)) / tan(phi) * cos(theta)) * u_th(:,i)
    end do

    u_th(:,par%nx+1) = u_th(:,par%nx)

  end subroutine compute_threshold_bedslope

  subroutine mix_toplayer(par, z, tide)

    type(parameters), intent(in) :: par
    real*8, dimension(:), intent(in) :: z
    real*8, intent(in) :: tide
    integer :: i, j, k, nmix
    real*8 :: th

    integer :: istat
    real(fp), dimension(:,:,:), pointer :: msed
    real(fp), dimension(:,:), pointer :: thlyr

    real*8, dimension(:,:), allocatable :: dist

    if (.not. par%mixtoplayer) return
       
    istat = bedcomp_getpointer_realfp(morlyr, 'layer_mass' , msed)
    istat = bedcomp_getpointer_realfp(morlyr, 'layer_thickness' , thlyr)
    if (istat/=0) call adderror(messages, message)

    allocate(dist(par%nfractions, par%nx+1))

    do k = 1,par%nx+1

       ! fix only if flooded
       if (z(k) <= tide) then

          ! determine mixing depth in terms of number of layers
          nmix = 0
          th = 0.d0
          do i = 1,size(thlyr(:,k))
             nmix = nmix + 1
             th = th + thlyr(i,k)
             if (th >= min(par%Hs, (tide - z(k)) * par%gamma) * par%facDOD) exit
          end do

          if (nmix == 0) continue

          ! determine average sediment distribution over mixing depth
          do i = 1,par%nfractions
             dist(i,k) = sum(msed(i,1:nmix,k))
          end do
          dist(:,k) = dist(:,k) / max(1e-10, sum(dist(:,k)))

          ! mix layer mass
          do j = 1,nmix
             do i = 1,par%nfractions
                msed(i,j,k) = thlyr(j,k) * par%rhom * dist(i,k)
             end do
          end do
       end if
       
    end do
        
  end subroutine mix_toplayer

  subroutine sweep_toplayer(par)

    type(parameters), intent(in) :: par
    integer :: i, k, n, m

    integer :: istat
    real(fp), dimension(:,:,:), pointer :: msed

    real*8, dimension(:), allocatable :: dist, mass

    if (.not. par%sweeptoplayer) return

    istat = bedcomp_getpointer_realfp(morlyr, 'layer_mass' , msed)
    if (istat/=0) call adderror(messages, message)

    n = par%nfractions
    m = par%nlayers
    
    allocate(dist(n))
    allocate(mass(n))

    do k = 1,par%nx+1

       ! determine sediment distribution in top layer
       dist = msed(:,1,k) / max(1e-10, sum(msed(:,1,k)))
       
       do i = 1,n
          if (dist(i) < par%minfrac) then ! if fraction is below threshold

             ! compute exchange of sediment between top two layers
             mass = 0.d0
             mass(i) = -msed(i,1,k)
             mass(i+1:n) = msed(i,1,k) * sum(msed(i+1:n,2:m,k), dim=2) / &
                  max(1e-10, sum(msed(i+1:n,2:m,k)))

             if (abs(sum(mass)) < 1e-10) then

                ! swap sediment
                msed(:,1,k) = msed(:,1,k) + mass
                msed(:,2,k) = msed(:,2,k) - mass

             end if
             
          end if
       end do
    end do

  end subroutine sweep_toplayer
  
end module bed_module

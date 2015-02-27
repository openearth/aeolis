module moist_module

  use constants_module
  use utils_module
  use input_module

  implicit none

  type meteorology
     real*8 :: solar_radiation = 1e4 ! [J/m2]
     real*8 :: air_temperature = 10.d0 ! [oC]
     real*8 :: relative_humidity = 0.4d0 ! [-]
     real*8 :: air_specific_heat = 1.0035e-3 ! [MJ/kg/K]
     real*8 :: atmospheric_pressure = 101.325 ! [kPa]
     real*8 :: latent_heat = 2.45 ! [MJ/kg]
  end type meteorology

contains

  subroutine generate_moist(par, z, m)

    type(parameters), intent(inout) :: par
    integer*4 :: fid, ierr, n, nz, i, j, nt, file_pos
    real*8, dimension(:,:), allocatable, intent(out) :: m
    real*8, dimension(:), allocatable, intent(out) :: z
    real*8, dimension(:,:), allocatable :: m0
    real*8, dimension(:), allocatable :: t
    character(2) :: fmt

    ! allocate arrays
    nt = nint(par%tstop / par%dt)
    allocate(m(nt,nz+2))

    if (trim(par%moist_file) /= '') then
       
       fid = 77

       ! read format
       open(fid, file=trim(par%moist_file))
       read(fid, *, iostat=ierr) fmt
       if (fmt .ne. 'TZ') then
          write(*,*) "ERROR: Unknown moist format: ", fmt
          stop 1
       end if
       
       ! read z values
       nz = 0
       read(fid, *, iostat=ierr) nz
       allocate(z(nz))
       do i = 1,nz
          read(fid, *, iostat=ierr) z(i)
       end do
       
       ! count lines
       n = -1
       ierr = 0
       do while (ierr == 0)
          read(fid, *, iostat=ierr)
          n = n + 1
       end do
       rewind(fid)
       
       ! allocate arrays
       allocate(t(n))
       allocate(m0(n,nz+2))
       
       ! read data
       ierr = 0
       do i = 1,nz+2
          read(fid, *, iostat=ierr)
       end do
       
       i = 1
       ierr = 0
       do while (ierr == 0)
          read(fid, *, iostat=ierr) t(i), m0(i,:)
          i = i + 1
       end do
       close(fid)

       ! checks
       if (maxval(t) < par%tstop) then
          write(*,*) "ERROR: moist definition file too short"
          stop 1
       end if
       
       ! interpolate time series
       do j = 1,nz+2
          do i = 1,nt
             m(i,j) = linear_interp(t, m0(:,j), i*par%dt)
          end do
       end do

    else

       m = 0.d0

    end if
    
  end subroutine generate_moist

  subroutine generate_tide(par, tide)

    type(parameters), intent(inout) :: par
    integer*4 :: fid, ierr, n, i, nt
    real*8, dimension(:), allocatable, intent(out) :: tide
    real*8, dimension(:), allocatable :: t, tmp

    ! allocate arrays
    nt = nint(par%tstop / par%dt)
    allocate(tide(nt))

    if (trim(par%tide_file) /= '') then

       fid = 99

       ! count lines
       n = -1
       ierr = 0
       open(fid, file=trim(par%tide_file))
       do while (ierr == 0)
          read(fid, *, iostat=ierr)
          n = n + 1
       end do
       rewind(fid)
       
       ! allocate arrays
       allocate(t(n))
       allocate(tmp(n))
       
       ! read data
       i = 1
       ierr = 0
       do while (ierr == 0)
          read(fid, *, iostat=ierr) t(i), tmp(i)
          i = i + 1
       end do
       close(fid)
       
       ! checks
       if (t(n) < par%tstop) then
          write(*,*) "ERROR: tide definition file too short"
          stop 1
       end if

       ! interpolate time series
       do i = 1,nt
          tide(i) = linear_interp(t, tmp, i*par%dt)
       end do

    else

       tide = 0.d0

    end if
    
  end subroutine generate_tide

  subroutine generate_meteo(par, meteo)

    type(parameters), intent(in) :: par
    type(meteorology), dimension(:), allocatable, intent(out) :: meteo

    allocate(meteo(1))

  end subroutine generate_meteo
  
  subroutine compute_threshold_moisture(par, moist, u_th)

    type(parameters), intent(in) :: par
    real*8, dimension(:), intent(in) :: moist
    real*8, dimension(par%nx+1) :: mg
    real*8, dimension(:,:), intent(inout) :: u_th
    real*8, dimension(par%nfractions, par%nx+1) :: u_th_m
    integer :: i, n

!    mg = map_moisture(par, zm, m, z)

    ! convert from volumetric content (percentage of volume) to
    ! geotechnical mass content (percentage of dry mass)
    mg = moist * par%rhow / (par%rhop * (1 - par%porosity))

    u_th_m = 0.d0
    do i = 1,par%nfractions
       if (par%method_moist .eq. 'belly_johnson') then
          u_th_m(i,:) = u_th(i,:) * max(1.d0, 1.8+0.6*log10(mg))
       else if (par%method_moist .eq. 'hotta') then
          u_th_m(i,:) = u_th(i,:) + 7.5 * mg
       else
          write(*, '(a,a)') "ERROR: Unknown moisture formulation: ", par%method_moist
          stop 1
       end if

       where (mg > 0.005)
          u_th(i,:) = u_th_m(i,:)
       end where

       where (mg > 0.064)
          ! should be .04 according to Pye and Tsoar
          ! should be .64 according to Delgado-Fernandez (10% vol.)
          u_th(i,:) = 999.d0
       end where
    end do

  end subroutine compute_threshold_moisture

  subroutine update_moisture(par, z, tide, meteo, u, moist)

    type(parameters), intent(in) :: par
    type(meteorology), intent(in) :: meteo
    real*8, dimension(:), intent(in) :: z
    real*8, dimension(:,:), intent(inout) :: moist
    real*8, intent(in) :: u, tide
    real*8 :: radiation, m, delta, gamma, evaporation
    integer :: i, t, nl, f

    nl = par%nlayers+2
    
    ! evaporation using Penman-Monteith
    radiation = meteo%solar_radiation / 1e6 / par%dt * 3600 * 24 ! conversion from J/m2 to MJ/m2/day
    m = vaporation_pressure_slope(meteo%air_temperature) ! [kPa/K]
    delta = saturation_pressure(meteo%air_temperature) * (1 - meteo%relative_humidity) ! [kPa]
    gamma = (meteo%air_specific_heat * meteo%atmospheric_pressure) / &
         (.622 * meteo%latent_heat) ! [kPa/K]
    evaporation = max(0.d0, (m * radiation + gamma * 6.43 * (1 + 0.536 * u) * delta) / &
         (meteo%latent_heat * (m + gamma)))
    evaporation = evaporation / 24 / 3600 / 1000 ! conversion from mm/day to m/s
    
    ! infiltration using Darcy
    do i = 1,par%nx
       if (tide >= z(i)) then
          moist(:,i) = par%porosity
       else
          moist(:,i) = moist(:,i) * exp(-par%F * par%dt)
!          moist(:,i) = max(0.d0, moist(:,i) - evaporation * par%dt / par%layer_thickness)
       end if
   end do

  end subroutine update_moisture
  
  function vaporation_pressure_slope(T) result (s)

    real*8, intent(in) :: T
    real*8 :: s

    ! Tetens, 1930; Murray, 1967
    s = 4098 * 0.6108 * exp((17.27 * T) / (T - 237.3)) / (T + 237.3)**2 ! [kPa/K]
    
  end function vaporation_pressure_slope

  function saturation_pressure(T) result (vp)

    real*8, intent(in) :: T
    real*8 :: TK, A, B, C, D, E, F, vp
    
    TK = T + 273.15
    A  = -1.88e4
    B  = -13.1
    C  = -1.5e-2
    D  =  8e-7
    E  = -1.69e-11
    F  =  6.456
    vp = exp(A/TK + B + C*TK + D*TK**2 + E*TK**3 + F*log(TK)) ! [kPa]

  end function saturation_pressure
    
  function map_moisture(par, zm, m, z) result (mg)
    
    type(parameters) :: par
    real*8, dimension(:) :: zm, m, z
    real*8, dimension(par%nx+1) :: mg
    integer :: i, n

    n = size(m)
    do i = 1,par%nx+1
       if (z(i) < minval(zm)) then
          mg(i) = m(1)
       elseif (z(i) > maxval(zm)) then
          mg(i) = m(n)
       else
          mg(i) = linear_interp(zm(2:n-1), m(2:n-1), z(i))
       end if
    end do

  end function map_moisture

end module moist_module

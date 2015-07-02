module moist_module

  use constants_module
  use utils_module
  use input_module

  implicit none

contains

  subroutine generate_moist(par, m)

    type(parameters), intent(inout) :: par
    type(moisture), dimension(:,:), allocatable, intent(out) :: m
    integer*4 :: fid, ierr, n, nz, i, j
    real*8 :: t, z
    character(2) :: fmt

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
       do i = 1,nz
          read(fid, *, iostat=ierr)
       end do
       
       ! count lines
       n = -1
       ierr = 0
       do while (ierr == 0)
          read(fid, *, iostat=ierr)
          n = n + 1
       end do
       rewind(fid)
       read(fid, *, iostat=ierr)
       read(fid, *, iostat=ierr)
       
       ! allocate arrays
       allocate(m(n,nz))
       
       ! read data
       ierr = 0
       do j = 1,nz
          read(fid, *, iostat=ierr) z
          m(:,j)%z = z
       end do
       
       i = 1
       ierr = 0
       do while (ierr == 0)
          read(fid, *, iostat=ierr) t, m(i,:)%moisture
          m(i,:)%t = t
          i = i + 1
       end do
       close(fid)

       ! checks
       if (maxval(m%t) < par%tstop) then
          write(*,*) "ERROR: moist definition file too short"
          stop 1
       end if
       
    else

       allocate(m(2,3))
       m(1,:)%t = 0.d0
       m(2,:)%t = par%tstop
       m%z = 0.d0
       m%moisture = 0.d0

    end if
    
  end subroutine generate_moist

  subroutine interpolate_moist(m, t, zb, moist)

    type(moisture), dimension(:,:), intent(in) :: m
    real*8, dimension(:), intent(in) :: zb
    real*8, intent(in) :: t
    real*8, dimension(:), allocatable :: m0
    real*8, dimension(:,:), intent(out) :: moist
    integer*4 :: j

    allocate(m0(size(m,2)))

    m0 = 0.d0
    do j = 1,size(m,2)
       m0(j) = linear_interp(m(:,j)%t, m(:,j)%moisture, t)
    end do

    moist = 0.d0
    do j = 1,size(moist,1)
       moist(j,:) = map_moisture(m(1,:)%z, m0, zb)
    end do

  end subroutine interpolate_moist

  subroutine generate_tide(par, zs)

    type(parameters), intent(inout) :: par
    type(tide), dimension(:), allocatable, intent(out) :: zs
    real*8, dimension(:), allocatable :: t, tmp
    integer*4 :: fid, ierr, n, i, nt

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
       allocate(zs(n))
       
       ! read data
       i = 1
       ierr = 0
       do while (ierr == 0)
          read(fid, *, iostat=ierr) zs(i)%t, zs(i)%level
          i = i + 1
       end do
       close(fid)
       
       ! checks
       if (maxval(zs%t) < par%tstop) then
          write(*,*) "ERROR: tide definition file too short"
          stop 1
       end if

    else

       ! allocate arrays
       allocate(zs(2))
       zs(1)%t = 0.d0
       zs(1)%level = 0.d0
       zs(2)%t = par%tstop
       zs(2)%level = 0.d0

    end if
    
  end subroutine generate_tide

  subroutine interpolate_tide(zs, t, field)

    type(tide), dimension(:), intent(in) :: zs
    real*8, intent(in) :: t
    real*8, dimension(:), intent(out) :: field
    real*8 :: level

    level = linear_interp(zs%t, zs%level, t)
    field(:) = level

  end subroutine interpolate_tide
  
  subroutine generate_meteo(par, meteo)

    type(parameters), intent(in) :: par
    type(meteorology), dimension(:), allocatable, intent(out) :: meteo
    real*8, dimension(7) :: tmp
    integer*4 :: fid, ierr, n, i, nt
    
    if (trim(par%meteo_file) /= '') then

       fid = 66

       ! count lines
       n = -1
       ierr = 0
       open(fid, file=trim(par%meteo_file))
       do while (ierr == 0)
          read(fid, *, iostat=ierr)
          n = n + 1
       end do
       rewind(fid)
       
       ! allocate arrays
       allocate(meteo(n))
       
       ! read data
       i = 1
       ierr = 0
       do while (ierr == 0)
          read(fid, *, iostat=ierr) tmp(:)
          meteo(i)%duration = tmp(1)
          meteo(i)%solar_radiation = tmp(2)
          meteo(i)%air_temperature = tmp(3)
          meteo(i)%relative_humidity = tmp(4)
          meteo(i)%air_specific_heat = tmp(5)
          meteo(i)%atmospheric_pressure = tmp(6)
          meteo(i)%latent_heat = tmp(7)
          i = i + 1
       end do
       close(fid)

       ! determine time axis
       meteo(2:n)%t = cumsum(meteo(1:n-1)%duration)

       ! checks
       if (meteo(n)%t < par%tstop) then
          write(*,*) "ERROR: meteo definition file too short"
          stop 1
       end if

    else

       allocate(meteo(2))
       meteo(1)%t = par%tstop
       meteo%solar_radiation = 0.d0
       meteo%air_temperature = 0.d0
       meteo%relative_humidity = 0.d0
       
    end if

  end subroutine generate_meteo

  subroutine interpolate_meteo(m, t, meteo)

    type(meteorology), dimension(:), intent(in) :: m
    type(meteorology), intent(out) :: meteo
    real*8, intent(in) :: t

    meteo%solar_radiation = linear_interp(m%t, m%solar_radiation, t)
    meteo%air_temperature = linear_interp(m%t, m%air_temperature, t)
    meteo%relative_humidity = linear_interp(m%t, m%relative_humidity, t)
    meteo%air_specific_heat = linear_interp(m%t, m%air_specific_heat, t)
    meteo%atmospheric_pressure = linear_interp(m%t, m%atmospheric_pressure, t)
    meteo%latent_heat = linear_interp(m%t, m%latent_heat, t)

  end subroutine interpolate_meteo
       
  subroutine compute_threshold_moisture(par, moist, u_th)

    type(parameters), intent(in) :: par
    real*8, dimension(:), intent(in) :: moist
    real*8, dimension(par%nc) :: mg
    real*8, dimension(:,:), intent(inout) :: u_th
    real*8, dimension(par%nfractions, par%nc) :: u_th_m
    integer :: i, j, n

    if (.not. par%th_moisture) return

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

  subroutine update_moisture(par, zb, zs, meteo, uw, moist)

    type(parameters), intent(in) :: par
    type(meteorology), intent(in) :: meteo
    real*8, dimension(:), intent(in) :: zb, zs, uw
    real*8, dimension(:,:), intent(inout) :: moist
    real*8 :: radiation, m, delta, gamma, evaporation
    integer :: i, j

    ! evaporation using Penman
    if (par%evaporation) then
       radiation = meteo%solar_radiation / 1e6 / par%dt * 3600 * 24 ! conversion from J/m2 to MJ/m2/day
       m = vaporation_pressure_slope(meteo%air_temperature) ! [kPa/K]
       delta = saturation_pressure(meteo%air_temperature) * (1 - meteo%relative_humidity) ! [kPa]
       gamma = (meteo%air_specific_heat * meteo%atmospheric_pressure) / &
            (.622 * meteo%latent_heat) ! [kPa/K]
    end if
    
    ! infiltration using Darcy
    do i = 1,par%nc
       if (zs(i) >= zb(i)) then
          moist(:,i) = par%porosity
       else
          moist(:,i) = moist(:,i) * exp(-par%F * par%dt)
          if (par%evaporation) then
             evaporation = max(0.d0, (m * radiation + gamma * 6.43 * (1 + 0.536 * uw(i)) * delta) / &
                  (meteo%latent_heat * (m + gamma)))
             evaporation = evaporation / 24 / 3600 / 1000 ! conversion from mm/day to m/s

             moist(:,i) = max(0.d0, moist(:,i) - evaporation * par%dt / par%layer_thickness)
          end if
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
    
  function map_moisture(zm, m, zb) result (mg)
    
    real*8, dimension(:) :: zm, m
    real*8, dimension(:) :: zb
    real*8, dimension(:), allocatable :: mg
    integer :: i, n

    allocate(mg(size(zb)))
    
    n = size(m)
    do i = 1,size(zb)
       if (zb(i) < minval(zm)) then
          mg(i) = m(1)
       elseif (zb(i) > maxval(zm)) then
          mg(i) = m(n)
       else
          mg(i) = linear_interp(zm(2:n-1), m(2:n-1), zb(i))
       end if
    end do

  end function map_moisture

end module moist_module

module run_module
! This is a docstring of run_module

  use init_module
  use output_module
  use input_module
  use wind_module
  use moist_module
  use bed_module
  use utils_module

  implicit none

contains

  subroutine step(par, s, sl, var)
    
    type(parameters), intent(inout) :: par
    type(spaceparams), intent(inout) :: s
    type(spaceparams_linear), intent(inout) :: sl
    type(variables), dimension(:), intent(inout) :: var
    integer*4 :: i, j, k, n
    real*8 :: err, alpha
    real*8, dimension(:,:,:), allocatable :: Ct2, Ct2p

    allocate(Ct2(par%nfractions, par%nx+1, par%ny+1))
    allocate(Ct2p(par%nfractions, par%nx+1, par%ny+1))
    Ct2 = 0.d0
    Ct2p = 0.d0

    ! interpolate wind
    call interpolate_wind(par%uw, par%t, s%uw, s%udir)
    s%uws = s%uw * cos(s%alfaz + s%udir)
    s%uwn = s%uw * sin(s%alfaz + s%udir)

    ! courant check
    if (trim(par%scheme) .eq. 'euler_forward') then
       if (par%CFL > 0.d0) then
          par%dt = par%CFL / (maxval(abs(s%uws) / s%dsz) + &
                              maxval(abs(s%uwn) / s%dnz))
       end if
    end if

    ! interpolate time series
    call interpolate_moist(par%moist, par%t, sl%zb, sl%moist)
    call interpolate_meteo(par%meteo, par%t, s%meteo)
    call interpolate_tide(par%zs, par%t, sl%zs)

    ! update moisture contents
    call update_moisture(par, sl%zb, sl%zs, s%meteo, sl%uw, sl%moist)
    call mix_toplayer(par, sl%zb, sl%zs)
       
    ! update threshold
    sl%uth = par%u_th
    call compute_threshold_grainsize(par, sl%uth)
    call compute_threshold_bedslope(par, s)
    call compute_threshold_moisture(par, sl%moist(1,:), sl%uth)

    ! get available mass
    sl%mass = get_layer_mass(par)
    sl%thlyr = get_layer_thickness(par)

    ! compute transport capacity by wind, including thresholds
    alpha = (0.174 / log10(par%z0/par%k))**3
    s%Cu = max(0.d0, alpha * par%Cb * par%rhoa / par%g * &
         (s%uw - s%uth)**3 / s%uw)

    ! determine first dry grid cell
    if (trim(par%scheme) .eq. 'euler_forward') then
       do i = 2,par%ny+1
          do j = 2,par%nx+1

             ! compute supply based on sediment availability
             call compute_supply(par, s%mass(:,1,j,i), &
                  par%accfac * s%Cu(:,j,i), s%Ct(:,j,i), s%supply(:,j,i), s%p(:,j,i))
          
             do k = 1,par%nfractions

                ! compute sediment advection by wind
                Ct2(k,j,i) = s%Ct(k,j,i) &
                     - s%uws(j,i) * par%dt / s%dsz(j,i) * s%dsdnzi(j,i) * (s%Ct(k,j,i) - s%Ct(k,j,i-1)) &
                     - s%uwn(j,i) * par%dt / s%dnz(j,i) * s%dsdnzi(j,i) * (s%Ct(k,j,i) - s%Ct(k,j-1,i)) &
                     + s%supply(k,j,i)
             
             end do
          end do
       end do
    elseif (trim(par%scheme) .eq. 'euler_backward') then
       do n = 1,par%max_iter

          Ct2p = Ct2

          do i = 2,par%ny+1
             do j = 2,par%nx+1
             
                ! compute supply based on sediment availability
                call compute_supply(par, s%mass(:,1,j,i), &
                     par%accfac * s%Cu(:,j,i), Ct2(:,j,i), s%supply(:,j,i), s%p(:,j,i))
             
                do k = 1,par%nfractions
                
                   ! compute sediment advection by wind
                   Ct2(k,j,i) = s%Ct(k,j,i) &
                     - s%uws(j,i) * par%dt / s%dsz(j,i) * s%dsdnzi(j,i) * (Ct2(k,j,i) - Ct2(k,j,i-1)) &
                     - s%uwn(j,i) * par%dt / s%dnz(j,i) * s%dsdnzi(j,i) * (Ct2(k,j,i) - Ct2(k,j-1,i)) &
                     + s%supply(k,j,i)
                   
                end do
             end do
          end do

          ! exit iteration if change is negligible
          err = sum(abs(Ct2 - Ct2p))
          if (err .le. par%max_error) exit

       end do

       if (err .gt. par%max_error) then
          write(0, '(a,i6,a,f10.4,a,e10.2,a)') &
               "WARNING: iteration not converged (i: ", par%nt, "; error: ", err, ")"
       end if

    end if

    s%Ct = Ct2
    s%Ct(:,1,:) = s%Ct(:,par%ny+1,:)
    
    ! add sediment deposit
    do i = 1,par%nfractions
       where (sl%zb < sl%zs)
          sl%supply(i,:) = sl%supply(i,:) - par%Cw * &
               min(par%w * par%dt, sl%zs - sl%zb) * &
               par%grain_dist(i) / max(1e-10, sum(par%grain_dist))
       end where
    end do
    s%supply(:,:,1) = 0.d0
    s%supply(:,1,:) = s%supply(:,par%ny+1,:)

    ! update bed elevation
    sl%zb = update_bed(par, sl%zb, -sl%supply, sl%rho)
    call sweep_toplayer(par)

    ! incremental output
    call output_update(var)

    par%t = par%t + par%dt
    par%nt = par%nt + 1
    
  end subroutine step
  
  subroutine compute_supply(par, mass, Cu, Ct, supply, p)

    type(parameters), intent(in) :: par
    real*8, dimension(:), intent(in) :: mass, Cu, Ct
    real*8, dimension(size(mass)), intent(inout) :: p, supply
    real*8, dimension(size(mass)) :: p_air, p_bed

    ! compute distribution in air
    p_air = Ct / max(1e-10, Cu)

    ! compute sediment distribution in bed
    p_bed = max(0.d0, mass) / max(1e-10, sum(mass))

    ! compute new sediment distributuion in the air
    p = (1 - par%bi) * p_air + p_bed * (1.d0 - min(1.d0, (1 - par%bi) * sum(p_air)))

    if (sum(p) == 0.d0) p = 1.d0
    p = p / sum(p)
    
    ! determine weighed supply
    supply = (Cu * p - Ct) / par%Tp * par%dt

    ! limit advection by available mass
    supply = min(mass, supply)
    
  end subroutine compute_supply
  
end module run_module

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
    integer*4 :: i, j, k, n, x1
    real*8 :: err, alpha
    real*8, dimension(:,:,:), allocatable :: Ct, Ct_prev

    allocate(Ct(par%nfractions, par%ny+1, par%nx+1))
    allocate(Ct_prev(par%nfractions, par%ny+1, par%nx+1))
    Ct = 0.d0
    Ct_prev = 0.d0

    ! interpolate wind
    call interpolate_wind(par%uw, par%t, s%uw)

    ! courant check
    par%dx = 5.d0
    if (trim(par%scheme) .eq. 'explicit') then
       if (par%CFL > 0.d0) then
          par%dt = par%CFL * par%dx / s%uw
          !write(0, '(a, f4.2)') "  adapted time step: ", par%dt
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
    call compute_threshold_bedslope(par, s%x, s%zb, s%uth)
    call compute_threshold_moisture(par, sl%moist(1,:), sl%uth)
    
    ! get available mass
    sl%mass = get_layer_mass(par)
    sl%thlyr = get_layer_thickness(par)

    ! compute transport capacity by wind, including thresholds
    alpha = (0.174 / log10(par%z0/par%k))**3
    s%Cu = max(0.d0, alpha * par%Cb * par%rhoa / par%g * &
         (s%uw - s%uth)**3 / (s%uw * par%VS))

    ! determine first dry grid cell
    x1 = max(2, first_exceedance(s%zb, s%zs))
    if (trim(par%scheme) .eq. 'explicit') then

       do i = x1,par%nx+1
          do j = 1,par%ny+1

             ! compute supply based on sediment availability
             s%supply(:,j,i) = compute_supply(par, s%mass(:,1,j,i), &
                  par%accfac * s%Cu(:,j,i), s%Ct(:,j,i))

             do k = 1,par%nfractions

                ! compute sediment advection by wind
                Ct(k,j,i) = -par%VS * s%uw * (s%Ct(k,j,i) - s%Ct(k,j,i-1)) * &
                     par%dt / par%dx + s%Ct(k,j,i) + s%supply(k,j,i)
             
             end do
          end do
       end do
    else

       Ct = s%Ct

       do n = 1,par%max_iter
          
          Ct_prev = Ct

          do i = x1,par%nx+1
             do j = 1,par%ny+1

                ! compute supply based on sediment availability
                s%supply(:,j,i) = compute_supply(par, s%mass(:,1,j,i), &
                     par%accfac * s%Cu(:,j,i), Ct(:,j,i))

                do k = 1,par%nfractions

                   ! compute sediment advection by wind
                   Ct(k,j,i) = (par%VS * s%uw * Ct(k,j,i-1) * &
                        par%dt / par%dx + s%Ct(k,j,i) + s%supply(k,j,i)) / &
                        (1 + par%VS * s%uw * par%dt / par%dx)
                   
                end do
             end do
          end do

          ! exit iteration if change is negligible
          err = sum(abs(Ct - Ct_prev)) / par%nc
          if (err .le. par%max_error) exit

       end do

       if (err .gt. par%max_error) then
          write(0, '(a,i6,a,f10.4,a,e10.2,a)') &
               "WARNING: iteration not converged (i: ", par%nt, "; error: ", err, ")"
       end if

    end if

    s%Ct = Ct
    
    ! add sediment deposit
    do i = 1,par%nfractions
       where (sl%zb < sl%zs)
          sl%supply(i,:) = sl%supply(i,:) - par%Cw * &
               min(par%w * par%dt, sl%zs - sl%zb) * &
               par%grain_dist(i) / max(1e-10, sum(par%grain_dist))
       end where
    end do
    s%supply(:,:,1) = 0.d0

    ! update bed elevation
    sl%zb = update_bed(par, sl%zb, -sl%supply, sl%rho)
    call sweep_toplayer(par)

    ! incremental output
    call output_update(var)

    par%t = par%t + par%dt
    par%nt = par%nt + 1
    
  end subroutine step
  
  function compute_supply(par, mass, Cu, Ct) result(supply)

    type(parameters), intent(in) :: par ! parameters structure
    real*8, dimension(:), intent(in) :: mass, Cu, Ct
    real*8, dimension(size(mass)) :: dist, dist2, supply

    dist = Ct / max(1e-10, Cu)

    if (sum(dist) < 1.d0) then ! erosion
    
       ! compute sediment distribution in bed
       dist2 = max(0.d0, mass) / max(1e-10, sum(mass))

       ! compute new sediment distributuion in the air
       dist = dist + dist2 * (1.d0 - sum(dist))

    end if

    ! compute distribution in air
    if (sum(dist) == 0.d0) dist = 1.d0
    dist = dist / sum(dist)

!    call assert(abs(sum(dist) - 1.d0) < 1e-3)

    ! determine weighed supply
    supply = (Cu * dist - Ct) / par%Tp * par%dt

    ! limit advection by available mass
    supply = min(mass, supply)

  end function compute_supply

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
  
end module run_module

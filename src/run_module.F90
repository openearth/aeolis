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

  subroutine step(par, s, var)
    
    type(parameters), intent(inout) :: par
    type(spaceparams), intent(inout) :: s
    type(variables), dimension(:), intent(inout) :: var
    integer*4 :: i, j, n, ti, x1
    real*8 :: err, alpha
    real*8, dimension(:,:), allocatable :: Ct2, Ct2_prev

    allocate(Ct2(par%nfractions, par%nx+1))
    allocate(Ct2_prev(par%nfractions, par%nx+1))
    Ct2 = 0.d0
    Ct2_prev = 0.d0

    ! update wind speed and meteo
    ti = nint(par%t / par%dt)
    s%uw = par%uw(ti)
    s%meteo = par%meteo(1)
    s%zs = par%zs(ti)

    ! update moisture contents
    call update_moisture(par, s%zb, s%zs, s%meteo, s%uw, s%moist)
    call mix_toplayer(par, s%zb, s%zs)
       
    ! update threshold
    s%uth = par%u_th
    call compute_threshold_grainsize(par, s%uth)
    call compute_threshold_bedslope(par, s%x, s%zb, s%uth)
    call compute_threshold_moisture(par, s%moist(1,:), s%uth)
    
    ! get available mass
    s%mass = get_layer_mass()
    s%thlyr = get_layer_thickness()

    ! compute transport capacity by wind, including thresholds
    alpha = (0.174 / log10(par%z0/par%k))**3
    s%Cu = max(0.d0, alpha * par%Cb * par%rhoa / par%g * &
         (s%uw - s%uth)**3 / (s%uw * par%VS))

    ! determine first dry grid cell
    x1 = 2
    do i=2,par%nx+1
       if (s%zb(i) >= s%zs) then
          x1 = i
          exit
       end if
    end do
    
    if (trim(par%scheme) .eq. 'explicit') then
       do j=x1,par%nx+1

          ! compute supply based on sediment availability
          s%supply(:,j) = compute_supply(par, s%mass(:,1,j), &
               par%accfac * s%Cu(:,j), s%Ct(:,j))

          do i=1,par%nfractions
             
             ! compute sediment advection by wind
             Ct2(i,j) = max(0.d0, -par%VS * s%uw * (s%Ct(i,j) - s%Ct(i,j-1)) * &
                  par%dt / par%dx + s%Ct(i,j) + s%supply(i,j))
             
          end do
       end do
    else
       do n=1,par%max_iter
          
          Ct2_prev = Ct2
          
          do j=x1,par%nx+1
             
             ! compute supply based on sediment availability
             s%supply(:,j) = compute_supply(par, s%mass(:,1,j), &
                  par%accfac * s%Cu(:,j), Ct2(:,j))
             
             do i=1,par%nfractions
                
                ! compute sediment advection by wind
                Ct2(i,j) = max(0.d0, (par%VS * s%uw * Ct2(i,j-1) * &
                     par%dt / par%dx + s%Ct(i,j) + s%supply(i,j)) / &
                     (1 + s%uw * par%dt / par%dx))

             end do
          end do

          ! exit iteration if change is negligible
          err = sum(abs(Ct2 - Ct2_prev))
          if (err .le. par%max_error) exit
          
       end do

       if (err .gt. par%max_error) then
          write(0, '(a,i6,a,f10.4,a,e10.2,a)') &
               "WARNING: iteration not converged (i: ", ti, "; error: ", err, ")"
       end if
       
    end if

    s%Ct = Ct2

    ! add sediment deposit
    do i = 1,par%nfractions
       where (s%zb < s%zs)
          s%supply(i,:) = s%supply(i,:) - par%Cw * &
               min(par%w * par%dt, s%zs - s%zb) * &
               par%grain_dist(i) / max(1e-10, sum(par%grain_dist))
       end where
    end do
    s%supply(:,1) = 0.d0

    ! update bed elevation
    s%zb = update_bed(par, s%zb, -s%supply, s%rho)
    call sweep_toplayer(par)

    ! incremental output
    call output_update(var)

  end subroutine step

  subroutine write_output(par, s, var)

    type(parameters), intent(inout) :: par
    type(spaceparams), intent(inout) :: s
    type(variables), dimension(:), intent(inout) :: var

    ! write output
    if (par%t .le. par%dt  .or. par%tout < par%dt .or. &
         mod(par%t, par%tout) < par%dt) then

       ! update derived variables
       if (is_output(var, 'mass')) s%mass = get_layer_mass()
       if (is_output(var, 'd10')) s%d10 = get_layer_percentile(par, 0.1d0)
       if (is_output(var, 'd50')) s%d50 = get_layer_percentile(par, 0.5d0)
       if (is_output(var, 'd90')) s%d90 = get_layer_percentile(par, 0.9d0)
       
       call output_write(var)
       call output_clear(var)
       
    end if

  end subroutine write_output
  
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

!    call assert(abs(sum(dist) - 1.d0) < 1e-10)

    ! determine weighed supply
    supply = (Cu * dist - Ct) / par%Tp * par%dt

    ! limit advection by available mass
    supply = min(mass, supply)

  end function compute_supply
  
end module run_module

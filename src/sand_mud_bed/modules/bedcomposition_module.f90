module bedcomposition_module
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2012.                                     
!                                                                               
!  This program is free software: you can redistribute it and/or modify         
!  it under the terms of the GNU General Public License as published by         
!  the Free Software Foundation version 3.                                      
!                                                                               
!  This program is distributed in the hope that it will be useful,              
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU General Public License for more details.                                 
!                                                                               
!  You should have received a copy of the GNU General Public License            
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.        
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D" and "Deltares"    
!  are registered trademarks of Stichting Deltares, and remain the property of  
!  Stichting Deltares. All rights reserved.                                     
!                                                                               
!-------------------------------------------------------------------------------
!  $Id: bedcomposition_module.f90 7697 2012-11-16 14:10:17Z boer_aj $
!  $HeadURL: https://svn.oss.deltares.nl/repos/openearthtools/trunk/programs/SandMudBedModule/03_Fortran/example/example/modules/bedcomposition_module.f90 $
!!--module description----------------------------------------------------------
!
! This module keeps track of the bed composition at one or more locations. The
! locations are indicated by indices nmlb to nmub. The bed composition consists
! of nfrac sediment fractions; their mass and volume fractions sum to 1. The
! bed may be schematized using one or more layers.
!
!!--module declarations---------------------------------------------------------

use precision
private
!
! public data types
!
public bedcomp_data
!
! public routines
!
public  bedcomposition_module_info
public  copybedcomp
public  updmorlyr
public  gettoplyr
public  detthcmud
public  getalluvthick
public  getporosity
public  getfrac
public  getmfrac
public  setmfrac
public  getvfrac
public  setvfrac
public  getsedthick
public  initmorlyr
public  setbedfracprop
public  allocmorlyr
public  clrmorlyr
public  bedcomp_use_bodsed
public  initpreload
!
public bedcomp_getpointer_integer
public bedcomp_getpointer_logical
public bedcomp_getpointer_realfp
public bedcomp_getpointer_realprec
!
! interfaces
!
interface bedcomp_getpointer_logical
   module procedure bedcomp_getpointer_logical_scalar
end interface
interface bedcomp_getpointer_integer
   module procedure bedcomp_getpointer_integer_scalar
end interface
interface bedcomp_getpointer_realfp
   module procedure bedcomp_getpointer_fp_scalar
   module procedure bedcomp_getpointer_fp_1darray
   module procedure bedcomp_getpointer_fp_2darray
   module procedure bedcomp_getpointer_fp_3darray
end interface
interface bedcomp_getpointer_realprec
   module procedure bedcomp_getpointer_prec_2darray
end interface
interface getsedthick
   module procedure getsedthick_1point
   module procedure getsedthick_allpoints
end interface
!
! morphology layers numerical settings
!
type morlyrnumericstype
    real(fp) :: MinMassShortWarning  ! minimum erosion thickness for a shortage warning
    integer  :: MaxNumShortWarning   ! maximum number of shortage warnings remaining
end type morlyrnumericstype

   
type bedcomp_settings
    !
    ! doubles
    !
    real(fp) :: thlalyr   !  thickness of Lagrangian underlayer layers
    real(fp) :: theulyr   !  thickness of Eulerian underlayer layers
    !
    ! reals
    !
    !
    ! integers
    !
    integer :: flufflyr   !  switch for fluff layer concept
                          !  0: no fluff layer
                          !  1: all mud to fluff layer, burial to bed layers
                          !  2: part mud to fluff layer, other part to bed layers (no burial)
    integer :: idiffusion !  switch for diffusion between layers
                          !  0: no diffusion
                          !  1: diffusion based on mass fractions
                          !  2: diffusion based on volume fractions
    integer :: iporosity  !  switch for porosity (simulate porosity if iporosity > 0)
                          !  0: porosity included in densities, set porosity to 0
                          !  1: ...
    integer :: iunderlyr  !  switch for underlayer concept
                          !  1: standard fully mixed concept
                          !  2: graded sediment concept
    integer :: keuler     !  index of first Eulerian (i.e. non-moving) layer
                          !  2   : standard Eulerian, only top layer moves with bed level
                          !  nlyr: fully Lagrangian (all layers move with bed level)
    integer :: nfrac      !  number of sediment fractions
    integer :: neulyr     !  number of Eulerian underlayers
    integer :: nlalyr     !  number of Lagrangian underlayers
    integer :: nlyr       !  number of layers (transport + exchange + under layers)
    integer :: ndiff      !  number of diffusion coefficients in vertical direction
    integer :: nmlb       !  start index of segments
    integer :: nmub       !  nm end index
    integer :: updbaselyr !  switch for computing composition of base layer
                          !  1: base layer is an independent layer
                          !  2: base layer composition is kept fixed
                          !  3: base layer composition is set equal to the
                          !     composition of layer above it
                          !  4: base layer composition and thickness constant
    !
    ! pointers
    !
    type (morlyrnumericstype) , pointer :: morlyrnum ! structure containing numerical settings
    real(fp) , dimension(:)   , pointer :: phi       ! D50 diameter expressed on phi scale, units : -
    real(fp) , dimension(:)   , pointer :: rhofrac   ! density of fraction (specific density or including pores), units : kg/m3
    integer  , dimension(:)   , pointer :: sedtyp    ! sediment type: 0=total/1=noncoh/2=coh
    real(fp) , dimension(:)   , pointer :: sigphi    ! standard deviation expressed on phi scale, units : -
    real(fp) , dimension(:)   , pointer :: thexlyr   ! thickness of exchange layer, units : m
    real(fp) , dimension(:)   , pointer :: thtrlyr   ! thickness of transport layer, units : m
    real(fp) , dimension(:)   , pointer :: zdiff     ! depth below bed level for which diffusion coefficients are defined, units : m
    real(fp) , dimension(:,:) , pointer :: bfluff0   ! burial parameter fluff layer (only when FluffLayer=1) [kg/m2/s]
    real(fp) , dimension(:,:) , pointer :: bfluff1   ! burial parameter fluff layer (only when FluffLayer=1) [1/s]
    real(fp) , dimension(:,:) , pointer :: kdiff     ! diffusion coefficients for mixing between layers, units : m2/s
    ! 
    ! logicals
    !
    logical :: exchlyr    !  flag for use of exchange layer (underlayer bookkeeping system)
    !
    ! characters
    !
end type bedcomp_settings
!
type bedcomp_state
    real(fp)   , dimension(:,:)  , pointer :: preload       ! historical largest load, units : kg
    real(prec) , dimension(:,:)  , pointer :: bodsed        ! Array with total sediment, units : kg /m2
    real(fp)   , dimension(:)    , pointer :: dpsed         ! Total depth sediment layer, units : m
    real(fp)   , dimension(:,:)  , pointer :: fluffshort    ! sediment shortage in fluff layer, units : kg /m2
    real(fp)   , dimension(:,:)  , pointer :: mfluff        ! composition of fluff layer: mass of mud fractions, units : kg /m2
    real(fp)   , dimension(:,:,:), pointer :: msed          ! composition of morphological layers: mass of sediment fractions, units : kg /m2
    real(fp)   , dimension(:,:)  , pointer :: sedshort      ! sediment shortage in transport layer, units : kg /m2
    real(fp)   , dimension(:,:)  , pointer :: svfrac        ! 1 - porosity coefficient, units : -
    real(fp)   , dimension(:,:)  , pointer :: thlyr         ! thickness of morphological layers, units : m
end type bedcomp_state
!
type bedcomp_work
    real(fp), dimension(:,:) , pointer :: msed2
    real(fp), dimension(:)   , pointer :: svfrac2
    real(fp), dimension(:)   , pointer :: thlyr2
end type bedcomp_work
!
type bedcomp_data
   private
   type (bedcomp_settings), pointer :: settings
   type (bedcomp_state)   , pointer :: state
end type bedcomp_data

contains
!
!
!
!==============================================================================
subroutine bedcomposition_module_info(messages)
    use message_module
    !
    type(message_stack), pointer :: messages
    !
    call addmessage(messages,'$Id: bedcomposition_module.f90 7697 2012-11-16 14:10:17Z boer_aj $')
    call addmessage(messages,'$URL: https://svn.oss.deltares.nl/')

end subroutine bedcomposition_module_info
!
!
!
!==============================================================================
function updmorlyr(this, dbodsd, dfluff, rhosol, dt, morfac, dz, messages) result (istat)
!!--description-----------------------------------------------------------------
!
!    Function: - Update underlayer bookkeeping system for erosion/sedimentation
!
!!--declarations----------------------------------------------------------------
    use precision
    use message_module
    implicit none
    !
    ! Call variables
    !
    type(bedcomp_data)                                                                           :: this 
    real(fp), dimension(this%settings%nfrac, this%settings%nmlb:this%settings%nmub), intent(in)  :: dfluff  !  change in sediment composition, units : kg/m2
    real(fp), dimension(this%settings%nfrac, this%settings%nmlb:this%settings%nmub), intent(in)  :: dbodsd  !  change in sediment composition, units : kg/m2
    real(fp)                                                                       , intent(in)  :: dt      !  time step, units : s
    real(fp)                                                                       , intent(in)  :: morfac  !  morphological scale factor, units : -
    real(fp), dimension(this%settings%nfrac)                                       , intent(in)  :: rhosol  !  density of sediment fractions, units : kg/m3
    real(fp), dimension(this%settings%nmlb:this%settings%nmub)                     , intent(out) :: dz      !  change in bed level, units : m
    type(message_stack)                                                                          :: messages
    integer                                                                                      :: istat
    !
    ! Local variables
    !
    integer                                 :: l
    integer                                 :: nm
    real(fp)                                :: dtmorfac
    real(fp)                                :: seddep0
    real(fp)                                :: seddep1
    real(fp)                                :: thdiff
    real(fp)                                :: fac
    real(fp)                                :: mfltot
    real(fp)                                :: temp
    real(fp)                                :: thick
    real(fp)                                :: thtrlyrnew
    real(fp), dimension(this%settings%nfrac):: dmi
    type (bedcomp_work)                     :: work
    !
    character(message_len)                  :: message
    type (morlyrnumericstype)     , pointer :: morlyrnum
    real(prec) , dimension(:,:)   , pointer :: bodsed
    real(fp)   , dimension(:,:)   , pointer :: bfluff0
    real(fp)   , dimension(:,:)   , pointer :: bfluff1
    real(fp)   , dimension(:)     , pointer :: dpsed
    real(fp)   , dimension(:,:)   , pointer :: fluffshort
    real(fp)   , dimension(:,:)   , pointer :: mfluff
    real(fp)   , dimension(:,:,:) , pointer :: msed
    real(fp)   , dimension(:)     , pointer :: rhofrac
    real(fp)   , dimension(:,:)   , pointer :: sedshort
    real(fp)   , dimension(:,:)   , pointer :: svfrac 
    real(fp)   , dimension(:,:)   , pointer :: thlyr
    real(fp)   , dimension(:)     , pointer :: thtrlyr   
!
!! executable statements -------------------------------------------------------
!
    bfluff0     => this%settings%bfluff0
    bfluff1     => this%settings%bfluff1
    rhofrac     => this%settings%rhofrac
    morlyrnum   => this%settings%morlyrnum
    rhofrac     => this%settings%rhofrac
    thtrlyr     => this%settings%thtrlyr
    bodsed      => this%state%bodsed
    dpsed       => this%state%dpsed
    fluffshort  => this%state%fluffshort
    mfluff      => this%state%mfluff
    msed        => this%state%msed
    sedshort    => this%state%sedshort
    svfrac      => this%state%svfrac
    thlyr       => this%state%thlyr  
    !
    istat = 0
    if (allocwork(this,work) /= 0) return
    ! 
    ! update fluff layer            
    !
    if (this%settings%flufflyr>0) then
        do nm = this%settings%nmlb,this%settings%nmub
            mfltot = 0.0_fp
            do l = 1, this%settings%nfrac
                temp  = mfluff(l, nm) + dfluff(l, nm)
                if (temp < 0.0_fp) then
                    if (temp < -morlyrnum%MinMassShortWarning .and. morlyrnum%MaxNumShortWarning>0) then
                        morlyrnum%MaxNumShortWarning = morlyrnum%MaxNumShortWarning - 1
                        write(message,'(a,i5,a,i3,a,e20.4,a,e20.4)') &
                            & 'Sediment erosion shortage in fluff layer at NM ', nm, ' Fraction: ', l, &
                            & ' Mass available   : ' ,mfluff(l, nm), &
                            & ' Mass to be eroded: ', dfluff(l, nm)
                        call addmessage(messages,message)
                        if (morlyrnum%MaxNumShortWarning == 0) then
                            message = 'Sediment erosion shortage messages suppressed'
                            call addmessage(messages,message)
                        endif
                    endif
                    fluffshort(l, nm) = fluffshort(l, nm) + temp
                    temp = 0.0_fp
                elseif ( fluffshort(l, nm) < 0.0 ) then
                    temp = temp + fluffshort(l, nm)
                    if ( temp < 0.0_fp ) then
                        fluffshort(l, nm) = temp
                        temp = 0.0_fp
                    else
                        fluffshort(l, nm) = 0.0_fp
                    endif
                endif
                mfluff(l, nm) = temp
                mfltot = mfltot + mfluff(l,nm)
            enddo               
            !
            ! Burial to bed layers
            !
            if (this%settings%flufflyr==1) then
                do l = 1, this%settings%nfrac
                    if (comparereal(mfltot,0.0_fp) == 1) then
                        fac = mfluff(l, nm)/mfltot
                    else
                        fac = 0.0_fp
                    endif
                    dmi(l)         = fac*min(mfltot*bfluff1(l, nm), bfluff0(l, nm))*dt
                    ! limit burial to available sediment in fluff layer
                    dmi(l)         = min(dmi(l),mfluff(l, nm))
                    mfluff(l, nm)  = mfluff(l, nm)  - dmi(l)                    
                    select case(this%settings%iunderlyr)
                    case(2)
                        msed(l, 1, nm) = msed(l, 1, nm) + dmi(l)*morfac ! klopt dit?
                    case default
                        bodsed(l, nm)  = bodsed(l, nm)  + dmi(l)*morfac ! klopt dit?
                    end select
                enddo
            endif
        enddo
    endif
    !
    select case(this%settings%iunderlyr)
    case(2)
       do nm = this%settings%nmlb,this%settings%nmub
          call getsedthick_1point(this, nm, seddep0)
          if (.not.this%settings%exchlyr) then
             !
             ! transport layer only
             !
             do l = 1, this%settings%nfrac
                temp  = msed(l, 1, nm) + dbodsd(l, nm)
                if (temp < 0.0_fp) then
                   if (temp < -morlyrnum%MinMassShortWarning .and. morlyrnum%MaxNumShortWarning>0) then
                      morlyrnum%MaxNumShortWarning = morlyrnum%MaxNumShortWarning - 1
                      write(message,'(a,i5,a,i3,a,e20.4,a,e20.4)') &
                         & 'Sediment erosion shortage at NM ', nm, ' Fraction: ', l, &
                         & ' Mass available   : ' ,msed(l, 1, nm), &
                         & ' Mass to be eroded: ', dbodsd(l, nm)
                      call addmessage(messages,message)
                      if (morlyrnum%MaxNumShortWarning == 0) then
                         message = 'Sediment erosion shortage messages suppressed'
                         call addmessage(messages,message)
                      endif
                   endif
                   sedshort(l, nm) = sedshort(l, nm) + temp
                   temp = 0.0_fp
                elseif ( sedshort(l, nm) < 0.0 ) then
                   temp = temp + sedshort(l, nm)
                   if ( temp < 0.0_fp ) then
                      sedshort(l, nm) = temp
                      temp = 0.0_fp
                   else
                      sedshort(l, nm) = 0.0_fp
                   endif
                endif
                msed(l, 1, nm) = temp
             enddo
             !
             ! get new requested transport layer thickness.
             !
             thtrlyrnew = thtrlyr(nm)
             !
             ! compute actual current thickness of top layer
             !
             call updateporosity(this, nm, 1)
             thick = 0.0_fp
             do l = 1, this%settings%nfrac
                 thick = thick + msed(l, 1, nm)/rhofrac(l)
             enddo
             thick = thick/svfrac(1, nm)
             !
             thdiff = thick-thtrlyrnew
             !
             ! get sediment from or put sediment into underlayers
             ! to get transport layer of requested thickness
             !
             if ( thdiff > 0.0_fp ) then
                !   
                ! sedimentation to underlayers
                ! determine surplus of mass per fraction
                ! 
                fac = thdiff/thick
                do l = 1, this%settings%nfrac
                    dmi(l) = msed(l, 1, nm)*fac
                    msed(l, 1, nm) = msed(l, 1, nm) - dmi(l)
                enddo
                !
                ! store surplus of mass in underlayers
                !
                call lyrsedimentation(this , nm, thdiff, dmi, svfrac(1, nm), work)
                !
             elseif ( thdiff < 0.0_fp ) then
                !
                ! erosion of underlayers
                ! total erosion thickness: thdiff
                ! associated mass returned in: dmi
                !
                thdiff = -thdiff
                !  
                call lyrerosion(this , nm, thdiff, dmi)
                !
                ! add to top layer
                ! 
                do l = 1, this%settings%nfrac 
                    msed(l, 1, nm)   = msed(l, 1, nm) + dmi(l)
                enddo
                !
                do l = 1, this%settings%nfrac
                    if (sedshort(l, nm) < 0.0_fp .and. msed(l, 1, nm) > 0.0_fp) then
                        sedshort(l, nm) = sedshort(l, nm) + msed(l, 1, nm)
                        if (sedshort(l, nm) > 0.0_fp) then
                            msed(l, 1, nm)  = sedshort(l, nm)
                            sedshort(l, nm) = 0.0_fp
                        else
                            msed(l, 1, nm) = 0.0_fp
                        endif
                    endif
                enddo
                !
                call updateporosity(this, nm, 1)
                thick = 0.0_fp
                do l = 1, this%settings%nfrac
                    thick = thick + msed(l, 1, nm)/rhofrac(l)
                enddo
                thick = thick/svfrac(1, nm)
                !
                ! if there is not enough sediment in the bed then the actual
                ! thickness thick of the top layer may not reach the desired
                ! thickness thtrlyrnew, so we should here use thick as the
                ! thickness instead of thtrlyrnew
                !
                thtrlyrnew = thick
             endif
             thlyr(1, nm) = thtrlyrnew  
          else
             !
             ! transport and exchange layer
             !
             message = 'Updating exchange layer not yet implemented'
             call adderror(messages,message)
             istat = -1
             exit
          endif
          call consolidate(this, nm)
          call getsedthick_1point(this, nm, seddep1)
          dz(nm) = seddep1-seddep0
          !
          ! mixing between layers
          !  
          dtmorfac = dt*morfac
          if (this%settings%idiffusion>0) call lyrdiffusion(this, nm, dtmorfac, rhosol)
          !
       enddo ! end loop over nm
    case default
       do nm = this%settings%nmlb,this%settings%nmub
          seddep0   = dpsed(nm)
          dpsed(nm) = 0.0_fp
          dz(nm)    = 0.0_fp
          do l = 1, this%settings%nfrac
             bodsed(l, nm) = bodsed(l, nm) + real(dbodsd(l, nm),prec)
             if (bodsed(l, nm) < 0.0_prec) then
                if (bodsed(l, nm) < real(-morlyrnum%MinMassShortWarning,prec) .and. morlyrnum%MaxNumShortWarning>0) then
                   morlyrnum%MaxNumShortWarning = morlyrnum%MaxNumShortWarning - 1
                   write(message,'(a,i0,a,i0,a,e20.4,a,e20.4)') &
                      & 'Sediment erosion shortage at NM ', nm, ' Fraction: ', l, &
                      & ' Mass available   : ' ,bodsed(l, nm), &
                      & ' Mass to be eroded: ', dbodsd(l, nm)
                   call addmessage(messages,message)
                   if (morlyrnum%MaxNumShortWarning == 0) then
                      message = 'Sediment erosion shortage messages suppressed'
                      call addmessage(messages,message)
                   endif
                endif
                bodsed(l, nm) = 0.0_prec
             endif
             dpsed(nm) = dpsed(nm) + real(bodsed(l, nm),fp)/rhofrac(l)
          enddo    
          dz(nm) = dpsed(nm) - seddep0
       enddo
    endselect
    istat = deallocwork(this,work)
end function updmorlyr
!
!
!
!==============================================================================
function gettoplyr(this, dz_eros, dbodsd, messages  ) result (istat)
!!--description-----------------------------------------------------------------
!
!    Function: - Retrieve the sediment in the top layer (thickness dz_eros
!                in m) from the underlayer bookkeeping system and update the
!                administration.
!
!!--declarations----------------------------------------------------------------
    use precision
    use message_module
    !
    implicit none
    !
    ! Function/routine arguments
    !
    type(bedcomp_data)                                                                           :: this
    real(fp), dimension(this%settings%nfrac, this%settings%nmlb:this%settings%nmub), intent(out) :: dbodsd  !  change in sediment composition, units : kg/m2
    real(fp), dimension(this%settings%nmlb:this%settings%nmub)                     , intent(in)  :: dz_eros !  change in sediment thickness, units : m
    type(message_stack)                                                                          :: messages
    integer                                                                                      :: istat
    !
    ! Local variables
    !
    integer                                 :: k
    integer                                 :: l
    integer                                 :: nm
    real(fp)                                :: dz
    real(fp)                                :: fac
    real(fp)                                :: thtrlyrnew
    real(fp)                                :: thick
    real(fp), dimension(this%settings%nfrac):: dmi 
    real(fp)                                :: dz_togo
    type (bedcomp_work)                     :: work
    !
    character(message_len)                  :: message
    real(prec) , dimension(:,:)   , pointer :: bodsed
    real(fp)   , dimension(:,:)   , pointer :: svfrac
    real(fp)   , dimension(:)     , pointer :: dpsed
    real(fp)   , dimension(:,:,:) , pointer :: msed
    real(fp)   , dimension(:,:)   , pointer :: sedshort 
    real(fp)   , dimension(:)     , pointer :: rhofrac
    real(fp)   , dimension(:,:)   , pointer :: thlyr
    real(fp)   , dimension(:)     , pointer :: thtrlyr
    !
    !! executable statements -------------------------------------------------------
    !
    thtrlyr     => this%settings%thtrlyr
    rhofrac     => this%settings%rhofrac
    svfrac      => this%state%svfrac
    bodsed      => this%state%bodsed
    dpsed       => this%state%dpsed
    msed        => this%state%msed
    sedshort    => this%state%sedshort
    thlyr       => this%state%thlyr
    !
    istat = allocwork(this,work)
    if (istat /= 0) return
    select case(this%settings%iunderlyr)
    case(2)
       do nm = this%settings%nmlb,this%settings%nmub
          if ( dz_eros(nm) < 0.0_fp ) then
             !
             ! erosion should not be negative
             !
             message = 'Negative dz_eros encountered'
             call adderror(messages,message)
             istat = -1
             return
          elseif (comparereal(dz_eros(nm),0.0_fp) == 0) then
             !
             ! skip it
             !
             do l = 1, this%settings%nfrac
                dbodsd (l, nm)    = 0.0_fp
             enddo
          elseif (.not.this%settings%exchlyr) then
             !
             ! transport layer only
             !  
             if (dz_eros(nm) <= thlyr(1, nm)) then
                !
                ! erosion less than transport layer thickness
                !
                thick = thlyr(1, nm) - dz_eros(nm)
                fac = dz_eros(nm)/thlyr(1, nm)
                do l = 1, this%settings%nfrac
                   dbodsd (l, nm) = msed(l, 1, nm)*fac
                   msed(l, 1, nm) = msed(l, 1, nm) - dbodsd(l, nm)
                enddo
             else
                !
                ! erosion more than transport layer thickness
                !
                do l = 1, this%settings%nfrac
                   dbodsd (l, nm) = msed(l, 1, nm)
                   msed(l, 1, nm) = 0.0_fp
                enddo
                dz_togo = dz_eros(nm)-thlyr(1, nm)
                thick  = 0.0_fp
                !
                ! get remaining dz_togo from underlayers
                ! get dmi from underlayers 
                !
                call lyrerosion(this , nm, dz_togo, dmi)
                !
                do l = 1, this%settings%nfrac
                   dbodsd(l, nm) = dbodsd(l, nm) + dmi(l)
                enddo
             endif
             !
             ! get new transport layer thickness.
             !
             thtrlyrnew = thtrlyr(nm)
             !
             ! compute new thickness of top layer
             ! get sediment from or put sediment into underlayers
             ! to get transport layer of requested thickness
             !
             dz = thick-thtrlyrnew
             !
             if ( dz > 0.0_fp ) then
                !   
                ! sedimentation to underlayers
                ! 
                fac = dz/thick
                do l = 1, this%settings%nfrac
                    dmi(l) = msed(l, 1, nm)*fac
                    msed(l, 1, nm) = msed(l, 1, nm) - dmi(l)
                enddo   
                !
                ! store surplus of mass in underlayers
                !
                call lyrsedimentation(this , nm, dz, dmi, svfrac(1, nm), work)
                !
             elseif ( dz < 0.0_fp ) then
                !
                ! erosion of underlayers
                !
                dz = -dz
                call lyrerosion(this , nm, dz, dmi)
                !
                ! add to top layer
                !  
                do l = 1, this%settings%nfrac 
                    msed(l, 1, nm)   = msed(l, 1, nm) + dmi(l) 
                enddo
                !
                do l = 1, this%settings%nfrac
                    if (sedshort(l, nm) < 0.0_fp .and. msed(l, 1, nm) > 0.0_fp) then
                        sedshort(l, nm) = sedshort(l, nm) + msed(l, 1, nm)
                        if (sedshort(l, nm) > 0.0_fp) then
                            msed(l, 1, nm)  = sedshort(l, nm)
                            sedshort(l, nm) = 0.0_fp
                        else
                            msed(l, 1, nm) = 0.0_fp
                        endif
                    endif
                enddo
                !
                call updateporosity(this, nm, 1)
                thick = 0.0_fp
                do l = 1, this%settings%nfrac
                    thick = thick + msed(l, 1, nm)/rhofrac(l)
                enddo
                thick = thick/svfrac(1, nm)
                !
                ! if there is not enough sediment in the bed then the actual
                ! thickness thick of the top layer may not reach the desired
                ! thickness thtrlyrnew, so we should here use thick as the
                ! thickness instead of thtrlyrnew
                !
                thtrlyrnew = thick
             endif
             thlyr(1, nm) = thtrlyrnew
          else
             !
             ! transport and exchange layer
             !
             message = 'Updating exchange layer not yet implemented'
             call adderror(messages,message)
             istat = -1
             exit
          endif
       enddo
    case default
       do nm = this%settings%nmlb,this%settings%nmub
          if (dz_eros(nm)<0.0_fp) then
             !
             ! erosion should not be negative
             !
             message = 'Negative dz_eros encountered'
             call adderror(messages,message)
             istat = -1
             return
          elseif (comparereal(dz_eros(nm),0.0_fp) == 0 .or. dpsed(nm)<=0.0_fp) then
             !
             ! skip it
             !
             do l = 1, this%settings%nfrac
                dbodsd (l, nm)    = 0.0_fp
             enddo
          elseif (dz_eros(nm) < dpsed(nm)) then
             !
             ! some sediment remaining
             !
             fac = dz_eros(nm)/dpsed(nm)
             dpsed(nm) = 0.0_fp
             do l = 1, this%settings%nfrac
                dbodsd(l, nm) = real(bodsed(l, nm),fp)*fac
                bodsed(l, nm) = bodsed(l, nm) - real(dbodsd(l, nm),prec)
                dpsed (nm)    = dpsed(nm) + real(bodsed(l, nm),fp)/rhofrac(l)
             enddo
          else
             !
             ! no sediment remaining
             !
             dpsed(nm) = 0.0_fp
             do l = 1, this%settings%nfrac
                dbodsd(l, nm) = real(bodsed(l, nm),fp)
                bodsed(l, nm) = 0.0_hp
             enddo
          endif
       enddo
    endselect
    istat = deallocwork(this,work)
end function gettoplyr
!
!
!
!==============================================================================
subroutine lyrerosion(this, nm, dzini, dmi)
!!--description-----------------------------------------------------------------
!
!    Function:
!     - lyrerosion implements the erosion of sediment from the layers below the
!       transport and exchange layers
!
!!--declarations----------------------------------------------------------------
    use precision
    !
    implicit none
    !
    ! Function/routine arguments
    !
    type(bedcomp_data)                                    :: this    
    integer                                  , intent(in) :: nm
    real(fp)                                 , intent(in) :: dzini   !  thickness of eroded layer, units : m
    real(fp), dimension(this%settings%nfrac) , intent(out):: dmi     !  density of sediment fractions, units : kg/m3
    !
    ! Local variables
    !
    logical                                            :: remove
    integer                                            :: k
    integer                                            :: kero1  ! top-most layer that has been (partially) eroded
    integer                                            :: l
    real(fp)                                           :: dz
    real(fp)                                           :: dm
    real(fp)                                           :: fac
    real(fp)                                           :: thick
    real(fp)                                           :: thbaselyr
    real(fp), dimension(this%settings%nfrac)           :: mbaselyr  
    real(fp)                                 , pointer :: thlalyr
    integer                                  , pointer :: keuler 
    integer                                  , pointer :: nlyr
    integer                                  , pointer :: updbaselyr
    real(fp), dimension(:,:)                 , pointer :: svfrac
    real(fp), dimension(:,:,:)               , pointer :: msed
    real(fp), dimension(:,:)                 , pointer :: thlyr
!
!! executable statements -------------------------------------------------------
!
    keuler      => this%settings%keuler
    nlyr        => this%settings%nlyr
    thlalyr     => this%settings%thlalyr
    updbaselyr  => this%settings%updbaselyr
    svfrac      => this%state%svfrac
    msed        => this%state%msed
    thlyr       => this%state%thlyr
    !
    k   = 2
    if (this%settings%exchlyr) k = 3
    !
    thbaselyr = thlyr(nlyr, nm)
    mbaselyr  = msed(:, nlyr, nm)
    dmi = 0.0_fp
    dz  = dzini
    !
    ! initially remove sediment irrespective of layer type
    ! then fill the Lagrangian layers again up to their
    ! original thickness
    ! kero1 represents the Lagrangian layer that was eroded and needs
    ! to be replenished
    ! remove indicates that sediment should be eroded (stored in dmi)
    ! rather than shifted to another Lagrangian layer
    !
    kero1 = k-1
    remove = .true.
    do while (dz>0.0_fp .and. k<=nlyr)
        if ( thlyr(k, nm) < dz ) then
            !
            ! more sediment is needed than there is available in layer
            ! k, so all sediment should be removed from this layer
            !          
            do l = 1, this%settings%nfrac
                if (remove) then
                   dmi(l) = dmi(l) + msed(l, k, nm)
                else
                   msed(l, kero1, nm) = msed(l, kero1, nm) + msed(l, k, nm)
                endif
                msed(l, k, nm) = 0.0_fp
            enddo
            dz          = dz - thlyr(k, nm)
            if (.not.remove) then
               svfrac(kero1, nm) = svfrac(kero1, nm)*thlyr(kero1, nm) + svfrac(k, nm)*thlyr(k, nm)
               thlyr(kero1, nm)  = thlyr(kero1, nm) + thlyr(k, nm)
               svfrac(kero1, nm) = svfrac(kero1, nm)/thlyr(kero1, nm)
            endif
            thlyr(k, nm) = 0.0_fp
            k           = k+1
        else ! thlyr(k)>=dz
            !
            ! layer k contains more sediment than is needed, so only part
            ! of the sediment has to be removed from the layer
            !
            fac = dz/thlyr(k, nm)
            do l = 1, this%settings%nfrac
                dm = msed(l, k, nm)*fac
                if (remove) then
                   dmi(l)             = dmi(l) + dm 
                else
                   msed(l, kero1, nm) = msed(l, kero1, nm) + dm
                endif
                msed(l, k, nm) = msed(l, k, nm) - dm
            enddo
            thlyr(k, nm)       = thlyr(k, nm)   - dz
            if (.not.remove) then
               svfrac(kero1, nm) = svfrac(kero1, nm)*thlyr(kero1, nm) + svfrac(k, nm)*dz
               thlyr(kero1, nm)  = thlyr(kero1, nm) + dz
               svfrac(kero1, nm) = svfrac(kero1, nm)/thlyr(kero1, nm)
            endif
            !
            ! erosion complete (dz=0) now continue to replenish the
            ! (partially) eroded Lagrangian layers as long as
            ! sediment is available in lower layers. Note that the
            ! Eulerian layers don't get replenished.
            !
            kero1 = kero1+1
            remove = .false.
            !
            ! do we have to fill again some of the Lagrangian layers?
            !
            if (kero1<keuler) then
                dz = max(thlalyr - thlyr(kero1, nm),0.0_fp)
                k = max(k,kero1+1)
            else
                dz = 0.0_fp
            endif
        endif
    enddo
    !
    ! update composition of base layer
    !
    select case (updbaselyr)
    case(1) ! compute separate composition for base layer
       !
       ! no change necessary
       !
    case(2) ! composition of base layer constant
       !
       ! compute new masses
       !
       fac = thlyr(nlyr, nm)/thbaselyr
       do l = 1, this%settings%nfrac
          msed(l, nlyr, nm) = mbaselyr(l)*fac
       enddo
    case(3) ! same as the (first non-empty) layer above it
       !
       ! find lowest non-empty layer
       !
       do k = nlyr-1,1,-1
          if ( thlyr(k, nm) > 0.0_fp ) exit
       enddo
       fac = thlyr(nlyr, nm)/thlyr(k, nm)
       do l = 1, this%settings%nfrac
          msed(l, nlyr, nm) = msed(l, k, nm)*fac
       enddo
    case(4) ! composition and thickness of base layer constant
       !
       ! reset thickness and masses
       !
       thlyr(nlyr, nm)  = thbaselyr
       msed(:, nlyr, nm) = mbaselyr
    case default
       !
       ! ERROR
       !
    endselect
end subroutine lyrerosion
!
!
!
!==============================================================================
subroutine lyrsedimentation(this, nm, dzini, dmi, svfracdep, work)
!!--description-----------------------------------------------------------------
!
!    Function:
!     - lyrsedimentation implements the deposition of sediment in the layers
!       below the transport and exchange layers
!
!!--declarations----------------------------------------------------------------
    use precision
    !
    implicit none
    !
!
! Function/routine arguments
!
    type(bedcomp_data)                                    :: this   
    integer                                  , intent(in) :: nm
    real(fp)                                 , intent(in) :: dzini
    real(fp)                                 , intent(in) :: svfracdep
    real(fp), dimension(this%settings%nfrac)              :: dmi
    type(bedcomp_work)                                    :: work
!
! Local variables
!
    integer                                     :: k
    integer                                     :: k2
    integer                                     :: kmin
    integer                                     :: kne
    integer                                     :: l
    real(fp)                                    :: depfrac
    real(fp)                                    :: dm
    real(fp)                                    :: dz
    real(fp)                                    :: dz2
    real(fp)                                    :: dzc
    real(fp)                                    :: dzk
    real(fp)                                    :: dzlmax
    real(fp)                                    :: fac
    real(fp)                                    :: newthlyr
    real(fp)                                    :: thick
    real(fp)                  , pointer         :: thlalyr
    integer                   , pointer         :: keuler
    integer                   , pointer         :: nlyr
    integer                   , pointer         :: updbaselyr
    real(fp), dimension(:,:)  , pointer         :: svfrac
    real(fp), dimension(:,:,:), pointer         :: msed
    real(fp), dimension(:,:)  , pointer         :: thlyr
    real(fp), dimension(this%settings%nfrac)    :: dmi2
!
!! executable statements -------------------------------------------------------
!
    thlalyr     => this%settings%thlalyr
    keuler      => this%settings%keuler
    nlyr        => this%settings%nlyr
    updbaselyr  => this%settings%updbaselyr
    svfrac      => this%state%svfrac
    msed        => this%state%msed
    thlyr       => this%state%thlyr
    !
    kmin = 2
    if (this%settings%exchlyr) kmin = 3
    dz = dzini
    !
    ! copy Lagrangian layer data to temporary array
    !
    do k = kmin,keuler-1
       do l = 1, this%settings%nfrac
          work%msed2(l, k) = msed(l, k, nm)
          msed(l, k, nm)   = 0.0_fp
       enddo
       work%svfrac2(k)  = svfrac(k, nm)
       work%thlyr2(k)   = thlyr(k, nm)
       thlyr(k, nm)     = 0.0_fp
    enddo
    !
    ! fill the Lagrangian layers from top to bottom
    !
    do k = kmin,keuler-1
       !
       ! while there is newly deposited sediment use that to fill the layers
       !
       if (dz>thlalyr) then
          !
          ! sufficient newly deposited sediment to fill whole layer
          !
          fac = thlalyr/dz
          do l = 1, this%settings%nfrac
             dm = dmi(l) * fac
             msed(l, k, nm) = dm
             dmi(l) = dmi(l) - dm
          enddo
          svfrac(k, nm) = svfracdep
          thlyr(k, nm)  = thlalyr
          dz            = dz - thlalyr
       elseif (dz>0.0_fp) then
          !
          ! last bit of newly deposited sediment to partially fill the layer
          !
          do l = 1, this%settings%nfrac
             msed(l, k, nm) = dmi(l)
             dmi(l) = 0.0_fp
          enddo
          svfrac(k, nm) = svfracdep
          thlyr(k, nm)  = dz
          dz            = 0.0_fp
       endif
       !
       ! as long as there is still space in this layer, fill it with sediment
       ! from the old Lagrangian layers
       !
       k2  = kmin
       dzc = thlalyr - thlyr(k, nm)
       do while (k2<keuler .and. dzc > 0.0_fp)
          !
          ! did this Lagrangian layer contain sediment?
          !
          if (work%thlyr2(k2)>0.0_fp) then
             if (dzc<work%thlyr2(k2)) then
                !
                ! sufficient sediment left in old layer k2 to fill the
                ! remainder of new layer k
                !
                fac = dzc/work%thlyr2(k2)
                do l = 1, this%settings%nfrac
                   dm = work%msed2(l, k2) * fac
                   msed(l, k, nm)    = msed(l, k, nm)    + dm
                   work%msed2(l, k2) = work%msed2(l, k2) - dm
                enddo
                svfrac(k, nm)   = svfrac(k, nm)*thlyr(k, nm) + work%svfrac2(k2)*dzc
                thlyr(k, nm)    = thlalyr
                svfrac(k, nm)   = svfrac(k, nm)/thlyr(k, nm)
                work%thlyr2(k2) = work%thlyr2(k2) - dzc
             else
                !
                ! all the sediment of old layer k2 fits into new layer k
                !
                do l = 1, this%settings%nfrac
                   msed(l, k, nm)    = msed(l, k, nm) + work%msed2(l, k2)
                   work%msed2(l, k2) = 0.0_fp
                enddo
                svfrac(k, nm)   = svfrac(k, nm)*thlyr(k, nm) + work%svfrac2(k2)*work%thlyr2(k2)
                thlyr(k, nm)    = thlyr(k, nm) + work%thlyr2(k2)
                svfrac(k, nm)   = svfrac(k, nm)/thlyr(k, nm)
                work%thlyr2(k2) = 0.0_fp
             endif
             dzc = thlalyr - thlyr(k, nm)
          endif
          k2 = k2+1
       enddo
    enddo
    !
    ! if there is any sediment left in the old Lagrangian layers then move it
    ! to the Eulerian layers below
    !
    do k2 = keuler-1,kmin,-1
       if (work%thlyr2(k2)>0.0_fp) then
          do l = 1, this%settings%nfrac
             dmi2(l) = work%msed2(l, k2)
          enddo
          call lyrsedimentation_eulerian(this, nm, work%thlyr2(k2), dmi2, work%svfrac2(k2))
       endif
    enddo
    !
    ! and finally if there is any sediment left in the original deposit that
    ! came from the active layer(s) then deposit that sediment
    !
    if (dz>0.0_fp) then
       call lyrsedimentation_eulerian(this, nm, dz, dmi, svfracdep)
    endif    
end subroutine lyrsedimentation
!
!
!
!==============================================================================
subroutine lyrsedimentation_eulerian(this, nm, dzini, dmi, svfracdep)
!!--description-----------------------------------------------------------------
!
!    Function:
!     - lyrsedimentation_eulerian implements the deposition of sediment in the 
!       Eulerian layers below the transport, exchange and other Lagrangian
!       layers
!
!!--declarations----------------------------------------------------------------
    use precision
    !
    implicit none
    !
!
! Function/routine arguments
!
    type(bedcomp_data)                                    :: this   
    integer                                  , intent(in) :: nm
    real(fp)                                 , intent(in) :: dzini
    real(fp)                                 , intent(in) :: svfracdep
    real(fp), dimension(this%settings%nfrac)              :: dmi
!
! Local variables
!
    integer                                     :: k
    integer                                     :: kmin
    integer                                     :: kne
    integer                                     :: l

    real(fp)                                    :: depfrac
    real(fp)                                    :: dm
    real(fp)                                    :: dz
    real(fp)                                    :: fac
    real(fp)                                    :: newthlyr
    real(fp)                                    :: thick
    real(fp)                  , pointer         :: theulyr
    integer                   , pointer         :: keuler
    integer                   , pointer         :: nlyr
    integer                   , pointer         :: updbaselyr
    real(fp), dimension(:,:)  , pointer         :: svfrac
    real(fp), dimension(:,:,:), pointer         :: msed
    real(fp), dimension(:,:)  , pointer         :: thlyr
!
!! executable statements -------------------------------------------------------
!
    theulyr     => this%settings%theulyr
    keuler      => this%settings%keuler
    nlyr        => this%settings%nlyr
    updbaselyr  => this%settings%updbaselyr
    svfrac      => this%state%svfrac
    msed        => this%state%msed
    thlyr       => this%state%thlyr
    !
    dz = dzini
    !
    ! find first (partially) filled underlayer
    !
    k = keuler
    do while (comparereal(thlyr(k, nm),0.0_fp)==0 .and. k<nlyr)
       k = k+1
    enddo
    !
    ! don't fill the last underlayer
    !
    if (k == nlyr) k = k-1
    !
    ! start filling upwards
    !
    do while ( k>=keuler .and. dz > 0.0_fp )
       if ( thlyr(k, nm) < theulyr ) then
          !
          ! sediment can be added to this layer
          !
          if ( dz > theulyr-thlyr(k, nm) ) then
             !
             ! not all sediment can be added to this layer
             !
             fac     = (theulyr-thlyr(k, nm))/dz
             dz      = dz*(1.0_fp-fac)
             do l = 1, this%settings%nfrac
                dm = dmi(l)*fac
                msed(l, k, nm) = msed(l, k, nm) + dm
                dmi(l)         = dmi(l)         - dm
             enddo
             svfrac(k, nm) = svfrac(k, nm)*thlyr(k, nm) + svfracdep*(theulyr-thlyr(k, nm))
             thlyr(k, nm)  = theulyr
             svfrac(k, nm) = svfrac(k, nm)/thlyr(k, nm)
          else
             !
             ! everything can be added to this layer
             !
             do l = 1, this%settings%nfrac            
                msed(l, k, nm) = msed(l, k, nm) + dmi(l)
                dmi(l) = 0.0_fp
             enddo
             svfrac(k, nm) = svfrac(k, nm)*thlyr(k, nm) + svfracdep*dz
             thlyr(k, nm)  = thlyr(k, nm) + dz
             svfrac(k, nm) = svfrac(k, nm)/thlyr(k, nm)
             dz            = 0.0_fp
          endif
       endif
       k = k-1
    enddo
    !
    ! the first Eulerian underlayer is completely filled or
    ! all sediment deposited
    !
    if ( dz > 0.0_fp ) then
       !
       ! still more sediment to be deposited
       !
       if (keuler == nlyr) then
          k = nlyr
          !
          ! no Eulerian underlayers, so put everything in
          ! the last (i.e. base) layer
          !
          select case (updbaselyr)
          case(1) ! compute separate composition for the base layer
             do l = 1, this%settings%nfrac           
                msed(l, k, nm) = msed(l, k, nm) + dmi(l)
             enddo
             svfrac(k, nm) = svfrac(k, nm)*thlyr(k, nm) + svfracdep*dz
             thlyr(k, nm)  = thlyr(k, nm) + dz
             svfrac(k, nm) = svfrac(k, nm)/thlyr(k, nm)
             dz            = 0.0_fp
          case(2) ! composition of base layer constant
             !
             ! composition of dz is lost, update thickness
             !
             fac = (thlyr(k, nm)+dz)/thlyr(k, nm)
             do l = 1, this%settings%nfrac
                msed(l, k, nm) = msed(l, k, nm)*fac
             enddo
             thlyr(k, nm) = thlyr(k, nm) + dz
          case(3) ! same as the (first non-empty) layer above it
             !
             ! composition of dz is lost, update thickness
             ! and set composition to that of layer nlyr-1
             !
             thlyr(k, nm) = thlyr(k, nm) + dz
             fac = thlyr(k, nm)/thlyr(k-1, nm)
             do l = 1, this%settings%nfrac
                msed(l, k, nm) = msed(l, kmin-1, nm)*fac
             enddo
          case default
             !
             ! ERROR
             !
          endselect
       else
          !
          ! thin underlayers filled shift administration and deposit
          ! remaining sediment
          !
          do while ( dz > 0.0_fp )
             !
             ! add last thin underlayer to the last (i.e. base) layer
             !
             newthlyr = thlyr(nlyr, nm)+thlyr(nlyr-1, nm)
             select case (updbaselyr)
             case(1) ! compute separate composition for the base layer
                if ( newthlyr > 0.0_fp ) then
                   do l = 1, this%settings%nfrac
                      msed(l, nlyr, nm) = msed(l, nlyr, nm) + msed(l, nlyr-1, nm) 
                   enddo
                   svfrac(nlyr, nm) = svfrac(nlyr, nm)*thlyr(nlyr, nm) + svfrac(nlyr-1, nm)*thlyr(nlyr-1, nm)
                   svfrac(nlyr, nm) = svfrac(nlyr, nm)/newthlyr
                endif
             case(2) ! composition of base layer constant
                !
                ! composition of layer nlyr-1 is lost; just the
                ! thickness of the base layer is updated
                !
                fac = newthlyr/thlyr(nlyr, nm)
                do l = 1, this%settings%nfrac
                    msed(l, nlyr, nm) = msed(l, nlyr, nm)*fac
                enddo
             case(3) ! same as the (first non-empty) layer above it
                !
                ! find lowest non-empty layer
                !
                do kne = nlyr-2,1,-1
                   if ( thlyr(kne, nm) > 0.0_fp ) exit
                enddo
                fac = newthlyr/thlyr(kne, nm)
                do l = 1, this%settings%nfrac
                   msed(l, nlyr, nm) = msed(l, kne, nm)*fac
                enddo
             case(4) ! composition and thickness of base layer constant
                !
                ! composition and sediment of layer nlyr-1 is lost
                ! make sure that the newthlyr of the base layer is equal
                ! to the old thickness
                !
                newthlyr = thlyr(nlyr, nm)
             case default
                !
                ! ERROR
                !
             end select
             thlyr(nlyr, nm) = newthlyr
             !
             ! shift layers down by one
             !
             do k = nlyr-1,keuler+1,-1
                do l = 1, this%settings%nfrac
                   msed(l, k, nm) = msed(l, k-1, nm)
                enddo
                thlyr(k, nm)  = thlyr(k-1, nm)
                svfrac(k, nm) = svfrac(k-1, nm)
             enddo
             !
             ! put all the sediment in one layer
             !
             k = keuler
             do l = 1, this%settings%nfrac 
                 msed(l, k, nm) = dmi(l)
             enddo
             thlyr(k, nm)  = dz
             svfrac(k, nm) = svfracdep
             dz           = 0.0_fp
          enddo
       endif
    endif
end subroutine lyrsedimentation_eulerian
!
!
!
!==============================================================================
subroutine lyrdiffusion(this, nm, dt, rhosol)
!!--description-----------------------------------------------------------------
!
!    Function:
!     - lyrdiffusion implements the mixing between the layers through diffusion
!       (porosity of layers is not updated yet!)
!
!!--declarations----------------------------------------------------------------
    use precision
    !
    implicit none
    !
    ! Function/routine arguments
    !
    type(bedcomp_data)                                    :: this    
    integer                                  , intent(in) :: nm  
    real(fp)                                 , intent(in) :: dt 
    real(fp), dimension(this%settings%nfrac) , intent(in) :: rhosol
    !
    ! Local variables
    !
    integer                                            :: k
    integer                                            :: km1
    integer                                            :: l
    integer                                            :: nd
    real(fp)                                           :: alpha
    real(fp)                                           :: Fdiff
    real(fp)                                           :: kd
    real(fp)                                           :: mtot
    real(fp)                                           :: mtot1
    real(fp)                                           :: pk
    real(fp)                                           :: pk_new
    real(fp)                                           :: pkm1_new
    real(fp)                                           :: thick
    real(fp)                                           :: thick1
    real(fp)                                           :: zd
    real(fp), dimension(this%settings%nfrac)           :: mkm1
    real(fp), dimension(this%settings%nfrac)           :: pkm1
    integer                                  , pointer :: idiffusion 
    integer                                  , pointer :: ndiff 
    integer                                  , pointer :: nlyr
    real(fp), dimension(:)                   , pointer :: rhofrac
    real(fp), dimension(:)                   , pointer :: zdiff
    real(fp), dimension(:,:)                 , pointer :: kdiff
    real(fp), dimension(:,:)                 , pointer :: svfrac
    real(fp), dimension(:,:)                 , pointer :: thlyr
    real(fp), dimension(:,:,:)               , pointer :: msed
!
!! executable statements -------------------------------------------------------
!
    idiffusion  => this%settings%idiffusion
    kdiff       => this%settings%kdiff
    ndiff       => this%settings%ndiff
    nlyr        => this%settings%nlyr
    rhofrac     => this%settings%rhofrac
    zdiff       => this%settings%zdiff
    svfrac      => this%state%svfrac
    msed        => this%state%msed
    thlyr       => this%state%thlyr  
    !
    km1     = 1
    !
    ! compute mass fractions in top layer
    !
    if (idiffusion==1) then
        ! diffusion based on mass fractions
        mtot   = 0.0_fp
        do l = 1, this%settings%nfrac
            mtot = mtot + msed(l, 1, nm)
        enddo 
        mtot = max(mtot,1e-20_fp)       
        do l = 1, this%settings%nfrac
            pkm1(l) = msed(l, 1, nm)/mtot
        enddo
        mtot1 = mtot
    elseif (idiffusion==2) then  
        ! diffusion based on volume fractions
        thick = svfrac(1, nm) * thlyr(1, nm)
        thick = max(thick,1e-20)
        do l = 1, this%settings%nfrac
            pkm1(l) = msed(l, 1, nm)/(rhofrac(l)*thick)
        enddo
        thick1 = thick
    endif    
    !     
    do l = 1, this%settings%nfrac
        mkm1(l) = msed(l, 1, nm)
    enddo
    !
    zd = thlyr(1, nm)   ! location of interface between the layers
    nd = 1              ! index of used diffusion coefficient
    !
    do k = 2, nlyr
        !
        ! apply mixing only between non-empty layers
        !
        if (comparereal(thlyr(k,nm),0.0_fp) == 0) cycle
        !
        ! Compute diffusion coefficient at interface between layers
        ! through linear interpolation
        ! Extrapolation is performed by assuming constant values
        !
        do while (nd < ndiff .and. zdiff(nd) < zd )
            nd = nd + 1
        enddo
        !
        if (nd == 1) then
            kd = kdiff(nd, nm)
        elseif (nd <= ndiff .and. zd < zdiff(ndiff)) then
            kd = kdiff(nd-1, nm) + (zd - zdiff(nd-1))/(zdiff(nd)-zdiff(nd-1)) * (kdiff(nd, nm) - kdiff(nd-1, nm))
        else
            kd = kdiff(nd, nm)
        endif
        !
        if (idiffusion==1) then
            mtot   = 0.0_fp
            do l = 1, this%settings%nfrac
                mtot = mtot + msed(l, k, nm)
            enddo
            mtot = max(mtot,1e-20_fp)
        elseif (idiffusion==2) then 
            thick = svfrac(k, nm) * thlyr(k, nm)
            thick = max(thick,1e-20)
        endif
        !
        !   Loop over all sediment fractions
        !        
        do l = 1, this%settings%nfrac
            !
            ! compute mass or volume fraction in layer k
            !
            if (idiffusion==1) then
                pk = msed(l, k, nm)/mtot
            elseif (idiffusion==2) then 
                pk = msed(l, k, nm)/(rhofrac(l)*thick)
            endif
            !
            ! Compute diffusion flux between layer k and km1
            ! based on difference in mass fraction
            !
            Fdiff = 0.5_fp*rhosol(l)*kd*(svfrac(km1, nm) + svfrac(k, nm))/min(thlyr(km1, nm),thlyr(k, nm)) * (pk - pkm1(l))
            !
            ! Diffusion flux is limited to half of the available sediment fraction in a layer
            ! to avoid sediment shortage
            !
            if ( Fdiff < 0.0_fp) then
                Fdiff = max( - msed(l, km1, nm)/(2.0*dt) , Fdiff)
            else
                Fdiff = min(   msed(l, k  , nm)/(2.0*dt) , Fdiff)
            endif
            !
            if (idiffusion==1) then
                pkm1_new = (msed(l, km1, nm) + Fdiff*dt)/mtot1
                pk_new   = (msed(l, k,   nm) - Fdiff*dt)/mtot
            elseif (idiffusion==2) then 
                pkm1_new = (msed(l, km1, nm) + Fdiff*dt)/(rhofrac(l)*thick1)
                pk_new   = (msed(l, k,   nm) - Fdiff*dt)/(rhofrac(l)*thick)
            endif
                
            if ( (pk-pkm1(l))*(pk_new-pkm1_new)<0 ) then
                if (idiffusion==1) then
                    alpha =            (pk-pkm1(l))*mtot* mtot1 /(Fdiff*dt*(mtot +mtot1 ))
                elseif (idiffusion==2) then 
                    alpha = rhofrac(l)*(pk-pkm1(l))*thick*thick1/(Fdiff*dt*(thick+thick1))
                endif    
                Fdiff = alpha * Fdiff
            endif
                        !
            pkm1(l) = pk
            mkm1(l) = msed(l, k, nm)
            !
            ! Update masses of sediment fractions
            !
            msed(l, km1, nm) = msed(l, km1, nm) + Fdiff*dt
            msed(l, k  , nm) = msed(l, k  , nm) - Fdiff*dt
            !
        enddo
        !
        mtot1   = mtot
        thick1  = thick
        km1 = k
        zd  = zd + thlyr(k, nm)
        !
    enddo
    !
    ! Update layer thicknesses
    !
    do k = 1, nlyr
        thick = 0.0_fp
        do l = 1, this%settings%nfrac
            thick = thick + msed(l, k, nm)/rhofrac(l)
        enddo
        thlyr(k, nm) = thick/svfrac(k, nm)
    enddo    
    !
end subroutine lyrdiffusion
!
!
!
!==============================================================================
subroutine detthcmud(this, thcmud)
!!--description-----------------------------------------------------------------
!
!    Function: - Determines the total thickness of mud fraction in the bed
!    DEPRECATED FUNCTIONALITY; use mudfrac (sum of frac over mud fractions) instead
!
!!--declarations----------------------------------------------------------------
    use precision
    !
    implicit none
    !
    include 'sedparams.inc'
!
! Function/routine arguments
!
    type(bedcomp_data)                                        , intent(in)  :: this
    real(fp), dimension(this%settings%nmlb:this%settings%nmub), intent(out) :: thcmud  !  Total thickness of the mud layers
!
! Local variables
!
    integer                             :: l
    integer                             :: nm
    real(prec), dimension(:,:), pointer :: bodsed
    real(fp)  , dimension(:)  , pointer :: rhofrac
!
!! executable statements -------------------------------------------------------
!
    bodsed      => this%state%bodsed
    rhofrac     => this%settings%rhofrac
    !
    do nm = this%settings%nmlb,this%settings%nmub
       thcmud (nm) = 0.0
        do l = 1, this%settings%nfrac
           if (this%settings%sedtyp(l) == SEDTYP_COHESIVE) then
              thcmud(nm) = thcmud(nm) + real(bodsed(l, nm),fp)/rhofrac(l)
           endif
        enddo
    enddo
end subroutine detthcmud
!
!
!
!==============================================================================
subroutine getalluvthick(this, seddep, nmfrom, nmto, nval)
!!--description-----------------------------------------------------------------
!
!    Function: - Determine sediment thickness optionally per sediment fraction
!    DEPRECATED FUNCTIONALITY; use getsedthick instead.
!
!!--declarations----------------------------------------------------------------
    use precision 
    !
    implicit none
    !
!
! Function/routine arguments
!
    type(bedcomp_data)                                                , intent(in)  :: this
    real(fp), dimension(nmfrom:nmto, nval)                            , intent(out) :: seddep
    integer                                                           , intent(in)  :: nmfrom
    integer                                                           , intent(in)  :: nmto
    integer                                                           , intent(in)  :: nval
!
! Local variables
!
    integer                              :: k
    integer                              :: l
    integer                              :: nm
    real(fp)                             :: thkl
    integer, pointer                     :: nlyr
    real(prec), dimension(:,:) , pointer :: bodsed
    real(fp)  , dimension(:,:) , pointer :: thlyr
    real(fp)                   , pointer :: thresh
    real(fp)  , dimension(:)   , pointer :: rhofrac
!
!! executable statements -------------------------------------------------------
!
    nlyr                => this%settings%nlyr
    rhofrac             => this%settings%rhofrac
    bodsed              => this%state%bodsed
    thlyr               => this%state%thlyr
    !
    select case(this%settings%iunderlyr)
    case(2)
       do nm = nmfrom,nmto
          thkl = 0.0_fp
          do k = 1, nlyr
             thkl = thkl + thlyr(k, nm)
          enddo
          do l = 1, nval
             seddep(nm, l) = thkl
          enddo
       enddo
    case default
       do nm = nmfrom,nmto
          if (nval==1) then
             seddep(nm, nval) = 0.0
             do l = 1, this%settings%nfrac
                seddep(nm, nval) = seddep(nm, nval) + real(bodsed(l, nm),fp)/rhofrac(l)
             enddo
          else
             do l = 1, this%settings%nfrac
                seddep(nm, l) = real(bodsed(l, nm),fp)/rhofrac(l)
             enddo
          endif
       enddo
    endselect
end subroutine getalluvthick
!
!
!
!==============================================================================
subroutine getfrac(this, frac, anymud, mudcnt, mudfrac, nmfrom, nmto)
!!--description-----------------------------------------------------------------
!
!    Function: Determines the (mass or volume) fractions
!
!!--declarations----------------------------------------------------------------
    use precision 
    implicit none
    include 'sedparams.inc'
    !
    ! Function/routine arguments
    !
    type(bedcomp_data)                                                , intent(in)  :: this    
    logical                                                           , intent(in)  :: anymud
    real(fp), dimension(nmfrom:nmto, this%settings%nfrac)             , intent(out) :: frac
    real(fp), dimension(nmfrom:nmto)                                  , intent(out) :: mudfrac
    real(fp), dimension(nmfrom:nmto)                                  , intent(in)  :: mudcnt
    integer                                                           , intent(in)  :: nmfrom
    integer                                                           , intent(in)  :: nmto
    !
    ! Local variables
    !
    integer  :: l
    integer  :: nm
    real(fp) :: nonmud
    !
    !! executable statements -------------------------------------------------------
    !
    !
    ! Calculate total bottom sediments and fractions
    !
    select case(this%settings%iunderlyr)
    case(1)
       call getmfrac(this ,frac, nmfrom, nmto)
    case default
       call getvfrac(this ,frac, nmfrom, nmto)
    endselect
    !
    ! Calculate mud fraction
    !
    if (anymud) then
       !
       ! Add simulated mud fractions.
       !
       mudfrac = 0.0
       do l = 1, this%settings%nfrac
          if (this%settings%sedtyp(l) == SEDTYP_COHESIVE) then
             do nm = nmfrom, nmto
                mudfrac(nm) = mudfrac(nm) + frac(nm,l)
             enddo
          endif
       enddo
    else
       !
       ! Take into account possible non-simulated
       ! mud fraction.
       !
       do nm = nmfrom, nmto
          mudfrac(nm) = mudcnt(nm)
          nonmud      = 1.0_fp - mudcnt(nm)
          do l = 1, this%settings%nfrac
             frac(nm,l) = frac(nm,l) * nonmud
          enddo
       enddo
    endif
end subroutine getfrac
!
!
!
!==============================================================================
subroutine getmfrac(this, frac, nmfrom, nmto)
!!--description-----------------------------------------------------------------
!
!    Function: Determines the mass fractions
!
!!--declarations----------------------------------------------------------------
    use precision 
    implicit none
    !
    ! Function/routine arguments
    !
    type(bedcomp_data)                                                , intent(in)  :: this
    real(fp), dimension(nmfrom:nmto, this%settings%nfrac)             , intent(out) :: frac
    integer                                                           , intent(in)  :: nmfrom
    integer                                                           , intent(in)  :: nmto
    !
    ! Local variables
    !
    integer                             :: l
    integer                             :: nm
    real(fp)                            :: sedtot
    real(prec), dimension(:,:), pointer :: bodsed
    real(fp), dimension(:,:,:), pointer :: msed
    !
    !! executable statements -------------------------------------------------------
    !
    bodsed      => this%state%bodsed
    msed        => this%state%msed
    !
    ! Calculate total bottom sediments and fractions
    !
    select case(this%settings%iunderlyr)
    case(2)
       do nm = nmfrom, nmto
          sedtot = 0.0_fp
          do l = 1, this%settings%nfrac
             sedtot = sedtot + msed(l, 1, nm)
          enddo
          if (comparereal(sedtot,0.0_fp) == 0) then
             frac(nm, :) = 1.0_fp/this%settings%nfrac
          else
            do l = 1, this%settings%nfrac
                frac(nm, l) = msed(l, 1, nm)/sedtot
            enddo
          endif
       enddo
    case default
       do nm = nmfrom, nmto
          sedtot = 0.0_fp
          do l = 1, this%settings%nfrac
             sedtot = sedtot + real(bodsed(l, nm),fp)
          enddo
          if (comparereal(sedtot,0.0_fp) == 0) then
             frac(nm, :) = 1.0_fp/this%settings%nfrac
          else
             do l = 1, this%settings%nfrac
                frac(nm, l) = real(bodsed(l, nm),fp)/sedtot
             enddo
          endif
       enddo
    endselect
end subroutine getmfrac
!
!
!
!==============================================================================
subroutine setmfrac(this, frac, nmfrom, nmto)
!!--description-----------------------------------------------------------------
!
!    Function: Update the bed composition given the mass fraction data
!
!!--declarations----------------------------------------------------------------
    use precision 
    implicit none
    !
    ! Function/routine arguments
    !
    type(bedcomp_data)                                                              :: this
    real(fp), dimension(nmfrom:nmto, this%settings%nfrac)             , intent(in)  :: frac
    integer                                                           , intent(in)  :: nmfrom
    integer                                                           , intent(in)  :: nmto
    !
    ! Local variables
    !
    integer                             :: l
    integer                             :: nm
    real(fp)                            :: sedtot
    real(prec), dimension(:,:), pointer :: bodsed
    real(fp), dimension(:,:,:), pointer :: msed
    !
    !! executable statements -------------------------------------------------------
    !
    bodsed      => this%state%bodsed
    msed        => this%state%msed
    !
    ! Update the sediment mass per fraction, but keep the total mass identical
    !
    select case(this%settings%iunderlyr)
    case(2)
       do nm = nmfrom, nmto
          sedtot = 0.0_fp
          do l = 1, this%settings%nfrac
             sedtot = sedtot + msed(l, 1, nm)
          enddo
          do l = 1, this%settings%nfrac
             msed(l, 1, nm) = frac(nm, l)*sedtot
          enddo
       enddo
    case default
       do nm = nmfrom, nmto
          sedtot = 0.0_fp
          do l = 1, this%settings%nfrac
             sedtot = sedtot + real(bodsed(l, nm),fp)
          enddo
          do l = 1, this%settings%nfrac
             bodsed(l, nm) = real(frac(nm, l)*sedtot,prec)
          enddo
       enddo
    endselect
end subroutine setmfrac
!
!
!
!==============================================================================
subroutine getvfrac(this, frac, nmfrom, nmto)
!!--description-----------------------------------------------------------------
!
!    Function: Determines the volume fractions
!
!!--declarations----------------------------------------------------------------
    use precision
    implicit none
    !
    ! Function/routine arguments
    !
    type(bedcomp_data)                                                , intent(in)  :: this
    real(fp), dimension(nmfrom:nmto, this%settings%nfrac)             , intent(out) :: frac
    integer                                                           , intent(in)  :: nmfrom
    integer                                                           , intent(in)  :: nmto
    !
    ! Local variables
    !
    integer                               :: l
    integer                               :: nm
    real(fp)                              :: thick
    real(prec), dimension(:,:)  , pointer :: bodsed
    real(fp)  , dimension(:)    , pointer :: dpsed
    real(fp)  , dimension(:,:)  , pointer :: svfrac
    real(fp)  , dimension(:,:,:), pointer :: msed
    real(fp)  , dimension(:,:)  , pointer :: thlyr
    real(fp)  , dimension(:)    , pointer :: rhofrac
!
!! executable statements -------------------------------------------------------
!
    rhofrac     => this%settings%rhofrac
    svfrac      => this%state%svfrac
    bodsed      => this%state%bodsed
    dpsed       => this%state%dpsed
    msed        => this%state%msed
    thlyr       => this%state%thlyr
    !
    ! Calculate total bottom sediments and fractions
    !
    select case(this%settings%iunderlyr)
    case(2)
       do nm = nmfrom, nmto
          if (comparereal(thlyr(1, nm),0.0_fp) == 0) then
             frac(nm, :) = 1.0_fp/this%settings%nfrac
          else
             thick = svfrac(1, nm) * thlyr(1, nm)
             do l = 1, this%settings%nfrac
                frac(nm, l) = msed(l, 1, nm)/(rhofrac(l)*thick)
             enddo
          endif
       enddo
    case default
       do nm = nmfrom, nmto
          if (comparereal(dpsed(nm),0.0_fp) == 0) then
             frac(nm, :) = 1.0_fp/this%settings%nfrac
          else
             do l = 1, this%settings%nfrac
                frac(nm, l) = real(bodsed(l, nm),fp)/(rhofrac(l)*dpsed(nm))
             enddo
          endif
       enddo
    endselect
end subroutine getvfrac
!
!
!
!==============================================================================
subroutine setvfrac(this, frac, nmfrom, nmto)
!!--description-----------------------------------------------------------------
!
!    Function: Update the bed composition given the volume fraction data
!
!!--declarations----------------------------------------------------------------
    use precision 
    implicit none
    !
    ! Function/routine arguments
    !
    type(bedcomp_data)                                                              :: this    
    real(fp), dimension(nmfrom:nmto, this%settings%nfrac)             , intent(in)  :: frac
    integer                                                           , intent(in)  :: nmfrom
    integer                                                           , intent(in)  :: nmto
    !
    ! Local variables
    !
    integer                                    :: l
    integer                                    :: nm
    real(fp)                                   :: sum
    real(fp)                                   :: sedtot
    real(fp)  , dimension(:,:)  , pointer      :: svfrac
    real(prec), dimension(:,:)  , pointer      :: bodsed
    real(fp)  , dimension(:,:,:), pointer      :: msed
    real(fp)  , dimension(:,:)  , pointer      :: thlyr
    real(fp)  , dimension(:)    , pointer      :: rhofrac
    !
    !! executable statements -------------------------------------------------------
    !
    rhofrac     => this%settings%rhofrac
    svfrac      => this%state%svfrac
    bodsed      => this%state%bodsed
    msed        => this%state%msed
    thlyr       => this%state%thlyr
    !
    ! Calculate total bottom sediments and fractions
    !
    select case(this%settings%iunderlyr)
    case(2)
       do nm = nmfrom, nmto
          sedtot = 0.0_fp
          do l = 1, this%settings%nfrac
             sedtot = sedtot + msed(l, 1, nm)
          enddo
          sum = 0.0_fp
          do l = 1, this%settings%nfrac
             sum = sum + frac(nm, l)*rhofrac(l)
          enddo
          do l = 1, this%settings%nfrac
             msed(l, 1, nm) = sedtot*(frac(nm, l)*rhofrac(l)/sum)
          enddo
       enddo
    case default
       do nm = nmfrom, nmto
          sedtot = 0.0_fp
          do l = 1, this%settings%nfrac
             sedtot = sedtot + real(bodsed(l, nm),fp)
          enddo
          sum = 0.0_fp
          do l = 1, this%settings%nfrac
             sum = sum + frac(nm, l)*rhofrac(l)
          enddo
          do l = 1, this%settings%nfrac
             bodsed(l, nm) = real(sedtot*(frac(nm, l)*rhofrac(l)/sum),prec)
          enddo
       enddo
    endselect
end subroutine setvfrac
!
!
!
!==============================================================================
subroutine getsedthick_allpoints(this, seddep)
!!--description-----------------------------------------------------------------
!
!    Function: Determines total thickness of sediment deposit at all points
!
!!--declarations----------------------------------------------------------------
    use precision 
    implicit none
    !
    ! Function/routine arguments
    !
    type(bedcomp_data)                                          , intent(in)  :: this 
    real(fp), dimension(this%settings%nmlb:this%settings%nmub)  , intent(out) :: seddep
    !
    ! Local variables
    !
    integer                           :: k
    integer                           :: nm
    real(fp), dimension(:), pointer   :: dpsed
    real(fp), dimension(:,:), pointer :: thlyr
    !
    !! executable statements -------------------------------------------------------
    !
    dpsed       => this%state%dpsed
    thlyr       => this%state%thlyr
    !
    ! Calculate total bottom sediments and fractions
    !
    select case(this%settings%iunderlyr)
    case(2)
        seddep = 0.0_fp
        do nm = this%settings%nmlb, this%settings%nmub
            do k = 1, this%settings%nlyr
                seddep(nm) = seddep(nm) + thlyr(k, nm)
            enddo
        enddo
    case default
       seddep = dpsed
    endselect
end subroutine getsedthick_allpoints
!
!
!
!==============================================================================
subroutine getsedthick_1point(this, nm, seddep)
!!--description-----------------------------------------------------------------
!
!    Function: Determines total thickness of sediment deposit at one point
!
!!--declarations----------------------------------------------------------------
    use precision 
    implicit none
    !
    ! Function/routine arguments
    !
    type(bedcomp_data)                      , intent(in)  :: this    
    integer                                 , intent(in)  :: nm
    real(fp)                                , intent(out) :: seddep
    !
    ! Local variables
    !
    integer                           :: k
    real(fp), dimension(:)  , pointer :: dpsed
    real(fp), dimension(:,:), pointer :: thlyr
    !
    !! executable statements -------------------------------------------------------
    !
    dpsed       => this%state%dpsed
    thlyr       => this%state%thlyr
    !
    ! Calculate total bottom sediments and fractions
    !
    select case(this%settings%iunderlyr)
    case(2)
       seddep = 0.0_fp
       do k = 1, this%settings%nlyr
          seddep = seddep + thlyr(k, nm)
       enddo
    case default
       seddep = dpsed(nm)
    endselect
end subroutine getsedthick_1point
!
!
!
!==============================================================================
function initmorlyr(this) result (istat)
!!--description-----------------------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    implicit none
    !
    real(fp), parameter :: rmissval = -999.0_fp
    !
    ! Function/routine arguments
    !
    type (bedcomp_data), intent(inout) :: this    
    integer                            :: istat
    !
    ! Local variables
    !
    type (bedcomp_settings), pointer   :: settings
    type (bedcomp_state   ), pointer   :: state
    !
    !! executable statements -------------------------------------------------------
    !
    istat = 0
    if (istat == 0) allocate (settings, stat = istat)
    if (istat == 0) allocate (state   , stat = istat)
    if (istat /= 0) then
       !error
       return
    endif
    !
    allocate (settings%morlyrnum , stat = istat)
    if (istat == 0) then
       settings%morlyrnum%MinMassShortWarning = 0.0_fp
       settings%morlyrnum%MaxNumShortWarning = 100
    endif
    !
    settings%keuler     = 2
    settings%ndiff      = 0
    settings%nfrac      = 0
    settings%nlyr       = 0
    settings%nmlb       = 0
    settings%nmub       = 0
    settings%idiffusion = 0
    settings%iunderlyr  = 1
    settings%iporosity  = 0
    settings%flufflyr   = 0
    settings%exchlyr    = .false.
    settings%neulyr     = 0
    settings%nlalyr     = 0
    settings%theulyr    = rmissval
    settings%thlalyr    = rmissval
    settings%updbaselyr = 1
    !
    nullify(settings%bfluff0)
    nullify(settings%bfluff1)
    nullify(settings%kdiff)
    nullify(settings%thexlyr)
    nullify(settings%thtrlyr)
    nullify(settings%zdiff)
    !
    nullify(state%preload)
    nullify(state%svfrac)
    nullify(state%bodsed)
    nullify(state%dpsed)
    nullify(state%fluffshort)
    nullify(state%mfluff)
    nullify(state%msed)
    nullify(state%thlyr)
    nullify(state%sedshort)
    !
    this%settings => settings
    this%state    => state
end function initmorlyr
!
!
!
!==============================================================================
function allocmorlyr(this) result (istat)
!!--description-----------------------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    implicit none
    !
    ! Function/routine arguments
    !
    type (bedcomp_data)             :: this    
    integer                         :: istat
    !
    ! Local variables
    !
    type (bedcomp_settings), pointer :: settings
    type (bedcomp_state), pointer    :: state
    integer                          :: nmlb
    integer                          :: nmub
    integer                          :: nfrac
    !
    !! executable statements -------------------------------------------------------
    !
    settings => this%settings
    state => this%state
    !
    ! Number of layers: 1       transport layer
    !                   1       exchange layer (optional)
    !                   nlalyr  lagrangian underlayers
    !                   neulyr  eulerian underlayers
    !                   1       persistent base layer
    !
    settings%nlyr   = 2 + settings%nlalyr + settings%neulyr
    settings%keuler = 2 + settings%nlalyr
    if (settings%exchlyr) then
       settings%nlyr   = settings%nlyr + 1
       settings%keuler = settings%keuler + 1
    endif
    !
    nmlb  = settings%nmlb
    nmub  = settings%nmub
    nfrac = settings%nfrac
    !
    istat = 0
    if (istat == 0) allocate (state%bodsed(nfrac,nmlb:nmub), stat = istat)
    if (istat == 0) state%bodsed = 0.0_fp
    if (istat == 0) allocate (state%dpsed(nmlb:nmub), stat = istat)
    if (istat == 0) state%dpsed = 0.0_fp
    !
    if (istat == 0) allocate (settings%sedtyp(nfrac), stat = istat)
    if (istat == 0) settings%sedtyp = 0
    if (istat == 0) allocate (settings%rhofrac(nfrac), stat = istat)
    if (istat == 0) settings%rhofrac = 0.0_fp
    if (istat == 0) allocate (settings%phi(nfrac), stat = istat)
    if (istat == 0) settings%phi = 0.0_fp
    if (istat == 0) allocate (settings%sigphi(nfrac), stat = istat)
    if (istat == 0) settings%sigphi = 0.0_fp
    !
    if (settings%iunderlyr==2) then
       if (istat == 0) allocate (settings%kdiff(settings%ndiff,nmlb:nmub), stat = istat)
       if (istat == 0) settings%kdiff = 0.0_fp
       if (istat == 0) allocate (settings%zdiff(settings%ndiff), stat = istat)
       if (istat == 0) settings%zdiff = 0.0_fp
       if (istat == 0) allocate (settings%thtrlyr(nmlb:nmub), stat = istat)
       if (istat == 0) settings%thtrlyr = 0.0_fp
       if (settings%exchlyr) then
          if (istat == 0) allocate (settings%thexlyr(nmlb:nmub), stat = istat)
          if (istat == 0) settings%thexlyr = 0.0_fp
       endif
       if (istat == 0) allocate (state%msed(nfrac,settings%nlyr,nmlb:nmub), stat = istat)
       if (istat == 0) state%msed = 0.0_fp
       if (istat == 0) allocate (state%thlyr(settings%nlyr,nmlb:nmub), stat = istat)
       if (istat == 0) state%thlyr = 0.0_fp
       if (istat == 0) allocate (state%sedshort(nfrac,nmlb:nmub), stat = istat)
       if (istat == 0) state%sedshort = 0.0_fp
       if (istat == 0) allocate (state%svfrac(settings%nlyr,nmlb:nmub), stat = istat)
       if (istat == 0) state%svfrac = 1.0_fp
       if (istat == 0) allocate (state%preload(settings%nlyr,nmlb:nmub), stat = istat)
       if (istat == 0) state%preload = 0.0_fp
    endif
    !
    if (settings%flufflyr>0) then
        if (istat == 0) allocate (state%fluffshort(nfrac,nmlb:nmub), stat = istat)
        if (istat == 0) state%fluffshort = 0.0_fp
        if (istat == 0) allocate (state%mfluff(nfrac,nmlb:nmub), stat = istat)
        if (istat == 0) state%mfluff = 0.0_fp
    endif
    !
    if (settings%flufflyr==1) then
      if (istat == 0) allocate (settings%bfluff0(nfrac,nmlb:nmub), stat = istat)
      if (istat == 0) settings%bfluff0 = 0.0_fp
      if (istat == 0) allocate (settings%bfluff1(nfrac,nmlb:nmub), stat = istat)
      if (istat == 0) settings%bfluff1 = 0.0_fp
    endif
    !
    !
end function allocmorlyr
!
!
!
!==============================================================================
function allocwork(this, work) result (istat)
!!--description-----------------------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    implicit none
    !
    ! Function/routine arguments
    !
    type (bedcomp_data), intent(in)  :: this    
    type (bedcomp_work), intent(out) :: work
    integer                          :: istat
    !
    ! Local variables
    !
    integer, pointer :: nfrac
    integer, pointer :: nlyr
    !
    !! executable statements -------------------------------------------------------
    !
    nfrac => this%settings%nfrac
    nlyr  => this%settings%nlyr
    !
    istat = 0
    if (istat == 0) allocate (work%msed2(nfrac, nlyr), stat = istat)
    if (istat == 0) allocate (work%thlyr2(nlyr)     , stat = istat)
    if (istat == 0) allocate (work%svfrac2(nlyr)    , stat = istat)
end function allocwork
!
!
!
!==============================================================================
function deallocwork(this, work) result (istat)
!!--description-----------------------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    implicit none
    !
    ! Function/routine arguments
    !
    type (bedcomp_data), intent(in)  :: this    
    type (bedcomp_work), intent(out) :: work
    integer                          :: istat
    !
    ! Local variables
    !
    !
    !! executable statements -------------------------------------------------------
    !
    istat = 0
    if (istat == 0) deallocate (work%msed2  , stat = istat)
    if (istat == 0) deallocate (work%thlyr2 , stat = istat)
    if (istat == 0) deallocate (work%svfrac2, stat = istat)
end function deallocwork
!
!
!
!==============================================================================
function clrmorlyr(this) result (istat)
!!--description-----------------------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    implicit none
    !
    ! Function/routine arguments
    !
    type (bedcomp_data)             :: this    
    integer                         :: istat
    !
    ! Local variables
    !
    type (bedcomp_settings), pointer :: settings
    type (bedcomp_state)   , pointer :: state
    !
    !! executable statements -------------------------------------------------------
    !
    if (associated(this%settings)) then
       settings => this%settings
       if (associated(settings%bfluff0))   deallocate(settings%bfluff0  , STAT = istat)
       if (associated(settings%bfluff1))   deallocate(settings%bfluff1  , STAT = istat)
       if (associated(settings%kdiff))     deallocate(settings%kdiff    , STAT = istat)
       if (associated(settings%morlyrnum)) deallocate(settings%morlyrnum, STAT = istat)
       if (associated(settings%thexlyr))   deallocate(settings%thexlyr  , STAT = istat)
       if (associated(settings%thtrlyr))   deallocate(settings%thtrlyr  , STAT = istat)
       if (associated(settings%zdiff))     deallocate(settings%zdiff    , STAT = istat)
       !
       if (associated(settings%sedtyp))    deallocate(settings%sedtyp , STAT = istat)
       if (associated(settings%phi))       deallocate(settings%phi    , STAT = istat)
       if (associated(settings%rhofrac))   deallocate(settings%rhofrac, STAT = istat)
       if (associated(settings%sigphi))    deallocate(settings%sigphi , STAT = istat)
       !
       deallocate(this%settings, STAT = istat)
       nullify(this%settings)
    endif
    !
    if (associated(this%state)) then
       state => this%state
       if (associated(state%svfrac))       deallocate(state%svfrac      , STAT = istat)
       if (associated(state%bodsed))       deallocate(state%bodsed      , STAT = istat)
       if (associated(state%dpsed))        deallocate(state%dpsed       , STAT = istat)
       if (associated(state%fluffshort))   deallocate(state%fluffshort  , STAT = istat)
       if (associated(state%mfluff))       deallocate(state%mfluff      , STAT = istat)
       if (associated(state%msed))         deallocate(state%msed        , STAT = istat)
       if (associated(state%thlyr))        deallocate(state%thlyr       , STAT = istat)
       if (associated(state%sedshort))     deallocate(state%sedshort    , STAT = istat)
       if (associated(state%preload))      deallocate(state%preload     , STAT = istat)
       !
       deallocate(this%state, STAT = istat)
       nullify(this%state)
    endif
end function clrmorlyr
!
!
!
!==============================================================================
subroutine setbedfracprop(this, sedtyp, sedd50, logsedsig, rhofrac)
!!--description-----------------------------------------------------------------
!
!    Function: - Set sediment fraction properties
!
!!--declarations----------------------------------------------------------------
    use string_module
    implicit none
    !
    ! Call variables
    !
    type(bedcomp_data)                  :: this    
    integer , dimension(:), intent(in)  :: sedtyp
    real(fp), dimension(:), intent(in)  :: sedd50
    real(fp), dimension(:), intent(in)  :: logsedsig
    real(fp), dimension(:), intent(in)  :: rhofrac
    !
    ! Local variables
    !
    integer :: l
    !
    !! executable statements -------------------------------------------------------
    !
    do l = 1, this%settings%nfrac
       this%settings%sedtyp(l) = sedtyp(l)
       if (sedd50(l)<=0.0_fp) then
          this%settings%phi(l) = 13.9_fp ! -log(65um)/log(2)
       else
          this%settings%phi(l)     = -log(sedd50(l))/log(2.0_fp)
       endif
       if (logsedsig(l)<=0.0_fp) then
          this%settings%sigphi(l)  = log(1.34) ! use default for "well sorted" sediment (see rdsed.f90)
       else
          this%settings%sigphi(l)  = logsedsig(l)/log(2.0_fp)
       endif
       this%settings%rhofrac(l) = rhofrac(l) ! either rhosol or cdryb
    enddo
end subroutine
!
!
!
!==============================================================================
function bedcomp_getpointer_logical_scalar(this, variable, val) result (istat)
!!--description-----------------------------------------------------------------
!
!    Function: - Get a scalar logical
!
!!--declarations----------------------------------------------------------------
    use string_module
    implicit none
    !
    ! Call variables
    !
    type(bedcomp_data)    , intent(in)  :: this    
    character(*)          , intent(in)  :: variable
    logical, pointer                    :: val
    integer                             :: istat
    !
    ! Local variables
    !
    character(len(variable))    :: localname
    !
    !! executable statements -------------------------------------------------------
    !
    istat = 0
    localname = variable
    call str_lower(localname)
    select case (localname)
    case ('exchange_layer','exchlyr')
       val => this%settings%exchlyr
    case default
       val => NULL()
    end select
    if (.not.associated(val)) istat = -1
end function
!
!
!
!==============================================================================
function bedcomp_getpointer_integer_scalar(this, variable, val) result (istat)
!!--description-----------------------------------------------------------------
!
!    Function: - Get a scalar integer
!
!!--declarations----------------------------------------------------------------
    use string_module
    implicit none
    !
    ! Call variables
    !
    type(bedcomp_data)    , intent(in)  :: this    
    character(*)          , intent(in)  :: variable
    integer, pointer                    :: val
    integer                             :: istat
    !
    ! Local variables
    !
    character(len(variable))    :: localname
    !
    !! executable statements -------------------------------------------------------
    !
    istat = 0
    localname = variable
    call str_lower(localname)
    select case (localname)
    case ('flufflayer_model_type','flufflyr')
       val => this%settings%flufflyr
    case ('diffusion_model_type','idiffusion')
       val => this%settings%idiffusion
    case ('bed_layering_type','iunderlyr')
       val => this%settings%iunderlyr
    case ('porosity_model_type','iporosity')
       val => this%settings%iporosity
    case ('keuler')
       val => this%settings%keuler
    case ('number_of_diffusion_values','ndiff')
       val => this%settings%ndiff
    case ('number_of_layers','nlyr')
       val => this%settings%nlyr
    case ('maxnumshortwarning')
       val => this%settings%morlyrnum%MaxNumShortWarning
    case ('number_of_eulerian_layers','neulyr')
       val => this%settings%neulyr
    case ('number_of_lagrangian_layers','nlalyr')
       val => this%settings%nlalyr
    case ('number_of_fractions','nfrac')
       val => this%settings%nfrac
    case ('first_column_number','nmlb')
       val => this%settings%nmlb
    case ('last_column_number','nmub')
       val => this%settings%nmub
    case ('base_layer_updating_type','updbaselyr')
       val => this%settings%updbaselyr
    case default
       val => NULL()
    end select
    if (.not.associated(val)) istat = -1
end function
!
!
!
!==============================================================================
function bedcomp_getpointer_fp_scalar(this, variable, val) result (istat)
!!--description-----------------------------------------------------------------
!
!    Function: - Get a scalar real(fp)
!
!!--declarations----------------------------------------------------------------
    use precision
    use string_module
    implicit none
    !
    ! Call variables
    !
    type(bedcomp_data)    , intent(in)  :: this    
    character(*)          , intent(in)  :: variable
    real(fp), pointer                   :: val
    integer                             :: istat
    !
    ! Local variables
    !
    character(len(variable))    :: localname
    !
    !! executable statements -------------------------------------------------------
    !
    istat = 0
    localname = variable
    call str_lower(localname)
    select case (localname)
    case ('thickness_of_eulerian_layers','theulyr')
       val => this%settings%theulyr
    case ('thickness_of_lagrangian_layers','thlalyr')
       val => this%settings%thlalyr
    case ('minmassshortwarning')
       val => this%settings%morlyrnum%MinMassShortWarning
    case default
       val => NULL()
    end select
    if (.not.associated(val)) istat = -1
end function
!
!
!
!==============================================================================
function bedcomp_getpointer_fp_1darray(this, variable, val) result (istat)
!!--description-----------------------------------------------------------------
!
!    Function: - Get a pointer to a 1D real(fp) array
!
!!--declarations----------------------------------------------------------------
    use precision
    use string_module
    implicit none
    !
    ! Call variables
    !
    type(bedcomp_data)             , intent(in)  :: this    
    character(*)                   , intent(in)  :: variable
    real(fp), dimension(:), pointer              :: val
    integer                                      :: istat
    !
    ! Local variables
    !
    character(len(variable))    :: localname
    !
    !! executable statements -------------------------------------------------------
    !
    istat = 0
    localname = variable
    call str_lower(localname)
    select case (localname)
    case ('total_sediment_thickness','dpsed')
       val => this%state%dpsed
    case ('sediment_density')
       val => this%settings%rhofrac
    case ('thickness_of_exchange_layer','thexlyr')
       val => this%settings%thexlyr
    case ('thickness_of_transport_layer','thtrlyr')
       val => this%settings%thtrlyr
    case ('diffusion_levels','zdiff')
       val => this%settings%zdiff
    case default
       val => NULL()
    end select
    if (.not.associated(val)) istat = -1    
end function
!
!
!
!==============================================================================
function bedcomp_getpointer_fp_2darray(this, variable, val) result (istat)
!!--description-----------------------------------------------------------------
!
!    Function: - Get a pointer to a 2D real(fp) array
!
!!--declarations----------------------------------------------------------------
    use precision
    use string_module
    implicit none
    !
    ! Call variables
    !
    type(bedcomp_data)               , intent(in)  :: this    
    character(*)                     , intent(in)  :: variable
    real(fp), dimension(:,:), pointer              :: val
    integer                                        :: istat
    !
    ! Local variables
    !
    character(len(variable))    :: localname
    !
    !! executable statements -------------------------------------------------------
    !
    istat = 0
    localname = variable
    call str_lower(localname)
    select case (localname)
    case ('bfluff0')
       val => this%settings%bfluff0
    case ('bfluff1')
       val => this%settings%bfluff1
    case ('diffusion_coefficients','kdiff')
       val => this%settings%kdiff
    case ('flufflayer_mass','mfluff')
       val => this%state%mfluff
    case ('solid_volume_fraction','svfrac')
       val => this%state%svfrac
    case ('layer_thickness','thlyr')
       val => this%state%thlyr
    case default
       val => NULL()
    end select
    if (.not.associated(val)) istat = -1
end function
!
!
!
!==============================================================================
function bedcomp_getpointer_fp_3darray(this, variable, val) result (istat)
!!--description-----------------------------------------------------------------
!
!    Function: - Get a pointer to a 3D real(fp) array
!
!!--declarations----------------------------------------------------------------
    use precision
    use string_module
    implicit none
    !
    ! Call variables
    !
    type(bedcomp_data)                 , intent(in)  :: this    
    character(*)                       , intent(in)  :: variable
    real(fp), dimension(:,:,:), pointer              :: val
    integer                                          :: istat
    !
    ! Local variables
    !
    character(len(variable))    :: localname
    !
    !! executable statements -------------------------------------------------------
    !
    istat = 0
    localname = variable
    call str_lower(localname)
    select case (localname)
    case ('layer_mass','msed')
       val => this%state%msed
    case default
       val => NULL()
    end select
    if (.not.associated(val)) istat = -1
end function
!
!
!
!==============================================================================
function bedcomp_getpointer_prec_2darray(this, variable, val) result (istat)
!!--description-----------------------------------------------------------------
!
!    Function: - Get a pointer to a 2D real(prec) array
!
!!--declarations----------------------------------------------------------------
    use precision
    use string_module
    implicit none
    !
    ! Call variables
    !
    type(bedcomp_data)                 , intent(in)  :: this    
    character(*)                       , intent(in)  :: variable
    real(prec), dimension(:,:), pointer              :: val
    integer                                          :: istat
    !
    ! Local variables
    !
    character(len(variable))    :: localname
    !
    !! executable statements -------------------------------------------------------
    !
    istat = 0
    localname = variable
    call str_lower(localname)
    select case (localname)
    case ('total_sediment_mass','bodsed')
       val => this%state%bodsed
    case default
       val => NULL()
    end select
    if (.not.associated(val)) istat = -1
end function
!
!
!
!==============================================================================
subroutine bedcomp_use_bodsed(this)
!!--description-----------------------------------------------------------------
!
!    Function: - Use the values of BODSED to compute other quantities
!
!!--declarations----------------------------------------------------------------
    use precision
    implicit none
    !
    ! Call variables
    !
    type(bedcomp_data)                         :: this
    !
    ! Local variables
    !
    integer                                    :: ised
    integer                                    :: k
    integer                                    :: kstart
    integer                                    :: nm
    real(fp)                                   :: fac
    real(fp)  , dimension(this%settings%nfrac) :: mfrac
    real(fp)                                   :: poros
    real(fp)                                   :: sedthick
    real(fp)                                   :: sedthicklim
    real(fp)                                   :: svf
    real(fp)                                   :: thsed
    real(fp)                                   :: totsed
    real(prec), dimension(:,:)       , pointer :: bodsed
    real(fp)  , dimension(:)         , pointer :: dpsed
    real(fp)  , dimension(:,:,:)     , pointer :: msed
    real(fp)  , dimension(:,:)       , pointer :: svfrac
    real(fp)  , dimension(:,:)       , pointer :: thlyr
    real(fp)  , dimension(:)         , pointer :: thtrlyr
    real(fp)  , dimension(:)         , pointer :: thexlyr
    real(fp)  , dimension(:)         , pointer :: rhofrac
    !
    !! executable statements -------------------------------------------------------
    thtrlyr    => this%settings%thtrlyr
    thexlyr    => this%settings%thexlyr
    rhofrac    => this%settings%rhofrac
    bodsed     => this%state%bodsed
    dpsed      => this%state%dpsed
    msed       => this%state%msed
    svfrac     => this%state%svfrac
    thlyr      => this%state%thlyr
    !
    ! Fill initial values of DPSED
    !
    do nm = this%settings%nmlb, this%settings%nmub
       !
       ! Compute thickness correctly in case of
       ! multiple fractions with different rhofrac
       !
       dpsed(nm) = 0.0_fp
       do ised = 1, this%settings%nfrac
          dpsed(nm) = dpsed(nm) + real(bodsed(ised, nm),fp)/rhofrac(ised)
       enddo
    enddo
    select case(this%settings%iunderlyr)
    case(2)
       !
       ! No file specified for initial bed composition: extract data from
       ! the BODSED data read above.
       !
       msed = 0.0_fp
       thlyr = 0.0_fp
       do nm = this%settings%nmlb, this%settings%nmub
          !if (kcs(nm)<1 .or. kcs(nm)>2) cycle  !TODO: find a solution for this line
          !
          ! nm = (m-1)*nmax + n
          !
          totsed = 0.0_fp
          do ised = 1, this%settings%nfrac
             totsed = totsed + real(bodsed(ised, nm),fp)
          enddo
          totsed         = max(totsed,1.0e-20_fp) ! avoid division by zero
          do ised = 1, this%settings%nfrac
             mfrac(ised) = real(bodsed(ised, nm),fp)/totsed
          enddo
          !
          call getporosity(this, mfrac, poros)
          svf = 1.0_fp - poros
          !
          thsed = 0.0_fp
          do ised = 1, this%settings%nfrac
             thsed = thsed + real(bodsed(ised, nm),fp)/rhofrac(ised)
          enddo
          thsed         = thsed/svf
          sedthick      = thsed
          thsed         = max(thsed,1.0e-20_fp) ! avoid division by zero
          !
          ! active/transport layer
          !
          thlyr(1, nm)  = min(thtrlyr(nm),sedthick)
          fac = thlyr(1, nm)/thsed
          do ised = 1, this%settings%nfrac
             msed(ised, 1, nm) = real(bodsed(ised, nm),fp)*fac
          enddo
          svfrac(1, nm) = svf
          sedthick      = sedthick - thlyr(1, nm)
          !
          ! exchange layer
          !
          kstart        = 1
          if (this%settings%exchlyr) then
             kstart       = 2
             thlyr(2, nm) = min(thexlyr(nm),sedthick)
             sedthick     = sedthick - thlyr(2, nm)
             fac = thlyr(2, nm)/thsed
             do ised = 1, this%settings%nfrac
                msed(ised, 2, nm) = real(bodsed(ised, nm),fp)*fac
             enddo
             svfrac(2, nm) = svf
          endif
          !
          ! Lagrangian layers
          !
          do k = kstart+1, this%settings%keuler-1
             thlyr(k, nm) = min(this%settings%thlalyr,sedthick)
             sedthick     = sedthick - thlyr(k, nm)
             fac = thlyr(k, nm)/thsed
             do ised = 1, this%settings%nfrac
                msed(ised, k, nm) = real(bodsed(ised, nm),fp)*fac
             enddo
             svfrac(k, nm) = svf
          enddo
          !
          ! Eulerian layers
          !
          do k = this%settings%keuler,this%settings%nlyr-1
             thlyr(k, nm) = min(this%settings%theulyr,sedthick)
             sedthick     = sedthick - thlyr(k, nm)
             fac = thlyr(k, nm)/thsed
             do ised = 1, this%settings%nfrac
                msed(ised, k, nm) = real(bodsed(ised, nm),fp)*fac
             enddo
             svfrac(k, nm) = svf
          enddo
          !
          ! base layer
          !
          thlyr(this%settings%nlyr, nm) = sedthick
          fac = thlyr(this%settings%nlyr, nm)/thsed
          do ised = 1, this%settings%nfrac
             msed(ised, this%settings%nlyr, nm) = real(bodsed(ised, nm),fp)*fac
          enddo
          svfrac(this%settings%nlyr, nm) = svf
       enddo
    case default
       !
       ! nothing to do, using bodsed as uniformly mixed sediment
       !
    endselect
end subroutine
!
!
!
!==============================================================================
subroutine copybedcomp(this, nmfrom, nmto)
!!--description-----------------------------------------------------------------
!
!    Function: - Copy the bed composition from nmfrom to nmto
!
!!--declarations----------------------------------------------------------------
    use precision
    implicit none
    !
    ! Call variables
    !
    type(bedcomp_data)                   :: this
    integer                , intent(in)  :: nmfrom
    integer                , intent(in)  :: nmto
    !
    ! Local variables
    !
    integer                               :: k
    integer                               :: l
    real(prec), dimension(:,:)  , pointer :: bodsed
    real(fp)  , dimension(:)    , pointer :: dpsed
    real(fp)  , dimension(:,:)  , pointer :: mfluff
    real(fp)  , dimension(:,:,:), pointer :: msed
    real(fp)  , dimension(:,:)  , pointer :: svfrac
    real(fp)  , dimension(:,:)  , pointer :: thlyr
    !
    !! executable statements -------------------------------------------------------
    bodsed     => this%state%bodsed
    dpsed      => this%state%dpsed
    mfluff     => this%state%mfluff
    msed       => this%state%msed
    svfrac     => this%state%svfrac
    thlyr      => this%state%thlyr
    !
    do l = 1, this%settings%nfrac
        mfluff(l, nmto) = mfluff(l, nmfrom)
    enddo
    !
    select case(this%settings%iunderlyr)
    case(2)
       do k = 1, this%settings%nlyr
          do l = 1, this%settings%nfrac
             msed(l, k, nmto) = msed(l, k, nmfrom)
          enddo
          thlyr(k, nmto)  = thlyr(k, nmfrom)
          svfrac(k, nmto) = svfrac(k, nmfrom)
       enddo
    case default
       do l = 1, this%settings%nfrac
          bodsed(l, nmto) = bodsed(l, nmfrom)
       enddo
       dpsed(nmto) = dpsed(nmfrom)
    end select
end subroutine
!
!
!
!==============================================================================
subroutine updateporosity(this, nm, k)
!!--description-----------------------------------------------------------------
!
!    Function: - Update the porosity for layer k in column nm
!
!!--declarations----------------------------------------------------------------
    use precision
    implicit none
    !
    ! Call variables
    !
    type(bedcomp_data)                   :: this
    integer                , intent(in)  :: nm
    integer                , intent(in)  :: k
    !
    ! Local variables
    !
    integer                                   :: l
    real(fp)                                  :: poros
    real(fp)                                  :: totmass
    real(fp) , dimension(:,:,:), pointer      :: msed
    real(fp) , dimension(:,:)  , pointer      :: svfrac
    real(fp) , dimension(:,:)  , pointer      :: thlyr
    real(fp) , dimension(this%settings%nfrac) :: mfrac
    !
    !! executable statements -------------------------------------------------------
    msed       => this%state%msed
    svfrac     => this%state%svfrac
    thlyr      => this%state%thlyr
    !
    select case(this%settings%iunderlyr)
    case(2)
       totmass       = 0.0_fp
       do l = 1, this%settings%nfrac
          totmass    = totmass + msed(l, k, nm)
       enddo
       if (totmass>0.0_fp) then
          do l = 1, this%settings%nfrac
            mfrac(l)   = msed(l, k, nm)/totmass
          enddo
          !
          call getporosity(this, mfrac, poros)
       else
          poros = 0.0_fp
       endif
       svfrac(k, nm) = 1.0_fp - poros
    case default
       ! option not available for this bed composition model
    end select
end subroutine
!
!
!
!==============================================================================
subroutine getporosity(this, mfrac, poros)
!!--description-----------------------------------------------------------------
!
!    Function: - Compute the porosity
!
!!--declarations----------------------------------------------------------------
    use precision
    implicit none
    !
    ! Call variables
    !
    type(bedcomp_data)                       , intent(in)  :: this
    real(fp) , dimension(this%settings%nfrac), intent(in)  :: mfrac
    real(fp)                                 , intent(out) :: poros
    !
    ! Local variables
    !
    integer                              :: l
    real(fp)                             :: a
    real(fp)                             :: b
    real(fp)                             :: phim
    real(fp)                             :: sigmix
    real(fp)                             :: x
    real(fp) , dimension(:)    , pointer :: phi
    real(fp) , dimension(:)    , pointer :: sigphi
    !
    !! executable statements -------------------------------------------------------
    phi        => this%settings%phi
    sigphi     => this%settings%sigphi
    !
    phim = 0.0_fp
    do l = 1, this%settings%nfrac
       phim   = phim + phi(l)*mfrac(l)
    enddo
    sigmix = 0.0_fp
    do l = 1, this%settings%nfrac
       sigmix = sigmix + mfrac(l)*((phi(l)-phim)**2 + sigphi(l)**2)
    enddo
    sigmix = sqrt(sigmix)
    !
    select case (this%settings%iporosity)
    case (1)
       !
       ! R. Frings (May 2009)
       !
       a = -0.06_fp
       b = 0.36_fp
       poros = max(0.0_fp,a*sigmix + b)
    case (2)
       !
       ! G.J. Weltje based on data by Beard & Weyl (AAPG Bull., 1973)
       !
       x             = 3.7632_fp * sigmix**(-0.7552_fp)
       poros         = 0.45_fp*x/(1+x)
    case default
       poros         = 0.0_fp
    end select
end subroutine

subroutine consolidate(this, nm)
!!--description-----------------------------------------------------------------
!
!    Function: - Consolidate the bed of column nm
!
!!--declarations----------------------------------------------------------------
    use precision
    implicit none
    !
    ! Call variables
    !
    type(bedcomp_data)                   :: this
    integer                , intent(in)  :: nm
    !
    ! Local variables
    !
    integer                                   :: k
    integer                                   :: l
    real(fp)                                  :: load
    real(fp)                                  :: thnew
    real(fp)                                  :: dzc
    real(fp) , dimension(:,:,:), pointer      :: msed
    real(fp) , dimension(:,:)  , pointer      :: preload
    real(fp) , dimension(:,:)  , pointer      :: svfrac
    real(fp) , dimension(:,:)  , pointer      :: thlyr
    !
    !! executable statements -------------------------------------------------------
    msed       => this%state%msed
    preload    => this%state%preload
    svfrac     => this%state%svfrac
    thlyr      => this%state%thlyr
    !
    select case(this%settings%iunderlyr)
    case(2)
       load = 0.0_fp
       do k = 2, this%settings%nlyr
          do l = 1, this%settings%nfrac
             load = load + msed(l, k-1, nm)
          enddo
          if (comparereal(thlyr(k, nm),0.0_fp)==0) then
             !
             ! layers with zero thickness don't consolidate
             !
          elseif (load > preload(k, nm)) then
             !
             ! compute consolidation
             !
             dzc = 0.0_fp ! function of load-preload(k,nm)
             !
             ! reduce layer thickness and porosity, i.e. increase svfrac
             !
             thnew = thlyr(k, nm) - dzc
             svfrac(k, nm) = svfrac(k, nm)*thlyr(k, nm)/thnew
             thlyr(k, nm) = thnew
             preload(k, nm) = load
          else
             !
             ! no consolidation in this layer and any layer below
             !
             return
          endif
       enddo
    case default
       ! option not available for this bed composition model
    end select
end subroutine consolidate
!
!
!
!==============================================================================
subroutine initpreload(this)
!!--description-----------------------------------------------------------------
!

! Initialize the preload array assuming that all sediment is fully consolidated.
!

!!--declarations----------------------------------------------------------------
    use precision
    implicit none
    !
    ! Function/routine arguments
    !
    type (bedcomp_data), intent(inout) :: this    
    !
    ! Local variables
    !
    integer                                   :: k
    integer                                   :: l
    integer                                   :: nm
    real(fp)                                  :: load
    real(fp) , dimension(:,:,:), pointer      :: msed
    real(fp) , dimension(:,:)  , pointer      :: preload
    !
    !! executable statements -------------------------------------------------------
    !
    msed       => this%state%msed
    preload    => this%state%preload
    !
    select case(this%settings%iunderlyr)
    case(2)
       do nm = this%settings%nmlb, this%settings%nmub
          preload(1, nm) = 0.0_fp
          load = 0.0_fp
          do k = 2, this%settings%nlyr
             do l = 1, this%settings%nfrac
                load = load + msed(l, k-1, nm)
             enddo
             preload(k, nm) = load
          enddo
       enddo
    case default
       ! option not available for this bed composition model
    end select
end subroutine initpreload

end module bedcomposition_module

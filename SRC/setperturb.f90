subroutine setperturb

!  Random noise

use vars
use params
use microphysics, only: micro_field, index_water_vapor, iqit, inci, ReffIce_P3, IceMassMixingRatio_P3, reffi, nCat_ice_P3
use micro_params, only: iqv
use sgs, only: setperturb_sgs

implicit none

integer i,j,k,ptype,it,jt
integer kk
real rrr,ranf_
real xxx,yyy,zzz
! parameters
real :: iwc_target_kgm3
real :: reff_top_um, reff_base_um
real :: nL_base, nL_top
real :: radius
integer :: ii

! --- variables used by case(30) and case(31) (must be declared up here)
real xc, yc, rxy, noise
real sum_noise, mean_noise, fracz, reff_um, nL, invrho, qice_mmr, n_mix, Tlocal
integer npts

! --- additional variables used by case(32)
real :: edge_fraction, rho_cloud_ice, T_tol_local
real :: Tv_target, T_cloud, T_increment_last, Tv_bias_last
real :: qtot_cloud, Tv_cloud_local, Tv_bias, T_increment, fac, qi_cloud

call ranset_(3*rank)

ptype = perturb_type

call setperturb_sgs(ptype)  ! set sgs fields

select case (ptype)

  case(0)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(k.le.5) then
            t(i,j,k)=t(i,j,k)+0.02*rrr*(6-k)
         endif
       end do
      end do
     end do

  case(1)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(q0(k).gt.6.e-3) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
         endif
       end do
      end do
     end do

  case(2) ! warm bubble

     if(masterproc) then
       print*, 'initialize with warm bubble:'
       print*, 'bubble_x0=',bubble_x0
       print*, 'bubble_y0=',bubble_y0
       print*, 'bubble_z0=',bubble_z0
       print*, 'bubble_radius_hor=',bubble_radius_hor
       print*, 'bubble_radius_ver=',bubble_radius_ver
       print*, 'bubble_dtemp=',bubble_dtemp
       print*, 'bubble_dq=',bubble_dq
     end if

     call task_rank_to_index(rank,it,jt)
     do k=1,nzm
       zzz = z(k)
       do j=1,ny
         yyy = dy*(j+jt)
         do i=1,nx
          xxx = dx*(i+it)
           if((xxx-bubble_x0)**2+YES3D*(yyy-bubble_y0)**2.lt.bubble_radius_hor**2 &
            .and.(zzz-bubble_z0)**2.lt.bubble_radius_ver**2) then
              rrr = cos(pi/2.*(xxx-bubble_x0)/bubble_radius_hor)**2 &
               *cos(pi/2.*(yyy-bubble_y0)/bubble_radius_hor)**2 &
               *cos(pi/2.*(zzz-bubble_z0)/bubble_radius_ver)**2
              t(i,j,k) = t(i,j,k) + bubble_dtemp*rrr
              micro_field(i,j,k,index_water_vapor) = &
                  micro_field(i,j,k,index_water_vapor) + bubble_dq*rrr
           end if
         end do
       end do
     end do

  case(3)   ! gcss wg1 smoke-cloud case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(q0(k).gt.0.5e-3) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
         endif
       end do
      end do
     end do

  case(4)  ! gcss wg1 arm case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(z(k).le.200.) then
            t(i,j,k)=t(i,j,k)+0.1*rrr*(1.-z(k)/200.)
         endif
       end do
      end do
     end do

  case(5)  ! gcss wg1 BOMEX case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(z(k).le.1600.) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
            micro_field(i,j,k,index_water_vapor)= &
                      micro_field(i,j,k,index_water_vapor)+0.025e-3*rrr
         endif
       end do
      end do
     end do

  case(6)  ! GCSS Lagragngian ASTEX

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(q0(k).gt.6.e-3) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
            micro_field(i,j,k,index_water_vapor)= &
                      micro_field(i,j,k,index_water_vapor)+2.5e-5*rrr
         endif
       end do
      end do
     end do

  case(22) !bloss: Try to make a general perturbation for boundary layer cloud simulations.
           ! Add noise everywhere that the water mass mixing ratio is more than half
           !   the value at the surface.  The noise has amplitude of 0.1K and 2% of the initial q0(k).

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(q0(k).gt.0.5*q0(1)) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
            micro_field(i,j,k,index_water_vapor)= &
                 (1. + 0.02*rrr)*micro_field(i,j,k,index_water_vapor)
         endif
       end do
      end do
     end do

  case(30)
    !============================================================
    ! Isolated anvil cloud perturbation (SAM-consistent)
    !
    ! Adds weak random temperature noise (±0.01 K) ONLY in a
    ! prescribed anvil region:
    !   base = 8 km, top = 13 km, radius = 30 km,
    ! centered in the 256-km square domain.
    !
    ! NOTE: This does NOT initialize ice condensate. In your SAM build,
    ! ice mass/number indices are not exposed here, so the cloud must
    ! already exist (typically via restart).
    !============================================================

    ! Domain center (global coordinates)
    xc = 0.5 * nx * dx
    yc = 0.5 * ny * dy

    call task_rank_to_index(rank,it,jt)

    do k = 1, nzm
      zzz = z(k)

      if (zzz .ge. 8000.0 .and. zzz .le. 13000.0) then

        do j = 1, ny
          yyy = dy * (j + jt)

          do i = 1, nx
            xxx = dx * (i + it)

            rxy = sqrt( (xxx - xc)**2 + YES3D*(yyy - yc)**2 )
            if (rxy .le. 30000.0) then
            ! Weak temperature perturbation (±0.01 K)
            noise = 0.01 * (2.0*ranf_() - 1.0)
            t(i,j,k) = t(i,j,k) + noise

            ! VERY weak water vapor perturbation (safe, SAM-consistent)
            micro_field(i,j,k,index_water_vapor) = &
                  micro_field(i,j,k,index_water_vapor) * (1.0 + 1.0e-6*noise)
            endif
          end do
        end do
      endif
    end do

  case(31)
    !============================================================
    ! Isolated anvil cloud initialization (Case 31)
    !
    ! Initializes an optically thick anvil cloud with:
    !  - Base = 8 km, Top = 13 km, Radius = 30 km (centered domain)
    !  - Ice water content (IWC) = 0.3 g/m3 (uniform target)
    !  - Effective radius: 20 um at top -> 40 um at base (linear)
    !  - Ice number concentration: 900 -> 3600 L^-1 (base->top),
    !      enforce >3000 L^-1 above 10.5 km
    !  - Random potential temperature perturbations ~ ±0.01 K (zero mean)
    !  - Attempts a simple buoyant compensation so cloud is near-neutral
    !============================================================

    if(masterproc) then
      print*, 'Initialize isolated anvil cloud (case 31): IWC=0.3 g/m3, reff 20->40um, N 900->3600 L^-1'
    end if

    xc = 0.5 * nx_gl * dx
    yc = 0.5 * ny_gl * dy

    iwc_target_kgm3 = 0.3e-3 ! 0.3 g/m3 -> kg/m3
    reff_top_um = 20.0
    reff_base_um = 40.0
    nL_base = 900.0
    nL_top = 3600.0
    radius = 30000.0
    ii = 1

    ! First pass: set microphysics fields and accumulate statistics for mean noise
    sum_noise = 0.0
    npts = 0

    call task_rank_to_index(rank,it,jt)

    do k = 1, nzm
      zzz = z(k)
      if (zzz .ge. zbase .and. zzz .le. ztop) then
        ! fractional height from base (0) to top (1)
        fracz = (zzz - zbase) / (ztop - zbase)

        ! determine target reff and number concentration at this level
        reff_um = reff_base_um + (reff_top_um - reff_base_um) * fracz
        nL = nL_base + (nL_top - nL_base) * fracz
        if (zzz .ge. 10500.0) nL = max(nL, 3000.0)

        do j = 1, ny
          yyy = dy * (j + jt)
          do i = 1, nx
            xxx = dx * (i + it)
            rxy = sqrt( (xxx - xc)**2 + YES3D*(yyy - yc)**2 )
            if (rxy .le. radius) then

              ! convert to mixing ratios (#/kg or kg/kg)
              invrho = 1.0 / rho(k)

              qice_mmr = iwc_target_kgm3 * invrho  ! kg/kg
              n_mix = (nL / 1000.0) * invrho       ! convert L^-1 -> m^-3 then /rho -> #/kg

              ! set total ice mass and number in micro_field (category ii)
              micro_field(i,j,k, iqit(ii)) = qice_mmr
              micro_field(i,j,k, inci(ii)) = n_mix

              ! set P3 arrays (per-category)
              IceMassMixingRatio_P3(i,j,k,ii) = qice_mmr
              ReffIce_P3(i,j,k,ii) = reff_um
              reffi(i,j,k) = reff_um

              ! random small temperature perturbation (±0.01 K)
              noise = 0.01 * (2.0*ranf_() - 1.0)
              t(i,j,k) = t(i,j,k) + noise
              sum_noise = sum_noise + noise
              npts = npts + 1

            end if
          end do
        end do

      end if
    end do

    ! remove mean of noise in the cloud so the perturbations have zero mean
    if (npts .gt. 0) then
      mean_noise = sum_noise / dble(npts)

      do k = 1, nzm
        zzz = z(k)
        if (zzz .ge. zbase .and. zzz .le. ztop) then
          do j = 1, ny
            do i = 1, nx
              xxx = dx * (i + it)
              yyy = dy * (j + jt)
              rxy = sqrt( (xxx - xc)**2 + YES3D*(yyy - yc)**2 )
              if (rxy .le. radius) then
                t(i,j,k) = t(i,j,k) - mean_noise
              end if
            end do
          end do
        end if
      end do
    end if

    ! Buoyant compensation: simple heuristic to offset virtual temperature decrease
    ! due to added condensate mass: add ~ T * qice (small correction)
    do k = 1, nzm
      zzz = z(k)
      if (zzz .ge. zbase .and. zzz .le. ztop) then
        do j = 1, ny
          do i = 1, nx
            xxx = dx * (i + it)
            yyy = dy * (j + jt)
            rxy = sqrt( (xxx - xc)**2 + YES3D*(yyy - yc)**2 )
            if (rxy .le. radius) then
              qice_mmr = micro_field(i,j,k, iqit(ii))
              Tlocal = tabs(i,j,k)
              ! add a small warming to offset the virtual temperature reduction
              t(i,j,k) = t(i,j,k) + Tlocal * qice_mmr
            end if
          end do
        end do
      end if
    end do

  case(32)
    !============================================================
    ! Isolated anvil cloud initialization
    ! Based on Blaz Gasparini case (8)
    !
    ! Cloud:
    !   Base   = zbase
    !   Top    = ztop
    !   Radius = 30 km
    !   IWC    = 0.3 g m^-3 in cloud interior
    !
    ! Effective radius:
    !   40 um at cloud base
    !   20 um at cloud top
    !
    ! Major differences from old Case 31:
    !   1. Smooth cosine taper at cloud edge
    !   2. Cloud vapor is initialized at ice saturation
    !   3. Cloud temperature is iteratively adjusted so that
    !      Tv_cloud = Tv_background
    !   4. Ice number concentration is diagnosed from qi and reff
    !============================================================

    if(masterproc) then
      print*, 'Initialize isolated anvil cloud (case 31)'
      print*, 'IWC = 0.3 g/m3, reff = 40 um (base) -> 20 um (top)'
      print*, 'Ice-saturated and neutrally buoyant initialization'
    end if

    !------------------------------------------------------------
    ! Cloud geometry
    !------------------------------------------------------------

    xc = 0.5 * real(nx_gl) * dx
    yc = 0.5 * real(ny_gl) * dy

    radius = 30000.0

    ! Fraction of cloud radius occupied by transition region.
    ! Example: 0.10 means outermost 10% of radius is tapered.
    edge_fraction = 0.10

    ! Target cloud-interior IWC
    iwc_target_kgm3 = 0.3e-3

    ! Effective radius
    reff_base_um = 40.0
    reff_top_um  = 20.0

    ! P3 ice category
    ii = 1

    ! Ice bulk density used in Blaz's number-concentration relation
    rho_cloud_ice = 917.0

    ! Tolerance for neutral-buoyancy iteration
    T_tol_local = 1.e-6

    call task_rank_to_index(rank,it,jt)

    !============================================================
    ! Initialize cloud
    !============================================================

    do k = 1,nzm

      if((z(k) .ge. zbase) .and. (z(k) .lt. ztop)) then

        !--------------------------------------------------------
        ! Environmental virtual temperature
        !
        ! The environmental sounding contains no cloud ice here:
        !
        ! Tv_env = T_env * (1 + epsv*qv_env)
        !--------------------------------------------------------

        Tv_target = tabs0(k) * (1.0 + epsv*q0(k))

        ! Fractional height:
        !     0 = cloud base
        !     1 = cloud top
        fracz = (z(k)-zbase)/(ztop-zbase)

        ! Effective radius:
        !     40 um at base
        !     20 um at top
        reff_um = reff_base_um +                         &
             (reff_top_um-reff_base_um)*fracz

        do j = 1,ny

          yyy = dy * (j+jt)

          do i = 1,nx

            xxx = dx * (i+it)

            ! Normalized horizontal radius
            rrr = sqrt((xxx-xc)**2 +                     &
                 YES3D*(yyy-yc)**2) / radius

            if(rrr .lt. 1.0) then

              !==================================================
              ! Smooth cloud edge
              !==================================================

              if(rrr .lt. 1.0-edge_fraction) then

                fac = 1.0

              else

                fac = 0.5 *                              &
                     (1.0-cos(pi*(1.0-rrr)/edge_fraction))

              end if


              !==================================================
              ! Ice mass mixing ratio
              !
              ! Constant IWC in cloud interior, smoothly reduced
              ! to zero at the cloud edge.
              !==================================================

              qi_cloud = fac * iwc_target_kgm3 / rho(k)


              !==================================================
              ! Find neutrally buoyant cloud temperature
              !==================================================

              ! Initial guess
              T_cloud = tabs(i,j,k)

              T_increment_last = 0.0
              Tv_bias_last = 0.0

              do kk = 1,5

                ! Ice-saturated water vapor + condensate
                qtot_cloud = qsati(T_cloud,pres(k))       &
                           + qi_cloud

                !------------------------------------------------
                ! Virtual temperature including condensate
                !
                ! qtot = qv + qi
                !
                ! Tv = T*[1 + epsv*qtot
                !          - (1+epsv)*qi]
                !
                ! which is equivalent to approximately
                !
                ! Tv = T*(1 + epsv*qv - qi)
                !------------------------------------------------

                Tv_cloud_local = T_cloud *                &
                     (1.0 + epsv*qtot_cloud               &
                     - (1.0+epsv)*qi_cloud)

                Tv_bias = Tv_target - Tv_cloud_local


                ! Converged
                if(abs(Tv_bias) .lt. T_tol_local) then

                  exit

                else

                  if(kk .eq. 1) then

                    ! First iteration: simple correction
                    T_increment = Tv_bias

                  else

                    ! Secant-method correction
                    if(abs(Tv_bias-Tv_bias_last)           &
                         .gt. 1.e-12) then

                      T_increment =                       &
                           -T_increment_last*Tv_bias       &
                           /(Tv_bias-Tv_bias_last)

                    else

                      T_increment = Tv_bias

                    end if

                  end if

                  T_cloud = T_cloud + T_increment

                  Tv_bias_last = Tv_bias
                  T_increment_last = T_increment

                end if

              end do


              !==================================================
              ! Update SAM thermodynamic fields
              !==================================================

              ! Vapor mixing ratio:
              ! qtot = qv + qi
              micro_field(i,j,k,iqv) =                   &
                   qtot_cloud - qi_cloud

              ! Ice mass mixing ratio
              micro_field(i,j,k,iqit(ii)) = qi_cloud

              ! SAM prognostic temperature variable
              t(i,j,k) = T_cloud + gamaz(k)              &
                   - fac_sub*micro_field(i,j,k,iqit(ii))

              ! Absolute temperature
              tabs(i,j,k) = T_cloud


              !==================================================
              ! P3 fields
              !==================================================

              IceMassMixingRatio_P3(i,j,k,ii) = qi_cloud

              ReffIce_P3(i,j,k,ii) = reff_um
              reffi(i,j,k) = reff_um


              !==================================================
              ! Diagnose ice number concentration
              !
              ! Same relation as Blaz Case (8).
              !
              ! reff_um is converted to meters with 1.e-6.
              !
              ! Result is number per kg of air.
              !==================================================

              micro_field(i,j,k,inci(ii)) =               &
                   micro_field(i,j,k,iqit(ii))             &
                   /(rho_cloud_ice*pi)                     &
                   *(3.0/(2.e-6*reff_um))**3


            end if

          end do

        end do

      end if

    end do  

  case(-1)
    !bloss: no perturbation

  case default

       if(masterproc) print*,'perturb_type is not defined in setperturb(). Exitting...'
       call task_abort()

end select

end subroutine setperturb

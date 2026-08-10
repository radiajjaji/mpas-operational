!************************************************************************
!
! * Author: Jiawang Feng and Chun Zhao from AEMOL USTC
! Last update: Dec 2019
!
! Contact:
! Chun Zhao, PhD
! Professor 
! University of Science and Technology of China 
! Hefei, Anhui, 230026 
! Email: ChunZhao@ustc.edu.cn
!
! Please report any bugs or problems to Chun Zhao
!
! References: 
!   Copied from WRF-Chem-v3.5:chem/module_dep_simple.F
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************
MODULE module_dep_simple
use mpas_kind_types
use mpas_pool_routines
use mpas_atmchem_configs
use mpas_atmchem_vars
use mpas_atmphys_vars
use mpas_atmchem_utilities
use module_data_binaer_asect
use module_data_binaer_other
use mpas_log


      IMPLICIT NONE

!--------------------------------------------------
! many of these parameters will depend on the RADM mechanism!
! if you change it, lets talk about it and get it done!!!
!--------------------------------------------------

      INTEGER, PARAMETER :: dep_seasons = 5
      INTEGER, PARAMETER :: nlu = 25
      REAL, parameter    :: small_value = 1.e-36
      REAL, parameter    :: large_value = 1.e36

!--------------------------------------------------
! following currently hardwired to USGS
!--------------------------------------------------
      integer, parameter :: isice_temp   = 24
      integer, parameter :: iswater_temp = 16
      integer, parameter :: wrf2mz_lt_map(nlu) = (/ 1, 2, 2, 2, 2, &
                                                    4, 3, 3, 3, 3, &
                                                    4, 5, 4, 5, 6, &
                                                    7, 9, 6, 8, 9, &
                                                    6, 6, 8, 0, 0 /)
      real, parameter    :: wh2o = 18.0153
      real, parameter    :: wpan = 121.04793
      character(len=4), parameter :: mminlu = 'USGS'

      INTEGER :: month = 0
      INTEGER :: ixxxlu(nlu)
!     include modis landuse 
      INTEGER, allocatable :: luse2usgs(:)
!--

      REAL    :: kpart(nlu)
      REAL    :: rac(nlu,dep_seasons), rclo(nlu,dep_seasons), rcls(nlu,dep_seasons)
      REAL    :: rgso(nlu,dep_seasons), rgss(nlu,dep_seasons)
      REAL    :: ri(nlu,dep_seasons), rlu(nlu,dep_seasons)
      REAL    :: ri_pan(5,11)
      real    :: c0_pan(11) = (/ 0.000, 0.006, 0.002, 0.009, 0.015, &
                                 0.006, 0.000, 0.000, 0.000, 0.002, 0.002 /)
      real    :: k_pan (11) = (/ 0.000, 0.010, 0.005, 0.004, 0.003, &
                                 0.005, 0.000, 0.000, 0.000, 0.075, 0.002 /)

!--------------------------------------------------
! NO MORE THAN 1000 SPECIES FOR DEPOSITION
!--------------------------------------------------
      REAL    :: dratio(1000), hstar(1000), hstar4(1000)
      REAL    :: f0(1000), dhr(1000), scpr23(1000)

      type wesely_pft
        integer          :: npft
        integer          :: months
        INTEGER, pointer :: seasonal_wes(:,:,:,:)
        logical          :: is_allocated
      end type wesely_pft

      type(wesely_pft), allocatable :: seasonal_pft(:)

!--------------------------------------------------
! .. Default Accessibility ..
!--------------------------------------------------
    PUBLIC

!   logical, allocatable :: is_aerosol(:) ! true if field is aerosol (any phase)

    CONTAINS

SUBROUTINE wesely_driver( id, ktau, dtstep, config_flags, current_month,  &
                          gmt, julday, t_phy,moist, p8w, t8w, raincv,     &
                          p_phy, chem, rho_phy, dz8w, ddvel, aer_res_def, &
                          aer_res_zcen, ivgtyp, tsk, gsw, vegfra, pbl,    &
                          rmol, ust, znt, xlat, xlong,                    &
                          z, z_at_w, snowh, numgas,                       &
                          ids,ide, jds,jde, kds,kde,                      &
                          ims,ime, jms,jme, kms,kme,                      &
                          its,ite, jts,jte, kts,kte                       )
!--------------------------------------------------
!  Wesely dry dposition driver
!--------------------------------------------------

! USE module_model_constants 
! USE module_configure
! USE module_state_description                       
! USE module_data_sorgam
! USE module_state_description, only:  param_first_scalar        
   INTEGER,      INTENT(IN   ) :: id,julday,                              &
                                  numgas, current_month,                  &
                                  ids,ide, jds,jde, kds,kde,              &
                                  ims,ime, jms,jme, kms,kme,              &
                                  its,ite, jts,jte, kts,kte     
   INTEGER,      INTENT(IN   ) :: ktau            
      REAL,      INTENT(IN   ) :: dtstep,gmt

!--------------------------------------------------
! advected moisture variables
!--------------------------------------------------
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ), INTENT(IN ) :: &
                                                      moist  
!--------------------------------------------------
! advected chemical species
!--------------------------------------------------
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ), INTENT(INOUT ) :: &
                                                      chem
!--------------------------------------------------
! deposition velocities
!--------------------------------------------------
   REAL, DIMENSION( its:ite, jts:jte, num_chem ), INTENT(INOUT ) ::      &
                                                      ddvel                     
!--------------------------------------------------
! input from met model
!--------------------------------------------------
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) ::      &
                                                      t_phy,              &
                                                      p_phy,              &
                                                      dz8w,               &
                                                      z,                  &
                                                      t8w,                &
                                                      p8w,                &
                                                      z_at_w,             &
                                                      rho_phy
   INTEGER,DIMENSION( ims:ime , jms:jme ), INTENT(IN   ) ::               &
                                                     ivgtyp
   REAL,  DIMENSION( ims:ime , jms:jme ), INTENT(INOUT   ) ::             &
                                                     tsk,                 &
                                                     gsw,                 &
                                                     vegfra,              &
                                                     pbl,                 &
                                                     rmol,                &
                                                     ust,                 &
                                                     xlat,                &
                                                     xlong,               &
                                                     raincv,              &
                                                     znt
   REAL, intent(inout) ::                            aer_res_def(its:ite,jts:jte)
   REAL, intent(inout) ::                            aer_res_zcen(its:ite,jts:jte)
   REAL,  optional, INTENT(IN)  ::                   snowh(ims:ime,jms:jme)
   TYPE(MPAS_chem_config),  INTENT(IN) ::            config_flags

!--------------------------------------------------
! .. Local Scalars
!--------------------------------------------------
      REAL    ::  clwchem, dvfog, dvpart, pa, rad, dep_vap
      REAL    ::  rhchem, ta, ustar, vegfrac, z1, zntt
      INTEGER :: i, iland, iprt, iseason, j, jce, jcs, n, nr, ipr,jpr,nvr
      LOGICAL :: highnh3, rainflag, vegflag, wetflag
!--------------------------------------------------
! .. Local Arrays
!--------------------------------------------------
      REAL :: p(kts:kte)
      REAL :: srfres(numgas)
      REAL :: ddvel0d(numgas)

!-----------------------------------------------------------
! necessary for aerosols (module dependent)         
!-----------------------------------------------------------
      real :: rcx(numgas)

!-----------------------------------------------------------
! .. Intrinsic Functions
!-----------------------------------------------------------
      INTRINSIC max, min                             

      !jwfeng: Just a random value to prevent compile error. Not a valid value!!!
      dep_vap= 1.0 !config_flags%depo_fact

!     CALL wrf_debug(15,'in dry_dep_wesely')

      if( config_flags%chem_opt /= MOZART_KPP .and. &
          config_flags%chem_opt /= MOZCART_KPP .and. &
          config_flags%chem_opt /= MOZART_MOSAIC_4BIN_VBS0_KPP ) then
         if( julday < 90 .or. julday > 270 ) then
            iseason = 2
!            CALL wrf_debug(15,'setting iseason to 2')
         else
            iseason = 1
         endif
      end if


tile_lat_loop : &
      do j = jts,jte 
tile_lon_loop : &
         do i = its,ite 
            iprt  = 0

            iland = luse2usgs( ivgtyp(i,j) )
!--

            if( config_flags%chem_opt == MOZART_KPP .or. &
                config_flags%chem_opt == MOZCART_KPP .or. &
                config_flags%chem_opt == MOZART_MOSAIC_4BIN_VBS0_KPP ) then
               if( snowh(i,j) < .01 ) then 
                  iseason = seasonal_pft(id)%seasonal_wes(i,j,iland,current_month)
               else
                  iseason = 4
               endif
            end if
            ta    = tsk(i,j)
            rad   = gsw(i,j)
            vegfrac = vegfra(i,j)
            pa      = .01*p_phy(i,kts,j)
            clwchem = moist(i,kts,j,p_qc)
            ustar = ust(i,j)
            zntt  = znt(i,j)
            z1    = z_at_w(i,kts+1,j) - z_at_w(i,kts,j)
!-----------------------------------------------------------
!     Set logical default values
!-----------------------------------------------------------
!           rainflag = .FALSE.
!           wetflag  = .FALSE.
!           highnh3  = .FALSE.
!           if(p_qr > 1) then
!              if(moist(i,kts,j,p_qr) > 1.e-18 .or. raincv(i,j) > 0.) then
!                 rainflag = .true.
!              endif
!           endif
!           rhchem = MIN( 100.,100. * moist(i,kts,j,p_qv) / &
!                    (3.80*exp(17.27*(t_phy(i,kts,j)-273.)/(t_phy(i,kts,j)-36.))/pa))
!           rhchem = MAX(5.,RHCHEM)
!           if (rhchem >= 95.) wetflag = .true.
!
!           if( p_nh3 > 1 .and. p_so2 > 1 ) then
!              if( chem(i,kts,j,p_nh3) > 2.*chem(i,kts,j,p_so2) ) then
!                 highnh3 = .true.
!              endif
!           endif

!-----------------------------------------------------------
!--- deposition
!-----------------------------------------------------------
!     if(snowc(i,j).gt.0.)iseason=4
!            CALL rc( rcx, ta, rad, rhchem, iland, &
!                     iseason, numgas, wetflag, rainflag, highnh3, &
!                     iprt, moist(i,kts,j,p_qv), p8w(i,kts,j), config_flags%chem_opt )
!            if( config_flags%chem_opt /= MOZART_KPP .and. &
!                config_flags%chem_opt /= MOZCART_KPP .and. &
!                config_flags%chem_opt /= MOZART_MOSAIC_4BIN_VBS0_KPP)  then
!               srfres(1:numgas-2) = rcx(1:numgas-2)
!               srfres(numgas-1:numgas) = 0.
!            else
!               srfres(1:numgas) = rcx(1:numgas)
!            end if
!            CALL deppart( rmol(i,j), ustar, rhchem, clwchem, iland, dvpart, dvfog )
            ddvel0d(1:numgas) = 0.
            aer_res_def(i,j)  = 0.
            aer_res_zcen(i,j) = 0.
            CALL landusevg( ddvel0d, ustar, rmol(i,j), zntt, z1, dvpart, iland,        &
                            numgas, srfres, aer_res_def(i,j), aer_res_zcen(i,j), p_sulf )

!-----------------------------------------------------------
!wig: CBMZ does not have HO and HO2 last so need to copy all species
!      ddvel(i,j,1:numgas-2)=ddvel0d(1:numgas-2)
!-----------------------------------------------------------
            ddvel(i,j,1:numgas) = ddvel0d(1:numgas)
!           if ( (config_flags%chem_opt == RADM2           ) .or.   &
!                (config_flags%chem_opt == RADM2SORG       ) .or.   &
!                (config_flags%chem_opt == RADM2SORG_AQ    ) .or.   &
!                (config_flags%chem_opt == RADM2SORG_AQCHEM) ) then
!                    ddvel(i,j,p_hcl)         = ddvel(i,j,p_hno3)
!           end if
         end do tile_lon_loop
      end do tile_lat_loop

END SUBROUTINE wesely_driver

!      SUBROUTINE rc( rcx, t, rad, rh, iland, &
!                     iseason, numgas, wetflag, rainflag, highnh3, &
!                     iprt, spec_hum, p_srf, chem_opt )
!!----------------------------------------------------------------------
!!     THIS SUBROUTINE CALCULATES SURFACE RESISTENCES ACCORDING
!!     TO THE MODEL OF
!!     M. L. WESELY,
!!     ATMOSPHERIC ENVIRONMENT 23 (1989), 1293-1304
!!     WITH SOME ADDITIONS ACCORDING TO
!!     J. W. ERISMAN, A. VAN PUL, AND P. WYERS,
!!     ATMOSPHERIC ENVIRONMENT 28 (1994), 2595-2607
!!     WRITTEN BY  WINFRIED SEIDL, APRIL 1997
!!     MODYFIED BY WINFRIED SEIDL, MARCH 2000
!!                    FOR MM5 VERSION 3
!!----------------------------------------------------------------------
!
!! USE module_state_description                       
!
!!----------------------------------------------------------------------
!!	... dummy arguments
!!----------------------------------------------------------------------
!        INTEGER, intent(in) :: iland, iseason, numgas
!        INTEGER, intent(in) :: iprt
!        INTEGER, intent(in) :: chem_opt
!        REAL, intent(in)    :: rad, rh
!        REAL, intent(in)    :: t                            ! surface temp (K)
!        REAL, intent(in)    :: p_srf                        ! surface pressure (Pa)
!        REAL, intent(in)    :: spec_hum                     ! surface specific humidity (kg/kg)
!        real, intent(inout)   :: rcx(numgas)
!        LOGICAL, intent(in) :: highnh3, rainflag, wetflag
!
!!----------------------------------------------------------------------
!! .. Local Scalars ..
!!----------------------------------------------------------------------
!        REAL, parameter :: t0    = 298.
!        REAL, parameter :: tmelt = 273.16
!        INTEGER :: lt, n
!        REAL    :: rclx, rdc, resice, rgsx, rluo1, rluo2
!        REAL    :: rlux, rmx, rs, rsmx, rdtheta, z, wrk
!        REAL    :: qs, es, ws, dewm, dv_pan, drat
!        REAL    :: crs, tc
!        REAL    :: rs_pan, tc_pan
!        LOGICAL :: has_dew
!!----------------------------------------------------------------------
!! .. Local Arrays ..
!!----------------------------------------------------------------------
!        REAL :: hstary(numgas)
!
!!----------------------------------------------------------------------
!! .. Intrinsic Functions ..
!!----------------------------------------------------------------------
!        INTRINSIC exp
!
!        rcx(1:numgas) = 1.
!
!        tc = t - 273.15
!        rdtheta = 0.
!
!        z = 200./(rad+0.1)
!
!!!!  HARDWIRE VALUES FOR TESTING
!!       z=0.4727409
!!       tc=22.76083
!!       t=tc+273.15
!!       rad = 412.8426
!!       rainflag=.false.
!!       wetflag=.false.
!
!        IF ( tc<=0. .OR. tc>=40. ) THEN
!          rs = 9999.
!        ELSE
!          rs = ri(iland,iseason)*(1+z*z)*(400./(tc*(40.-tc)))
!        END IF
!        rdc   = 100.*(1. + 1000./(rad + 10.))/(1. + 1000.*rdtheta)
!        rluo1 = 1./(1./3000. + 3./rlu(iland,iseason))
!        rluo2 = 1./(1./1000. + 3./rlu(iland,iseason))
!        resice = 1000.*exp( -(tc + 4.) )
!        wrk    = (t0 - t)/(t0*t)
!
!
!        DO n = 1, numgas
!          IF( hstar(n) /= 0. ) then
!             hstary(n) = hstar(n)*exp( dhr(n)*wrk )
!!----------------------------------------------------------------------
!!     SPECIAL TREATMENT FOR HNO3, HNO4, H2O2, PAA
!!----------------------------------------------------------------------
!is_mozart :  if( chem_opt == MOZART_KPP .or. chem_opt == MOZCART_KPP ) then 
!                 if( n == p_hno3 ) then
!                    hstary(n) = 2.6e6*exp( 8700.*wrk )*1.e5
!                 else if( n == p_hno4 ) then
!                    hstary(n) = 32.e5
!                 else if( n == p_h2o2 ) then
!                    hstary(n) = hstary(n)*(1. + 2.2e-12*exp( -3730.*wrk )*1.e5)
!                 else if( n == p_paa ) then
!                    hstary(n) = hstary(n)*(1. + 1.8e-4*exp( -1510.*wrk )*1.e5)
!                 end if
!             endif is_mozart
!             rmx = 1./(hstary(n)/3000. + 100.*f0(n))
!             rsmx = rs*dratio(n) + rmx
!             if( (chem_opt == MOZART_KPP .or. chem_opt == MOZCART_KPP .or. chem_opt == MOZART_MOSAIC_4BIN_VBS0_KPP ) .and. &
!                 iseason /= 4 .and. p_pan > 1 ) then
!                if( iland /= iswater_temp .and. n == p_pan ) then
!!-------------------------------------------------------------------------------------
!! 	saturation vapor pressure (Pascals)
!! 	saturation mixing ratio
!! 	saturation specific humidity
!!-------------------------------------------------------------------------------------
!                    es = 611.*exp( 5414.77*(t - tmelt)/(tmelt*t) )
!                    ws = .622*es/(p_srf - es)
!                    qs = ws/(1. + ws)
!                    has_dew = .false.
!                    if( qs <= spec_hum ) then
!                      has_dew = .true.
!                    end if
!                    if( t < tmelt ) then
!                       has_dew = .false.
!                    end if
!                    if( has_dew .or. rainflag ) then
!                       dewm = 3.
!                    else
!                       dewm = 1.
!                    endif
!                    drat   = wpan/wh2o
!                    tc_pan = t - tmelt
!                    if( t > tmelt .and. t < 313.15 ) then
!                       crs = (1. + (200./(rad + .1))**2) * (400./(tc_pan*(40. - tc_pan)))
!                    else
!                       crs = large_value
!                    end if
!                    lt     = wrf2mz_lt_map(iland)
!                    rs_pan = ri_pan(iseason,lt)*crs
!                    dv_pan = c0_pan(lt) * (1. - exp( -k_pan(lt)*(dewm*rs_pan*drat)*1.e-2 ))
!                    if( dv_pan > 0. ) then
!                       rsmx = 1./dv_pan
!                    end if
!                endif
!             endif
!             rclx = 1./(1.e-5*hstary(n)/rcls(iland,iseason) &
!                        + f0(n)/rclo(iland,iseason)) + resice
!             rgsx = 1./(1.e-5*hstary(n)/rgss(iland,iseason) &
!                        + f0(n)/rgso(iland,iseason)) + resice
!             rlux = rlu(iland,iseason)/(1.e-5*hstary(n) + f0(n)) + resice
!             IF( wetflag ) THEN
!               rlux = 1./(1./(3.*rlu(iland,iseason)) + 1.e-7*hstary(n) + f0(n)/rluo1)
!             END IF
!             IF( rainflag ) THEN
!               rlux = 1./(1./(3.*rlu(iland,iseason)) + 1.e-7*hstary(n) + f0(n)/rluo2)
!             END IF
!             rcx(n) = 1./(1./rsmx + 1./rlux + 1./(rdc + rclx) + 1./(rac(iland,iseason) + rgsx))
!             rcx(n) = max( 1.,rcx(n) )
!          end IF
!        END DO
!
!!--------------------------------------------------
!!     SPECIAL TREATMENT FOR OZONE
!!--------------------------------------------------
!        if(p_o3 > 1)then
!           hstary(p_o3) = hstar(p_o3)*exp( dhr(p_o3)*(298. - t)/(298.*t) )
!           rmx = 1./(hstary(p_o3)/3000.+100.*f0(p_o3))
!           rsmx = rs*dratio(p_o3) + rmx
!           rlux = rlu(iland,iseason)/(1.E-5*hstary(p_o3)+f0(p_o3)) + resice
!           rclx = rclo(iland,iseason) + resice
!           rgsx = rgso(iland,iseason) + resice
!           IF (wetflag) rlux = rluo1
!           IF (rainflag) rlux = rluo2
!           rcx(p_o3) = 1. &
!                       /(1./rsmx+1./rlux+1./(rdc+rclx)+1./(min(100.,rac(iland,iseason))+rgsx))
!           rcx(p_o3) = max( 1.,rcx(p_o3) )
!        endif
!
!!     SPECIAL TREATMENT FOR SO2 (Wesely)
!!       HSTARY(P_SO2)=HSTAR(P_SO2)*EXP(DHR(P_SO2)*(1./T-1./298.))
!!       RMX=1./(HSTARY(P_SO2)/3000.+100.*F0(P_SO2))
!!       RSMX=RS*DRATIO(P_SO2)+RMX
!!       RLUX=RLU(ILAND,ISEASON)/(1.E-5*HSTARY(P_SO2)+F0(P_SO2))
!!    &       +RESICE
!!       RCLX=RCLS(ILAND,ISEASON)+RESICE
!!       RGSX=RGSS(ILAND,ISEASON)+RESICE
!!       IF ((wetflag).OR.(RAINFLAG)) THEN
!!         IF (ILAND.EQ.1) THEN
!!           RLUX=50.
!!         ELSE
!!           RLUX=100.
!!         END IF
!!       END IF
!!       RCX(P_SO2)=1./(1./RSMX+1./RLUX+1./(RDC+RCLX)
!!    &                +1./(RAC(ILAND,ISEASON)+RGSX))
!!       IF (RCX(P_SO2).LT.1.) RCX(P_SO2)=1.
!
!!--------------------------------------------------
!!     SO2 according to Erisman et al. 1994
!!       R_STOM
!!--------------------------------------------------
!is_so2 : &
!     if( p_so2 > 1 ) then
!        rsmx = rs*dratio(p_so2)
!!--------------------------------------------------
!!       R_EXT
!!--------------------------------------------------
!        IF (tc> -1. ) THEN
!          IF (rh<81.3) THEN
!            rlux = 25000.*exp(-0.0693*rh)
!          ELSE
!            rlux = 0.58E12*exp(-0.278*rh)
!          END IF
!        END IF
!        IF (((wetflag) .OR. (rainflag)) .AND. (tc> -1. )) THEN
!          rlux = 1.
!        END IF
!        IF ((tc>= -5. ) .AND. (tc<= -1. )) THEN
!          rlux = 200.
!        END IF
!        IF (tc< -5. ) THEN
!          rlux = 500.
!        END IF
!!--------------------------------------------------
!!       INSTEAD OF R_INC R_CL and R_DC of Wesely are used
!!--------------------------------------------------
!        rclx = rcls(iland,iseason)
!!--------------------------------------------------
!!       DRY SURFACE
!!--------------------------------------------------
!        rgsx = 1000.
!!--------------------------------------------------
!!       WET SURFACE
!!--------------------------------------------------
!        IF ((wetflag) .OR. (rainflag)) THEN
!          IF (highnh3) THEN
!            rgsx = 0.
!          ELSE
!            rgsx = 500.
!          END IF
!        END IF
!!--------------------------------------------------
!!       WATER
!!--------------------------------------------------
!        IF (iland==iswater_temp) THEN
!          rgsx = 0.
!        END IF
!!--------------------------------------------------
!!       SNOW
!!--------------------------------------------------
!        IF( iseason==4 .OR. iland==isice_temp ) THEN
!          IF( tc > 2. ) THEN
!            rgsx = 0.
!          else IF ( tc >= -1. .AND. tc <= 2. ) THEN
!            rgsx = 70.*(2. - tc)
!          else IF ( tc < -1. ) THEN
!            rgsx = 500.
!          END IF
!        END IF
!!--------------------------------------------------
!!       TOTAL SURFACE RESISTENCE
!!--------------------------------------------------
!        IF ((iseason/=4) .AND. (ixxxlu(iland)/=1) .AND. (iland/=iswater_temp) .AND. &
!            (iland/=isice_temp)) THEN
!          rcx(p_so2) = 1./(1./rsmx+1./rlux+1./(rclx+rdc+rgsx))
!        ELSE
!          rcx(p_so2) = rgsx
!        END IF
!        rcx(p_so2) = max( 1.,rcx(p_so2) )
!     end if is_so2
!!--------------------------------------------------
!!     NH3 according to Erisman et al. 1994
!!       R_STOM
!!--------------------------------------------------
!is_nh3: if( p_nh3 > 1 ) then
!        rsmx = rs*dratio(p_nh3)
!!--------------------------------------------------
!!       GRASSLAND (PASTURE DURING GRAZING)
!!--------------------------------------------------
!        IF (ixxxlu(iland)==3) THEN
!          IF (iseason==1) THEN
!!--------------------------------------------------
!!           SUMMER
!!--------------------------------------------------
!            rcx(p_nh3) = 1000.
!          END IF
!          IF ((iseason==2) .OR. (iseason==3) .OR. (iseason==5)) THEN
!!--------------------------------------------------
!!           WINTER, NO SNOW
!!--------------------------------------------------
!            IF (tc>-1.) THEN
!              IF (rad/=0.) THEN
!                rcx(p_nh3) = 50.
!              ELSE
!                rcx(p_nh3) = 100.
!              END IF
!              IF ((wetflag) .OR. (rainflag)) THEN
!                rcx(p_nh3) = 20.
!              END IF
!            END IF
!            IF ((tc>=(-5.)) .AND. (tc<=-1.)) THEN
!              rcx(p_nh3) = 200.
!            END IF
!            IF (tc<(-5.)) THEN
!              rcx(p_nh3) = 500.
!            END IF
!          END IF
!        END IF
!!--------------------------------------------------
!!       AGRICULTURAL LAND (CROPS AND UNGRAZED PASTURE)
!!--------------------------------------------------
!        IF (ixxxlu(iland)==2) THEN
!          IF (iseason==1) THEN
!!--------------------------------------------------
!!           SUMMER
!!--------------------------------------------------
!            IF (rad/=0.) THEN
!              rcx(p_nh3) = rsmx
!            ELSE
!              rcx(p_nh3) = 200.
!            END IF
!            IF ((wetflag) .OR. (rainflag)) THEN
!              rcx(p_nh3) = 50.
!            END IF
!          END IF
!          IF ((iseason==2) .OR. (iseason==3) .OR. (iseason==5)) THEN
!!--------------------------------------------------
!!           WINTER, NO SNOW
!!--------------------------------------------------
!            IF (tc>-1.) THEN
!              IF (rad/=0.) THEN
!                rcx(p_nh3) = rsmx
!              ELSE
!                rcx(p_nh3) = 300.
!              END IF
!              IF ((wetflag) .OR. (rainflag)) THEN
!                rcx(p_nh3) = 100.
!              END IF
!            END IF
!            IF ((tc>=(-5.)) .AND. (tc<=-1.)) THEN
!              rcx(p_nh3) = 200.
!            END IF
!            IF (tc<(-5.)) THEN
!              rcx(p_nh3) = 500.
!            END IF
!          END IF
!        END IF
!!--------------------------------------------------
!!       SEMI-NATURAL ECOSYSTEMS AND FORESTS
!!--------------------------------------------------
!        IF ((ixxxlu(iland)==4) .OR. (ixxxlu(iland)==5) .OR. (ixxxlu( &
!            iland)==6)) THEN
!          IF (rad/=0.) THEN
!            rcx(p_nh3) = 500.
!          ELSE
!            rcx(p_nh3) = 1000.
!          END IF
!          IF ((wetflag) .OR. (rainflag)) THEN
!            IF (highnh3) THEN
!              rcx(p_nh3) = 100.
!            ELSE
!              rcx(p_nh3) = 0.
!            END IF
!          END IF
!          IF ((iseason==2) .OR. (iseason==3) .OR. (iseason==5)) THEN
!!--------------------------------------------------
!!           WINTER, NO SNOW
!!--------------------------------------------------
!            IF ((tc>=(-5.)) .AND. (tc<=-1.)) THEN
!              rcx(p_nh3) = 200.
!            END IF
!            IF (tc<(-5.)) THEN
!              rcx(p_nh3) = 500.
!            END IF
!          END IF
!        END IF
!!--------------------------------------------------
!!       WATER
!!--------------------------------------------------
!        IF (iland==iswater_temp) THEN
!          rcx(p_nh3) = 0.
!        END IF
!!--------------------------------------------------
!!       URBAN AND DESERT (SOIL SURFACES)
!!--------------------------------------------------
!        IF (ixxxlu(iland)==1) THEN
!          IF ( .NOT. wetflag) THEN
!            rcx(p_nh3) = 50.
!          ELSE
!            rcx(p_nh3) = 0.
!          END IF
!        END IF
!!--------------------------------------------------
!!       SNOW COVERED SURFACES OR PERMANENT ICE
!!--------------------------------------------------
!        IF ((iseason==4) .OR. (iland==isice_temp)) THEN
!          IF (tc>2.) THEN
!            rcx(p_nh3) = 0.
!          END IF
!          IF ((tc>=(-1.)) .AND. (tc<=2.)) THEN
!            rcx(p_nh3) = 70.*(2.-tc)
!          END IF
!          IF (tc<(-1.)) THEN
!            rcx(p_nh3) = 500.
!          END IF
!        END IF
!        rcx(p_nh3) = max( 1.,rcx(p_nh3) )
!        endif is_nh3
!      END SUBROUTINE rc

      SUBROUTINE deppart( rmol, ustar, rh, clw, iland, &
                          dvpart, dvfog )
!--------------------------------------------------
!     THIS SUBROUTINE CALCULATES SURFACE DEPOSITION VELOCITIES
!     FOR FINE AEROSOL PARTICLES ACCORDING TO THE MODEL OF
!     J. W. ERISMAN, A. VAN PUL, AND P. WYERS,
!     ATMOSPHERIC ENVIRONMENT 28 (1994), 2595-2607
!     WRITTEN BY WINFRIED SEIDL, APRIL 1997
!     MODIFIED BY WINFRIED SEIDL, MARCH 2000
!            FOR MM5 VERSION 3
!--------------------------------------------------

!--------------------------------------------------
! .. Scalar Arguments ..
!--------------------------------------------------
        INTEGER, intent(in) :: iland
        REAL, intent(in)    :: clw, rh, rmol, ustar
        REAL, intent(inout)   :: dvfog, dvpart

!--------------------------------------------------
! .. Intrinsic Functions ..
!--------------------------------------------------
        INTRINSIC exp

        dvpart = ustar/kpart(iland)
        IF (rmol<0.) THEN
!--------------------------------------------------
!         UNSTABLE LAYERING CORRECTION
!--------------------------------------------------
          dvpart = dvpart*(1.+(-300.*rmol)**0.66667)
        END IF
        IF (rh>80.) THEN
!--------------------------------------------------
!         HIGH RELATIVE HUMIDITY CORRECTION
!         ACCORDING TO J. W. ERISMAN ET AL.
!         ATMOSPHERIC ENVIRONMENT 31 (1997), 321-332
!--------------------------------------------------
          dvpart = dvpart*(1.+0.37*exp((rh-80.)/20.))
        END IF

!--------------------------------------------------
!       SEDIMENTATION VELOCITY OF FOG WATER ACCORDING TO
!       R. FORKEL, W. SEIDL, R. DLUGI AND E. DEIGELE
!       J. GEOPHYS. RES. 95D (1990), 18501-18515
!--------------------------------------------------
        dvfog = 0.06*clw
        IF (ixxxlu(iland)==5) THEN
!--------------------------------------------------
!         TURBULENT DEPOSITION OF FOG WATER IN CONIFEROUS FOREST ACCORDI
!         A. T. VERMEULEN ET AL.
!         ATMOSPHERIC ENVIRONMENT 31 (1997), 375-386
!--------------------------------------------------
          dvfog = dvfog + 0.195*ustar*ustar
        END IF

      END SUBROUTINE deppart

      SUBROUTINE landusevg( vgs, ustar, rmol, z0, zz, &
                            dvparx, iland, numgas, srfres, aer_res_def, &
                            aer_res_zcen, p_sulf )
!--------------------------------------------------
!     This subroutine calculates the species specific deposition velocit
!     as a function of the local meteorology and land use.  The depositi
!     Velocity is also landuse specific.
!     Reference: Hsieh, C.M., Wesely, M.L. and Walcek, C.J. (1986)
!                A Dry Deposition Module for Regional Acid Deposition
!                EPA report under agreement DW89930060-01
!     Revised version by Darrell Winner (January 1991)
!        Environmental Engineering Science 138-78
!           California Institute of Technology
!              Pasadena, CA  91125
!     Modified by Winfried Seidl (August 1997)
!       Fraunhofer-Institut fuer Atmosphaerische Umweltforschung
!                    Garmisch-Partenkirchen, D-82467
!          for use of Wesely and Erisman surface resistances
!     Inputs:
!        Ustar  : The grid average friction velocity (m/s)
!        Rmol   : Reciprocal of the Monin-Obukhov length (1/m)
!        Z0     : Surface roughness height for the grid square (m)
!        SrfRes : Array of landuse/atmospheric/species resistances (s/m)
!        Slist  : Array of chemical species codes
!        Dvparx : Array of surface deposition velocity of fine aerosol p
!     Outputs:
!        Vgs    : Array of species and landuse specific deposition
!                 velocities (m/s)
!        Vg     : Cell-average deposition velocity by species (m/s)
!     Variables used:
!        SCPR23  : (Schmidt #/Prandtl #)**(2/3) Diffusion correction fac
!        Zr      : Reference Height (m)
!        Iatmo   : Parameter specifying the stabilty class (Function of
!        Z0      : Surface roughness height (m)
!        karman  : Von Karman constant (from module_model_constants)
!--------------------------------------------------

        USE mpas_atmchem_constants, only: karman

!--------------------------------------------------
! .. Scalar Arguments ..
!--------------------------------------------------
        INTEGER, intent(in) :: iland, numgas, p_sulf
        REAL, intent(in)    :: dvparx, ustar, z0, zz
        REAL, intent(inout) :: rmol
        REAL, intent(inout) :: aer_res_def
        REAL, intent(inout) :: aer_res_zcen
!--------------------------------------------------
! .. Array Arguments ..
!--------------------------------------------------
        REAL, intent(in)  :: srfres(numgas)
        REAL, intent(inout) :: vgs(numgas)

!--------------------------------------------------
! .. Local Scalars ..
!--------------------------------------------------
        INTEGER :: jspec
        REAL    :: vgp, vgpart, zr
        REAL    :: rmol_tmp
!--------------------------------------------------
! .. Local Arrays ..
!--------------------------------------------------
        REAL :: vgspec(numgas)

        !czhao initialize vgs
        vgs(:)=0.0

!--------------------------------------------------
!   Calculate aerodynamic resistance for reference
!   height = layer center
!--------------------------------------------------
        zr = zz*.5
        rmol_tmp = rmol
        CALL depvel( numgas, rmol_tmp, zr, z0, ustar, &
                     vgspec, vgpart, aer_res_zcen )
!--------------------------------------------------
!   Set the reference height (2.0 m)
!--------------------------------------------------
!       zr = 10.0
        zr = 2.0

!--------------------------------------------------
!   CALCULATE THE DEPOSITION VELOCITY without any surface
!   resistance term, i.e. 1 / (ra + rb)
!--------------------------------------------------
        CALL depvel( numgas, rmol, zr, z0, ustar, &
                     vgspec, vgpart, aer_res_def )

!--------------------------------------------------
!   Calculate the deposition velocity for each species
!   and grid cell by looping through all the possibile combinations
!   of the two
!--------------------------------------------------
        vgp = 1.0/((1.0/vgpart)+(1.0/dvparx))
!--------------------------------------------------
!   Loop through the various species
!--------------------------------------------------
        DO jspec = 1, numgas
!--------------------------------------------------
!   Add in the surface resistance term, rc (SrfRes)
!--------------------------------------------------
          vgs(jspec) = 1.0/(1.0/vgspec(jspec) + srfres(jspec))
        END DO
        vgs(p_sulf) = vgp

      END SUBROUTINE landusevg

      SUBROUTINE depvel( numgas, rmol, zr1, z01, ustar, &
                         depv, vgpart, aer_res )
!--------------------------------------------------
!     THIS FUNCTION HAS BEEN DESIGNED TO EVALUATE AN UPPER LIMIT
!     FOR THE POLLUTANT DEPOSITION VELOCITY AS A FUNCTION OF THE
!     SURFACE ROUGHNESS AND METEOROLOGICAL CONDITIONS.
!     PROGRAM WRITTEN BY GREGORY J.MCRAE (NOVEMBER 1977)
!         Modified by Darrell A. Winner  (Feb. 1991)
!                  by Winfried Seidl     (Aug. 1997)
!.....PROGRAM VARIABLES...
!     RMOL     - RECIPROCAL OF THE MONIN-OBUKHOV LENGTH
!     ZR       - REFERENCE HEIGHT
!     Z0       - SURFACE ROUGHNESS HEIGHT
!     SCPR23   - (Schmidt #/Prandtl #)**(2/3) Diffusion correction fact
!     UBAR     - ABSOLUTE VALUE OF SURFACE WIND SPEED
!     DEPVEL   - POLLUTANT DEPOSITION VELOCITY
!     Vk       - VON KARMAN CONSTANT
!     USTAR    - FRICTION VELOCITY U*
!     POLINT   - POLLUTANT INTEGRAL
!     AER_RES  - AERODYNAMIC RESISTANCE
!.....REFERENCES...
!     MCRAE, G.J. ET AL. (1983) MATHEMATICAL MODELING OF PHOTOCHEMICAL
!       AIR POLLUTION, ENVIRONMENTAL QUALITY LABORATORY REPORT 18,
!       CALIFORNIA INSTITUTE OF TECHNOLOGY, PASADENA, CALIFORNIA.
!.....RESTRICTIONS...
!     1. THE MODEL EDDY DIFFUSIVITIES ARE BASED ON MONIN-OBUKHOV
!        SIMILARITY THEORY AND SO ARE ONLY APPLICABLE IN THE
!        SURFACE LAYER, A HEIGHT OF O(30M).
!     2. ALL INPUT UNITS MUST BE CONSISTENT
!     3. THE PHI FUNCTIONS USED TO CALCULATE THE FRICTION
!        VELOCITY U* AND THE POLLUTANT INTEGRALS ARE BASED
!        ON THE WORK OF BUSINGER ET AL.(1971).
!     4. THE MOMENTUM AND POLLUTANT DIFFUSIVITIES ARE NOT
!        THE SAME FOR THE CASES L<0 AND L>0.
!--------------------------------------------------

       USE mpas_atmchem_constants, only: karman

!--------------------------------------------------
! .. Scalar Arguments ..
!--------------------------------------------------
        INTEGER, intent(in) :: numgas
        REAL, intent(in)    :: ustar, z01, zr1
        REAL, intent(inout)   :: vgpart, aer_res
        REAL, intent(inout) :: rmol
!--------------------------------------------------
! .. Array Arguments ..
!--------------------------------------------------
        REAL, intent(inout) :: depv(numgas)
!--------------------------------------------------
! .. Local Scalars ..
!--------------------------------------------------
        REAL   :: z0, zr
        INTEGER :: l
        REAL    :: ao, ar, polint, vk
!--------------------------------------------------
! .. Intrinsic Functions ..
!--------------------------------------------------
        INTRINSIC alog
!--------------------------------------------------
!     Set the von Karman constant
!--------------------------------------------------
        vk = karman
        z0 = z01
        zr = zr1
!--------------------------------------------------
!     Calculate the diffusion correction factor
!     SCPR23 is calculated as (Sc/Pr)**(2/3) using Sc= 1.15 and Pr= 1.0
!     SCPR23 = 1.10
!--------------------------------------------------
!     DETERMINE THE STABILITY BASED ON THE CONDITIONS
!             1/L < 0 UNSTABLE
!             1/L = 0 NEUTRAL
!             1/L > 0 STABLE
!--------------------------------------------------
        !czhao +++++++++++++ 
        if (((zr/z0).lt.1.05).and.((zr/z0).gt.0.95)) then 
          if (zr .gt. z0) zr=z0*1.05
          if (zr .lt. z0) zr=z0*0.95
        endif   
        !czhao ------------- 

        if(abs(rmol) < 1.E-6 ) rmol = 0.

        IF (rmol<0) THEN
          ar = ((1.0-9.0*zr*rmol)**(0.25)+0.001)**2
          ao = ((1.0-9.0*z0*rmol)**(0.25)+0.001)**2
          polint = 0.74*(alog((ar-1.0)/(ar+1.0))-alog((ao-1.0)/(ao+1.0)))
        ELSE IF (rmol==0.) THEN
          polint = 0.74*alog(zr/z0)
        ELSE
          polint = 0.74*alog(zr/z0) + 4.7*rmol*(zr-z0)
        END IF

!--------------------------------------------------
!     CALCULATE THE Maximum DEPOSITION VELOCITY
!--------------------------------------------------
        DO l = 1, numgas
          depv(l) = ustar*vk/(2.0*scpr23(l)+polint)
        END DO
        vgpart = ustar*vk/polint
        aer_res = polint/(karman*max(ustar,1.0e-4))

      END SUBROUTINE depvel
    END MODULE module_dep_simple

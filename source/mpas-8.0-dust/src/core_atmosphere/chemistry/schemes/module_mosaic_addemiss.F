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
!   Copied from WRF-Chem-v3.5:chem/module_mosaic_addemiss.F
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************
module module_mosaic_addemiss

   USE mpas_atmchem_vars
   USE mpas_atmphys_vars
   USE mpas_kind_types
   USE mpas_derived_types
   USE mpas_pool_routines
   USE module_data_binaer_asect
   use mpas_atmchem_gocart_dust_init, only: has_done_gocart_init
   use mpas_atmchem_gocart_dust_init, only: gocart_dust_init_forMPAS

   USE mpas_atmchem_utilities
   USE mpas_log, only : mpas_log_write
 
contains

! jwfeng:chem,emis_ant are useless in MPAS. Use tend_emis instead
!        Change tend_emis() instead of chem()?
subroutine mosaic_addemiss(debug_on,configs,dminfo,                         &
            slai,ust,smois,ivgtyp,isltyp,                                   &
            dtstep, u10, v10, alt, dz8w, xland, sstemp,                     &
            ktau,num_soil_layers,u_phy,v_phy,rho_phy,p8w,erod,              &
            xlat,xlong,dx,g,emis_dust_factor,emis_seas_factor,              &
            tend_emis, e_dust, e_seas,                                      &
            dust_emiss_active, seasalt_emiss_active,                        &
            ids,ide, jds,jde, kds,kde,                                      &
            ims,ime, jms,jme, kms,kme,                                      &
            its,ite, jts,jte, kts,kte                                       )

       implicit none

       logical,intent(in) :: debug_on
       type(mpas_pool_type),intent(in):: configs
       type(dm_info),intent(in):: dminfo

       character(len=strkind), intent(in   ) :: dust_emiss_active
       integer, intent(in   ) :: seasalt_emiss_active
       integer, intent(in   ) :: ktau, num_soil_layers
       integer, intent(in   ) :: ids,ide, jds,jde, kds,kde,           &
                                 ims,ime, jms,jme, kms,kme,           &
                                 its,ite, jts,jte, kts,kte

       real, intent(in   ) :: dtstep,emis_dust_factor,emis_seas_factor,g

    ! 10-m wind speed components (m/s)
       real,  dimension( ims:ime , jms:jme ),                              &
              intent(in   ) ::   u10,v10,xland,slai,ust,dx,xlat,xlong

       integer,  dimension( ims:ime , jms:jme ),                           &
              intent(in   ) ::   ivgtyp, isltyp

       real,  dimension( ims:ime , kms:kme , jms:jme ),                    &
              intent(in   ) ::   alt, dz8w,p8w,u_phy,v_phy,rho_phy 

       real, dimension( ims:ime, num_soil_layers, jms:jme ) ,              &
              intent(inout) ::   smois
       real,  dimension( ims:ime , jms:jme, 3 ) ,                          &
          intent(in   ) ::    erod
       real,  dimension( ims:ime , jms:jme) ,                              &
          intent(in   ) ::    sstemp
       real, dimension( ims:ime, kms:kme, jms:jme, num_chem ) ,            &
              intent(inout) ::   tend_emis
       real, dimension( ims:ime, jms:jme),                                 &
              intent(inout) ::   e_dust, e_seas

! do seasalt emissions
if( lptr_seas_aer(1,1,1) > param_first_scalar ) then
    if( seasalt_emiss_active == 3 ) then
          call mosaic_seasalt_emiss2(                               &
          dtstep, u10, v10, u_phy,v_phy, alt, rho_phy, dz8w,        &
          xland, dx,g,tend_emis,sstemp,e_seas,emis_seas_factor,     &
          ids,ide, jds,jde, kds,kde,                                &
          ims,ime, jms,jme, kms,kme,                                &
          its,ite, jts,jte, kts,kte                                 )
    else if( seasalt_emiss_active > 0) then
        call mosaic_seasalt_emiss(                                  &
        dtstep, u10, v10, alt, dz8w, xland,                         &
        tend_emis, e_seas,emis_seas_factor,                         &
        ids,ide, jds,jde, kds,kde,                                  &
        ims,ime, jms,jme, kms,kme,                                  &
        its,ite, jts,jte, kts,kte, seasalt_emiss_active             )
    else if( seasalt_emiss_active > 4) then
        write(mpas_err_message,'(A,I3)') 'illegal value for config_seas_opt, must chosen from 1~4:', seasalt_emiss_active
        call chemistry_error_fatal(mpas_err_message)
    endif
else
    call chemistry_message("Don't do seasalt emission.")
endif

    if (dust_emiss_active == 'mosaic') then
        call mosaic_dust_emiss (                                            &
                 slai = slai,      ust = ust,        smois = smois        , &
                 ivgtyp = ivgtyp,  isltyp = isltyp                        , &
                 dtstep = dtstep,   u10 = u10,        v10 = v10           , &
                 alt = alt,        dz8w = dz8w,     xland = xland         , &
                 tend_emis = tend_emis, e_dust = e_dust                   , &
                 ids = ids, ide = ide, kds = kds, kde = kde, jds = jds, jde = jde , &
                 ims = ims, ime = ime, kms = kms, kme = kme, jms = jms, jme = jme , &
                 its = its, ite = ite, kts = kts, kte = kte, jts = jts, jte = jte   &
                )
    elseif (dust_emiss_active == 'gocart') then

!$OMP MASTER
        if (.not. has_done_gocart_init) then
            if(debug_on) call mpas_log_write("Reading from file: SOILPARM.TBL")
            call gocart_dust_init_forMPAS(dminfo, configs)
            has_done_gocart_init = .true.
        endif
!$OMP END MASTER
        call mosaic_dust_gocartemis (debug_on                             , &
                 ktau = 0 ,dt = dtstep ,num_soil_layers = num_soils       , &
                 alt  = alt,        u_phy = u_phy, v_phy = v_phy          , &
                 rho_phy = rho_phy, dz8w = dz8w, smois = smois            , &
                 u10 = u10, v10 = v10, p8w = p8w, erod = erod             , &
                 ivgtyp = ivgtyp, isltyp = isltyp, xland = xland          , &
                 xlat = xlat, xlong = xlong, dx = dx                      , &
                 g = g, emis_dust_factor = emis_dust_factor               , &
                 tend_emis = tend_emis, e_dust = e_dust                   , &
                 ids = ids, ide = ide, kds = kds, kde = kde, jds = jds, jde = jde , &
                 ims = ims, ime = ime, kms = kms, kme = kme, jms = jms, jme = jme , &
                 its = its, ite = ite, kts = kts, kte = kte, jts = jts, jte = jte &
                 )
    endif

end subroutine mosaic_addemiss  

!----------------------------------------------------------------------

subroutine mosaic_dust_emiss(        slai,ust, smois, ivgtyp, isltyp  , &
            dtstep, u10, v10, alt, dz8w, xland                        , &
            tend_emis, e_dust                                         , &
            ids,ide, jds,jde, kds,kde                                 , &
            ims,ime, jms,jme, kms,kme                                 , &
            its,ite, jts,jte, kts,kte                                   )

!
! adds dust emissions for mosaic aerosol species (i.e. emission tendencies
! over time dtstep are applied to the aerosol mixing ratios)
!
! This is a simple dust scheme based on Shaw et al. (2008) to appear in
! Atmospheric Environment, recoded by Jerome Fast
!
! NOTE: 
! 1) This version only works with the 8-bin version of MOSAIC.
! 2) Dust added to MOSAIC's other inorganic specie, OIN.  If Ca and CO3 are 
!    activated in the Registry, a small fraction also added to Ca and CO3.
! 3) The main departure from Shaw et al., is now alphamask is computed since
!    the land-use categories in that paper and in WRF differ.  WRF currently 
!    does not have that many land-use categories and adhoc assumptions had to
!    be made. This version was tested for Mexico in the dry season.  The main
!    land-use categories in WRF that are likely dust sources are grass, shrub,
!    and savannna (that WRF has in the desert regions of NW Mexico).  Having
!    dust emitted from these types for other locations and other times of the
!    year is not likely to be valid.
! 4) An upper bound on ustar was placed because the surface parameterizations
!    in WRF can produce unrealistically high values that lead to very high
!    dust emission rates.
! 5) Other departures' from Shaw et al. noted below, but are probably not as
!    important as 2) and 3).
! 6) Now the dust added into dust species

!  USE module_configure, only:  grid_config_rec_type
   IMPLICIT NONE


   INTEGER,      INTENT(IN   ) ::                                          &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   ) ::    dtstep

! 10-m wind speed components (m/s)
   REAL,  DIMENSION( ims:ime , jms:jme ),                                  &
          INTENT(IN   ) ::   u10, v10, xland, slai, ust
   INTEGER,  DIMENSION( ims:ime , jms:jme ),                               &
          INTENT(IN   ) ::   ivgtyp, isltyp

! alt  = 1.0/(dry air density) in (m3/kg)
! dz8w = layer thickness in (m)
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )        , &
          INTENT(IN   ) ::   alt, dz8w
   REAL, DIMENSION( ims:ime, num_soils, jms:jme )         , &
          INTENT(INOUT) ::   smois
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ) , &
          INTENT(INOUT) ::   tend_emis
   REAL, DIMENSION( ims:ime, jms:jme)                     , &
          INTENT(INOUT) ::   e_dust
         
! local variables
        integer i, j, k, l, l_oin, l_ca, l_co3, n, ii
        integer iphase, itype, izob
        integer p1st

        real dum, dumdlo, dumdhi, dumlandfrac, dumspd10
        real factaa, factbb, fracoin, fracca, fracco3, fractot
        real ustart, ustar1, ustart0
        real alphamask, f8, f50, f51, f52, wetfactor, sumdelta, ftot
        real smois_grav, wp, pclay
        real :: beta(4,7)
        real :: gamma(4), delta(4)
        real :: sz(8)
        real :: dustflux, densdust, mass1part
        real :: c_const

        p1st = param_first_scalar
!
! from Nickovic et al., JGR, 2001 and Shaw et al. 2007
! beta: fraction of clay, small silt, large silt, and sand correcsponding to Zobler class (7)
! beta (1,*) for 0.5-1 um
! beta (2,*) for 1-10 um
! beta (3,*) for 10-25 um
! beta (4,*) for 25-50 um
!
        beta(1,1)=0.12
        beta(2,1)=0.04
        beta(3,1)=0.04
        beta(4,1)=0.80
        beta(1,2)=0.34
        beta(2,2)=0.28
        beta(3,2)=0.28
        beta(4,2)=0.10
        beta(1,3)=0.45
        beta(2,3)=0.15
        beta(3,3)=0.15
        beta(4,3)=0.25
        beta(1,4)=0.12
        beta(2,4)=0.09
        beta(3,4)=0.09
        beta(4,4)=0.70
        beta(1,5)=0.40
        beta(2,5)=0.05
        beta(3,5)=0.05
        beta(4,5)=0.50
        beta(1,6)=0.34
        beta(2,6)=0.18
        beta(3,6)=0.18
        beta(4,6)=0.30
        beta(1,7)=0.22
        beta(2,7)=0.09
        beta(3,7)=0.09
        beta(4,7)=0.60
        gamma(1)=0.08
        gamma(2)=1.00
        gamma(3)=1.00
        gamma(4)=0.12

        itype=1
  if (nsize_aer(itype) .eq. 8) then
!jwfeng:04/2020
!Modified to 10 bins, up to 40um
          sz(1)=0
          sz(2)=1.78751e-06
          sz(3)=0.000273786
          sz(4)=0.00847978
          sz(5)=0.056055
          sz(6)=0.0951896
          sz(7)=0.17
          sz(8)=0.67
  else if (nsize_aer(itype) .eq. 4) then
          !BSINGH(12/11/2013): Based on the suggestions by Manish Shrivastva, sz variable is modified below.
          !Original values are commented out
          !sz(1)=0.0
          !sz(2)=0.01
          !sz(3)=0.07
          !sz(4)=0.92
          sz(1)=1.78751e-06
          sz(2)=0.00875357
          sz(3)=0.1512446
          sz(4)=0.84
          sz(5)=0.0
          sz(6)=0.0
          sz(7)=0.0
          sz(8)=0.0
        endif

!   for now just do itype=1
        itype = 1
        iphase = ai_phase

!   loop over i,j and apply dust emissions
        k = kts
        do 1830 j = jts, jte
        do 1820 i = its, ite

    if( xland(i,j) > 1.5 ) cycle

! compute wind speed anyway, even though ustar is used below

        dumlandfrac = 1.
        dumspd10=(u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j))**(0.5)
        if(dumspd10 >= 5.0) then
           dumspd10 = dumlandfrac* &
         ( dumspd10*dumspd10*(dumspd10-5.0))
         else
            dumspd10=0.
         endif

! part1 - compute vegetation mask
!
! * f8, f50, f51, f52 refer to vegetation classes from the Olsen categories
!   for desert, sand desert, grass semi-desert, and shrub semi-desert
! * in WRF, vegetation types 7, 8 and 10 are grassland, shrubland, and savanna
!   that are dominate types in Mexico and probably have some erodable surface
!   during the dry season
! * currently modified these values so that only a small fraction of cell
!   area is erodable
! * these values are highly tuneable!

         alphamask=0.001
         if (ivgtyp(i,j) .eq. 7) then
           f8=0.005
           f50=0.00
           f51=0.10
           f51=0.066
           f52=0.00
           alphamask=(f8+f50)*1.0+(f51+f52)*0.5
         endif
         if (ivgtyp(i,j) .eq. 8) then
           f8=0.010
           f50=0.00
           f51=0.00
           f52=0.15
           f52=0.10
           alphamask=(f8+f50)*1.0+(f51+f52)*0.5
         endif
         if (ivgtyp(i,j) .eq. 10) then
           f8=0.00
           f50=0.00
           f51=0.01
           f52=0.00
           alphamask=(f8+f50)*1.0+(f51+f52)*0.5
         endif

! part2 - zobler
! 
! * in Shaw's paper, dust is computed for 4 size ranges:
!   0.5-1 um 
!    1-10 um  
!   10-25 um  
!   25-50 um
! * Shaw's paper also accounts for sub-grid variability in soil
!   texture, but here we just assume the same soil texture for each
!   grid cell
! * since MOSAIC is currently has a maximum size range up to 10 um,
!   neglect upper 2 size ranges and lowest size range (assume small)
! * map WRF soil classes arbitrarily to Zolber soil textural classes
! * skip dust computations for WRF soil classes greater than 13, i.e. 
!   do not compute dust over water, bedrock, and other surfaces
! * should be skipping for water surface at this point anyway
!
         izob=0
         if(isltyp(i,j).eq.1) izob=1
         if(isltyp(i,j).eq.2) izob=1
         if(isltyp(i,j).eq.3) izob=4
         if(isltyp(i,j).eq.4) izob=2
         if(isltyp(i,j).eq.5) izob=2
         if(isltyp(i,j).eq.6) izob=2
         if(isltyp(i,j).eq.7) izob=7
         if(isltyp(i,j).eq.8) izob=2
         if(isltyp(i,j).eq.9) izob=6
         if(isltyp(i,j).eq.10) izob=5
         if(isltyp(i,j).eq.11) izob=2
         if(isltyp(i,j).eq.12) izob=3
         if(isltyp(i,j).ge.13) izob=0
         if(izob.eq.0) goto 1840
!
! part3 - dustprod
!
         do ii=1,4
           delta(ii)=0.0
         enddo
         sumdelta=0.0
         do ii=1,4
           delta(ii)=beta(ii,izob)*gamma(ii)
           if(ii.lt.4) then
             sumdelta=sumdelta+delta(ii)
           endif
         enddo
         do ii=1,4
           delta(ii)=delta(ii)/sumdelta
         enddo

! part4 - wetness
!
! * assume dry for now, have passed in soil moisture to this routine
!   but needs to be included here
! * wetfactor less than 1 would reduce dustflux
! * convert model soil moisture (m3/m3) to gravimetric soil moisture
!   (mass of water / mass of soil in %) assuming a constant density 
!   for soil
         pclay=beta(1,izob)*100.
         wp=0.0014*pclay*pclay+0.17*pclay
         smois_grav=(smois(i,1,j)/2.6)*100.
         if(smois_grav.gt.wp) then
           wetfactor=sqrt(1.0+1.21*(smois_grav-wp)**0.68)
         else
           wetfactor=1.0
         endif
!        wetfactor=1.0

! part5 - dustflux
! lower bound on ustar = 20 cm/s as in Shaw et al, but set upper
! bound to 100 cm/s
         c_const=1.e-14  ! default

         ustar1=ust(i,j)*100.0
         if(ustar1.gt.100.0) ustar1=100.0
         ustart0=20.0
         ustart=ustart0*wetfactor
         if(ustar1.le.ustart) then
           dustflux=0.0
         else
           dustflux=c_const*(ustar1**4)*(1.0-(ustart/ustar1))
         endif
         dustflux=dustflux*10.0
! units kg m-2 s-1
         ftot=0.0
         do ii=1,2
           ftot=ftot+dustflux*alphamask*delta(ii)
         enddo
! convert to ug m-2 s-1  
         ftot=ftot*1.0e+09
         e_dust(i,j) = ftot

!   apportion other inorganics only
         factaa = alt(i,k,j)/dz8w(i,k,j)
         factbb = factaa * ftot
         fracoin = 0.97
         fracca = 0.03*0.4
         fracco3 = 0.03*0.6
         fractot = fracoin + fracca + fracco3
         
         do n = 1,num_chem
            if (mod(n,8) < 1) then
                tend_emis(i,k,j,n) = factbb * sz(8) * fractot 
            end if
            if (mod(n,8) >= 1) then
                tend_emis(i,k,j,n) = factbb * sz(mod(n,8)) * fractot 
            end if
         end do
!   if (ivgtyp(i,j) .eq. 8) print*,'jdf',i,j,ustar1,ustart0,factaa,ftot
!         do 1810 n = 1, nsize_aer(itype)
!            l_oin = lptr_oin_aer(n,itype,iphase)
!            l_ca = lptr_ca_aer(n,itype,iphase)
!            l_co3 = lptr_co3_aer(n,itype,iphase)
!            if (l_oin >= p1st) then
!               chem(i,k,j,l_oin) = chem(i,k,j,l_oin) +      &
!               factbb * sz(n) * fracoin
!               if (l_ca >= p1st)                            &
!                  chem(i,k,j,l_ca) = chem(i,k,j,l_ca) +     &
!                  factbb * sz(n) * fracca
!               if (l_co3 >= p1st)                           &
!                  chem(i,k,j,l_co3) = chem(i,k,j,l_co3) +   &
!                  factbb * sz(n) * fracco3
!! mass1part is mass of a single particle in ug, density of dust ~2.5 g cm-3
!               densdust=2.5
!               mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
!               l = numptr_aer(n,itype,iphase)
!               if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) +   &
!                        factbb * sz(n) * fractot / mass1part
!            end if
!1810    continue

1840    continue

1820    continue
1830    continue

        return

   END subroutine mosaic_dust_emiss

!====================================================================================
!add another dust emission scheme following GOCART mechanism  --czhao  09/17/2009
!====================================================================================
  subroutine mosaic_dust_gocartemis (debug_on,ktau,dt,num_soil_layers,alt,u_phy,    &
         v_phy,rho_phy,dz8w,smois,u10,v10,p8w,erod,                   &
         ivgtyp,isltyp,xland,xlat,xlong,dx,g,emis_dust_factor,             &
         tend_emis,e_dust,                                                 &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )

  USE module_data_gocart_dust

  IMPLICIT NONE

   logical :: debug_on
   INTEGER,      INTENT(IN   ) :: ktau, num_soil_layers,           &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,               &
          INTENT(IN   ) ::                                                 &
                                                     ivgtyp,               &
                                                     isltyp

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ) ,                      &
          INTENT(INOUT) ::   tend_emis
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
          INTENT(INOUT   ) ::    e_dust
   !czhao

   REAL, INTENT(IN) :: emis_dust_factor 

  REAL, DIMENSION( ims:ime, num_soil_layers, jms:jme ) ,      &
      INTENT(INOUT) ::                               smois
   REAL,  DIMENSION( ims:ime , jms:jme, 3 )                   ,               &
          INTENT(IN   ) ::    erod
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
          INTENT(IN   ) ::                                                 &
                                                     u10,                  &
                                                     v10,                  &
                                                     xland,                &
                                                     dx,                   &
                                                     xlat,                 &
                                                     xlong                 
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::                                                 &
                                                        alt,               &
                                                     dz8w,p8w,             &
                                              u_phy,v_phy,rho_phy

   REAL, INTENT(IN   ) :: dt,g

!
! local variables
!
  integer :: nmx,i,j,k,ndt,imx,jmx,lmx
  integer ilwi, start_month
  real*8, DIMENSION (3) :: erodin
  real*8, DIMENSION (5) :: bems
  real*8, DIMENSION( ims:ime , jms:jme ) :: dxy
  real*8  w10m,gwet,airden,airmas
  real*8  cdustemis,jdustemis,cdustcon,jdustcon
  real*8  cdustdens,jdustdens,mass1part,jdustdiam,cdustdiam,dp_meanvol_tmp
  real factaa, factbb, fracoin, fracca, fracco3, fractot
  real*8  conver,converi
  real dttt
  real soilfacj,rhosoilj,rhosoilc
  integer p1st,l_dust,l_ca,l_co3,l,n,l_tr1r1,l_tr1r2,l_tr1r3,l_tr1r4
  real :: densdust
  real sz(10),ftot,frac_10um
  real totalemis,accfrac,corfrac,rscale
  integer iphase, itype

  rscale=1.01  !  to account for the emitted dust with radius larger than 10 um
  accfrac=dust_emass_jfrac  ! assign 7% to accumulation mode
  corfrac=dust_emass_cfrac  ! assign 93% to coarse mode
  accfrac=accfrac*rscale
  corfrac=corfrac*rscale
  p1st = param_first_scalar

!czhao give an example of MOSAIC 8 size bins here (unit: um) from Jerome 2005 paper
! 0.039-0.078; 0.078-0.156; 0.156-0.312; 0.312-0.625; 0.625-1.25; 1.25-2.5; 2.5-5.0; 5.0-10.0
!in the model the size bound for 8 bins is calculated in the code using dlo_sec and dhi_sec
!size3
!       sz(1)=0.0
!       sz(2)=0.0
!       sz(3)=0.0003
!       sz(4)=0.0048
!       sz(5)=0.0130
!       sz(6)=0.0250
!       sz(7)=0.1400
!       sz(8)=0.8169
!the size is following the best match2 of modal distribution
!we have ~3.6% within 1 um radius and ~7.4% within 1.25 um radius
!versus ~6.0% within 1 um radius following GOCART original and ~9.3% within 1.25 um radius
! totally ~82% are within the MOSAIC size range 0.03-10 um diameter
! now the dust distribution follows the parameter defined in module_data_sorgam.F 
        sz(1)=1.0e-8
        sz(2)=1.0e-6
        sz(3)=3.0e-4
        sz(4)=3.5e-3
        sz(5)=0.018
        sz(6)=0.070
        sz(7)=0.2595
        sz(8) =0.42
        sz(9) =0.20
        sz(10)=0.0287
!!based on Kok et al., 2011, we can get another size distribution with lambda=12 and Cv=12.62  -czhao
!       sz(1)=7.7e-7
!       sz(2)=1.5e-5
!       sz(3)=1.9e-4
!       sz(4)=1.8e-3
!       sz(5)=0.0119
!       sz(6)=0.0574
!       sz(7)=0.2018
!       sz(8)=0.4564
!       !sz(9)=0.2702   !10um-20um in diameter


!!based on Kok et al., 2011, we can get another size distribution with lambda=8 and Cv=6.83  -czhao
!       sz(1)=1.43e-6
!       sz(2)=2.71e-5
!       sz(3)=3.57e-4
!       sz(4)=3.34e-3
!       sz(5)=0.0220
!       sz(6)=0.1048
!       sz(7)=0.3425
!       sz(8)=0.4782
!       !sz(9)=0.0484   !10um-20um in diameter


!        frac_10um=sum(sz(1:8))
!        sz(:)=sz(:)/sum(sz(1:8))  ! relative ratio in diameter range of 0-10 um

!  for now just do itype=1
        itype = 1
        iphase = ai_phase

  conver=1.e-9
  converi=1.e9
!
! number of dust bins
  nmx=5
  k=kts
  do j=jts,jte
  do i=its,ite


!
! don't do dust over water!!!
     if(xland(i,j).lt.1.5)then

     ilwi=1
     start_month = 3   ! it doesn't matter, ch_dust is not a month dependent now, a constant
     w10m=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
!    airmas=-(p8w(i,kts+1,j)-p8w(i,kts,j))*dx*dx/g   ! kg 
     airmas=1.0  ! czhao initialize here but not used any more

! we don't trust the u10,v10 values, if model layers are very thin near surface
     if(dz8w(i,kts,j).lt.12.)w10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
    !erodin(1)=erod(i,j,1)/dx/dx   ! czhao erod shouldn't be scaled to the area, because it's a fraction
    !erodin(2)=erod(i,j,2)/dx/dx
    !erodin(3)=erod(i,j,3)/dx/dx
     erodin(1)=erod(i,j,1)
     erodin(2)=erod(i,j,2)
     erodin(3)=erod(i,j,3)
!
!  volumetric soil moisture over porosity

!jwfeng : Add porosity through lsm_noah_init.
    porosity(isltyp(i,j))=maxsmc(isltyp(i,j))

     gwet=smois(i,1,j)/porosity(isltyp(i,j))
     ndt=ifix(dt)
     airden=rho_phy(i,kts,j)
     dxy(i,j)=dx(i,j)*dx(i,j)

    call mosaic_source_du( nmx, dt, &
                     erodin, ilwi, dxy(i,j), w10m, gwet, airden, airmas, &
                     bems,start_month,g)

!bems: kg/timestep/cell
    !sum up the dust emission from 0.1-10 um in radius 
    ! unit change from kg/timestep/cell to ug/m2/s
    totalemis=(sum(bems(1:5))/dt)*converi/dxy(i,j)
    !tune total dust emission, previously it's done for ch_dust
    totalemis=totalemis*emis_dust_factor
    !totalemis=totalemis*rscale !to account for the particles larger than 10 um 
                              ! based on assumed size distribution
    jdustemis = totalemis*accfrac   ! accumulation mode
    cdustemis = totalemis*corfrac   ! coarse mode 

    e_dust(i,j)=cdustemis+jdustemis
    ftot=jdustemis+cdustemis  ! total dust emission in ug/m2/s 
!    ftot=ftot*frac_10um  ! only in diameter range of 0-10um

!   apportion other inorganics only
! jwfeng: Using tend_emis instead of directly changing chem,
!         Dont multi dt in factaa !!!Important
!        factaa = (dt/dz8w(i,k,j))*alt(i,k,j) 
         factaa = alt(i,k,j)/dz8w(i,k,j)
         factbb = factaa * ftot

         fracoin = 0.99
         fracca = 0.01*0.4
         fracco3 = 0.01*0.6
         fractot = fracoin

         do 1810 n = 1, nsize_aer(itype)

            l_dust = lptr_dust_aer(n,itype,iphase)
!jwfeng: Only dust for now
!l_ca = lptr_ca_aer(n,itype,iphase)
!l_co3 = lptr_co3_aer(n,itype,iphase)
!l_tr1r1 = lptr_tr1r1_aer(n,itype,iphase)
!l_tr1r2 = lptr_tr1r2_aer(n,itype,iphase)
!l_tr1r3 = lptr_tr1r3_aer(n,itype,iphase)
!l_tr1r4 = lptr_tr1r4_aer(n,itype,iphase)

!           if  (l_ca > p1st) then
!            fracoin = 0.99
!            fracca = 0.01*0.4
!            fracco3 = 0.01*0.6
!            fractot = fracoin
!           else 
             fracoin = 1.0
             fracca = 0.0
             fracco3 = 0.0
             fractot = fracoin
!           endif

            if (l_dust > p1st) then
! jwfeng:04/2020 Add bin 5,9,10
                if (nsize_aer(itype).eq.4 .or. nsize_aer(itype).eq.5) then
                   tend_emis(i,k,j,l_dust) = factbb * (sz(2*n-1)+sz(2*n)) * fracoin
                   if (n <= 3) densdust=2.5
                   if (n > 3 ) densdust=2.65
                   ! mass1part is mass of a single particle in ug
                   mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
                   l = numptr_aer(n,itype,iphase)
                   if (l > p1st) &
                       tend_emis(i,k,j,l) = tend_emis(i,k,j,l) + factbb * (sz(2*n-1)+sz(2*n)) * fractot/mass1part

                else if (nsize_aer(itype).eq.8 .or. nsize_aer(itype).eq.9 .or. nsize_aer(itype).eq.10) then
                    tend_emis(i,k,j,l_dust) = factbb * sz(n) * fracoin
                    if (n <= 5) densdust=2.5
                    if (n > 5 ) densdust=2.65
                 !  mass1part is mass of a single particle in ug
                    mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
                    l = numptr_aer(n,itype,iphase)
                    if (l > p1st)  &
                        tend_emis(i,k,j,l) = tend_emis(i,k,j,l) + factbb * sz(n) * fractot / mass1part
                endif !numptr_aer
            endif !l_dust

!JWFeng : No source tracking for now !
            !for source tracking
!           if (l_tr1r1 > p1st) then
!            if (nsize_aer(itype).eq.4) then
!             if ((xlong(i,j).ge.-140.and.xlong(i,j).le.-80).and. & !region 1 for US 
!                 (xlat(i,j).ge.15.and.xlat(i,j).le.50)) then
!             chem(i,k,j,l_tr1r1) = chem(i,k,j,l_tr1r1) +      &
!             factbb * (sz(2*n-1)+sz(2*n)) * fracoin
!             else if ((xlong(i,j).ge.75.and.xlong(i,j).le.150).and. & !region 2 for East Asia
!                      (xlat(i,j).ge.25.and.xlat(i,j).le.50)) then 
!             chem(i,k,j,l_tr1r2) = chem(i,k,j,l_tr1r2) +      &
!             factbb * (sz(2*n-1)+sz(2*n)) * fracoin
!             else if ((xlong(i,j).ge.-20.and.xlong(i,j).le.35).and. & !region 3 for North Africa 
!                      (xlat(i,j).ge.0.and.xlat(i,j).le.40)) then
!             chem(i,k,j,l_tr1r3) = chem(i,k,j,l_tr1r3) +      &
!             factbb * (sz(2*n-1)+sz(2*n)) * fracoin
!             else !region 4 for all others 
!             chem(i,k,j,l_tr1r4) = chem(i,k,j,l_tr1r4) +      &
!             factbb * (sz(2*n-1)+sz(2*n)) * fracoin
!             endif
!             if (l_ca > p1st) chem(i,k,j,l_ca)=chem(i,k,j,l_ca)+ &
!             factbb*(sz(2*n-1)+sz(2*n))*fracca
!             if (l_co3 > p1st) chem(i,k,j,l_co3)=chem(i,k,j,l_co3)+ &
!             factbb*(sz(2*n-1)+sz(2*n))*fracco3
!             if (n <= 3) densdust=2.5
!             if (n > 3 ) densdust=2.65
!             ! mass1part is mass of a single particle in ug    
!             mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
!             l = numptr_aer(n,itype,iphase)
!             if (l > p1st) chem(i,k,j,l) = chem(i,k,j,l) +   &
!                        factbb*(sz(2*n-1)+sz(2*n))*fractot/mass1part
!            else if (nsize_aer(itype).eq.8) then
!             if ((xlong(i,j).ge.-140.and.xlong(i,j).le.-80).and. & !region 1 for US 
!                 (xlat(i,j).ge.15.and.xlat(i,j).le.50)) then
!             chem(i,k,j,l_tr1r1) = chem(i,k,j,l_tr1r1) +      &
!             factbb * sz(n) * fracoin
!             else if ((xlong(i,j).ge.75.and.xlong(i,j).le.150).and. & !region 2 for East Asia
!                      (xlat(i,j).ge.25.and.xlat(i,j).le.50)) then
!             chem(i,k,j,l_tr1r2) = chem(i,k,j,l_tr1r2) +      &
!             factbb * sz(n) * fracoin
!             else if ((xlong(i,j).ge.-20.and.xlong(i,j).le.35).and. & !region 3 for North Africa 
!                      (xlat(i,j).ge.0.and.xlat(i,j).le.40)) then
!             chem(i,k,j,l_tr1r3) = chem(i,k,j,l_tr1r3) +      &
!             factbb * sz(n) * fracoin
!             else    ! region 4 for all others
!             chem(i,k,j,l_tr1r4) = chem(i,k,j,l_tr1r4) +      &
!             factbb * sz(n) * fracoin
!             endif
!             if (l_ca > p1st) chem(i,k,j,l_ca)=chem(i,k,j,l_ca)+ &
!             factbb*sz(n)*fracca
!             if (l_co3 > p1st) chem(i,k,j,l_co3)=chem(i,k,j,l_co3)+ &
!             factbb*sz(n)*fracco3
!             if (n <= 5) densdust=2.5
!             if (n > 5 ) densdust=2.65
!             ! mass1part is mass of a single particle in ug
!             mass1part=0.523598*(dcen_sect(n,itype)**3)*densdust*1.0e06
!             l = numptr_aer(n,itype,iphase)
!             if (l > p1st) chem(i,k,j,l) = chem(i,k,j,l) +   &
!                       factbb * sz(n) * fractot / mass1part
!            endif !nsize_aer
!           endif !l_tr1r1

1810    continue

     endif ! xland
  enddo ! i
  enddo ! j

end subroutine mosaic_dust_gocartemis

  SUBROUTINE mosaic_source_du( nmx, dt1, &
                     erod, ilwi, dxy, w10m, gwet, airden, airmas, &
                     bems,month,g0)

! ****************************************************************************
! *  Evaluate the source of each dust particles size classes  (kg/m3)        
! *  by soil emission.
! *  Input:
! *         EROD      Fraction of erodible grid cell                (-)
! *                   for 1: Sand, 2: Silt, 3: Clay
! *         DUSTDEN   Dust density                                  (kg/m3)
! *         DXY       Surface of each grid cell                     (m2)
! *         AIRVOL    Volume occupy by each grid boxes              (m3)
! *         NDT1      Time step                                     (s)
! *         W10m      Velocity at the anemometer level (10meters)   (m/s)
! *         u_tresh   Threshold velocity for particule uplifting    (m/s)
! *         CH_dust   Constant to fudge the total emission of dust  (s2/m2)
! *      
! *  Output:
! *         DSRC      Source of each dust type           (kg/timestep/cell) 
! *
! *  Working:
! *         SRC       Potential source                   (kg/m/timestep/cell)
! *
! ****************************************************************************

 USE module_data_gocart_dust

  INTEGER, INTENT(IN)    :: nmx
  REAL*8,    INTENT(IN)  :: erod(ndcls)
  INTEGER, INTENT(IN)    :: ilwi,month

  REAL*8,    INTENT(IN)    :: w10m, gwet
  REAL*8,    INTENT(IN)    :: dxy
  REAL*8,    INTENT(IN)    :: airden, airmas
  REAL*8,    INTENT(OUT)   :: bems(nmx)

  REAL*8    :: den(nmx), diam(nmx)
  REAL*8    :: tsrc, u_ts0, cw, u_ts, dsrc, srce
  REAL, intent(in)    :: g0
  REAL    :: rhoa, g,dt1
  INTEGER :: i, j, n, m, k

  ! default is 1 ug s2 m-5 == 1.e-9 kg s2 m-5
  !ch_dust(:,:)=0.8D-9   ! ch_dust is defined here instead of in the chemics_ini.F if with SORGAM  -czhao
  !ch_dust(:,:)=1.0D-9  ! default 
!  ch_dust(:,:)=0.50D-9  ! ch_dust is scaled for Sahara 
   ch_dust(:,:)=1.0D-9   ! now totalemis is tuned by emis_dust_factor 


  ! executable statemenst
  DO n = 1, nmx
     ! jwfeng:
     ! Use Marticorena and Bergametti (1995). to calculate threshold velocity instead of the original Bagnold(1941)
     den(n) = den_dust(n)*1.0D-3
     diam(n) = 2.0*reff_dust(n)*1.0D2
     g = g0*1.0E2
     ! Pointer to the 3 classes considered in the source data files
     m = ipoint(n)
     tsrc = 0.0
              rhoa = airden*1.0D-3
              u_ts0 = 0.13*1.0D-2*SQRT(den(n)*g*diam(n)/rhoa)* &
                   SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
                   SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0)

              ! Case of surface dry enough to erode
!             IF (gwet < 0.5) THEN  !  Pete's modified value
              IF (gwet < 0.2) THEN
                 u_ts = 1.20D0 * MAX(0.0D0, u_ts0*(1.2D0 + 2.0D-1*LOG10(MAX(1.0D-3, gwet))))
              ELSE
                 ! Case of wet surface, no erosion
                 u_ts = 100.0
              END IF
              srce = frac_s(n)*erod(m)*dxy  ! (m2)
              IF (ilwi == 1 ) THEN
                 dsrc = ch_dust(n,month)*srce*w10m**2 &
                      * (w10m - u_ts)*dt1  ! (kg)
              ELSE
                 dsrc = 0.0
              END IF
              IF (dsrc < 0.0) dsrc = 0.0

              ! Update dust mixing ratio at first model level.
              !tc(n) = tc(n) + dsrc / airmas    !kg/kg-dryair -czhao
              bems(n) = dsrc     ! kg/timestep/cell

  ENDDO

END SUBROUTINE mosaic_source_du

!----------------------------------------------------------------------
subroutine mosaic_seasalt_emiss(                                       &
           dtstep, u10, v10, alt, dz8w, xland,                         &
           tend_emis, e_seas, emis_seas_factor,                        &
           ids,ide, jds,jde, kds,kde,                                  &
           ims,ime, jms,jme, kms,kme,                                  &
           its,ite, jts,jte, kts,kte, seasalt_emiss_active             )
!
! adds seasalt emissions for mosaic aerosol species
! (i.e., seasalt emissions tendencies over time dtstep are applied 
! to the aerosol mixing ratios)
!

   USE mpas_atmchem_vars, only:  num_chem, param_first_scalar

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   ) ::    dtstep, emis_seas_factor
   INTEGER, INTENT(IN) ::    seasalt_emiss_active
! 10-m wind speed components (m/s)
   REAL,  DIMENSION( ims:ime , jms:jme ),                                  &
          INTENT(IN   ) ::   u10, v10, xland
! trace species mixing ratios (aerosol mass = ug/kg; number = #/kg)
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::   tend_emis
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT  ) :: e_seas
! alt  = 1.0/(dry air density) in (m3/kg)
! dz8w = layer thickness in (m)
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )                       , &
          INTENT(IN   ) ::   alt, dz8w

! local variables
    integer i, j, k, l, n, l_oc, l_seas
    integer iphase, itype
    integer p1st

    real dum, dumdlo, dumdhi, dumoceanfrac, dumspd10
    real factaa, factbb, fracseas, fracorg

    real :: ssemfact_numb( maxd_asize, maxd_atype )
    real :: ssemfact_mass( maxd_asize, maxd_atype )


    ! Factors for Fuentes et al (ACP, 2011, doi:10.5194/acp-11-2585-2011) 
    !                                              seasalt emission scheme.
    ! These are for calculating the seawater_OC content dependent factors.
    real, save :: alpha_f1(4) =   &
                 (/ 12.328, 38.077, 102.31, 281.65 /)
    real, save :: alpha_f2(4) =   &
                 (/ 2.2958, 8.0935, 25.251, 46.80 /)
    real, save :: alpha_f3(4) =   &
                 (/ 0.00452, 0.00705, 0.00080, 0.000761 /)
    real, save :: beta_f1(4) =   &
                 (/ 0.0311, -0.031633, 0.013154, -0.0017762 /)
    real, save :: beta_f2(4) =   &
                 (/ -13.916, 35.73, -9.7651, 1.1665 /)
    real, save :: beta_f3(4) =   &
                 (/ 4747.8, 12920.0, 7313.4, 6610.0 /)
    real :: nti(4), dp0gi(4)
    ! seawater OC<0.2um concentration (in uM) - first is for low activity, the 2nd for high activity
    !    High activity conc is from the average for RHaMBLe cruise in high activity regions
    real, save :: oc02um(2) = (/ 0.0, 280.0 /)        
    ! Organic fraction for seasalt emissions (by size bin).
    ! These are estimated from Figure 5 of Fuentes et al (ACP, 2011), assuming a value of 
    !        2 ug/l for Chl-a for the high activity
    real, save :: org_frac_low_activity(8) =   &
       (/ 0.05,  0.05,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 /)
    real, save :: org_frac_high_activity(8) =   &
       (/ 0.2,  0.2,  0.1,  0.01,  0.01,  0.01,  0.01,  0.01 /)
    

    p1st = PARAM_FIRST_SCALAR

!   for now just do itype=1
    itype = 1
    iphase = ai_phase

    ! if using Feuntes et al (2011) then calculate the seawater OC dependent factors
    if(seasalt_emiss_active .eq. 3 .or. seasalt_emiss_active .eq. 4)then
        do i=1,4
            nti(i)   = beta_f1(i) * oc02um(seasalt_emiss_active-2)**2.0 + beta_f2(i) &
                            * oc02um(seasalt_emiss_active-2) + beta_f3(i)
            dp0gi(i) = alpha_f1(i) + alpha_f2(i) &
                            * exp(-alpha_f3(i)*oc02um(seasalt_emiss_active-2))
        end do
    end if


!   compute emissions factors for each size bin
    do n = 1, nsize_aer(itype)
!     changed the lower bound of the emission size limit to match lowest section edge, Qing.Yang@pnl.gov
        ! DL - 30/3/2012 - select emission scheme (1=original, 3&4=Feuntes et al, 2011)
        if(seasalt_emiss_active == 1)then
!             changed the lower bound of the emission size limit to match lowest section edge, Qing.Yang@pnl.gov
!            dumdlo = max( dlo_sect(n,itype), 0.1e-4 )
!            dumdhi = max( dhi_sect(n,itype), 0.1e-4 )
            dumdlo = max( dlo_sect(n,itype), 0.02e-4 )
            dumdhi = max( dhi_sect(n,itype), 0.02e-4 )
            call seasalt_emitfactors_1bin( 1, dumdlo, dumdhi,   &
            ssemfact_numb(n,itype), dum, ssemfact_mass(n,itype) )
        elseif(seasalt_emiss_active .eq. 3 .or. seasalt_emiss_active .eq. 4)then
            call seasalt_emit_feuntes_1bin( dlo_sect(n,itype), dhi_sect(n,itype),  &
                ssemfact_numb(n,itype), dum, ssemfact_mass(n,itype), nti, dp0gi, oc02um(seasalt_emiss_active-2) )
        endif

!   convert mass emissions factor from (g/m2/s) to (ug/m2/s)
        ssemfact_mass(n,itype) = ssemfact_mass(n,itype)*1.0e6
        !jwfeng
        ssemfact_mass(n,itype) = ssemfact_mass(n,itype)*emis_seas_factor
        ssemfact_numb(n,itype) = ssemfact_numb(n,itype)*emis_seas_factor
    end do


!   loop over i,j and apply seasalt emissions
    k = kts
    do j = jts, jte
    do i = its, ite
    
        ! Skip this point if over land. xland=1 for land and 2 for water.
        ! Also, there is no way to differentiate fresh from salt water.
        ! Currently, this assumes all water is salty.
        if( xland(i,j) < 1.5 ) cycle
    
        !wig: As far as I can tell, only real.exe knows the fractional breakdown
        !     of land use. So, in wrf.exe, dumoceanfrac will always be 1.
        dumoceanfrac = 1. !fraction of grid i,j that is salt water
        dumspd10 = dumoceanfrac* &
             ( (u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j))**(0.5*3.41) )

!store the emissions  (ug/m2/s)
        e_seas(i,j)=0.0
        do n = 1, nsize_aer(itype)
            e_seas(i,j) = e_seas(i,j) + &
                          dumspd10*ssemfact_mass(n,itype)
        enddo
        ! factaa is (s*m2/kg-air)
        ! factaa*ssemfact_mass(n) is (s*m2/kg-air)*(ug/m2/s) = ug/kg-air
        ! factaa*ssemfact_numb(n) is (s*m2/kg-air)*( #/m2/s) =  #/kg-air
! jwfeng: Using tend_emis instead of directly changing chem,
!         Dont multi dt in factaa !!!Important
!       factaa = (dt/dz8w(i,k,j))*alt(i,k,j) 
        factaa = alt(i,k,j)/dz8w(i,k,j)
        factbb = factaa * dumspd10

        if(seasalt_emiss_active == 1)then
            ! apportion seasalt mass emissions assumming that seasalt is pure nacl
            fracseas = 1.0
            do n = 1, nsize_aer(itype)
                ! only apply emissions if bin has seas species
                l_seas = lptr_seas_aer(n,itype,iphase)
                if ( l_seas >= p1st ) then
                    tend_emis(i,k,j,l_seas) = factbb * ssemfact_mass(n,itype) * fracseas
                    l = numptr_aer(n,itype,iphase)
                    if (l >= p1st) tend_emis(i,k,j,l) = tend_emis(i,k,j,l) + &
                                                        factbb * ssemfact_numb(n,itype)
                end if
            end do ! n = 1, nsize_aer(itype)

        elseif (seasalt_emiss_active.eq.3 .or. seasalt_emiss_active.eq.4) then
            do n = 1, nsize_aer(itype)
                ! apportion seasalt mass emissions assumming that seasalt is a 
                ! simple mix of pure nacl and a single primary organic 
                if(seasalt_emiss_active.eq.3)then
                    fracorg = org_frac_low_activity(n)
                elseif(seasalt_emiss_active.eq.4)then
                    fracorg = org_frac_high_activity(n)
                endif
                fracseas = 1.0 - fracorg
                ! only apply emissions if bin has both na and cl species
                l_oc   = lptr_oc_aer(n,itype,iphase)
                l_seas = lptr_seas_aer(n,itype,iphase)
                if ((l_seas >= p1st) .and. (l_oc >= p1st)) then
                    tend_emis(i,k,j,l_seas) = factbb * ssemfact_mass(n,itype) * fracseas
                    tend_emis(i,k,j,l_oc) = factbb * ssemfact_mass(n,itype) * fracorg

                    l = numptr_aer(n,itype,iphase)
                    if (l >= p1st) tend_emis(i,k,j,l) = tend_emis(i,k,j,l) + &
                                                        factbb * ssemfact_numb(n,itype)
                end if
            end do ! n = 1, nsize_aer(itype)
        endif

    end do ! i = its, ite
    end do ! j = jts, jte

    return

END subroutine mosaic_seasalt_emiss

!c----------------------------------------------------------------------
!c   following is from gong06b.f in
!c    /net/cirrus/files1/home/rce/oldfiles1/box/seasaltg
!c----------------------------------------------------------------------
subroutine seasalt_emitfactors_1bin( ireduce_smallr_emit,    &
           dpdrylo_cm, dpdryhi_cm,      &
           emitfact_numb, emitfact_surf, emitfact_mass )
!c
!c   computes seasalt emissions factors for a specifed 
!c   dry particle size range
!c    dpdrylo_cm  = lower dry diameter (cm)
!c    dpdryhi_cm  = upper dry diameter (cm)
!c
!c   number and mass emissions are then computed as
!c    number   emissions (#/m2/s)   == emitfact_numb * (spd10*3.41)
!c    dry-sfc  emissions (cm2/m2/s) == emitfact_surf * (spd10*3.41)
!c    dry-mass emissions (g/m2/s)   == emitfact_mass * (spd10*3.41)
!c
!c   where spd10 = 10 m windspeed in m/s
!c
!c   uses bubble emissions formula (eqn 5a) from 
!c    Gong et al. [JGR, 1997, p 3805-3818]
!c
!c   *** for rdry < rdry_star, this formula overpredicts emissions.
!c    A strictly ad hoc correction is applied to the formula,
!c    based on sea-salt size measurements of
!c    O'Dowd et al. [Atmos Environ, 1997, p 73-80]
!c
!c   *** the correction is only applied when ireduce_smallr_emit > 0
!c
    implicit none

!c   subr arguments
    integer ireduce_smallr_emit
    real dpdrylo_cm, dpdryhi_cm,                &
                emitfact_numb, emitfact_surf, emitfact_mass

!c   local variables
    integer isub_bin, nsub_bin

    real alnrdrylo
    real drydens, drydens_43pi_em12, x_4pi_em8
    real dum, dumadjust, dumb, dumexpb
    real dumsum_na, dumsum_ma, dumsum_sa
    real drwet, dlnrdry
    real df0drwet, df0dlnrdry, df0dlnrdry_star
    real relhum
    real rdry, rdrylo, rdryhi, rdryaa, rdrybb
    real rdrylowermost, rdryuppermost, rdry_star
    real rwet, rwetaa, rwetbb
    real rdry_cm, rwet_cm
    real sigmag_star
    real xmdry, xsdry

    real pi
    parameter (pi = 3.1415936536)

!c   c1-c4 are constants for seasalt hygroscopic growth parameterization
!c   in Eqn 3 and Table 2 of Gong et al. [1997]
    real c1, c2, c3, c4, onethird
    parameter (c1 = 0.7674)
    parameter (c2 = 3.079)
    parameter (c3 = 2.573e-11)
    parameter (c4 = -1.424)
    parameter (onethird = 1.0/3.0)


!c   dry particle density (g/cm3)
    drydens = 2.165
!c   factor for radius (micrometers) to mass (g)
    drydens_43pi_em12 = drydens*(4.0/3.0)*pi*1.0e-12
!c   factor for radius (micrometers) to surface (cm2)
    x_4pi_em8 = 4.0*pi*1.0e-8
!c   bubble emissions formula assume 80% RH
    relhum = 0.80

!c   rdry_star = dry radius (micrometers) below which the
!c   dF0/dr emission formula is adjusted downwards
    rdry_star = 0.1
    if (ireduce_smallr_emit .le. 0) rdry_star = -1.0e20
!c   sigmag_star = geometric standard deviation used for
!c   rdry < rdry_star
    sigmag_star = 1.9

!c   initialize sums
    dumsum_na = 0.0
    dumsum_sa = 0.0
    dumsum_ma = 0.0

!c   rdrylowermost, rdryuppermost = lower and upper 
!c   dry radii (micrometers) for overall integration
        rdrylowermost = dpdrylo_cm*0.5e4
        rdryuppermost = dpdryhi_cm*0.5e4

!c
!c   "section 1"
!c   integrate over rdry > rdry_star, where the dF0/dr emissions 
!c   formula is applicable
!c   (when ireduce_smallr_emit <= 0, rdry_star = -1.0e20,
!c   and the entire integration is done here)
!c
    if (rdryuppermost .le. rdry_star) goto 2000

!c   rdrylo, rdryhi = lower and upper dry radii (micrometers) 
!c   for this part of the integration
        rdrylo = max( rdrylowermost, rdry_star )
        rdryhi = rdryuppermost

    nsub_bin = 1000

    alnrdrylo = log( rdrylo )
    dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

!c   compute rdry, rwet (micrometers) at lowest size
    rdrybb = exp( alnrdrylo )
    rdry_cm = rdrybb*1.0e-4
    rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/        &
              ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
    rwetbb = rwet_cm*1.0e4

    do 1900 isub_bin = 1, nsub_bin

!c   rdry, rwet at sub_bin lower boundary are those
!c   at upper boundary of previous sub_bin
    rdryaa = rdrybb
    rwetaa = rwetbb

!c   compute rdry, rwet (micrometers) at sub_bin upper boundary
    dum = alnrdrylo + isub_bin*dlnrdry
    rdrybb = exp( dum )

    rdry_cm = rdrybb*1.0e-4
    rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/        &
              ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
    rwetbb = rwet_cm*1.0e4

!c   geometric mean rdry, rwet (micrometers) for sub_bin
    rdry = sqrt(rdryaa * rdrybb)
    rwet = sqrt(rwetaa * rwetbb)
    drwet = rwetbb - rwetaa

!c   xmdry is dry mass in g
    xmdry = drydens_43pi_em12 * (rdry**3.0)

!c   xsdry is dry surface in cm2
    xsdry = x_4pi_em8 * (rdry**2.0)

!c   dumb is "B" in Gong's Eqn 5a
!c   df0drwet is "dF0/dr" in Gong's Eqn 5a
    dumb = ( 0.380 - log10(rwet) ) / 0.650
    dumexpb = exp( -dumb*dumb)
    df0drwet = 1.373 * (rwet**(-3.0)) *             &
              (1.0 + 0.057*(rwet**1.05)) *             &
              (10.0**(1.19*dumexpb))

    dumsum_na = dumsum_na + drwet*df0drwet
    dumsum_ma = dumsum_ma + drwet*df0drwet*xmdry
    dumsum_sa = dumsum_sa + drwet*df0drwet*xsdry

1900    continue


!c
!c   "section 2"
!c   integrate over rdry < rdry_star, where the dF0/dr emissions 
!c   formula is just an extrapolation and predicts too many emissions
!c
!c   1.  compute dF0/dln(rdry) = (dF0/drwet)*(drwet/dlnrdry) 
!c    at rdry_star
!c   2.  for rdry < rdry_star, assume dF0/dln(rdry) is lognormal,
!c    with the same lognormal parameters observed in 
!c    O'Dowd et al. [1997]
!c
2000    if (rdrylowermost .ge. rdry_star) goto 3000

!c   compute dF0/dln(rdry) at rdry_star
    rdryaa = 0.99*rdry_star
    rdry_cm = rdryaa*1.0e-4
    rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/        &
              ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
    rwetaa = rwet_cm*1.0e4

    rdrybb = 1.01*rdry_star
    rdry_cm = rdrybb*1.0e-4
    rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/        &
              ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
    rwetbb = rwet_cm*1.0e4

    rwet = 0.5*(rwetaa + rwetbb)
    dumb = ( 0.380 - log10(rwet) ) / 0.650
    dumexpb = exp( -dumb*dumb)
    df0drwet = 1.373 * (rwet**(-3.0)) *             &
              (1.0 + 0.057*(rwet**1.05)) *             &
              (10.0**(1.19*dumexpb))

    drwet = rwetbb - rwetaa
    dlnrdry = log( rdrybb/rdryaa )
    df0dlnrdry_star = df0drwet * (drwet/dlnrdry)


!c   rdrylo, rdryhi = lower and upper dry radii (micrometers) 
!c   for this part of the integration
        rdrylo = rdrylowermost
        rdryhi = min( rdryuppermost, rdry_star )

    nsub_bin = 1000

    alnrdrylo = log( rdrylo )
    dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

    do 2900 isub_bin = 1, nsub_bin

!c   geometric mean rdry (micrometers) for sub_bin
    dum = alnrdrylo + (isub_bin-0.5)*dlnrdry
    rdry = exp( dum )

!c   xmdry is dry mass in g
    xmdry = drydens_43pi_em12 * (rdry**3.0)

!c   xsdry is dry surface in cm2
    xsdry = x_4pi_em8 * (rdry**2.0)

!c   dumadjust is adjustment factor to reduce dF0/dr
    dum = log( rdry/rdry_star ) / log( sigmag_star )
    dumadjust = exp( -0.5*dum*dum )

    df0dlnrdry = df0dlnrdry_star * dumadjust

    dumsum_na = dumsum_na + dlnrdry*df0dlnrdry
    dumsum_ma = dumsum_ma + dlnrdry*df0dlnrdry*xmdry
    dumsum_sa = dumsum_sa + dlnrdry*df0dlnrdry*xsdry

2900    continue


!c
!c  all done
!c
3000    emitfact_numb = dumsum_na
    emitfact_mass = dumsum_ma
    emitfact_surf = dumsum_sa

    return
end subroutine seasalt_emitfactors_1bin

!----------------------------------------------------------------------
!----------------------------------------------------------------------

!c----------------------------------------------------------------------
!c   Following is an adaption of the subroutine above.
!c   It has been modified to include flux from Fuentes et al, ACP, 2010
!c   for dealing with smaller particle emissions.
!c----------------------------------------------------------------------
subroutine seasalt_emit_feuntes_1bin( &
           dpdrylo_cm, dpdryhi_cm,      &
           emitfact_numb, emitfact_surf, emitfact_mass, nti, dp0gi, oc02um  )
!c
!c   computes seasalt emissions factors for a specifed 
!c   dry particle size range
!c    dpdrylo_cm  = lower dry diameter (cm)
!c    dpdryhi_cm  = upper dry diameter (cm)
!c
!c   number and mass emissions are then computed as
!c    number   emissions (#/m2/s)   == emitfact_numb * (spd10*3.41)
!c    dry-sfc  emissions (cm2/m2/s) == emitfact_surf * (spd10*3.41)
!c    dry-mass emissions (g/m2/s)   == emitfact_mass * (spd10*3.41)
!c
!c   where spd10 = 10 m windspeed in m/s
!c
!c  Uses bubble emissions formula (eqn 5a) from 
!c    Gong et al. [JGR, 1997, p 3805-3818] for particles 
!c  with dry diameters of 0.45 um or more
!c
!c  For smaller particles we use the parameterisation of 
!c  Fuentes et al, ACP, 2010.
!c


    implicit none

    ! subr arguments
    !integer ireduce_smallr_emit
    real dpdrylo_cm, dpdryhi_cm,                &
                emitfact_numb, emitfact_surf, emitfact_mass
    real, intent(in) :: nti(4), dp0gi(4), oc02um

    ! local variables
    integer isub_bin, nsub_bin, jd, nsub_lower_bin, nsub_upper_bin

    real alnrdrylo
    real drydens, drydens_43pi_em12, x_4pi_em8, drydens_16pi_em12, x_pi_em8
    real dum, dumadjust, dumb, dumexpb
    real dumsum_na, dumsum_ma, dumsum_sa
    real drwet, dlnrdry, ddwet, ddry, dwet
    real df0drwet, df0dlnrdry, df0dlnrdry_star
    real relhum
    real rdry, rdrylo, rdryhi, rdryaa, rdrybb
    real rdrylowermost, rdryuppermost, rdry_star
    real rwet, rwetaa, rwetbb
    real rdry_cm, rwet_cm
    real sigmag_star
    real xmdry, xsdry
    
    real ddrylo, ddryhi, alogddrylo
    real ddrybb, ddry_cm, dwet_cm, dwetbb
    real ddryaa, dwetaa, dlogddry, df0dlogddry

    real pi
    parameter (pi = 3.1415936536)

    ! c1-c4 are constants for seasalt hygroscopic growth parameterization
    ! in Eqn 3 and Table 2 of Gong et al. [1997]
    real c1, c2, c3, c4, onethird
    parameter (c1 = 0.7674)
    parameter (c2 = 3.079)
    parameter (c3 = 2.573e-11)
    parameter (c4 = -1.424)
    parameter (onethird = 1.0/3.0)


    !  constants for seasalt/organic aerosol distribution parameterisation
    !    from Fuentes et al, ACP, 2010
    real, save :: width_ssinorg_f(4) = &
                 (/ 1.55, 1.7, 1.5, 1.7 /)
    real, save :: width_ssorg_f(4) = &
                 (/ 1.55, 1.9, 1.5, 1.7 /)
    real :: width_ss_f(4), frac
    ! scale_factor is a combination of:
    !   1) (Q)  sweep air flow (58.3 cm3/s)  
    !   2) (Ab) total surface area covered by bubbles (0.0146 m2)
    !   4) scaling factor of 3.84e-6 from the whitecap coverage formulation of 
    !          Monahan and O'Muircheartaigh (1980) - as the Gong et al formulation
    !          includes a whitecap coverage factor too
    ! this is for converting from dNt to dFp for Fuentes et al (ACP, 2010)
    real, parameter :: scale_factor = 58.3/(0.0146)*3.84e-6 !0.01533


    ! select which distribution width to use for Fuentes et al,
    if(oc02um.gt.0e0)then
        width_ss_f = width_ssorg_f
    else
        width_ss_f = width_ssinorg_f
    end if


    ! dry particle density (g/cm3)
    drydens = 2.165
    ! factor for radius (micrometers) to mass (g)
    drydens_43pi_em12 = drydens*(4.0/3.0)*pi*1.0e-12
    ! factor for diameter (micrometers) to mass (g)
    drydens_16pi_em12 = drydens*(1.0/6.0)*pi*1.0e-12
    ! factor for radius (micrometers) to surface (cm2)
    x_4pi_em8 = 4.0*pi*1.0e-8
    ! factor for diameter (micrometers) to surface (cm2)
    x_pi_em8 = pi*1.0e-8
    ! bubble emissions formula assume 80% RH
    relhum = 0.80

    ! rdry_star = dry radius (micrometers) below which the
    ! dF0/dr emission formula from Fuentes et al, ACP, 2010
    ! will be used.
    rdry_star = 0.45 / 2.0
    !if (ireduce_smallr_emit .le. 0) rdry_star = -1.0e20
    ! sigmag_star = geometric standard deviation used for
    ! rdry < rdry_star
    sigmag_star = 1.9

    ! initialize sums
    dumsum_na = 0.0
    dumsum_sa = 0.0
    dumsum_ma = 0.0

    ! rdrylowermost, rdryuppermost = lower and upper 
    ! dry radii (micrometers) for overall integration
    rdrylowermost = dpdrylo_cm*0.5e4
    rdryuppermost = dpdryhi_cm*0.5e4


    ! "section 1"
    ! integrate over rdry > rdry_star, where the dF0/dr emissions 
    ! formula from the original subroutine is applicable
    !  (so using Gong et al for all emission profile)
    if (rdrylowermost .ge. rdry_star) then

        ! rdrylo, rdryhi = lower and upper dry radii (micrometers) 
        ! for this part of the integration
        rdrylo = rdrylowermost
        rdryhi = rdryuppermost

        nsub_bin = 1000

        alnrdrylo = log( rdrylo )
        dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

        ! compute rdry, rwet (micrometers) at lowest size
        rdrybb = exp( alnrdrylo )
        rdry_cm = rdrybb*1.0e-4
        rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/        &
                  ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
        rwetbb = rwet_cm*1.0e4

        do isub_bin = 1, nsub_bin

            ! rdry, rwet at sub_bin lower boundary are those
            ! at upper boundary of previous sub_bin
            rdryaa = rdrybb
            rwetaa = rwetbb

            ! compute rdry, rwet (micrometers) at sub_bin upper boundary
            dum = alnrdrylo + isub_bin*dlnrdry
            rdrybb = exp( dum )

            rdry_cm = rdrybb*1.0e-4
            rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/        &
                      ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
            rwetbb = rwet_cm*1.0e4

            ! geometric mean rdry, rwet (micrometers) for sub_bin
            rdry = sqrt(rdryaa * rdrybb)
            rwet = sqrt(rwetaa * rwetbb)
            drwet = rwetbb - rwetaa

            ! xmdry is dry mass in g
            xmdry = drydens_43pi_em12 * (rdry**3.0)

            ! xsdry is dry surface in cm2
            xsdry = x_4pi_em8 * (rdry**2.0)

            ! dumb is "B" in Gong's Eqn 5a
            ! df0drwet is "dF0/dr" in Gong's Eqn 5a
            dumb = ( 0.380 - log10(rwet) ) / 0.650
            dumexpb = exp( -dumb*dumb)
            df0drwet = 1.373 * (rwet**(-3.0)) *             &
                      (1.0 + 0.057*(rwet**1.05)) *             &
                      (10.0**(1.19*dumexpb))

            dumsum_na = dumsum_na + drwet*df0drwet
            dumsum_ma = dumsum_ma + drwet*df0drwet*xmdry
            dumsum_sa = dumsum_sa + drwet*df0drwet*xsdry

        end do


    ! "section 2"
    ! integrate over rdry < rdry_star for old parameterisation, and use
    ! Fuentes et al below that
    !
    ! 1.  for rdry < rdry_star, use Fuentes et al parameterisation
    ! 2.  for rdry > rdry_star, use Gong et al parameterisation, as above
    elseif (rdryuppermost .gt. rdry_star) then

        ! determine the fraction of size bin below rdry_star
        frac = (log(rdry_star)-log(rdrylowermost)) / (log(rdryuppermost)-log(rdrylowermost))
        nsub_lower_bin = floor(frac*1000.0)     ! calc number of size bins to use for section below rdry_star
        nsub_upper_bin = 1000-nsub_lower_bin ! use remaining size bins for section above rdry_star


    ! 1.................
        ! ddrylo, ddryhi = lower and upper dry diameter (micrometers) 
        ! for this part of the integration
        ddrylo = rdrylowermost*2.0
        ddryhi = rdry_star*2.0

        nsub_bin = nsub_lower_bin

        alogddrylo = log10( ddrylo )
        dlogddry = (log10( ddryhi ) - alogddrylo)/nsub_bin

        ! compute ddry, dwet (micrometers) at lowest size
        ddrybb = 10.0**( alogddrylo )
        ddry_cm = ddrybb*1.0e-4
        dwet_cm = ( ddry_cm**3 + (c1*(ddry_cm**c2))/        &
                  ( (c3*(ddry_cm**c4)) - log10(relhum) ) )**onethird
        dwetbb = dwet_cm*1.0e4

        do isub_bin = 1, nsub_bin

            ! ddry, dwet at sub_bin lower boundary are those
            ! at upper boundary of previous sub_bin
            ddryaa = ddrybb
            !dwetaa = dwetbb

            ! compute ddry, dwet (micrometers) at sub_bin upper boundary
            dum = alogddrylo + isub_bin*dlogddry
            ddrybb = 10.0**( dum )

            ddry_cm = ddrybb*1.0e-4
            !dwet_cm = ( ddry_cm**3 + (c1*(ddry_cm**c2))/        &
              !        ( (c3*(ddry_cm**c4)) - log10(relhum) ) )**onethird
            !dwetbb = dwet_cm*1.0e4

            ! geometric mean rdry, rwet (micrometers) for sub_bin
            ddry = sqrt(ddryaa * ddrybb)
            !dwet = sqrt(dwetaa * dwetbb)

            ! xmdry is dry mass in g
            xmdry = drydens_16pi_em12 * (ddry**3.0)

            ! xsdry is dry surface in cm2
            xsdry = x_pi_em8 * (ddry**2.0)

            ! Equation 2 from Fuentes et al (ACP, 2011). Wet diameter is scaled by a factor
            ! of 1e3 to convert from um to nm for this calculation
            df0dlogddry = 0.0
            do jd = 1, 4
                df0dlogddry = df0dlogddry + nti(jd)/( (2.0*pi)**0.5 * log10(width_ss_f(jd)) ) &
                            * exp(-1.0/2.0 * (log10(ddry*1e3/dp0gi(jd))/log10(width_ss_f(jd)))**2.0 )
            end do

            dumsum_na = dumsum_na + dlogddry*df0dlogddry*scale_factor
            dumsum_ma = dumsum_ma + dlogddry*df0dlogddry*scale_factor*xmdry
            dumsum_sa = dumsum_sa + dlogddry*df0dlogddry*scale_factor*xsdry

        end do

    
    ! 2................
        ! rdrylo, rdryhi = lower and upper dry radii (micrometers) 
        ! for this part of the integration
        rdrylo = rdry_star
        rdryhi = rdryuppermost

        nsub_bin = nsub_upper_bin

        alnrdrylo = log( rdrylo )
        dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

        ! compute rdry, rwet (micrometers) at lowest size
        rdrybb = exp( alnrdrylo )
        rdry_cm = rdrybb*1.0e-4
        rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/        &
                  ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
        rwetbb = rwet_cm*1.0e4

        do isub_bin = 1, nsub_bin

            ! rdry, rwet at sub_bin lower boundary are those
            ! at upper boundary of previous sub_bin
            rdryaa = rdrybb
            rwetaa = rwetbb

            ! compute rdry, rwet (micrometers) at sub_bin upper boundary
            dum = alnrdrylo + isub_bin*dlnrdry
            rdrybb = exp( dum )

            rdry_cm = rdrybb*1.0e-4
            rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/        &
                      ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
            rwetbb = rwet_cm*1.0e4

            ! geometric mean rdry, rwet (micrometers) for sub_bin
            rdry = sqrt(rdryaa * rdrybb)
            rwet = sqrt(rwetaa * rwetbb)
            drwet = rwetbb - rwetaa

            ! xmdry is dry mass in g
            xmdry = drydens_43pi_em12 * (rdry**3.0)

            ! xsdry is dry surface in cm2
            xsdry = x_4pi_em8 * (rdry**2.0)

            ! dumb is "B" in Gong's Eqn 5a
            ! df0drwet is "dF0/dr" in Gong's Eqn 5a
            dumb = ( 0.380 - log10(rwet) ) / 0.650
            dumexpb = exp( -dumb*dumb)
            df0drwet = 1.373 * (rwet**(-3.0)) *             &
                      (1.0 + 0.057*(rwet**1.05)) *             &
                      (10.0**(1.19*dumexpb))

            dumsum_na = dumsum_na + drwet*df0drwet
            dumsum_ma = dumsum_ma + drwet*df0drwet*xmdry
            dumsum_sa = dumsum_sa + drwet*df0drwet*xsdry

        end do



    ! "section 3"
    ! where rdry < rdry_star for the whole bin we use Fuentes et al 
    ! for all emissions
    !
    else
    
        ! ddrylo, ddryhi = lower and upper dry diameter (micrometers) 
        ! for this part of the integration
        ddrylo = rdrylowermost*2.0
        ddryhi = rdryuppermost*2.0

        nsub_bin = 1000

        alogddrylo = log10( ddrylo )
        dlogddry = (log10( ddryhi ) - alogddrylo)/nsub_bin

        ! compute ddry, dwet (micrometers) at lowest size
        ddrybb = 10.0**( alogddrylo )
        ddry_cm = ddrybb*1.0e-4
        dwet_cm = ( ddry_cm**3 + (c1*(ddry_cm**c2))/        &
                  ( (c3*(ddry_cm**c4)) - log10(relhum) ) )**onethird
        dwetbb = dwet_cm*1.0e4

        do isub_bin = 1, nsub_bin

            ! ddry, dwet at sub_bin lower boundary are those
            ! at upper boundary of previous sub_bin
            ddryaa = ddrybb
            !dwetaa = dwetbb

            ! compute ddry, dwet (micrometers) at sub_bin upper boundary
            dum = alogddrylo + isub_bin*dlogddry
            ddrybb = 10.0**( dum )

            ddry_cm = ddrybb*1.0e-4
            !dwet_cm = ( ddry_cm**3 + (c1*(ddry_cm**c2))/        &
              !        ( (c3*(ddry_cm**c4)) - log10(relhum) ) )**onethird
            !dwetbb = dwet_cm*1.0e4

            ! geometric mean rdry, rwet (micrometers) for sub_bin
            ddry = sqrt(ddryaa * ddrybb)
            !dwet = sqrt(dwetaa * dwetbb)

            ! xmdry is dry mass in g
            xmdry = drydens_16pi_em12 * (ddry**3.0)

            ! xsdry is dry surface in cm2
            xsdry = x_pi_em8 * (ddry**2.0)

            ! Equation 2 from Fuentes et al (ACP, 2011). Wet diameter is scaled by a factor
            ! of 1e3 to convert from um to nm for this calculation
            df0dlogddry = 0.0
            do jd = 1, 4
                df0dlogddry = df0dlogddry + nti(jd)/( (2.0*pi)**0.5 * log10(width_ss_f(jd)) ) &
                            * exp(-1.0/2.0 * (log10(ddry*1e3/dp0gi(jd))/log10(width_ss_f(jd)))**2.0 )
            end do

            dumsum_na = dumsum_na + dlogddry*df0dlogddry*scale_factor
            dumsum_ma = dumsum_ma + dlogddry*df0dlogddry*scale_factor*xmdry
            dumsum_sa = dumsum_sa + dlogddry*df0dlogddry*scale_factor*xsdry

        end do
    


    end if
    !c
    !c  all done
    !c    
    emitfact_numb = dumsum_na
    emitfact_mass = dumsum_ma
    emitfact_surf = dumsum_sa


    return
end subroutine seasalt_emit_feuntes_1bin

subroutine mosaic_seasalt_emiss2(                            &
           dtstep, u10, v10,u_phy,v_phy,alt,rho_phy,dz8w,    &
           xland,dx,g,tend_emis,sstemp,e_seas,               &
           emis_seas_factor,                                 &
           ids,ide, jds,jde, kds,kde,                        &
           ims,ime, jms,jme, kms,kme,                        &
           its,ite, jts,jte, kts,kte    )

! Routine to calculate seasalt emissions for MOSAIC over the time dtstep... 
! Chun.Zhao@pnnl.gov; 10-Aug-2011
! New formula including correction of particle radius < 0.2 um and SST correction
! Reference:
! Gong 2003 and Jeagle et al. 2011
!------------------------------------------------------------------------

    use mpas_atmchem_vars, only:  num_chem, param_first_scalar


    implicit none

    integer, intent(in   ) :: ids,ide, jds,jde, kds,kde,            &
                              ims,ime, jms,jme, kms,kme,            &
                              its,ite, jts,jte, kts,kte
    real, intent(in   ) :: dtstep, g, emis_seas_factor
    real, dimension(ims:ime,        jms:jme), intent(inout) :: e_seas
    real, dimension(ims:ime,        jms:jme), intent(in   ) :: sstemp
    real, dimension(ims:ime,        jms:jme), intent(in   ) :: dx
    real, dimension(ims:ime,        jms:jme), intent(in   ) :: u10       , &
                                                               v10       , &
                                                               xland
    real, dimension(ims:ime,kms:kme,jms:jme), intent(in   ) :: alt       , &
                                                               dz8w      , &
                                                               u_phy     , &
                                                               v_phy     , &
                                                               rho_phy
    real, dimension(ims:ime,kms:kme,jms:jme,num_chem ), intent(inout) :: tend_emis

! local variables
    integer :: i,j,k,n,l
    real dumdlo,dumdhi
    !jwfeng: There is no need to use 1000 sub-bins, use 5 instead
    integer, parameter :: nsub_bin = 5 ! sub-bin for integration
    real bems_mass(nsub_bin),bems_num(nsub_bin)  ! sea-salt emission for bins 
    real total_mass,total_num
    real w10m,gwet,airden,airmas
    real*8, DIMENSION( ims:ime , jms:jme ) :: dxy
    real fracseas, factaa
    real conver,converi
    real ssemfact_numb( maxd_asize, maxd_atype )
    real ssemfact_mass( maxd_asize, maxd_atype )
    real sst_c
    integer l_seas,iphase,itype,p1st
    character(len=strkind) :: debugString

    p1st    = PARAM_FIRST_SCALAR
    itype   = 1
    iphase  = ai_phase
    conver  = 1.d-9
    converi = 1.d9

    k=kts
    do j=jts,jte
    do i=its,ite
        e_seas(i,j)=0.0
        !sea-salt only over water
        if ( xland(i,j).gt.1.5 ) then
        !-----------------------------------------------------
        !1. Pre-setting inputs 
        !-----------------------------------------------------
            !1.1 10-m wind speed
            w10m = sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
            ! we don't trust the u10,v10 values, for thin model surface layer
            if (dz8w(i,kts,j).le.10.) then
                w10m = sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
            end if
            !1.2 Grid area 
            dxy = dx(i,j)*dx(i,j)
            !1.3 Sea-surface temperature in degree C
            sst_c=sstemp(i,j)-273.15
        !-----------------------------------------------------
        !2. Emissions 
        !-----------------------------------------------------
            !2.1 Emission for bins 
            do n = 1, nsize_aer(itype)
                dumdlo = max(dlo_sect(n,itype), 0.02e-4 )
                dumdhi = max(dhi_sect(n,itype), 0.02e-4 )
                dumdlo = dumdlo * 1.e4 ! cm to um
                dumdhi = dumdhi * 1.e4 ! cm to um 
                call seasalt_source(dtstep,w10m,sst_c,dumdlo,dumdhi,nsub_bin,bems_mass,bems_num)
                !bems_mass: particles kg m-2 s-1
                !bems_num: particles # m-2 s-1
                ssemfact_mass(n,itype) = sum(bems_mass(1:nsub_bin))*converi !kg to ug
                ssemfact_numb(n,itype) = sum(bems_num(1:nsub_bin))
                !jwfeng: add seasalt emis factor
                ssemfact_mass(n,itype) = ssemfact_mass(n,itype)*emis_seas_factor
                ssemfact_numb(n,itype) = ssemfact_numb(n,itype)*emis_seas_factor
            enddo

            !2.2 Sum the total emission
            total_mass = 0.0
            total_num  = 0.0
            do n = 1, nsize_aer(itype)
                total_mass = total_mass + ssemfact_mass(n,itype)
                total_num  = total_num  + ssemfact_numb(n,itype)
            enddo
            !ug m-2 s-1
            e_seas(i,j) = total_mass

            !2.3 Apportion seasalt mass emissions assumming that seasalt is pure NaCl
            fracseas = 1.0

         !-----------------------------------------------------
         !3. Update chemistry tendency (Accumulation and Coarse modes) 
         !-----------------------------------------------------
            !3.1 Factor converting emission to concentration
            !factaa is (s*m2/kg-air)
            !factaa*ssemfact_mass is (s*m2/kg-air)*(ug/m2/s) = ug/kg-air
            !factaa*ssemfact_numb is (s*m2/kg-air)*( #/m2/s) =  #/kg-air

! jwfeng: Using tend_emis instead of directly changing chem,
!         Dont multi dt in factaa !!!Important
!         factaa = (dtstep/dz8w(i,k,j))*alt(i,k,j) 
            factaa = alt(i,k,j)/dz8w(i,k,j)

            !3.2 Update 
            do n = 1, nsize_aer(itype)
                l_seas = lptr_seas_aer(n,itype,iphase)
                if ( l_seas >= p1st ) then
                    tend_emis(i,k,j,l_seas) = factaa * ssemfact_mass(n,itype) * fracseas

                    l = numptr_aer(n,itype,iphase)
                    if (l >= p1st) tend_emis(i,k,j,l) = tend_emis(i,k,j,l) +   &
                                                        factaa * ssemfact_numb(n,itype)
                endif
            enddo

        endif !xland 

    enddo
    enddo
!    do n = 1, nsize_aer(itype)
!        l_seas = lptr_seas_aer(n,itype,iphase)
!        write(debugString,*) "l_seas = ", l_seas, "max/min of tend_emis = ",maxval(tend_emis(its:ite,kts:kte,jts:jte,l_seas)), minval(tend_emis(its:ite,kts:kte,jts:jte,l_seas))
!        call mpas_log_write(debugString)
!    enddo

end subroutine mosaic_seasalt_emiss2

subroutine seasalt_source(dt,w10m,sst_c,diamlo,diamhi,nsub_bin,bems_mass,bems_num)
! ****************************************************************************
! * Number flux density: Original formula by Monahan et al. (1986) adapted
! * by Sunling Gong (JGR 1997 (old) and GBC 2003 (new)).  The new version is
! * to better represent emission of sub-micron sea salt particles.
!
! * dFn/dr = c1*u10**c2/(r**A) * (1+c3*r**c4)*10**(c5*exp(-B**2))
! * where B = (b1 -log(r))/b2
! * see c_old, c_new, b_old, b_new below for the constants.
! * number fluxes are at 80% RH.
! *
! * To calculate the flux:
! * 1) Calculate dFn based on Monahan et al. (1986) and Gong (2003)
! * 2) The wet radius r is estimated at 80% RH
! * 3) Convert particles flux to mass flux :
! *    dFM/dr_d = 4/3*pi*rho_d*r_d^3 *(dr_w/dr_d) * dFn/dr_w
! *    or 
! *    dFM = 4/3*pi*rho_d*r_d^3 * dFn
! *               where rho_p is particle density [kg/m3]
! *    The factor 1.e-18 is to convert in micro-meter r_d^3
! ****************************************************************************
  IMPLICIT NONE

  !input & output
  REAL,    INTENT(IN)    :: dt
  REAL,    INTENT(IN)    :: w10m
  REAL,    INTENT(IN)    :: sst_c
  REAL,    INTENT(IN)    :: diamlo,diamhi
  INTEGER, INTENT(IN)    :: nsub_bin
  REAL,    INTENT(INOUT) :: bems_mass(nsub_bin),bems_num(nsub_bin)

  !local variables
  !values used in Gong's paper (2003)
  REAL, PARAMETER :: c_old(5)=(/1.373, 3.41, 0.057, 1.05, 1.190/)
  REAL, PARAMETER :: c_new(5)=(/1.373, 3.41, 0.057, 3.45, 1.607/)
  !Change suggested by MC
  !REAL, PARAMETER :: c_old(5)=(/1.373, 3.2, 0.057, 1.05, 1.190/)
  !REAL, PARAMETER :: c_new(5)=(/1.373, 3.2, 0.057, 3.45, 1.607/)
  REAL, PARAMETER :: b_old(2)=(/0.380, 0.650/)
  REAL, PARAMETER :: b_new(2)=(/0.433, 0.433/)
  REAL :: c0(5), b0(2)
 !REAL*8, PARAMETER :: dr=5.0D-2 ! um, radius interval for integration 
  REAL, PARAMETER :: theta=30.0
  REAL, PARAMETER :: den_seas=2.2e3 !sea-salt density, kg/m3
  REAL :: r_to_mass,xmdry
 !Swelling coefficient frh (d rwet / d rd)
 !REAL*8,    PARAMETER :: frh = 1.65
 !REAL*8,    PARAMETER :: frh = 2.d0
  LOGICAL, PARAMETER :: new=.TRUE.  ! new from Gong (2003) or old (1997)
  LOGICAL, PARAMETER :: sstadjust=.true. ! SST adjustment 
  REAL SSTFACTOR
  REAL :: rho_d, rdrylo, rdryhi
  real :: rd1, rd2, rw1,rw2, rd1_cm,rd2_cm,rw1_cm,rw2_cm
  real :: rdry,rwet,drdry,drwet
  real,parameter :: c1 = 0.7674
  real,parameter :: c2 = 3.079
  real,parameter :: c3 = 2.573e-11
  real,parameter :: c4 = -1.424
  real,parameter :: onethird = 1.0/3.0
  real,parameter :: pi = 3.1415936536
  real,parameter :: relhum = 0.8
  real :: a, b, dfndr,dfn,dfm
  integer isub_bin

 !-----------------------------------------------------
 !1. Initializing
 !-----------------------------------------------------
     bems_mass = 0.0
     bems_num = 0.0
     rho_d = den_seas
     !factor for radius (micrometers) to mass (kg)
     r_to_mass = rho_d*(4.0/3.0)*pi*1.0e-18
     rdrylo = diamlo/2.0
     rdryhi = diamhi/2.0
     drdry = (rdryhi - rdrylo)/nsub_bin

 !-----------------------------------------------------
 !2. Emission source 
 !-----------------------------------------------------
     !2.1 compute rdry,rwet (micrometers) at lowest bound 
     rd2 = rdrylo
     rd2_cm = rdrylo*1.e-4
     rw2_cm = ( rd2_cm**3 +  &
               (c1*(rd2_cm**c2))/((c3*(rd2_cm**c4))-log10(relhum)) &
              )**onethird
     rw2 = rw2_cm*1.0e4

     do 2000 isub_bin = 1, nsub_bin

       !2.2 rdry&rwet at the sub_bin lower and upper bound 
       rd1 = rd2
       rw1 = rw2
       rd2 = rdrylo + isub_bin*drdry
       rd2_cm = rd2*1.e-4
       rw2_cm = ( rd2_cm**3 +  &
                 (c1*(rd2_cm**c2))/((c3*(rd2_cm**c4))-log10(relhum)) &
                )**onethird
       rw2 = rw2_cm*1.0e4

       !2.3 geometric mean rdry&rwet (micrometers) for sub_bin
       rdry = 0.5*(rd1+rd2)
       rwet = 0.5*(rw1+rw2)
       drwet= rw2 - rw1

       !2.4 convert factor (xmdry),from wet number to dry mass in kg
       !xmdry = r_to_mass * (drwet/drdry) * (rdry**3.0)
       xmdry = r_to_mass * (rdry**3.0)

       !2.5 Coefficients for source formula 
       IF (new) THEN
          a = -1.0*(4.7*(1.0 + theta*rwet)**(-0.017*rwet**(-1.44)))
          c0 = c_new
          b0 = b_new
       ELSE
          a = -3.0
          c0 = c_old
          b0 = b_old
       END IF
       b = (b0(1) - LOG10(rwet))/b0(2)

       !2.6 Source function 
       !particles # m-2 s-1 um-1
       dfndr = (c0(1)*(rwet**a))*(w10m**c0(2))*(1.0+c0(3)*rwet**c0(4))* &
              10**(c0(5)*EXP(-(b**2)))
       if (sstadjust) then ! SST adjustment 
        if(sst_c.ge.0.and.sst_c.le.30) then
        sstfactor = (0.3+0.1*sst_c-0.0076*sst_c**2+0.00021*sst_c**3)
        dfndr = dfndr*sstfactor
        endif
       endif

       !particles # m-2 s-1
       dfn = dfndr*drwet
       !particles kg m-2 s-1
       dfm = xmdry*dfn
       bems_mass(isub_bin) = dfm
       bems_num(isub_bin) = dfn

 2000  CONTINUE

       return

END SUBROUTINE seasalt_source


end module module_mosaic_addemiss

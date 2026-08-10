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
!   Copied from WRF-Chem-v3.5:chem/module_aer_drydep.F
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

MODULE module_aer_drydep

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

CONTAINS


!-----------------------------------------------------------------------
    subroutine aer_drydep_driver(debug_on,                             &
            id, ktau, dtstep, config_flags, aer_mech_id,               &
            gmt, julday,                                               &
            t_phy, rho_phy, p_phy,                                     &
            alt, p8w, t8w, dz8w, z, z_at_w,                            &
            ust, aer_res, ivgtyp, vegfra, pbl, rmol, znt,              &
            moist, chem, ddvel,                                        &
            h2oai, h2oaj, numgas,                                      &
            ids,ide, jds,jde, kds,kde,                                 &
            ims,ime, jms,jme, kms,kme,                                 &
            its,ite, jts,jte, kts,kte                                  )
!
! This routines allow different aerosol dry-deposition velocity parameterizations
! to be used with various MOSAIC and SORGAM chem_opt packages.  
!
! This routine is called when config_flags%aer_drydep_opt >= 100.
! However, values other than 101,111,201,301,302,311,312 cause an error.
!
! For MOSAIC packages, calls are made to mosaic_aer_drydep_prep (to compute
! wet/ambient size and density), then to one of the dry-deposition velocity
! routines, then to mosaic_aer_drydep_load_ddvel.
! For SORGAM packages, calls are made to sorgam_aer_drydep_prep (to compute
! wet/ambient size and density), then to one of the dry-deposition velocity
! routines, then to sorgam_aer_drydep_load_ddvel.
!
! config_flags%aer_drydep_opt == 101 -- call sorgam_aer_drydepvel_2, which
!    uses the parameterization of subr vdvg_2 of module_aerosols_sorgam.F
!
! config_flags%aer_drydep_opt == 111 -- call sorgam_aer_drydepvel_1, which
!    uses the parameterization of subr vdvg of module_aerosols_sorgam.F
!
! config_flags%aer_drydep_opt == 201 -- call mosaic_aer_drydepvel_1, which
!    uses the parameterization of subr aerosol_depvel_2 of module_mosaic_drydep.F
!
! config_flags%aer_drydep_opt == 301,302,311,312 -- call zhang2001_aer_drydepvel_1,
!    which uses the parameterization of zhang et al. (2001, atmos. environ., p. 549-560).
!    301 & 302 - uses the parameterization as described in the 2001 paper
!    311 & 312 - has some modifications that increase the interception term (thus 
!        increasing depvel for particles in the accumulation size range)
!        and reduce the brownian diffusion and impaction terms (thus reducing 
!        depvel for particles in the Aitken and coarse size ranges) 
!    301 & 311 - iseason is based on julian date, as in module_dep_simple.F
!    302 & 312 - iseason is always 1 (summer)
!

!   use module_configure, only:  MPAS_chem_config, num_moist, num_chem
!   use module_state_description, only:  param_first_scalar

!   use modal_aero_data, only:  ntot_amode_cam_mam => ntot_amode
    implicit none

    logical, intent(in) :: debug_on
!   subr arguments
    integer, intent(in) ::   &
        id, ktau, julday, aer_mech_id,   &
        numgas

    integer, intent(in) ::   &
        ids, ide, jds, jde, kds, kde,   &
        ims, ime, jms, jme, kms, kme,   &
        its, ite, jts, jte, kts, kte

    real, intent(in) :: dtstep, gmt

    real, intent(in),   &
        dimension( ims:ime, kms:kme, jms:jme ) :: &
        t_phy, rho_phy, p_phy, &
        alt, p8w, t8w, dz8w, z, z_at_w

    integer, intent(in),   &
        dimension( ims:ime, jms:jme ) :: &
        ivgtyp

    real, intent(in),   &
        dimension( ims:ime, jms:jme ) :: &
        ust, vegfra, pbl, rmol, znt

    real, intent(in),   &
        dimension( its:ite, jts:jte ) :: &
        aer_res

    real, intent(in),   &
        dimension( ims:ime, kms:kme, jms:jme, 1:num_moist ) :: &
        moist
 
    real, intent(inout),   &
        dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
        chem

    real, intent(inout),   &
        dimension( its:ite, jts:jte, 1:num_chem ) :: &
        ddvel

    real, intent(inout),   &
        dimension( ims:ime, kms:kme, jms:jme ) :: &
        h2oai, h2oaj

    type(MPAS_chem_config), intent(in) :: config_flags


!   local variables
    integer :: i, ibin, iflag_coarse, iflag_nac, &
               iok_aer_drydep_opt, ioptaa, iphase, iseason, isurftype, itype, iwetsurf
    integer :: itmpa
    integer :: j
    integer :: n, nbin_add, numaer

    real :: airdensbox, airkinvisc
    real :: depresist_a, depresist_d0, depresist_d3
    real :: depresist_unstabpblfact
    real :: freepath
    real :: presbox
    real :: tempbox
    real :: tmpa, tmp_lnsg, tmp_num
    real :: ustar
    real :: vsettl_0, vsettl_3
    real :: wetdgnum, wetdens, wstar
    real :: wetdgnum_si, wetdens_si

    real, allocatable :: alnsg_add(:)
    real, allocatable :: depvel0_add(:), depvel3_add(:)
    real, allocatable :: wetdgnum_add(:), wetdens_add(:)
    real, allocatable :: yextra_add(:,:)

    character(len=160) :: txtaa


    if(debug_on) call mpas_log_write( &
        '     *** aer_drydep_driver - aer_mech_id, aer_drydep_opt = $i,$i', &
        intArgs = (/aer_mech_id, config_flags%aer_drydep_opt/))

!   allocate working arrays
    if ((aer_mech_id == 1) .or.   &
        (aer_mech_id == 2)) then
        nbin_add = 3

    else if (aer_mech_id == 3) then
        nbin_add = sum( nsize_aer(1:ntype_aer) )
        else if (aer_mech_id == 4) then
            nbin_add = ntot_amode_cam_mam
    else
        if(debug_on) write(txtaa,'(a,2(1x,i10))') &
        '     *** aer_drydep_driver - bad aer_mech_id, aer_drydep_opt = ', &
        aer_mech_id, config_flags%aer_drydep_opt
        if(debug_on) call mpas_log_write( trim(txtaa) )
        call chemistry_error_fatal( txtaa )
    end if

    allocate( alnsg_add(    nbin_add ) )
    allocate( depvel0_add(  nbin_add ) )
    allocate( depvel3_add(  nbin_add ) )
    allocate( wetdgnum_add( nbin_add ) )
    allocate( wetdens_add(  nbin_add ) )
    allocate( yextra_add(   nbin_add, 10 ) )

main_j_loop:   &
    do j = jts, jte
main_i_loop:   &
    do i = its, ite


!   get aerosol properties at lowest model layer
    if ((aer_mech_id == 1) .or.   &
        (aer_mech_id == 2)) then
        call sorgam_aer_drydep_prep(                                   &
        id, ktau, dtstep, config_flags, aer_mech_id,               &
        gmt, julday,                                               &
        t_phy, rho_phy, p_phy, alt, p8w, t8w,                      &
        moist, chem,                                               &
        h2oai, h2oaj,                                              &
        nbin_add, alnsg_add, wetdgnum_add, wetdens_add,            &
        yextra_add,                                                &
        i, j,                                                      &
        ids,ide, jds,jde, kds,kde,                                 &
        ims,ime, jms,jme, kms,kme,                                 &
        its,ite, jts,jte, kts,kte                                  )

    else if (aer_mech_id == 3) then
        call mosaic_aer_drydep_prep(debug_on,                      &
        id, ktau, dtstep, config_flags, aer_mech_id,               &
        gmt, julday,                                               &
        t_phy, rho_phy, p_phy, alt, p8w, t8w,                      &
        moist, chem,                                               &
        nbin_add, alnsg_add, wetdgnum_add, wetdens_add,            &
        yextra_add,                                                &
        i, j,                                                      &
        ids,ide, jds,jde, kds,kde,                                 &
        ims,ime, jms,jme, kms,kme,                                 &
        its,ite, jts,jte, kts,kte                                  )
    else if (aer_mech_id == 4) then
        call cam_mam_aer_drydep_prep(                              &
        id, ktau, dtstep, config_flags, aer_mech_id,               &
        gmt, julday,                                               &
        t_phy, rho_phy, p_phy, alt, p8w, t8w,                      &
        moist, chem,                                               &
        nbin_add, alnsg_add, wetdgnum_add, wetdens_add,            &
        yextra_add,                                                &
        i, j,                                                      &
        ids,ide, jds,jde, kds,kde,                                 &
        ims,ime, jms,jme, kms,kme,                                 &
        its,ite, jts,jte, kts,kte                                  )

    else
        txtaa = '*** aer_drydep_driver - fatal error 200'
        call chemistry_error_fatal( txtaa )
    end if


!   calculate deposition velocities
    do ibin = 1, nbin_add

    ustar = ust(i,j)
    depresist_a = aer_res(i,j)
    tempbox = t_phy(i,kts,j)
    presbox = p_phy(i,kts,j)

    wetdgnum = wetdgnum_add(ibin)
    wetdens  = wetdens_add(ibin)
    tmp_lnsg = alnsg_add(ibin)

    iok_aer_drydep_opt = 0

    if ((config_flags%aer_drydep_opt >= 100) .and.   &
        (config_flags%aer_drydep_opt <= 199)) then
        airkinvisc = 0.0
        freepath = 0.0
        ustar = max( 1.e-1, ust(i,j) )   ! this done in sorgam_depdriver
        tmpa = .001*p_phy(i,kts,j)
        presbox = 1.e3*tmpa

        airdensbox = 1.0/alt(i,kts,j)

        iflag_coarse = 0
        if (aer_mech_id <= 2) then
        if (ibin == 3) iflag_coarse = 1
        else
        ! yextra_add(ibin,9) = "reference" diameter (m)
        if (yextra_add(ibin,9) >= 1.0e-6) iflag_coarse = 1
        end if

        if (config_flags%aer_drydep_opt == 101) then
        iok_aer_drydep_opt = 1

        call sorgam_aer_drydepvel_2(   &
            iflag_coarse,   &
            wetdgnum, tmp_lnsg, wetdens,   &
            tempbox, presbox, airdensbox,   &
            ustar, depresist_a,   &
            pbl(i,j), znt(i,j), rmol(i,j),   &
            airkinvisc, freepath,   &
            depvel0_add(ibin), depvel3_add(ibin),   &
            depresist_d0, depresist_d3,   &
            vsettl_0, vsettl_3 )

        else if (config_flags%aer_drydep_opt == 111) then
        iok_aer_drydep_opt = 1
        wstar = 0.0

        call sorgam_aer_drydepvel_1(   &
            iflag_coarse,   &
            wetdgnum, tmp_lnsg, wetdens,   &
            tempbox, presbox, airdensbox,   &
            ustar, wstar, depresist_a,   &
            airkinvisc, freepath,   &
            depvel0_add(ibin), depvel3_add(ibin),   &
            depresist_d0, depresist_d3,   &
            vsettl_0, vsettl_3 )

        end if
    
        depresist_unstabpblfact = 0.0
!       depresist_d0 = 0.0
!       depresist_d3 = 0.0

    else if ((config_flags%aer_drydep_opt >= 200) .and.   &
             (config_flags%aer_drydep_opt <= 299)) then
        airkinvisc = 0.0
        freepath = 0.0
        depresist_unstabpblfact = 0.0

        if ((config_flags%aer_drydep_opt == 201) .or.   &
            (config_flags%aer_drydep_opt == 211)) then
        iok_aer_drydep_opt = 1
        airdensbox = rho_phy(i,kts,j)
        ioptaa = 1
        if (config_flags%aer_drydep_opt == 211) ioptaa = 2

        call mosaic_aer_drydepvel_1(   &
            ioptaa, wetdgnum, tmp_lnsg, wetdens,   &
            tempbox, airdensbox, ustar, depresist_a,   &
            airkinvisc, freepath, depresist_unstabpblfact,   &
            depvel0_add(ibin), depvel3_add(ibin),   &
            depresist_d0, depresist_d3,   &
            vsettl_0, vsettl_3 )

        end if
    
    else if ((config_flags%aer_drydep_opt >= 300) .and.   &
             (config_flags%aer_drydep_opt <= 399)) then
        airkinvisc = 0.0
        freepath = 0.0
        depresist_unstabpblfact = 0.0

        if ((config_flags%aer_drydep_opt == 301) .or. &
            (config_flags%aer_drydep_opt == 302) .or. &
            (config_flags%aer_drydep_opt == 311) .or. &
            (config_flags%aer_drydep_opt == 312)) then
        iok_aer_drydep_opt = 1
        airdensbox = rho_phy(i,kts,j)
        isurftype = ivgtyp(i,j)
        ! should set iwetsurf=1 when rainfall>0 or relhum>~95%
        iwetsurf = 0
        ! this is how iseason is set in module_dep_simple
        ! could do a better job using leaf area index
        iseason = 1
        if (julday < 90 .or. julday > 270) iseason = 2
        if ((config_flags%aer_drydep_opt == 302) .or. &
            (config_flags%aer_drydep_opt == 312)) iseason = 1
        ioptaa = 1
        if ((config_flags%aer_drydep_opt == 311) .or. &
            (config_flags%aer_drydep_opt == 312)) ioptaa = 2

        call zhang2001_aer_drydepvel_1(   &
            wetdgnum, tmp_lnsg, wetdens,   &
            tempbox, airdensbox, ustar, depresist_a,   &
            ioptaa, iseason, isurftype, 2, iwetsurf,   &
            airkinvisc, freepath,   &
            1, 1,   &
            depvel0_add(ibin), depvel3_add(ibin),   &
            depresist_d0, depresist_d3,   &
            vsettl_0, vsettl_3 )

        end if
    
    end if

    if (iok_aer_drydep_opt <= 0) then
        if(debug_on) write(txtaa,'(a,2(1x,i10))') &
        '     *** aer_drydep_driver - bad aer_mech_id, aer_drydep_opt = ', &
        aer_mech_id, config_flags%aer_drydep_opt
        if(debug_on) write(0,'(a)') trim(txtaa)
        call chemistry_error_fatal( txtaa )
    end if


    end do   ! ibin


!   load deposition velocities into ddvel array
    if ((aer_mech_id == 1) .or.   &
        (aer_mech_id == 2)) then
        call sorgam_aer_drydep_load_ddvel(                             &
        id, ktau, dtstep, config_flags, aer_mech_id,               &
        numgas, ddvel,                                             &
        nbin_add, depvel0_add, depvel3_add,                        &
        i, j,                                                      &
        ids,ide, jds,jde, kds,kde,                                 &
        ims,ime, jms,jme, kms,kme,                                 &
        its,ite, jts,jte, kts,kte                                  )

    else if (aer_mech_id == 3) then
        call mosaic_aer_drydep_load_ddvel(debug_on,                &
        id, ktau, dtstep, config_flags, aer_mech_id,               &
        ddvel,                                                     &
        nbin_add, depvel0_add, depvel3_add,                        &
        i, j,                                                      &
        ids,ide, jds,jde, kds,kde,                                 &
        ims,ime, jms,jme, kms,kme,                                 &
        its,ite, jts,jte, kts,kte                                  )

        else if (aer_mech_id == 4) then
            call cam_mam_aer_drydep_load_ddvel(                            &
                id, ktau, dtstep, config_flags, aer_mech_id,               &
                ddvel,                                                     &
                nbin_add, depvel0_add, depvel3_add,                        &
                i, j,                                                      &
                ids,ide, jds,jde, kds,kde,                                 &
                ims,ime, jms,jme, kms,kme,                                 &
                its,ite, jts,jte, kts,kte                                  )

    else
        txtaa = '*** aer_drydep_driver - fatal error 400'
        call chemistry_error_fatal( txtaa )
    end if


    end do main_i_loop
    end do main_j_loop


!   deallocate working arrays
    deallocate( alnsg_add    )
    deallocate( depvel0_add  )
    deallocate( depvel3_add  )
    deallocate( wetdgnum_add )
    deallocate( wetdens_add  )
    deallocate( yextra_add   )


    return
    end subroutine aer_drydep_driver



!-----------------------------------------------------------------------
    subroutine mosaic_aer_drydep_prep(debug_on,                        &
            id, ktau, dtstep, config_flags, aer_mech_id,               &
            gmt, julday,                                               &
            t_phy, rho_phy, p_phy, alt, p8w, t8w,                      &
            moist, chem,                                               &
            nbin_add, alnsg_add, wetdgnum_add, wetdens_add,            &
            yextra_add,                                                &
            i, j,                                                      &
            ids,ide, jds,jde, kds,kde,                                 &
            ims,ime, jms,jme, kms,kme,                                 &
            its,ite, jts,jte, kts,kte                                  )

!   use module_configure, only:  num_moist, num_chem, &
!       MPAS_chem_config
!   use module_state_description, only:  param_first_scalar

    implicit none

!   subr arguments
    logical, intent(in) :: debug_on
    integer, intent(in) ::   &
        id, ktau, julday, aer_mech_id, nbin_add,   &
        i, j

    integer, intent(in) ::   &
        ids, ide, jds, jde, kds, kde,   &
        ims, ime, jms, jme, kms, kme,   &
        its, ite, jts, jte, kts, kte

    real, intent(in) :: dtstep, gmt

    real, intent(in),   &
        dimension( ims:ime, kms:kme, jms:jme ) :: &
        t_phy, rho_phy, p_phy, alt, p8w, t8w

    real, intent(in),   &
        dimension( ims:ime, kms:kme, jms:jme, 1:num_moist ) :: &
        moist

    real, intent(in),   &
        dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
        chem

    real, intent(inout),   &
        dimension( 1:nbin_add ) :: &
        alnsg_add, wetdgnum_add, wetdens_add

    real, intent(inout),   &
        dimension( 1:nbin_add, 10 ) :: &
        yextra_add

    type(MPAS_chem_config), intent(in) :: config_flags


!   local variables
    integer :: ibin, ijcount, iok, iphase, itype
    integer :: itmpa
    integer :: k
    integer :: l, ll, l1, n
    integer :: p1st

    real, parameter :: densdefault = 2.0
    real, parameter :: smallmassaa = 1.0e-20
    real, parameter :: smallmassbb = 1.0e-30
    real, parameter :: piover6 = pi/6.0
    real, parameter :: onethird = 1.0/3.0

    real :: drydens, drydp, drymass, dryvol
    real :: rnum
    real :: tmpa, tmp_lnsg
    real :: wetdgnum, wetdens, wetdp, wetmass, wetvol
    real :: wetdgnum_si, wetdens_si

    character(len=160) :: txtaa


    n = sum( nsize_aer(1:ntype_aer) )
    if (n /= nbin_add) then
        if(debug_on) write(txtaa,'(a,2(1x,i10))') &
        '     *** mosaic_aer_drydep_prep - bad nbin_add, sum(nsize)', &
        nbin_add, n
        if(debug_on) write(0,'(a)') trim(txtaa)
        call chemistry_error_fatal( txtaa )
        stop
    end if


!   set some variables to their wrf-chem "standard" values
    p1st = param_first_scalar
    lunerr = -1
    lunout = -1

    ijcount = 0

    k = kts

    ijcount = ijcount + 1

!   *** for now, just calc vdep_aer for iphase = ai_phase
    iphase = ai_phase

!   calculate vdep_aer for each type and size
    ibin = 0
    do 2900 itype = 1, ntype_aer
    do 2800 n = 1, nsize_aer(itype)
    ibin = ibin + 1

!   calculate dry and wet mass, volume, density, and mass-mean diameters
    dryvol = 0.0
    drymass = 0.0
    do ll = 1, ncomp_aer(itype)
        l1 = massptr_aer(ll,n,itype,iphase)
        tmpa = max( chem(i,k,j,l1), 0.0 )
        drymass = drymass + tmpa
        dryvol = dryvol + tmpa/dens_aer(ll,itype)
    end do

    l1 = waterptr_aer(n,itype)
    tmpa = max( chem(i,k,j,l1), 0.0 )
    wetmass = drymass + tmpa
    wetvol = dryvol + tmpa/dens_water_aer

    l1 = numptr_aer(n,itype,iphase)
    rnum = max( chem(i,k,j,l1), 0.0 )

    drymass = drymass*28.966e-9 ! ug/kg-air --> g/mol-air
    wetmass = wetmass*28.966e-9
    dryvol  = dryvol *28.966e-9 ! "u-cm3"/kg-air --> cm3/mol-air
    wetvol  = wetvol *28.966e-9
    rnum    = rnum*28.966e-3    ! #/kg-air  --> #/mol-air

        if (drymass <= smallmassbb) then
            drydp = dcen_sect(n,itype)
            drydens = densdefault
            wetdp = drydp
            wetdens = drydens

        else
            if (drymass <= smallmassaa) then
                wetmass = drymass
                wetvol = dryvol
            end if
            drydens = drymass/dryvol
            wetdens = wetmass/wetvol

            if (rnum >= dryvol/volumlo_sect(n,itype)) then
                drydp = dlo_sect(n,itype)
            else if (rnum <= dryvol/volumhi_sect(n,itype)) then
                drydp = dhi_sect(n,itype)
            else
                drydp = (dryvol/(piover6*rnum))**onethird
            end if

            if (abs(wetvol) > (1000.*abs(dryvol))) then
              tmpa=10.0
            else
              tmpa=abs(wetvol/dryvol)**onethird
              tmpa=max(1.0,min(tmpa,10.0))
            endif
            wetdp = drydp*tmpa

        end if


!   "within section" size distribution currently assumed to be
!   monodisperse, but could assume otherwise
    tmp_lnsg = log( 1.0 )
    wetdgnum = wetdp * exp( -1.5*tmp_lnsg*tmp_lnsg )

    wetdgnum_add(ibin) = wetdgnum * 1.0e-2 ! cm --> m
    wetdens_add(ibin) = wetdens * 1.0e3    ! g/cm3 --> kg/m3
    alnsg_add(ibin) = tmp_lnsg


    yextra_add(ibin,1) = drymass/28.966e-9 ! g/mol-air --> ug/kg-air = "n-kg"/kg-air
    yextra_add(ibin,5) = wetmass/28.966e-9
    yextra_add(ibin,2) = dryvol/28.966e-6  ! cm3/mol-air --> "n-m3"/kg-air
    yextra_add(ibin,6) = wetvol/28.966e-6
    yextra_add(ibin,3) = drydens*1.0e3     ! g/cm3 --> kg/m3
    yextra_add(ibin,7) = wetdens*1.0e3
    yextra_add(ibin,4) = drydp*1.0e-2      ! cm --> m
    yextra_add(ibin,8) = wetdp*1.0e-2

    yextra_add(ibin,9) = dcen_sect(n,itype)*1.0e-2   ! cm --> m
    yextra_add(ibin,10) = rnum/28.966e-3    ! #/mol-air  --> #/kg-air


2800    continue
2900    continue


    return
    end subroutine mosaic_aer_drydep_prep



!-----------------------------------------------------------------------
    subroutine mosaic_aer_drydep_load_ddvel(debug_on,              &
        id, ktau, dtstep, config_flags, aer_mech_id,               &
        ddvel,                                                     &
        nbin_add, depvel0_add, depvel3_add,                        &
        i, j,                                                      &
        ids,ide, jds,jde, kds,kde,                                 &
        ims,ime, jms,jme, kms,kme,                                 &
        its,ite, jts,jte, kts,kte                                  )

!use module_configure, only:  num_moist, num_chem, &
!   MPAS_chem_config
!use module_state_description, only:  param_first_scalar

    implicit none

!   subr arguments
    logical, intent(in) :: debug_on
    integer, intent(in) ::   &
        id, ktau, aer_mech_id, nbin_add,   &
        i, j

    integer, intent(in) ::   &
        ids, ide, jds, jde, kds, kde,   &
        ims, ime, jms, jme, kms, kme,   &
        its, ite, jts, jte, kts, kte

    real, intent(in) :: dtstep

    real, intent(in),   &
        dimension( 1:nbin_add ) :: &
        depvel0_add, depvel3_add

    real, intent(inout),   &
        dimension( its:ite, jts:jte, 1:num_chem ) :: &
        ddvel

    type(MPAS_chem_config), intent(in) :: config_flags


!   local variables
    integer :: ibin, iphase, itype
    integer :: l, ll, n

    character(len=160) :: txtaa


    n = sum( nsize_aer(1:ntype_aer) )
    if (n /= nbin_add) then
        if(debug_on) write(txtaa,'(a,2(1x,i10))') &
        '     *** mosaic_aer_drydep_load_ddvel - bad nbin_add, sum(nsize)', &
        nbin_add, n
        if(debug_on) write(0,'(a)') trim(txtaa)
        call chemistry_error_fatal( txtaa )
        stop
    end if


! copy from vdep_aer array to ddvel array
    do iphase = 1, nphase_aer
    ibin = 0
    do itype = 1, ntype_aer
    do n = 1, nsize_aer(itype)
        ibin = ibin + 1
        do ll = -2, ncomp_plustracer_aer(itype)
        if (ll == -2) then
            l = numptr_aer(n,itype,iphase)
        else if (ll == -1) then
            l = -1
            if (iphase .eq. ai_phase) l = waterptr_aer(n,itype)
        else if (ll == 0) then
            l = -1
            if (iphase == ai_phase) l = hyswptr_aer(n,itype)
        else
            l = massptr_aer(ll,n,itype,iphase)
        end if
        if (iphase > 1) l = -1

        if ((l >= param_first_scalar) .and. (l <= num_chem)) then
            ddvel(i,j,l) = depvel3_add(ibin)
        end if
        end do
    end do
    end do
    end do


    return
    end subroutine mosaic_aer_drydep_load_ddvel



!------------------------------------------------------------------------
    subroutine mosaic_aer_drydepvel_1(   &
        ioptaa, dgnum, alnsg, aerodens,   &
        temp, airdens, ustar, depresist_a,   &
        airkinvisc, freepath, depresist_unstabpblfact,   &
        depvel_0, depvel_3,   &
        depresist_d0, depresist_d3,   &
        vsettl_0, vsettl_3 )
!
!   computes the surface layer resistance term and the
!   gravitational settling velocity term for the 3rd moment
!   of a log-normal aerosol mode
!
!   input parameters
!   ioptaa = 2 - apply a sticking fraction for larger particles
!              = other - sticking fraction = 1 always
!   dgnum - geometric mean diameter for aerosol number (m)
!   alnsg - natural logarithm of the geometric standard deviation
!       for aerosol number
!   aerodens - aerosol density (dgnum and aerodens are for the
!       actual wet distribution)
!   temp - temperature (K)
!   airdens - air density (kg/m^3)
!   airkinvisc - air kinematic viscosity (m^2/s)
!   freepath - air molecular freepath (m)
!   ustar - friction velocity (m/s)
!   depresist_a = aerodynamic resistance (s/m)
!   depresist_unstabpblfact = weseley et al. 1985 factor for increasing
!       depvel in unstable pbl -- either
!       1. + (-0.3*zi/L)**0.667 OR
!       1. + 0.24*((wstar/ustar)**2)
!   output parameters
!   depvel_0 - deposition velocity for 0th (number) moment (m/s)
!   depvel_3 - deposition velocity for 3rd (mass  ) moment (m/s)
!   depresist_d0/_d3 - surface layer resistance for 0th/3rd moment (s/m)
!   vsettl_0/_3 - gravitational settling velocity for 0th/3rd moment (m/s)
!   

    implicit none

    integer, intent(in) ::   &
        ioptaa
    real, intent(in) ::   &
        dgnum, alnsg, aerodens,   &
        temp, airdens, ustar, depresist_a
    real, intent(inout) ::   &
        airkinvisc, freepath, depresist_unstabpblfact
    real, intent(inout) ::   &
        depvel_0, depvel_3,   &
        depresist_d0, depresist_d3,   &
        vsettl_0, vsettl_3

    real :: aerodiffus_0, schmidt_0, stokes_0, facdepresist_d0
    real :: aerodiffus_3, schmidt_3, stokes_3, facdepresist_d3
    common / aerosol_depvel_cmn01 /   &
            aerodiffus_0, schmidt_0, stokes_0, facdepresist_d0,   &
            aerodiffus_3, schmidt_3, stokes_3, facdepresist_d3

    real :: aerodiffus_dgnum, alnsg2,   &
        stickfrac,   &
        tmpa, tmpb, tmp_dvm,   &
        vsettl_dgnum, xknudsen, xknudsenfact

!   gravity = gravitational acceleration in m/s^2
    real, parameter :: gravity = 9.80616
!   boltzmann constant in erg/deg-K
    real, parameter :: boltzmann = 1.3807e-23


    if (airkinvisc <= 0.0) then
!   air kinematic viscosity (m^2/s)
        airkinvisc = ( 1.8325e-5 * (416.16/(temp+120.0)) *   &
                ((temp/296.16)**1.5) ) / airdens
    end if
    if (freepath <= 0.0) then
!   air molecular freepath (m)
        freepath = 7.39758e-2 * airkinvisc / sqrt(temp)
    end if
    if (depresist_unstabpblfact <= 0) then
!   enhancement factor for unstable pbl
        depresist_unstabpblfact = 1.0
    end if

    xknudsen = 2.*freepath/dgnum
    xknudsenfact = xknudsen*1.246
    alnsg2 = alnsg*alnsg

    vsettl_dgnum = (gravity*aerodens*dgnum*dgnum)/   &
            (18.*airkinvisc*airdens)
    vsettl_0 = vsettl_dgnum *   &
            ( exp(2.*alnsg2) + xknudsenfact*exp(0.5*alnsg2) )
    vsettl_3 = vsettl_dgnum *   &
            ( exp(8.*alnsg2) + xknudsenfact*exp(3.5*alnsg2) )

    aerodiffus_dgnum = (boltzmann*temp)/   &
            (3.*pi*airkinvisc*airdens*dgnum)
    aerodiffus_0 = aerodiffus_dgnum *   &
            ( exp(+0.5*alnsg2) + xknudsenfact*exp(+2.*alnsg2) )
    aerodiffus_3 = aerodiffus_dgnum *   &
            ( exp(-2.5*alnsg2) + xknudsenfact*exp(-4.*alnsg2) )

    schmidt_0 = airkinvisc/aerodiffus_0
    schmidt_3 = airkinvisc/aerodiffus_3

    stokes_0 = ustar*ustar*vsettl_0/(gravity*airkinvisc)
    stokes_3 = ustar*ustar*vsettl_3/(gravity*airkinvisc)
    
    tmp_dvm = dgnum * exp(1.5*alnsg2) ! volume mean diameter

    tmpa = (schmidt_0**(-0.66666666)) +   &
            (10.**(-3./max(0.03,stokes_0)))
!   tmpb = tmpa*ustar*(1. + 0.24*wstaroverustar*wstaroverustar)
    tmpb = tmpa*ustar*depresist_unstabpblfact
    if (ioptaa == 2) then
        if (tmp_dvm >= 5.0e-6) then
        stickfrac = exp( -2.0*sqrt(max(0.0,stokes_0)) )
        tmpb = tmpb*stickfrac
        end if
    end if
    depresist_d0 = 1./max( tmpb, 1.e-22 )
    facdepresist_d0 = tmpa

    tmpa = (schmidt_3**(-0.66666666)) +   &
            (10.**(-3./max(0.03,stokes_3)))
!   tmpb = tmpa*ustar*(1. + 0.24*wstaroverustar*wstaroverustar)
    tmpb = tmpa*ustar*depresist_unstabpblfact
    if (ioptaa == 2) then
        if (tmp_dvm >= 5.0e-6) then
        stickfrac = exp( -2.0*sqrt(max(0.0,stokes_3)) )
        tmpb = tmpb*stickfrac
        end if
    end if
    depresist_d3 = 1./max( tmpb, 1.e-22 )
    facdepresist_d3 = tmpa

!   compute overall deposition velocity (binkowski/shankar eqn a33)
    tmpa = depresist_a + depresist_d3 +   &
            depresist_a*depresist_d3*vsettl_3
    depvel_3 = vsettl_3 + (1./tmpa)

    tmpa = depresist_a + depresist_d0 +   &
            depresist_a*depresist_d0*vsettl_0
    depvel_0 = vsettl_0 + (1./tmpa)

    return
    end subroutine mosaic_aer_drydepvel_1



!-----------------------------------------------------------------------
    subroutine sorgam_aer_drydep_prep(                                 &
        id, ktau, dtstep, config_flags, aer_mech_id,               &
        gmt, julday,                                               &
        t_phy, rho_phy, p_phy, alt, p8w, t8w,                      &
        moist, chem,                                               &
        h2oai, h2oaj,                                              &
        nbin_add, alnsg_add, wetdgnum_add, wetdens_add,            &
        yextra_add,                                                &
        i, j,                                                      &
        ids,ide, jds,jde, kds,kde,                                 &
        ims,ime, jms,jme, kms,kme,                                 &
        its,ite, jts,jte, kts,kte                                  )
    implicit none
!   subr arguments
    integer, intent(in) ::   &
        id, ktau, julday, aer_mech_id,   &
        nbin_add, i, j
    integer, intent(in) ::   &
        ids, ide, jds, jde, kds, kde,   &
        ims, ime, jms, jme, kms, kme,   &
        its, ite, jts, jte, kts, kte
    real, intent(in) :: dtstep, gmt
    real, intent(in),dimension( ims:ime, kms:kme, jms:jme ) :: &
        t_phy, rho_phy, p_phy, alt, p8w, t8w
    real, intent(in),   &
        dimension( ims:ime, kms:kme, jms:jme, 1:num_moist ) :: &
        moist
    real, intent(in),   &
        dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
        chem
    real, intent(in),   &
        dimension( ims:ime, kms:kme, jms:jme ) :: &
        h2oai, h2oaj
    real, intent(inout),   &
        dimension( 1:nbin_add ) :: &
        alnsg_add, wetdgnum_add, wetdens_add
    real, intent(inout),   &
        dimension( 1:nbin_add, 10 ) :: &
        yextra_add
    type(MPAS_chem_config), intent(in) :: config_flags
    end subroutine sorgam_aer_drydep_prep

!-----------------------------------------------------------------------
    subroutine sorgam_aer_drydep_load_ddvel(                           &
        id, ktau, dtstep, config_flags, aer_mech_id,               &
        numgas, ddvel,                                             &
        nbin_add, depvel0_add, depvel3_add,                        &
        i, j,                                                      &
        ids,ide, jds,jde, kds,kde,                                 &
        ims,ime, jms,jme, kms,kme,                                 &
        its,ite, jts,jte, kts,kte                                  )

!!!!use module_configure, only:  num_moist, num_chem, &
!!!!    MPAS_chem_config
!!!!use module_state_description, only:  param_first_scalar, &
!!!!    p_naai, p_naaj, p_clai, p_claj

!use module_data_sorgam

    implicit none

!   subr arguments
    integer, intent(in) ::   &
        id, ktau, aer_mech_id, &
        numgas, nbin_add, i, j

    integer, intent(in) ::   &
        ids, ide, jds, jde, kds, kde,   &
        ims, ime, jms, jme, kms, kme,   &
        its, ite, jts, jte, kts, kte

    real, intent(in) :: dtstep

    real, intent(in),   &
        dimension( 1:nbin_add ) :: &
        depvel0_add, depvel3_add

    real, intent(inout),   &
        dimension( its:ite, jts:jte, 1:num_chem ) :: &
        ddvel

    type(MPAS_chem_config), intent(in) :: config_flags
    end subroutine sorgam_aer_drydep_load_ddvel



!------------------------------------------------------------------------
    subroutine sorgam_aer_drydepvel_1(   &
        iflag_coarse,   &
        dgnum, alnsg, aerodens,   &
        temp, pres, airdens, ustar, wstar, depresist_a,   &
        airkinvisc, freepath,   &
        depvel_0, depvel_3,   &
        depresist_d0, depresist_d3,   &
        vsettl_0, vsettl_3 )
    implicit none

    integer, intent(in) ::   &
        iflag_coarse
    real, intent(in) ::   &
        dgnum, alnsg, aerodens,   &
        temp, pres, airdens, ustar, wstar, depresist_a
    real, intent(inout) ::   &
        airkinvisc, freepath
    real, intent(inout) ::   &
        depvel_0, depvel_3,   &
        depresist_d0, depresist_d3,   &
        vsettl_0, vsettl_3


    end subroutine sorgam_aer_drydepvel_1



!------------------------------------------------------------------------
    subroutine sorgam_aer_drydepvel_2(   &
        iflag_coarse,   &
        dgnum, alnsg, aerodens,   &
        temp, pres, airdens, ustar, depresist_a,   &
        pblh, zntt, rmolm,   &
        airkinvisc, freepath,   &
        depvel_0, depvel_3,   &
        depresist_d0, depresist_d3,   &
        vsettl_0, vsettl_3 )
!
!   computes the surface layer resistance term and the
!   gravitational settling velocity term for the 0th and 3rd moment
!   of a log-normal aerosol mode
!
!   input parameters
!   dgnum - geometric mean diameter for aerosol number (m)
!   alnsg - natural logarithm of the geometric standard deviation
!       for aerosol number
!   aerodens - aerosol density (dgnum and aerodens are for the
!       actual wet distribution)
!   temp - temperature (K)
!   pres - pressure (Pa)
!   airdens - air density (kg/m^3)
!   ustar - friction velocity (m/s)
!   wstar - convective scale velocity (m/s)
!   depresist_a = aerodynamic resistance (s/m)
!   airkinvisc - air kinematic viscosity (m^2/s)
!   freepath - air molecular freepath (m)
!   output parameters
!   depvel_0 - deposition velocity for 0th (number) moment (m/s)
!   depvel_3 - deposition velocity for 3rd (mass  ) moment (m/s)
!   depresist_d0/_d3 - surface layer resistance for 0th/3rd moment (s/m)
!   vsettl_0/_3 - gravitational settling velocity for 0th/3rd moment (m/s)
!   

!   use module_data_sorgam, only: boltz, grav, pss0, sqrt2, sqrtpi, tss0, threepi, two3
!   use module_data_sorgam, only: ngausdv, wgaus, y_gq

    implicit none

    integer, intent(in) ::   &
        iflag_coarse
    real, intent(in) ::   &
        dgnum, alnsg, aerodens,   &
        temp, pres, airdens, ustar, depresist_a, pblh, zntt, rmolm
    real, intent(inout) ::   &
        airkinvisc, freepath
    real, intent(inout) ::   &
        depvel_0, depvel_3,   &
        depresist_d0, depresist_d3,   &
        vsettl_0, vsettl_3


    ! idowescor = 1 means do wesley (85) convective correction to pm dry dep velocities; 0 means dont do correction
    integer, parameter :: idowescor = 0

! following are from module_data_sorgam
! avogadro's constant [ 1/mol ]
        real, parameter :: avo=6.0221367e23
! universal gas constant [ j/mol-k ]
        real, parameter :: rgasuniv=8.314510
! boltzmann's constant [ j / k]
        real, parameter :: boltz=rgasuniv/avo
! mean gravitational acceleration [ m/sec**2 ]
        real, parameter :: grav=9.80622
!  starting standard surface pressure [ pa ]
        real, parameter :: pss0=101325.0
!  starting standard surface temperature [ k ]
        real, parameter :: tss0=288.15
        real*8, parameter :: pirs=3.14159265358979324
    real, parameter :: sqrt2=1.4142135623731
    real, parameter :: sqrtpi=1.7724539
        real, parameter :: two3=2.0/3.0
        real, parameter :: threepi=3.0*pirs


    real, parameter :: bhat = 1.246   ! constant from cunningham slip correction.
    real, parameter :: colctr_bigd=2.e-3,colctr_smald=20.e-6  
                     ! collector diameters in stokes number and interception efficiency (needleleaf forest)

    integer :: n

    real :: alnsgy2, amu
    real :: blkta, blkprs, blkdens
    real :: cunq, czh 
    real :: dconst1, dconst2, dconst3, dconst3y
    real :: dgyyy, dq
    real :: eff_dif, eff_imp, eff_int
    real :: ey1, esy04, esy08, esy16, esy20, esy28, esy32, esy64
    real :: knyyy, knq
    real :: nu
    real :: pdensy
    real :: ra, rbcor, rsurfq
    real :: scq, stq, sum0, sum0b, sum3, sum3b
    real :: tmpa
    real :: utscale
    real :: vdplim, vsedq
    real :: xlm
    real :: yextra_vdvg(100,3)

    integer, parameter :: ngausdv = 7
    real, save :: y_gq(ngausdv), wgaus(ngausdv)


    ! define these locally because those in module_data_sorgam
    ! are not always initialized
    y_gq = (/ -2.651961356835233, &
    -1.673551628767471, -0.816287882858965, -0.0, &
    0.816287882858965, 1.673551628767471, 2.651961356835233 /)
    wgaus = (/ 0.0009717812450995, &
    0.05451558281913, 0.4256072526101, 0.8102646175568, &
    0.4256072526101, 0.05451558281913, 0.0009717812450995 /)

    blkta = temp
    blkprs = pres
    blkdens = airdens

!   amu = air dynamic viscosity (kg/m/s)
    if (airkinvisc <= 0.0) then
        amu = 1.458e-6*blkta*sqrt(blkta)/(blkta+110.4)
        airkinvisc = amu/airdens
    else
        amu = airkinvisc*airdens
    end if

!   air molecular freepath (m)
    if (freepath <= 0.0) then
        freepath = 6.6328e-8*pss0*blkta/(tss0*blkprs)
    end if
    xlm = freepath

    ra = depresist_a

    dgyyy = dgnum
    pdensy = aerodens

        alnsgy2 = alnsg**2
        ey1 = exp(0.125*alnsgy2)
        esy04 = ey1**4
        esy08 = esy04*esy04
        esy16 = esy08*esy08
        esy20 = esy16*esy04
        esy28 = esy20*esy08
        esy32 = esy16*esy16
        esy64 = esy32*esy32

    dconst1 = boltz * blkta / ( threepi * amu )
    dconst2 = grav / ( 18.0 * amu )
    dconst3 =  ustar/(9.*amu*colctr_bigd)

    nu = amu / blkdens
                                                                                                                                 
    utscale =  1.
    if (idowescor.eq.1) then
    ! wesley (1985) monin-obukov dependence for convective conditions (sam 10/08)
        if (rmolm.lt.0.) then
        czh = -1.*pblh*rmolm
        if (czh.gt.30.0) then
            utscale=0.45*(czh**0.6667)
        else
            utscale=1.+((-300.*rmolm)**0.6667)
        endif
        endif
    endif   ! end of (idowescor.eq.1) test
    utscale = ustar*utscale

        sum0 = 0. ; sum0b = 0.
        sum3 = 0. ; sum3b = 0.

        do n = 1,ngausdv
        dq = dgyyy*exp(y_gq(n)*sqrt2*alnsg)  ! diameter (m) at quadrature point
        knq = 2.*xlm/dq  ! knudsen number at quadrature point
        cunq = 1.+knq*(1.257+.4*exp(-1.1/knq))  ! cunningham correction factor; pruppacher and klett (1980) eq (12-16)
        vsedq = pdensy*dconst2*cunq*dq*dq  ! gravitational sedimentation velocity m/s
        scq = nu*dq/dconst1/cunq  ! schmidt number, brownian diffusion parameter - same as binkowski and shankar
        eff_dif = scq**(-two3)    ! efficiency term for diffusion - same as binkowski and shankar
        stq = dconst3*pdensy*dq**2  ! stokes number, peters and eiden (1992)
        eff_imp = (stq/(0.8+stq))**2   ! efficiency term for impaction - peters and eiden (1992)
!       eff_int = 0.3*dq/(colctr_smald+dq) ! slinn (1982) interception term, 0.3 prefac insures .2 cm/s at .2 micron diam.
        eff_int = (0.00116+0.0061*zntt)*dq/1.414e-7 ! mckeen(2008) intercptn trm, val of .00421 @ ustr=0.475, diam=.1414 micrn, stable, needleleaf evergreen
        if (iflag_coarse > 0) then
           eff_int = min(1.,eff_int)
           rbcor=exp(-2.0*(stq**0.5)) ! rebound correction factor used in slinn (1982)
           !czhao add limit
           rbcor=max(rbcor,1.e-10)
        else
           rbcor = 1. ! rebound correction factor
        endif
        vdplim = utscale*(eff_dif+eff_imp+eff_int)*rbcor
!       vdplim = .002*utscale
        vdplim = min(vdplim,.02)
!PMA
        vdplim = max(vdplim,1.e-20)
        rsurfq = ra+1./vdplim
        sum0 = sum0 + wgaus(n)*(vsedq + 1./rsurfq)  ! quadrature sum for 0 moment
        sum3 = sum3 + wgaus(n)*(vsedq + 1./rsurfq)*dq**3  ! quadrature sum for 3rd moment
        sum0b = sum0b + wgaus(n)*(1./rsurfq)  ! quadrature sum for 0 moment
        sum3b = sum3b + wgaus(n)*(1./rsurfq)*dq**3  ! quadrature sum for 3rd moment
        end do ! n

    ! normalize 0 moment vdep quadrature sum to sqrt(pi) (and number =1 per unit volume)
        depvel_0 = sum0/sqrtpi  
    !normalize 3 moment quad. sum to sqrt(pi) and 3rd moment analytic sum
        depvel_3 = sum3/(sqrtpi*exp((1.5*sqrt2*alnsg)**2)*dgyyy**3) 

        sum0b = sum0b/sqrtpi  
        sum3b = sum3b/(sqrtpi*exp((1.5*sqrt2*alnsg)**2)*dgyyy**3) 

    tmpa = 1.0/max( sum0b, 1.0e-10 )
    depresist_d0 = max( (tmpa - ra), 1.0e-10 )
    tmpa = 1.0/max( sum3b, 1.0e-10 )
    depresist_d3 = max( (tmpa - ra), 1.0e-10 )


! calculate settling velocity 
!   knyyy = 2.0*xlm/dgyyy
!   dconst3y = dconst2 * pdensy * dgyyy**2
!   vsettl_0 = dconst3y * ( esy16  + bhat * knyyy * esy04 )
!   vsettl_3 = dconst3y * ( esy64  + bhat * knyyy * esy28 )

! following are consistent with the quadrature calcs
    vsettl_0 = max( (depvel_0 - sum0b), 0.0 )
    vsettl_3 = max( (depvel_3 - sum3b), 0.0 )


    return
    end subroutine sorgam_aer_drydepvel_2



!-----------------------------------------------------------------------
     subroutine zhang2001_aer_drydepvel_1( &
         dg_num_in, lnsig_in, dens_part_in, &
         airtemp_in, airdens_in, ustar_in, resist_a_in, &
         zhang_optaa, iseason_in, &
         isurftype_in, modesurftype, iwetsurf, &
         airkinvisc, freepath,   &
         ido_mom0, ido_mom3,   &
         depvel_0, depvel_3,   &
         resist_s0, resist_s3,   &
         vsettl_0, vsettl_3 )
!    calculates surface deposition velocity of particles
!    L. Zhang, S. Gong, J. Padro, and L. Barrie
!    A size-seggregated particle dry deposition scheme for an atmospheric aerosol module
!    Atmospheric Environment, 35, 549-560, 2001.
! use precision
! use physconst, only: gravit, tmelt, rhodair

      implicit none

!     in
      real     :: dg_num_in    ! geometric mean diameter of number distribution (m)
      real     :: lnsig_in     ! nat. log. of geometric standard deviation of particles
      real     :: dens_part_in ! density of particle material (kg/m3)
      real     :: airdens_in   ! air density (kg/m3)
      real     :: airtemp_in   ! air temperature (K)
      integer  :: zhang_optaa  ! 1=as published; 2=modified
      integer  :: iseason_in   ! seasonal category (same as wesely)
      integer  :: isurftype_in ! surface type id
      integer  :: modesurftype ! if <= 1, isurftype is zhang et al surface type
                               ! if >= 2, isurftype is usgs surface type
      integer  :: iwetsurf     ! 1 if surface is wet, 0 if dry
      real     :: resist_a_in  ! aerodynamic resistance (s/m)
      real     :: ustar_in     ! friction velocity (m/s)
      real     :: airkinvisc   ! kinematic viscosity of air (m2/s)
      real     :: freepath     ! mean free path of air molecules (m)
      integer  :: ido_mom0, ido_mom3 ! flag for do/skip moment 0 & 3 calculations

!     out
      real     :: depvel_0, depvel_3    ! surface deposition velocity for 
                                        ! 0th (number) and 3rd (mass) moments (m/s)
      real     :: resist_s0, resist_s3  ! surface resistances (s/m)
      real     :: vsettl_0, vsettl_3    ! gravitational settling velocities (m/s)

      end subroutine zhang2001_aer_drydepvel_1

!-----------------------------------------------------------------------
        subroutine cam_mam_aer_drydep_prep(                                 &
                id, ktau, dtstep, config_flags, aer_mech_id,               &
                gmt, julday,                                               &
                t_phy, rho_phy, p_phy, alt, p8w, t8w,                      &
                moist, chem,                                               &
                nbin_add, alnsg_add, wetdgnum_add, wetdens_add,            &
                yextra_add,                                                &
                i, j,                                                      &
                ids,ide, jds,jde, kds,kde,                                 &
                ims,ime, jms,jme, kms,kme,                                 &
                its,ite, jts,jte, kts,kte                                  )
        implicit none

!   subr arguments
        integer, intent(in) ::   &
                id, ktau, julday, aer_mech_id, nbin_add,   &
                i, j

        integer, intent(in) ::   &
                ids, ide, jds, jde, kds, kde,   &
                ims, ime, jms, jme, kms, kme,   &
                its, ite, jts, jte, kts, kte

        real, intent(in) :: dtstep, gmt

        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme ) :: &
                t_phy, rho_phy, p_phy, alt, p8w, t8w

        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme, 1:num_moist ) :: &
                moist

        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
                chem

        real, intent(inout),   &
                dimension( 1:nbin_add ) :: &
                alnsg_add, wetdgnum_add, wetdens_add

        real, intent(inout),   &
                dimension( 1:nbin_add, 10 ) :: &
                yextra_add

        type(MPAS_chem_config), intent(in) :: config_flags
        end subroutine cam_mam_aer_drydep_prep

!-----------------------------------------------------------------------
        subroutine cam_mam_aer_drydep_load_ddvel(                           &
                id, ktau, dtstep, config_flags, aer_mech_id,               &
                ddvel,                                                     &
                nbin_add, depvel0_add, depvel3_add,                        &
                i, j,                                                      &
                ids,ide, jds,jde, kds,kde,                                 &
                ims,ime, jms,jme, kms,kme,                                 &
                its,ite, jts,jte, kts,kte                                  )
        implicit none
!   subr arguments
        integer, intent(in) ::   &
                id, ktau, aer_mech_id, nbin_add,   &
                i, j

        integer, intent(in) ::   &
                ids, ide, jds, jde, kds, kde,   &
                ims, ime, jms, jme, kms, kme,   &
                its, ite, jts, jte, kts, kte

        real, intent(in) :: dtstep

        real, intent(in),   &
                dimension( 1:nbin_add ) :: &
                depvel0_add, depvel3_add

        real, intent(inout),   &
                dimension( its:ite, jts:jte, 1:num_chem ) :: &
                ddvel

        type(MPAS_chem_config), intent(in) :: config_flags

        end subroutine cam_mam_aer_drydep_load_ddvel



!-----------------------------------------------------------------------
END MODULE module_aer_drydep

!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC)
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to declare or initialize variables needed by chemistry processes.
! * Author: Jiawang Feng and Chun Zhao
! Last update: Nov 2019
!
! Contact:
! Chun Zhao, PhD
! Professor
! University of Science and Technology of China
! Hefei, Anhui, 230026
! Email: ChunZhao@ustc.edu.cn
!
! Please report any bugs or problems to Chun Zhao

! References:
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

 module mpas_atmchem_vars
 use mpas_kind_types
 use mpas_timekeeping
 use mpas_pool_routines
 use mpas_atmchem_configs
 use mpas_atmphys_vars, only: num_moist
 use mpas_constants , only: pii
 use mpas_log,only : mpas_log_write
 implicit none
 public :: atmchem_vars_init
 save

!----------------------------------------------------------------------------!
!---------------------------Wrf vars-----------------------------------------!
!----------------------------------------------------------------------------!

 type(MPAS_chem_config),public :: config_flags

! Number of different fields.
 integer,public :: num_species
 integer,public :: num_chem
 integer,public :: num_scalars
 integer,public :: num_scalar
 integer,public :: num_vert_mix
 integer,public :: num_gas
!  integer,public :: nswbands, nlwbands

 integer,public :: ntot_amode_cam_mam = 3
 integer,public :: kemit
!  integer,public :: nbin_a

 real(kind=RKIND),public :: emis_dust_factor,dust_emass_jfrac,dust_emass_cfrac

 integer,public :: curr_year            !Current year.
 integer,public :: curr_month           !Current month.
 integer,public :: curr_secs            !Current seconds.
 integer,public :: julday               !Initial Julian day.
 logical,public :: LeapYear             !Is leap year or not.
 logical,dimension(:),allocatable,public :: is_aerosol

 real(kind=RKIND),public :: curr_julday !Current Julian day (= 0.0 at 0Z on January 1st).
 real(kind=RKIND),public :: gmt         !Greenwich mean time hour of model start (hr)
! Species indents used in wrf-chem.

! How to get p_p25???
 integer,parameter :: p_p25              = 1
 integer,public,parameter :: param_first_scalar = 0
 real   ,parameter :: pi                 = pii
 integer,public,pointer :: &
        p_qndrop, p_qniceaer, p_qnicemix                      , &
        p_qv, p_qc,  p_qi, p_qs                               , &
        p_seas_a01,  p_seas_a02,  p_seas_a03,  p_seas_a04                , &
        p_seas_a05,  p_seas_a06,  p_seas_a07,  p_seas_a08                , &
        p_oc_a01,    p_oc_a02,    p_oc_a03,    p_oc_a04                  , &
        p_oc_a05,    p_oc_a06,    p_oc_a07,    p_oc_a08                  , &
        p_dust_a01,  p_dust_a02,  p_dust_a03,  p_dust_a04,  p_dust_a05   , &
        p_dust_a06,  p_dust_a07,  p_dust_a08,  p_dust_a09,  p_dust_a10   , &
        p_hysw_a01,  p_hysw_a02,  p_hysw_a03,  p_hysw_a04,  p_hysw_a05   , &
        p_hysw_a06,  p_hysw_a07,  p_hysw_a08,  p_hysw_a09,  p_hysw_a10   , &
        p_water_a01, p_water_a02, p_water_a03, p_water_a04, p_water_a05  , &
        p_water_a06, p_water_a07, p_water_a08, p_water_a09, p_water_a10  , &
        p_num_a01,   p_num_a02,   p_num_a03,   p_num_a04,   p_num_a05    , &
        p_num_a06,   p_num_a07,   p_num_a08,   p_num_a09,   p_num_a10    , &
        p_seas_cw01, p_seas_cw02, p_seas_cw03, p_seas_cw04               , &
        p_seas_cw05, p_seas_cw06, p_seas_cw07, p_seas_cw08               , &
        p_oc_cw01,   p_oc_cw02,   p_oc_cw03,   p_oc_cw04                 , &
        p_oc_cw05,   p_oc_cw06,   p_oc_cw07,   p_oc_cw08                 , &
        p_dust_cw01, p_dust_cw02, p_dust_cw03, p_dust_cw04, p_dust_cw05  , &
        p_dust_cw06, p_dust_cw07, p_dust_cw08, p_dust_cw09, p_dust_cw10  , &
        p_num_cw01,  p_num_cw02,  p_num_cw03,  p_num_cw04,  p_num_cw05   , &
        p_num_cw06,  p_num_cw07,  p_num_cw08,  p_num_cw09,  p_num_cw10
!jwfeng from AERMOL USTC 10/2021
!Not included species, but will be included in future.
!Set to -1 for now, remember to change them into pointer vars.
 integer,public :: &
        p_co=-1,     p_so2=-1,    p_hno3=-1,   p_hno4=-1                 , &
        p_n2o5=-1,   p_nh3=-1,    p_h2o2=-1,   p_o3=-1,      p_ora1=-1   , &
        p_op1=-1,    p_paa=-1,    p_sulf=-1,   p_so4aj=-1,   p_nh4aj=-1  , &
        p_no3aj=-1,  p_bc1=-1,    p_bc2=-1,    p_oc1=-1,     p_oc2=-1    , &
        p_seas_1=-1, p_seas_2=-1, p_seas_3=-1, p_seas_4=-1 , p_dms=-1
!jwfeng from AEMOL USTC 12/2019
!----------------------------------------------------------------------------!
! Vars from MPAS to Wrf------------------------------------------------------!
! In MPAS:dimension i&j are replaced by nCells-----------------------------!
!----------------------------------------------------------------------------!
! WRF Dimensions:longitude,latitude (i,j)
 real(kind=RKIND),dimension(:,:),allocatable ::     &
    u10_c,                                          &!u at 10 m                [m/s]
    v10_c,                                          &!v at 10 m                [m/s]
    ust_c,                                          &!u* in similarity theory  [m/s]
    xland_c,                                        &
    dms_0_c,                                        &
    sstemp_c,                                       &
    xlat_c,                                         &
    xlong_c,                                        &
    dx_c,                                           &
    hpbl_c,                                         &
    znt_c,                                          &
    psfc_c,                                         &
    psfc_hyd_c,                                     &
    psfc_hydd_c,                                    &
    raincv_c,                                       &
    snowh_c,                                        &
    tsk_c,                                          &
    gsw_c,                                          &
    hfx_c,                                          &
    vegfra_c,                                       &
    rmol_c,                                         &
    slai_c,                                         &
    e_seas_c,                                       &
    e_dust_c

 integer,dimension(:,:),allocatable ::              &
    ivgtyp_c, isltyp_c

!----------------------------------------------------------------------------!
! WRF Dimensions:Longitude,Vert,Latitude (i,k,j)
 real(kind=RKIND),dimension(:,:,:),allocatable ::   &
    w_c,                                            &
    z_c,                                            &
    qv_c,                                           &
    alt_c,                                          &!alt = 1.0/rho in (m3/kg)
    rho_phy_c,                                      &
    t_phy_c,                                        &
    th_c,                                           &
    t8w_c,                                          &
    t2_c,                                           &
    dz_c,                                           &!dz = layer thickness in(m)
    pres_c,                                         &
    pres2_c,                                        &
    pres2_hyd_c,                                    &
    pres2_hydd_c,                                   &
    pres_hyd_c,                                     &
    pres_hydd_c,                                    &
    znu_hyd_c,                                      &
    p8w_c,                                          &
    erod_c,                                         &
    exch_h_c,                                       &
    u_phy_c,                                        &
    v_phy_c,                                        &
    ddvel_c,                                        &!deposition velocities(only on surface)
    smois_c,                                        &
    qlsink_c,                                       &
    ilsink_c,                                       &
    precr_c,                                        &
    preci_c,                                        &
    precs_c,                                        &
    precg_c,                                        &
!jwfeng : Start vars in module_mixactivate.F
      cldfra_c,                                     &!horizontal cloud fraction
  cldfra_old_c,                                     &!horizontal cloud fraction
       f_ice_c,                                     &!horizontal cloud fraction
      f_rain_c,                                     &!horizontal cloud fraction
     tke_pbl_c,                                     &!turbulent kinetic energy from PBL
    relhum3d_c,                                     &!Relative humidity respect to ice
      nihf3d_c,                                     &
     niimm3d_c,                                     &
     nidep3d_c,                                     &
      wsub3d_c,                                     &
        ccn1_c,                                     &
        ccn2_c,                                     &
        ccn3_c,                                     &
        ccn4_c,                                     &
        ccn5_c,                                     &
        ccn6_c,                                     &
     nsource_c,                                     &
      rndst1_c,                                     &
      rndst2_c,                                     &
      rndst3_c,                                     &
      rndst4_c,                                     &
      nacon1_c,                                     &
      nacon2_c,                                     &
      nacon3_c,                                     &
      nacon4_c,                                     &
aq_gas_ratio_c

 real(kind=RKIND),dimension(:,:,:),allocatable :: &
    tend_wetconv_c
!jwfeng : Start vars in module_mixactivate.F
!----------------------------------------------------------------------------!
! WRF Dimensions:Longitude,Vert,Latitude,num_*** (i,k,j,l)
 real(kind=RKIND),dimension(:,:,:,:),allocatable :: &
    tend_emis_c,                                    &
    tend_ctrans_c,                                  &
    moist_c,                                        &
    scalars_c,                                      &
    chem_c,                                         &
    tauaersw_c,                                     &
    absaersw_c,                                     &
    gaersw_c,                                       &
    waersw_c,                                       &
    tauaersw_dust_c,                                &
    absaersw_dust_c,                                &
    gaersw_dust_c,                                  &
    waersw_dust_c,                                  &
    tauaersw_seas_c,                                &
    absaersw_seas_c,                                &
    gaersw_seas_c,                                  &
    waersw_seas_c,                                  &
    extaersw_c,                                     &
    bscoefsw_c,                                     &
    tauaerlw_c,                                     &
    tauaerlw_dust_c,                                &
    tauaerlw_seas_c,                                &
    extaerlw_c,                                     &
    l2aer_c,                                        &
    l3aer_c,                                        &
    l4aer_c,                                        &
    l5aer_c,                                        &
    l6aer_c,                                        &
    l7aer_c

!... jwfeng:2021-06-03, add for kain fritsch
 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    umf_c,            &
    uer_c,            &
    udr_c,            &
    dmf_c,            &
    der_c,            &
    ddr_c,            &
     wu_c

 real(kind=RKIND),dimension(:,:),allocatable:: &
    nca_c,            &
    cutop_c,          &
    cubot_c,          &
    xmb_shalc_c

 integer,dimension(:,:),allocatable:: &
    ldb_c,            &
    let_c,            &
    ksl_c,            &
    lfs_c,            &
    klc_c,            &
    ishall_c,         &
    check_kfeta_c,    &
    k22_shalc_c,      &
    kbcon_shalc_c,    &
    ktop_shalc_c,     &
    kbcon_deepc_c,    &
    ktop_deepc_c

 contains
!=================================================================================================================
subroutine atmchem_vars_init(configs,mesh,state,diag_chemistry,clock)
!=================================================================================================================
  type(mpas_pool_type),intent(in) :: mesh
  type(mpas_pool_type),intent(in) :: configs
  type(mpas_pool_type),intent(in) :: state
  type(mpas_pool_type),intent(inout) :: diag_chemistry

  integer,pointer :: nCellsSolve,nVertLevels,nSoils,nChem,nScalars
  integer,pointer :: species_start,species_end
  integer,pointer :: moist_start,moist_end
  integer         :: nSpecies,nMoist
  real(kind=RKIND),pointer :: config_dt

  type(MPAS_Time_Type) :: currTime
  type(MPAS_Time_Type) :: startTime
  type(MPAS_Clock_type),intent(in) :: clock

  character(len=StrKIND) :: timeStamp
  integer :: ierr
  real(kind=RKIND) :: utc_s,utc_h
  integer :: h,m,s,s_n,s_d,DoY,yr,month
  integer :: n,nhi,itype,isize
  real    :: dum

!local config variables
  integer,pointer :: gas_drydep_opt,aer_drydep_opt,kdepvel,settling_option
!   integer,pointer :: nswbands_pointer, nlwbands_pointer
  character(len=StrKIND),pointer :: chem_opt

!jwfeng: 2020/10/09 add for aerosol optical variables
!local pointers:
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaersw, absaersw, extaersw, gaersw, waersw, bscoefsw
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaerlw, extaerlw
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaersw_dust, absaersw_dust, gaersw_dust, waersw_dust
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaersw_seas, absaersw_seas, gaersw_seas, waersw_seas
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaerlw_dust, tauaerlw_seas
 real(kind=RKIND),dimension(:,:,:),pointer :: l2aer, l3aer, l4aer, l5aer, l6aer, l7aer

  call mpas_pool_get_config(configs  ,     'config_dt'     ,      config_dt      )
  call mpas_pool_get_dimension(state ,     'species_start' ,      species_start  )
  call mpas_pool_get_dimension(state ,     'species_end'   ,      species_end    )
  call mpas_pool_get_dimension(state ,     'moist_start'   ,      moist_start    )
  call mpas_pool_get_dimension(state ,     'moist_end'     ,      moist_end      )
  call mpas_pool_get_dimension(state ,     'num_chem'      ,      nChem          )
  call mpas_pool_get_dimension(state ,     'num_scalars'   ,      nScalars       )
  call mpas_pool_get_dimension(mesh  ,     'nSoilLevels'   ,      nSoils         )
  call mpas_pool_get_dimension(mesh  ,     'nCellsSolve'   ,      nCellsSolve    )
  call mpas_pool_get_dimension(mesh  ,     'nVertLevels'   ,      nVertLevels    )

  call mpas_pool_get_config(configs,'config_chem_opt'        ,chem_opt           )
  call mpas_pool_get_config(configs,'config_gas_drydep_opt'  ,gas_drydep_opt     )
  call mpas_pool_get_config(configs,'config_aer_drydep_opt'  ,aer_drydep_opt     )
  call mpas_pool_get_config(configs,'config_kdepvel'         ,kdepvel            )
  call mpas_pool_get_config(configs,'config_settling_option' ,settling_option    )
!jwfeng: 2020/10/09 add for aerosol optical variables
  call mpas_pool_get_array(diag_chemistry,'tauaersw'  , tauaersw )
  call mpas_pool_get_array(diag_chemistry,'absaersw'  , absaersw )
  call mpas_pool_get_array(diag_chemistry,'gaersw'    , gaersw   )
  call mpas_pool_get_array(diag_chemistry,'waersw'    , waersw   )
  call mpas_pool_get_array(diag_chemistry,'tauaersw_dust', tauaersw_dust )
  call mpas_pool_get_array(diag_chemistry,'absaersw_dust', absaersw_dust )
  call mpas_pool_get_array(diag_chemistry,'gaersw_dust'  , gaersw_dust   )
  call mpas_pool_get_array(diag_chemistry,'waersw_dust'  , waersw_dust   )
  call mpas_pool_get_array(diag_chemistry,'tauaersw_seas', tauaersw_seas )
  call mpas_pool_get_array(diag_chemistry,'absaersw_seas', absaersw_seas )
  call mpas_pool_get_array(diag_chemistry,'gaersw_seas'  , gaersw_seas   )
  call mpas_pool_get_array(diag_chemistry,'waersw_seas'  , waersw_seas   )
  call mpas_pool_get_array(diag_chemistry,'extaersw'  , extaersw )
  call mpas_pool_get_array(diag_chemistry,'bscoefsw'  , bscoefsw )
  call mpas_pool_get_array(diag_chemistry,'tauaerlw'  , tauaerlw )
  call mpas_pool_get_array(diag_chemistry,'tauaerlw_dust'  , tauaerlw_dust )
  call mpas_pool_get_array(diag_chemistry,'tauaerlw_seas'  , tauaerlw_seas )
  call mpas_pool_get_array(diag_chemistry,'extaerlw'  , extaerlw )
  call mpas_pool_get_array(diag_chemistry,'l2aer'     , l2aer    )
  call mpas_pool_get_array(diag_chemistry,'l3aer'     , l3aer    )
  call mpas_pool_get_array(diag_chemistry,'l4aer'     , l4aer    )
  call mpas_pool_get_array(diag_chemistry,'l5aer'     , l5aer    )
  call mpas_pool_get_array(diag_chemistry,'l6aer'     , l6aer    )
  call mpas_pool_get_array(diag_chemistry,'l7aer'     , l7aer    )

  config_flags%chem_opt         = chem_opt
  config_flags%chem_in_opt      = 1
  config_flags%gas_drydep_opt   = gas_drydep_opt
  config_flags%aer_drydep_opt   = aer_drydep_opt

  config_flags%kdepvel          = kdepvel
  config_flags%settling_option  = settling_option

!!! All of these will be passed from Registry.xml
!!! and are belonged to phys in WRF-chem

! if(radt_sw_scheme == ...)
  config_flags%ra_sw_physics    = 1
! end if
  config_flags%bl_pbl_physics   = 1
  config_flags%wsubscheme       = 0
  config_flags%ice_para         = 1
  config_flags%mixcloud         = 0

  nSpecies = species_end - species_start + 1
  nMoist   =   moist_end -   moist_start + 1

!   nswbands = nswbands_pointer
!   nlwbands = nlwbands_pointer

! Get start gmt and julday:
  startTime   = mpas_get_clock_time(clock,MPAS_START_TIME, ierr)
  call mpas_get_time(curr_time=startTime,YYYY=yr,MM=month,H=h,M=m,S=s,S_n=s_n,S_d=s_d,DoY=DoY,ierr=ierr)
  gmt         = real(h) + real(m) / 60.0 + real(s + s_n / s_d) / 3600.0
  julday      = DoY

! Update the current julian day and current year:
  currTime    = mpas_get_clock_time(clock,MPAS_NOW,ierr)
  call mpas_get_time(curr_time=currTime,dateTimeString=timeStamp,YYYY=yr,MM=month,H=h,M=m,&
                     S=s,S_n=s_n,S_d=s_d,DoY=DoY,ierr=ierr)
  utc_h       = real(h) + real(m) / 60.0 + real(s + s_n / s_d) / 3600.0
  utc_s       = real(s + s_n + s_d)
  curr_year   = yr
  curr_month  = month
  curr_secs   = utc_s
  curr_julday = real(DoY-1) + utc_h / 24.0
  LeapYear    = isLeapYear(yr)

  num_species      = nSpecies
  num_chem         = nChem
  num_scalars      = nScalars
  num_scalar       = num_scalars
  num_moist        = nMoist ! num_moist is declared in mpas_atmphys_vars.F
  num_gas          = 0
  num_vert_mix     = num_gas

  if(.not.allocated(is_aerosol)) allocate(is_aerosol(1:num_chem))
  is_aerosol(:)    = .false.

  dust_emass_jfrac = 0.07
  dust_emass_cfrac = 0.93
!jwfeng: 2020/10/09 add for aerosol optical variables
  tauaersw(:,:,:) = 0.0
  absaersw(:,:,:) = 0.0
  extaersw(:,:,:) = 0.0
    gaersw(:,:,:) = 0.0
    waersw(:,:,:) = 0.0
  bscoefsw(:,:,:) = 0.0
  extaerlw(:,:,:) = 0.0
  tauaerlw(:,:,:) = 0.0
     l2aer(:,:,:) = 0.0
     l3aer(:,:,:) = 0.0
     l4aer(:,:,:) = 0.0
     l5aer(:,:,:) = 0.0
     l6aer(:,:,:) = 0.0
     l7aer(:,:,:) = 0.0
  call mpas_log_write( ' ' )
  call mpas_log_write( '--- Basic settings of chemistry species' )
  call mpas_log_write( '       nSpecies = $i',intArgs = (/num_species/) )
  call mpas_log_write( '         nMoist = $i',intArgs = (/num_moist  /) )
  call mpas_log_write( '       nScalars = $i',intArgs = (/num_scalars/) )
  call mpas_log_write( '           nGas = $i',intArgs = (/num_gas    /) )
  call mpas_log_write( '          nChem = $i',intArgs = (/num_chem   /) )
  call mpas_log_write( '       chem_opt = '//trim(config_flags%chem_opt) )

end subroutine atmchem_vars_init

!=================================================================================================================
 end module mpas_atmchem_vars
!=================================================================================================================

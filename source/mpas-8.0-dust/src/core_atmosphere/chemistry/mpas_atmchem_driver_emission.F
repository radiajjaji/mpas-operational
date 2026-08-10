!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to compute dry deposition tendencies for chem(:,:,:,:). 
! * Author: Jiawang Feng and Chun Zhao from AEMOL USTC
! Last update: jwfeng 11/2019
!
! Contact:
! Chun Zhao, PhD
! Professor 
! University of Science and Technology of China 
! Hefei, Anhui, 230026 
! Email: ChunZhao@ustc.edu.cn
!
! References: 
!     The code structure is similar to the microphysics parameterization drivers
!     in MPAS-v5.2:src/core_atmosphere/physics/physics_wrf/mpas_atmphys_driver_*.F
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

 module mpas_atmchem_driver_emission
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines
 use mpas_atmchem_vars
 use mpas_atmchem_constants
 use mpas_atmphys_vars
 use mpas_atmchem_configs
 use mpas_atmchem_interface
 use mpas_log, only : mpas_log_write
!chemistry schemes:
 use module_mosaic_addemiss
 use module_data_binaer_asect, only : lptr_seas_aer

 implicit none
 private
 public:: allocate_emission,   &
          deallocate_emission, &
          driver_emission

!MPAS driver for chemistry emission.

! subroutines in mpas_atmchem_driver_emission:
! ---------------------------------------
! allocate_emission  : allocate local arrays for emission.
! deallocate_emission: deallocate local arrays for emission.
! driver_emission    : main driver (called from subroutine chemistry_driver).
! emission_from_MPAS : initialize local arrays.
! emission_to_MPAS   : copy local arrays to MPAS arrays.
!
! chemistry schemes called from driver_emission:
! -----------------------------------
! * module_mosaic_addemiss : mosaic emission scheme.

 contains

!=================================================================================================================
 subroutine allocate_emission
!=================================================================================================================
    if(.not.allocated(    dms_0_c))    allocate(    dms_0_c(ims:ime,           jms:jme          ) )
    if(.not.allocated(   sstemp_c))    allocate(   sstemp_c(ims:ime,           jms:jme          ) )
    if(.not.allocated(     erod_c))    allocate(     erod_c(ims:ime,           jms:jme, 3       ) )
    if(.not.allocated(    smois_c))    allocate(    smois_c(ims:ime,num_soils, jms:jme          ) )
    if(.not.allocated(tend_emis_c))    allocate(tend_emis_c(ims:ime,kms:kme,   jms:jme, num_chem) )
    if(.not.allocated(   e_dust_c))    allocate(   e_dust_c(ims:ime,           jms:jme          ) )
    if(.not.allocated(   e_seas_c))    allocate(   e_seas_c(ims:ime,           jms:jme          ) )
 end subroutine allocate_emission

!=================================================================================================================
 subroutine deallocate_emission
!=================================================================================================================
    if(     allocated(    dms_0_c))    deallocate(    dms_0_c)
    if(     allocated(   sstemp_c))    deallocate(   sstemp_c)
    if(     allocated(     erod_c))    deallocate(     erod_c)
    if(     allocated(    smois_c))    deallocate(    smois_c)
    if(     allocated(tend_emis_c))    deallocate(tend_emis_c)
    if(     allocated(   e_dust_c))    deallocate(   e_dust_c)
    if(     allocated(   e_seas_c))    deallocate(   e_seas_c)
 end subroutine deallocate_emission

!=================================================================================================================
 subroutine emission_from_MPAS(state,sfc_input,diag_chemistry,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(in):: diag_chemistry
 integer,intent(in):: its,ite

!  integer,pointer:: index_dust_a01
 real(kind=RKIND),dimension(:),pointer    :: dms_0,sst,e_dust,e_seas
 real(kind=RKIND),dimension(:,:),pointer  :: smois,erod

!local variables:
 integer:: i,k,j,h

!-----------------------------------------------------------------------------------------------------------------
!  call mpas_pool_get_dimension(state,'index_dust_a01',index_dust_a01)
 call mpas_pool_get_array(sfc_input,'erod' ,erod )
 call mpas_pool_get_array(sfc_input,'dms_0',dms_0)
 call mpas_pool_get_array(sfc_input,'smois',smois)
 call mpas_pool_get_array(sfc_input,'sst'  ,sst  )
 call mpas_pool_get_array(diag_chemistry,'e_dust',e_dust)
 
 if ( lptr_seas_aer(1,1,1)>param_first_scalar ) then
    call mpas_pool_get_array(diag_chemistry,'e_seas',e_seas)
 endif

 do h = 1,num_chem
 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    tend_emis_c(i,k,j,h) = 0._RKIND
 enddo
 enddo
 enddo
 enddo

 do j = jts,jte
 do k = 1,num_soils
 do i = its,ite
    smois_c(i,k,j) = smois(k,i)
 enddo
 enddo
 enddo

 do j = jts,jte
 do i = its,ite
    dms_0_c(i,j)  = dms_0(i)
    sstemp_c(i,j) = sst(i)
    e_dust_c(i,j) = 0._RKIND
    e_seas_c(i,j) = 0._RKIND
    do h = 1,3
        erod_c(i,j,h) = erod(h,i)
    enddo
 enddo
 enddo

 end subroutine emission_from_MPAS

!=================================================================================================================
 subroutine emission_to_MPAS(debug_on,state,tend_chemistry,diag_chemistry,its,ite)
!=================================================================================================================

 type(mpas_pool_type),intent(in):: state
!inout arguments:
 type(mpas_pool_type),intent(inout):: tend_chemistry
 type(mpas_pool_type),intent(inout):: diag_chemistry

 integer,intent(in):: its,ite
!local variables:
 integer:: i,k,j,h
 logical:: debug_on

!local pointers:
 real(kind=RKIND),dimension(:,:,:),pointer:: tend_emis
 real(kind=RKIND),dimension(:),pointer:: e_dust
 real(kind=RKIND),dimension(:),pointer:: e_seas
 character(len=StrKIND):: debugString
!-----------------------------------------------------------------------------------------------------------------
!  call mpas_pool_get_dimension(state,'index_dust_a01' ,index_dust_a01)
 call mpas_pool_get_array(tend_chemistry,'tend_emis' ,tend_emis)
 call mpas_pool_get_array(diag_chemistry,'e_dust'    ,e_dust)
 if ( lptr_seas_aer(1,1,1)>param_first_scalar ) then
    call mpas_pool_get_array(diag_chemistry,'e_seas'    ,e_seas)
 endif

 do h = 1,num_chem 
 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    tend_emis(h,k,i)  = tend_emis_c(i,k,j,h)
 enddo
 enddo
 enddo
 enddo

 do j = jts,jte
 do i = its,ite
    e_dust(i)  = e_dust_c(i,j)
    if ( lptr_seas_aer(1,1,1)>param_first_scalar ) then
        e_seas(i)  = e_seas_c(i,j)
    endif
 enddo
 enddo

 if(debug_on) then
    write(debugString,*) 'max/min of tend_emis         = ',maxval(tend_emis(:,:,its:ite)),minval(tend_emis(:,:,its:ite))
    call mpas_log_write(debugString)
    write(debugString,*) 'max/min of e_dust            = ',maxval(e_dust(its:ite)),minval(e_dust(its:ite))
    call mpas_log_write(debugString)
    if ( lptr_seas_aer(1,1,1)>param_first_scalar ) then
        write(debugString,*) 'max/min of e_seas            = ',maxval(e_seas(its:ite)),minval(e_seas(its:ite))
        call mpas_log_write(debugString)
    endif
 end if
 end subroutine emission_to_MPAS

!=================================================================================================================
 subroutine driver_emission(itimestep,dminfo,configs,mesh,state,sfc_input,diag,diag_physics,diag_chemistry,tend_chemistry,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(in):: diag_physics

 integer,intent(in):: its,ite
 integer,intent(in):: itimestep
!inout arguments:
 type(mpas_pool_type),intent(inout):: tend_chemistry
 type(mpas_pool_type),intent(inout):: diag_chemistry
!local pointers:
 character(len=StrKIND),pointer:: config_dust_opt
 integer,pointer:: config_seas_opt
 real,pointer:: config_dust_factor, config_seas_factor
 logical,pointer:: debug_on

!local variables:
 integer:: i,k,j

 call mpas_pool_get_config(configs,'config_seas_opt'   ,config_seas_opt)
 call mpas_pool_get_config(configs,'config_dust_opt'   ,config_dust_opt)
 call mpas_pool_get_config(configs,'config_dust_factor',config_dust_factor)
 call mpas_pool_get_config(configs,'config_seas_factor',config_seas_factor)
 call mpas_pool_get_config(configs,'config_print_debug_info' ,debug_on )

!copy MPAS arrays to local arrays:
 call emission_from_MPAS(state,sfc_input,diag_chemistry,its,ite)
 if(debug_on) call mpas_log_write('emission from MPAS success')

 call mosaic_addemiss(debug_on,configs,dminfo                                               , &
         slai = slai_c     ,    ust = ust_c  ,                 smois = smois_c              , &
       ivgtyp = ivgtyp_c   ,                                  isltyp = isltyp_c             , &
       dtstep = dt_dyn     ,    u10 = u10_c  ,                   v10 = v10_c                , &
          alt = alt_c      ,   dz8w = dz_c   ,                 xland = xland_c              , & 
       sstemp = sstemp_c   ,   ktau = 1      ,       num_soil_layers = num_soils            , &
        u_phy = u_phy_c    ,  v_phy = v_phy_c,               rho_phy = rho_phy_c            , &
          p8w = pres2_hyd_c,                                    erod = erod_c               , &
         xlat = xlat_c     ,  xlong = xlong_c,                    dx = dx_c                 , &
            g = gravity    ,                        emis_dust_factor = config_dust_factor   , &
                                                    emis_seas_factor = config_seas_factor   , &
    tend_emis = tend_emis_c, e_dust = e_dust_c,               e_seas = e_seas_c             , &
       dust_emiss_active = config_dust_opt,     seasalt_emiss_active = config_seas_opt      , &
          ids = ids, ide = ide, kds = kds, kde = kde, jds = jds, jde = jde                  , &
          ims = ims, ime = ime, kms = kms, kme = kme, jms = jms, jme = jme                  , &
          its = its, ite = ite, kts = kts, kte = kte, jts = jts, jte = jte                    &
        )

 if(debug_on) call mpas_log_write('emission scheme success')
!copy local arrays to MPAS grid:
 call emission_to_MPAS(debug_on,state,tend_chemistry,diag_chemistry,its,ite)
 if(debug_on) call mpas_log_write('emission to MPAS success')

 end subroutine driver_emission

!=================================================================================================================
 end module mpas_atmchem_driver_emission
!=================================================================================================================

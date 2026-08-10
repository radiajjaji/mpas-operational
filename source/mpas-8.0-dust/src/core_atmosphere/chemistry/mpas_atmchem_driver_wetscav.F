!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to compute wet scavenging tendencies for chem(:,:,:,:). 
! * Author: Jiawang Feng from AEMOL USTC
! Last update: jwfeng 03/2020
!
! Contact:
! Chun Zhao, PhD
! Professor 
! University of Science and Technology of China 
! Hefei, Anhui, 230026 
! Email: ChunZhao@ustc.edu.cn
!
! References: 
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

module mpas_atmchem_driver_wetscav
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines
 use mpas_atmchem_vars
 use mpas_atmchem_configs
 use mpas_atmchem_utilities
 use mpas_atmphys_vars

!chemistry schemes:
 use module_mosaic_wetscav

 implicit none
 private
 public:: allocate_wetscav,   &
          deallocate_wetscav, &
          driver_wetscav


REAL, PARAMETER :: mwso4 = 96.00   ! Molecular mass of SO4-- (g/mol)
REAL, PARAMETER :: mwno3 = 62.0    ! Molecular mass of NO3- (g/mol)

CONTAINS

!=================================================================================================================
 subroutine allocate_wetscav
!=================================================================================================================
    if(.not.allocated( qlsink_c) )    allocate( qlsink_c(ims:ime, kms:kme, jms:jme ) )
    if(.not.allocated( ilsink_c) )    allocate( ilsink_c(ims:ime, kms:kme, jms:jme ) )
    if(.not.allocated(  precr_c) )    allocate(  precr_c(ims:ime, kms:kme, jms:jme ) )
    if(.not.allocated(  preci_c) )    allocate(  preci_c(ims:ime, kms:kme, jms:jme ) )
    if(.not.allocated(  precs_c) )    allocate(  precs_c(ims:ime, kms:kme, jms:jme ) )
    if(.not.allocated(  precg_c) )    allocate(  precg_c(ims:ime, kms:kme, jms:jme ) )
 end subroutine allocate_wetscav

!=================================================================================================================
 subroutine deallocate_wetscav
!=================================================================================================================
    if(allocated( qlsink_c) )    deallocate( qlsink_c )
    if(allocated( ilsink_c) )    deallocate( ilsink_c )
    if(allocated(  precr_c) )    deallocate(  precr_c )
    if(allocated(  preci_c) )    deallocate(  preci_c )
    if(allocated(  precs_c) )    deallocate(  precs_c )
    if(allocated(  precg_c) )    deallocate(  precg_c )
 end subroutine deallocate_wetscav

!=================================================================================================================
 subroutine wetscav_from_MPAS(diag_physics, its, ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: diag_physics
 integer,intent(in):: its,ite

 real(kind=RKIND),dimension(:,:),pointer :: qlsink, ilsink
 real(kind=RKIND),dimension(:,:),pointer :: precr, preci, precs, precg

!local variables:
 integer:: i,k,j,h
 character*200, msg

!-----------------------------------------------------------------------------------------------------------------
 call mpas_pool_get_array(diag_physics,'qlsink' ,qlsink )
 call mpas_pool_get_array(diag_physics,'ilsink' ,ilsink )
 call mpas_pool_get_array(diag_physics,'precr'  , precr )
 call mpas_pool_get_array(diag_physics,'preci'  , preci )
 call mpas_pool_get_array(diag_physics,'precs'  , precs )
 call mpas_pool_get_array(diag_physics,'precg'  , precg )

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    qlsink_c (i,k,j) = qlsink(k,i)
    ilsink_c (i,k,j) = ilsink(k,i)
     precr_c (i,k,j) =  precr(k,i)
     preci_c (i,k,j) =  preci(k,i)
     precs_c (i,k,j) =  precs(k,i)
     precg_c (i,k,j) =  precg(k,i)
     if (precr(k,i) .ge. 1.e20 .or. precr(k,i) .le. -1.e20) then
        write(msg,*) 'i,k,j,precr/preci/precs =  ',i,k,j,precr(k,i),preci(k,i),precs(k,i)
        call chemistry_error_fatal(msg)
     end if
 enddo
 enddo
 enddo

 end subroutine wetscav_from_MPAS

!=================================================================================================================
 subroutine wetscav_to_MPAS(debug_on, state, time_lev, tend_chemistry, its, ite)
!=================================================================================================================

 logical:: debug_on
!inout arguments:
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(inout):: tend_chemistry
 integer,intent(in):: its,ite,time_lev

!local variables:
 integer:: i,k,j,h

!local pointers:
 real(kind=RKIND),dimension(:,:,:),pointer :: chem
 real(kind=RKIND),dimension(:,:,:),pointer :: tend_wetls
!-----------------------------------------------------------------------------------------------------------------
 call mpas_pool_get_array(state,'chem',chem,time_lev)
 call mpas_pool_get_array(tend_chemistry,'tend_wetls',tend_wetls)

 do h = 1,num_chem 
 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    !jwfeng:04/2020 NaN.gt.0.0 return .false. 
    if( .not. chem_c(i,k,j,h).gt.0.0_RKIND ) chem_c(i,k,j,h) = 0.0_RKIND
    tend_wetls(h,k,i) = (chem_c(i,k,j,h) - chem(h,k,i))/dt_dyn
    chem(h,k,i)  = chem_c(i,k,j,h)
 enddo
 enddo
 enddo
 enddo

 if(maxval(tend_wetls) .gt. 0.0_RKIND) then
    write(0,*) 'warning! max/min of tend_wetls = ',maxval(tend_wetls),minval(tend_wetls)
 end if

 end subroutine wetscav_to_MPAS

!===========================================================================
subroutine driver_wetscav(state, diag_physics, time_lev, configs, tend_chemistry, its, ite)

 IMPLICIT NONE
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(inout):: tend_chemistry
 integer,intent(in):: its,ite,time_lev

!----------------------------------------------------------------------
!
!  wet removal by grid-resolved precipitation
!  scavenging of cloud-phase aerosols and gases by collection, freezing, ...
!  scavenging of interstitial-phase aerosols by impaction
!  scavenging of gas-phase gases by mass transfer and reaction
!
!  This driver calls subroutines for wet scavenging.

!       Only mosaic for now
!
!----------------------------------------------------------------------

!local vars
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ) :: r2evap, i2evap
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme,1 ) :: gas_aqfrac
   integer :: numgas_aqfrac = 0 !jwfeng: no gas in MPAS-chem for now
   logical,pointer :: debug_on

   call mpas_pool_get_config(configs,'config_print_debug_info' ,debug_on )
   r2evap = 0.0
   i2evap = 0.0
   gas_aqfrac = 0.0

!copy MPAS arrays to local arrays:
 call wetscav_from_MPAS(diag_physics, its, ite)
 if(debug_on) call chemistry_message('wetscav from MPAS success')

    if(debug_on) CALL chemistry_message('WETSCAV_DRIVER CALLING MOSAIC_WETSCAV_DRIVER')
    call wetscav_cbmz_mosaic (  debug_on = debug_on,                                    &
         id = 1,ktau = 1,dtstep = dt_dyn,ktauc = 1,config_flags = config_flags,         &
         dtstepc = dt_dyn, alt = alt_c, t_phy = t_phy_c, p8w = pres2_hyd_c, t8w = t2_c, &
         p_phy = pres_hyd_c,chem = chem_c,rho_phy = rho_phy_c, &
         cldfra = cldfra_c, qlsink = qlsink_c, ilsink = ilsink_c , &
         precr = precr_c,preci = preci_c,precs = precs_c,precg = precg_c , &
!r2evap,i2evap are not important,not used further.
         r2evap = r2evap,i2evap = i2evap,    & 
         gas_aqfrac = gas_aqfrac, numgas_aqfrac = numgas_aqfrac,             &
         qv = moist_c(ims:ime,kms:kme,jms:jme,P_QV),                         &
         qc = moist_c(ims:ime,kms:kme,jms:jme,P_QC),                         &
         qi = moist_c(ims:ime,kms:kme,jms:jme,P_QI),                         &
         qs = moist_c(ims:ime,kms:kme,jms:jme,P_QS),                         &
         ids = ids, ide = ide, kds = kds, kde = kde, jds = jds, jde = jde,   &
         ims = ims, ime = ime, kms = kms, kme = kme, jms = jms, jme = jme,   &
         its = its, ite = ite, kts = kts, kte = kte, jts = jts, jte = jte    )

    ! Calculate wet deposition of the species by accumulating their (negative)
    ! column change. It is assumed that all precipitation reaches the ground.

 if(debug_on) call chemistry_message( 'wetscav scheme success')
!copy local arrays to MPAS grid:
 call wetscav_to_MPAS(debug_on, state, time_lev, tend_chemistry, its, ite)
 if(debug_on) call chemistry_message( 'wetscav to MPAS success')
end subroutine driver_wetscav


END MODULE mpas_atmchem_driver_wetscav

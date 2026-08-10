!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Interface between MPAS and chemistry schemes for all common varaiables in
! different chemistry processes.
! * Author: Jiawang Feng from AEMOL USTC
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
!     MPAS-v5.2:src/core_atmosphere/physics/mpas_atmphys_interface.F
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

MODULE mpas_atmchem_interface
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines
 use mpas_atmchem_vars
 use mpas_atmchem_constants
 use mpas_atmphys_vars
 use mpas_atmchem_utilities

 implicit none
 private
 public:: allocate_forall_chemistry,   &
          deallocate_forall_chemistry, &
          chemistry_from_MPAS

! MPAS interface for all chemistry.

! subroutines:
! ---------------------------------------
! allocate_forall_chemistry  : allocate local arrays for all chemistry.
! deallocate_forall_chemistry: deallocate local arrays for all chemistry.
! chemistry_from_MPAS        : initialize local arrays.
! chemistry_to_MPAS          : copy local arrays to MPAS arrays.
! --------------------------------------

 contains

!=================================================================================================================
 subroutine allocate_forall_chemistry
!=================================================================================================================
! Wind
    if(.not.allocated(        ust_c) )    allocate(        ust_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(        u10_c) )    allocate(        u10_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(        v10_c) )    allocate(        v10_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(          w_c) )    allocate(          w_c(ims:ime, kms:kme,   jms:jme            ) )
    if(.not.allocated(      u_phy_c) )    allocate(      u_phy_c(ims:ime, kms:kme,   jms:jme            ) )
    if(.not.allocated(      v_phy_c) )    allocate(      v_phy_c(ims:ime, kms:kme,   jms:jme            ) )
! Rho                   
    if(.not.allocated(        alt_c) )    allocate(        alt_c(ims:ime, kms:kme,   jms:jme            ) )
    if(.not.allocated(    rho_phy_c) )    allocate(    rho_phy_c(ims:ime, kms:kme,   jms:jme            ) )
! Temperature           
    if(.not.allocated(      t_phy_c) )    allocate(      t_phy_c(ims:ime, kms:kme,   jms:jme            ) )
    if(.not.allocated(         th_c) )    allocate(         th_c(ims:ime, kms:kme,   jms:jme            ) )
    if(.not.allocated(        t8w_c) )    allocate(        t8w_c(ims:ime, kms:kme,   jms:jme            ) )
    if(.not.allocated(         t2_c) )    allocate(         t2_c(ims:ime, kms:kme,   jms:jme            ) )
! Pressure              
    if(.not.allocated(       pres_c) )    allocate(       pres_c(ims:ime, kms:kme,   jms:jme            ) )
    if(.not.allocated(      pres2_c) )    allocate(      pres2_c(ims:ime, kms:kme,   jms:jme            ) )
    if(.not.allocated(  pres2_hyd_c) )    allocate(  pres2_hyd_c(ims:ime, kms:kme,   jms:jme            ) )
    if(.not.allocated( pres2_hydd_c) )    allocate( pres2_hydd_c(ims:ime, kms:kme,   jms:jme            ) )
    if(.not.allocated(   pres_hyd_c) )    allocate(   pres_hyd_c(ims:ime, kms:kme,   jms:jme            ) )
    if(.not.allocated(  pres_hydd_c) )    allocate(  pres_hydd_c(ims:ime, kms:kme,   jms:jme            ) )
    if(.not.allocated(       psfc_c) )    allocate(       psfc_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(   psfc_hyd_c) )    allocate(   psfc_hyd_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(  psfc_hydd_c) )    allocate(  psfc_hydd_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(    znu_hyd_c) )    allocate(    znu_hyd_c(ims:ime, kms:kme,   jms:jme            ) )
! Cell Information      
    if(.not.allocated(         dz_c) )    allocate(         dz_c(ims:ime, kms:kme,   jms:jme            ) )
    if(.not.allocated(       xlat_c) )    allocate(       xlat_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(      xlong_c) )    allocate(      xlong_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(         dx_c) )    allocate(         dx_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(          z_c) )    allocate(          z_c(ims:ime, kms:kme,   jms:jme            ) )
    if(.not.allocated(       hpbl_c) )    allocate(       hpbl_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(        znt_c) )    allocate(        znt_c(ims:ime,            jms:jme            ) )
! Field                 
    if(.not.allocated(     cldfra_c) )    allocate(     cldfra_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated( cldfra_old_c) )    allocate( cldfra_old_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(     exch_h_c) )    allocate(     exch_h_c(ims:ime, kms:kme,   jms:jme            ) )
    if(.not.allocated(     raincv_c) )    allocate(     raincv_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(        tsk_c) )    allocate(        tsk_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(        gsw_c) )    allocate(        gsw_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(     vegfra_c) )    allocate(     vegfra_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(       rmol_c) )    allocate(       rmol_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(      snowh_c) )    allocate(      snowh_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(        hfx_c) )    allocate(        hfx_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(      xland_c) )    allocate(      xland_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(       slai_c) )    allocate(       slai_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(     ivgtyp_c) )    allocate(     ivgtyp_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(     isltyp_c) )    allocate(     isltyp_c(ims:ime,            jms:jme            ) )
    if(.not.allocated(         qv_c) )    allocate(         qv_c(ims:ime, kms:kme,   jms:jme            ) )
    if(.not.allocated(      moist_c) )    allocate(      moist_c(ims:ime, kms:kme,   jms:jme, num_moist ) )
    if(.not.allocated(    scalars_c) )    allocate(    scalars_c(ims:ime, kms:kme,   jms:jme, num_scalars ) )
    if(.not.allocated(       chem_c) )    allocate(       chem_c(ims:ime, kms:kme,   jms:jme, num_chem  ) )


 end subroutine allocate_forall_chemistry

!=================================================================================================================
 subroutine deallocate_forall_chemistry
!=================================================================================================================
    if(allocated(        ust_c) )   deallocate(        ust_c)
    if(allocated(        u10_c) )   deallocate(        u10_c)
    if(allocated(        v10_c) )   deallocate(        v10_c)
    if(allocated(          w_c) )   deallocate(          w_c)
    if(allocated(      u_phy_c) )   deallocate(      u_phy_c)
    if(allocated(      v_phy_c) )   deallocate(      v_phy_c)
! Rho              
    if(allocated(        alt_c) )   deallocate(        alt_c)
    if(allocated(    rho_phy_c) )   deallocate(    rho_phy_c)
! Temperature           
    if(allocated(      t_phy_c) )   deallocate(      t_phy_c)
    if(allocated(         th_c) )   deallocate(         th_c)
    if(allocated(        t8w_c) )   deallocate(        t8w_c)
! Pressure              
    if(allocated(       pres_c) )   deallocate(       pres_c)
    if(allocated(      pres2_c) )   deallocate(      pres2_c)
    if(allocated(  pres2_hyd_c) )   deallocate(  pres2_hyd_c)
    if(allocated( pres2_hydd_c) )   deallocate( pres2_hydd_c)
    if(allocated(   pres_hyd_c) )   deallocate(   pres_hyd_c)
    if(allocated(  pres_hydd_c) )   deallocate(  pres_hydd_c)
    if(allocated(       psfc_c) )   deallocate(       psfc_c)
    if(allocated(   psfc_hyd_c) )   deallocate(   psfc_hyd_c)
    if(allocated(  psfc_hydd_c) )   deallocate(  psfc_hydd_c)
    if(allocated(    znu_hyd_c) )   deallocate(    znu_hyd_c)
! Cell Information      
    if(allocated(         dz_c) )   deallocate(         dz_c)
    if(allocated(       xlat_c) )   deallocate(       xlat_c)
    if(allocated(      xlong_c) )   deallocate(      xlong_c)
    if(allocated(         dx_c) )   deallocate(         dx_c)
    if(allocated(          z_c) )   deallocate(          z_c)
    if(allocated(       hpbl_c) )   deallocate(       hpbl_c)
    if(allocated(        znt_c) )   deallocate(        znt_c)
! Field            
    if(allocated(     cldfra_c) )   deallocate(     cldfra_c)
    if(allocated( cldfra_old_c) )   deallocate( cldfra_old_c)
    if(allocated(     exch_h_c) )   deallocate(     exch_h_c)
    if(allocated(     raincv_c) )   deallocate(     raincv_c)
    if(allocated(        tsk_c) )   deallocate(        tsk_c)
    if(allocated(        gsw_c) )   deallocate(        gsw_c)
    if(allocated(     vegfra_c) )   deallocate(     vegfra_c)
    if(allocated(       rmol_c) )   deallocate(       rmol_c)
    if(allocated(      snowh_c) )   deallocate(      snowh_c)
    if(allocated(        hfx_c) )   deallocate(        hfx_c)
    if(allocated(      xland_c) )   deallocate(      xland_c)
    if(allocated(       slai_c) )   deallocate(       slai_c)
    if(allocated(     ivgtyp_c) )   deallocate(     ivgtyp_c)
    if(allocated(     isltyp_c) )   deallocate(     isltyp_c)
    if(allocated(         qv_c) )   deallocate(         qv_c)
    if(allocated(      moist_c) )   deallocate(      moist_c)
    if(allocated(    scalars_c) )   deallocate(    scalars_c)
    if(allocated(       chem_c) )   deallocate(       chem_c)

 end subroutine deallocate_forall_chemistry

!=================================================================================================================
 subroutine chemistry_from_MPAS(mesh,state,time_lev,sfc_input,diag,diag_physics,its,ite)
!=================================================================================================================

 use mpas_atmchem_utilities, only: chemistry_message
!input variables:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: diag

 integer,intent(in):: its,ite,time_lev

 integer,pointer:: index_dust_a01,index_qv
 integer,dimension(:),pointer:: ivgtyp,isltyp

 real(kind=RKIND),dimension(:)    ,pointer :: areaCell,xlat,xlong,u10,v10,ust,xland,dms_0
 real(kind=RKIND),dimension(:)    ,pointer :: rdzw,surface_pressure
 real(kind=RKIND),dimension(:)    ,pointer :: raincv,skintemp,gsw,vegfra,rmol
 real(kind=RKIND),dimension(:)    ,pointer :: snowh,hfx,hpbl,znt
 real(kind=RKIND),dimension(:,:)  ,pointer :: pressure,smois,zgrid,erod
 real(kind=RKIND),dimension(:,:)  ,pointer :: u_phy,v_phy
 real(kind=RKIND),dimension(:,:)  ,pointer :: rho_zz,zz,w,exner,pressure_p,pressure_b
 real(kind=RKIND),dimension(:,:)  ,pointer :: qv,theta_m,kzh
!real(kind=RKIND),dimension(:,:)  ,pointer :: cldfrac_old
 real(kind=RKIND),dimension(:,:)  ,pointer :: cldfrac
 real(kind=RKIND),dimension(:,:,:),pointer :: scalars,moist,chem


!local variables:
 integer :: i,k,j,h
 real(kind=RKIND) :: z0,z1,z2,w1,w2,fzm,fzp
 real(kind=RKIND) :: rho_a,rho1,rho2,tem1,tem2


!-----------------------------------------------------------------------------------------------------------------
!Jiawang Feng:
!Gocart use a constantly dxy as area.Can we use array areaCell to reach the same results?
 call mpas_pool_get_array(mesh,'areaCell',areaCell)
 call mpas_pool_get_array(mesh,'latCell' ,xlat    )
 call mpas_pool_get_array(mesh,'lonCell' ,xlong   )
 call mpas_pool_get_array(mesh,'zgrid'   ,zgrid   )
 call mpas_pool_get_array(mesh,'zz'      ,zz      )
 call mpas_pool_get_array(mesh,'rdzw'    ,rdzw    )

 call mpas_pool_get_array(diag_physics,'u10'         ,u10           )
 call mpas_pool_get_array(diag_physics,'v10'         ,v10           )
 call mpas_pool_get_array(diag_physics,'ust'         ,ust           )
 call mpas_pool_get_array(diag_physics,'gsw'         ,gsw           )
 call mpas_pool_get_array(diag_physics,'hfx'         ,hfx           )
 call mpas_pool_get_array(diag_physics,'znt'         ,znt           )
 call mpas_pool_get_array(diag_physics,'rmol'        ,rmol          )
 call mpas_pool_get_array(diag_physics,'hpbl'        ,hpbl          )
 call mpas_pool_get_array(diag_physics,'raincv'      ,raincv        )
 call mpas_pool_get_array(diag_physics,'kzh'         ,kzh           )
 call mpas_pool_get_array(diag_physics,'cldfrac'     ,cldfrac       )
!call mpas_pool_get_array(diag_physics,'cldfrac_old' ,cldfrac_old   )

 call mpas_pool_get_array(diag,'surface_pressure'      ,surface_pressure)
 call mpas_pool_get_array(diag,'exner'                 ,exner           )
 call mpas_pool_get_array(diag,'uReconstructZonal'     ,u_phy           )
 call mpas_pool_get_array(diag,'uReconstructMeridional',v_phy           )
 call mpas_pool_get_array(diag,'pressure_base'         ,pressure_b      )
 call mpas_pool_get_array(diag,'pressure_p'            ,pressure_p      )

 call mpas_pool_get_array(state,'rho_zz' ,rho_zz ,time_lev)
 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)
 call mpas_pool_get_array(state,'w      ',w      ,time_lev)
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 call mpas_pool_get_array(state,'chem'   ,chem   ,time_lev)

 call mpas_pool_get_array(sfc_input,'xland'    ,xland    )
 call mpas_pool_get_array(sfc_input,'ivgtyp'   ,ivgtyp   )
 call mpas_pool_get_array(sfc_input,'isltyp'   ,isltyp   )
 call mpas_pool_get_array(sfc_input,'skintemp' ,skintemp )
 call mpas_pool_get_array(sfc_input,'vegfra'   ,vegfra   )
 call mpas_pool_get_array(sfc_input,'snowh'    ,snowh    )

 call mpas_pool_get_dimension(state,'index_dust_a01',index_dust_a01)
 call mpas_pool_get_dimension(state,'index_qv      ',index_qv      )
 qv    => scalars(index_qv,:,:)
 moist => scalars(1:num_moist,:,:)

 if(associated(areaCell         )  .and. associated(xlat            )  .and. &
    associated(xlong            )  .and. associated(zgrid           )  .and. &
    associated(zz               )  .and. associated(rdzw            )  .and. &
    associated(u10              )  .and. associated(v10             )  .and. &
    associated(ust              )  .and. associated(gsw             )  .and. &
    associated(hfx              )  .and. associated(znt             )  .and. &   
    associated(rmol             )  .and. associated(hpbl            )  .and. &  
    associated(raincv           )  .and. associated(surface_pressure)  .and. & 
    associated(exner            )  .and. associated(u_phy           )  .and. &
    associated(v_phy            )  .and. associated(pressure_b      )  .and. &
    associated(pressure_p       )  .and. associated(rho_zz          )  .and. &
    associated(theta_m          )  .and. associated(w               )  .and. &
    associated(scalars          )  .and. associated(chem            )  .and. &
    associated(xland            )  .and. associated(ivgtyp          )  .and. &
    associated(isltyp           )  .and. associated(skintemp        )  .and. &
    associated(vegfra           )  .and. associated(snowh           )  .and. &
!   associated(cldfrac          )  .and. associated(cldfrac_old     )  .and. &
    associated(kzh              )                                            &
    ) then
 else
    call mpas_log_write( '    areaCell       = $l', logicArgs=(/associated(areaCell)     /) )
    call mpas_log_write( '    xlat       = $l', logicArgs=(/associated(xlat)     /) )
    call mpas_log_write( '    xlong       = $l', logicArgs=(/associated(xlong)     /) )
    call mpas_log_write( '    zgrid       = $l', logicArgs=(/associated(zgrid)     /) )
    call mpas_log_write( '    zz       = $l', logicArgs=(/associated(zz)     /) )
    call mpas_log_write( '    rdzw       = $l', logicArgs=(/associated(rdzw)     /) )
    call mpas_log_write( '    u10       = $l', logicArgs=(/associated(u10)     /) )
    call mpas_log_write( '    v10       = $l', logicArgs=(/associated(v10)     /) )
    call mpas_log_write( '    ust       = $l', logicArgs=(/associated(ust)     /) )
    call mpas_log_write( '    gsw       = $l', logicArgs=(/associated(gsw)     /) )
    call mpas_log_write( '    hfx       = $l', logicArgs=(/associated(hfx)     /) )
    call mpas_log_write( '    znt       = $l', logicArgs=(/associated(znt)     /) )
    call mpas_log_write( '    rmol       = $l', logicArgs=(/associated(rmol)     /) )
    call mpas_log_write( '    hpbl       = $l', logicArgs=(/associated(hpbl)     /) )
    call mpas_log_write( '    raincv       = $l', logicArgs=(/associated(raincv)     /) )
    call mpas_log_write( '    surface_pressure       = $l', logicArgs=(/associated(surface_pressure)     /) )
    call mpas_log_write( '    exner       = $l', logicArgs=(/associated(exner)     /) )
    call mpas_log_write( '    u_phy       = $l', logicArgs=(/associated(u_phy)     /) )
    call mpas_log_write( '    v_phy       = $l', logicArgs=(/associated(v_phy)     /) )
    call mpas_log_write( '    pressure_b       = $l', logicArgs=(/associated(pressure_b)     /) )
    call mpas_log_write( '    pressure_p       = $l', logicArgs=(/associated(pressure_p)     /) )
    call mpas_log_write( '    rho_zz       = $l', logicArgs=(/associated(rho_zz)     /) )
    call mpas_log_write( '    theta_m       = $l', logicArgs=(/associated(theta_m)     /) )
    call mpas_log_write( '    w       = $l', logicArgs=(/associated(w)     /) )
    call mpas_log_write( '    scalars       = $l', logicArgs=(/associated(scalars)     /) )
    call mpas_log_write( '    chem       = $l', logicArgs=(/associated(chem)     /) )
    call mpas_log_write( '    xland       = $l', logicArgs=(/associated(xland)     /) )
    call mpas_log_write( '    ivgtyp       = $l', logicArgs=(/associated(ivgtyp)     /) )
    call mpas_log_write( '    isltyp       = $l', logicArgs=(/associated(isltyp)     /) )
    call mpas_log_write( '    skintemp       = $l', logicArgs=(/associated(skintemp)     /) )
    call mpas_log_write( '    vegfra       = $l', logicArgs=(/associated(vegfra)     /) )
    call mpas_log_write( '    snowh       = $l', logicArgs=(/associated(snowh)     /) )
    call mpas_log_write( '    cldfrac       = $l', logicArgs=(/associated(cldfrac)     /) )
!   call mpas_log_write( '    cldfrac_old       = $l', logicArgs=(/associated(cldfrac_old)     /) )
    call mpas_log_write( '    kzh       = $l', logicArgs=(/associated(kzh)     /) )
    call chemistry_error_fatal('Something wrong in chemistry interface:some arguments are not correctly associated!!!!!')
 endif

 do h = 1,num_moist
 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    moist_c(i,k,j,h) = max(0.,moist(h,k,i))
 enddo
 enddo
 enddo
 enddo

 do h = 1,num_scalars
 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    scalars_c(i,k,j,h) = max(0.,scalars(h,k,i))
 enddo
 enddo
 enddo
 enddo


 do h = 1,num_chem
 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    chem_c(i,k,j,h) = max(epsilc,chem(h,k,i))
 enddo
 enddo
 enddo
 enddo

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
         qv_c(i,k,j) = max(0.,qv(k,i))
    rho_phy_c(i,k,j) = zz(k,i) * rho_zz(k,i)
    rho_phy_c(i,k,j) = rho_phy_c(i,k,j)*(1._RKIND + qv_c(i,k,j))
         th_c(i,k,j) = theta_m(k,i)/(1._RKIND + R_v/R_d * qv_c(i,k,j))
      t_phy_c(i,k,j) = th_c(i,k,j)*exner(k,i)
        alt_c(i,k,j) = 1.0/rho_phy_c(i,k,j)
       pres_c(i,k,j) = pressure_p(k,i) + pressure_b(k,i)
      u_phy_c(i,k,j) = u_phy(k,i)
      v_phy_c(i,k,j) = v_phy(k,i)
         dz_c(i,k,j) = zgrid(k+1,i)-zgrid(k,i)
     cldfra_c(i,k,j) = cldfrac(k,i)
!cldfra_old_c(i,k,j) = cldfrac_old(k,i)
 enddo
 enddo
 enddo

!  exch_h_c = 0.0
 do j = jts,jte
 do i = its,ite
 do k = kts,kte
    exch_h_c(i,k,j) = kzh(k,i)
 enddo
 enddo
 enddo

 !calculation of the surface pressure using hydrostatic assumption down to the surface::
 do j = jts,jte
 do i = its,ite
                   tem1 = zgrid(2,i)-zgrid(1,i)
                   tem2 = zgrid(3,i)-zgrid(2,i)
                   rho1 = rho_zz(1,i) * zz(1,i) * (1. + qv_c(i,1,j))
                   rho2 = rho_zz(2,i) * zz(2,i) * (1. + qv_c(i,2,j))
    surface_pressure(i) = 0.5*gravity*(zgrid(2,i)-zgrid(1,i)) &
                             *(rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
    surface_pressure(i) = surface_pressure(i) + pressure_p(1,i) + pressure_b(1,i)
 enddo
 enddo

 !arrays located at w points:
 do j = jts, jte
 do k = kts,kte+1
 do i = its,ite
    w_c(i,k,j) = w(k,i)
    z_c(i,k,j) = zgrid(k,i)
 enddo
 enddo
 enddo

!interpolation of pressure and temperature from theta points to w points:
 do j = jts,jte
 do k = kts+1,kte
 do i = its,ite
              tem1 = 1./(zgrid(k+1,i)-zgrid(k-1,i))
               fzm = (zgrid(k,i)-zgrid(k-1,i)) * tem1
               fzp = (zgrid(k+1,i)-zgrid(k,i)) * tem1
       t2_c(i,k,j) = fzm*t_phy_c(i,k,j) + fzp*t_phy_c(i,k-1,j)
    pres2_c(i,k,j) = fzm*pres_c(i,k,j)  + fzp*pres_c(i,k-1,j)
 enddo
 enddo
 enddo

!interpolation of pressure and temperature from theta points to the top-of-the-model: follows
!the calculation of the top-of-the-model pressure and temperature in WRF (subroutine phy_prep
!in ./dyn_em/module_big_step_utilities.F):
 k = kte+1
 do j = jts,jte
 do i = its,ite
                z0 = zgrid(k,i)
                z1 = 0.5*(zgrid(k,i)+zgrid(k-1,i)) 
                z2 = 0.5*(zgrid(k-1,i)+zgrid(k-2,i))
                w1 = (z0-z2)/(z1-z2)
                w2 = 1.-w1
       t2_c(i,k,j) = w1*t_phy_c(i,k-1,j) + w2*t_phy_c(i,k-2,j)
    !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
    pres2_c(i,k,j) = exp(w1*log(pres_c(i,k-1,j))+w2*log(pres_c(i,k-2,j)))
 enddo
 enddo
 
 !ldf (2012-06-22): recalculates the pressure at the surface as an extrapolation of the
 !pressures in the 2 layers above the surface, as was originally done:
 k = kts
 do j = jts,jte
 do i = its,ite
                z0 = zgrid(k,i)
                z1 = 0.5*(zgrid(k,i)+zgrid(k+1,i)) 
                z2 = 0.5*(zgrid(k+1,i)+zgrid(k+2,i))
                w1 = (z0-z2)/(z1-z2)
                w2 = 1.-w1
       t2_c(i,k,j) = w1*t_phy_c(i,k,j) + w2*t_phy_c(i,k+1,j)
    pres2_c(i,k,j) = w1* pres_c(i,k,j) + w2* pres_c(i,k+1,j)
       psfc_c(i,j) = surface_pressure(i)
 enddo
 enddo

 !calculation of the hydrostatic pressure:
 do j = jts,jte
 do i = its,ite
    !pressure at w-points:
    k = kte+1
    pres2_hyd_c(i,k,j)  = pres2_c(i,k,j)
    pres2_hydd_c(i,k,j) = pres2_c(i,k,j)
    do k = kte,1,-1
       rho_a = rho_phy_c(i,k,j) / (1.+qv_c(i,k,j))
       pres2_hyd_c(i,k,j)  = pres2_hyd_c(i,k+1,j)  + gravity*rho_phy_c(i,k,j)*dz_c(i,k,j)
       pres2_hydd_c(i,k,j) = pres2_hydd_c(i,k+1,j) + gravity*rho_a*dz_c(i,k,j)
    enddo
    !pressure at theta-points:
    do k = kte,1,-1
       pres_hyd_c(i,k,j)  = 0.5*(pres2_hyd_c(i,k+1,j)+pres2_hyd_c(i,k,j))
       pres_hydd_c(i,k,j) = 0.5*(pres2_hydd_c(i,k+1,j)+pres2_hydd_c(i,k,j))
    enddo
    !surface pressure:
     psfc_hyd_c(i,j) = pres2_hyd_c(i,1,j)
    psfc_hydd_c(i,j) = pres2_hydd_c(i,1,j)
    !znu:
    do k = kte,1,-1
       znu_hyd_c(i,k,j) = pres_hyd_c(i,k,j) / psfc_hyd_c(i,j) 
    enddo
 enddo
 enddo

 do j = jts,jte
 do i = its,ite
       u10_c(i,j) = u10(i)
       v10_c(i,j) = v10(i)
       ust_c(i,j) = ust(i)
    ivgtyp_c(i,j) = ivgtyp(i)
    isltyp_c(i,j) = isltyp(i)
     xland_c(i,j) = xland(i)
        dx_c(i,j) = sqrt(areaCell(i))
      xlat_c(i,j) = xlat(i)*57.2957795
     xlong_c(i,j) = xlong(i)*57.2957795
    raincv_c(i,j) = raincv(i)
     snowh_c(i,j) = snowh(i)
    vegfra_c(i,j) = vegfra(i)
       tsk_c(i,j) = skintemp(i)
       gsw_c(i,j) = gsw(i)
       hfx_c(i,j) = hfx(i)
      rmol_c(i,j) = rmol(i)
      hpbl_c(i,j) = hpbl(i)
       znt_c(i,j) = znt(i)
 enddo
 enddo

 end subroutine chemistry_from_MPAS

END MODULE mpas_atmchem_interface

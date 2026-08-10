!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to compute dry deposition tendencies for chem(:,:,:,:). 
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
!
! Terms of Use:
!  1) Users are requested to consult the primary author prior to
!     modifying this module or incorporating it or its submodules in
!     another code. This is meant to ensure that the any linkages and/or
!     assumptions will not adversely affect the operation of this module.
!  2) The source code in this module is intended for research and
!     educational purposes. Users are requested to contact the primary
!     author regarding the use of the code for any commercial application.
!  3) Users preparing publications resulting from the usage of this code
!     are requested to cite one or more of the references below
!     (depending on the application) for proper acknowledgement.
!
! References: 
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

!===========================================================================
!IN THIS SUBROUTINE, IT INCLUDES FOLLOWING PARTS:
!  + Dry deposition calculation. 
!  + Tendency calculation for budget analysis 
!                                         ---jwfeng   12/2019
!=============================================================================

MODULE mpas_atmchem_driver_drydep
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines
 use mpas_atmchem_vars
 use mpas_atmchem_configs
 use mpas_atmchem_utilities
 use mpas_atmphys_vars

 use module_dep_simple
 use module_gocart_drydep
 use module_aer_drydep
 use module_vertmx_wrf
 use module_mosaic_settling
 use module_data_binaer_asect
 use module_mixactivate_wrappers

 implicit none
 private
 public:: allocate_drydep,   &
          deallocate_drydep, &
          driver_drydep

! subroutines in mpas_atmchem_driver_drydep:
! ---------------------------------------
! allocate_drydep  : allocate local arrays for drydep.
! deallocate_drydep: deallocate local arrays for emission.
! driver_drydep    : main driver (called from subroutine chemistry_driver).

! drydep_from_MPAS : initialize local arrays.
! drydep_to_MPAS   : copy local arrays to MPAS arrays.

 contains
!=================================================================================================================
 subroutine allocate_drydep
!=================================================================================================================
! mixactivate
    if(.not.allocated(      f_ice_c) )    allocate(      f_ice_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(     f_rain_c) )    allocate(     f_rain_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(    tke_pbl_c) )    allocate(    tke_pbl_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(   relhum3d_c) )    allocate(   relhum3d_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(     nihf3d_c) )    allocate(     nihf3d_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(    niimm3d_c) )    allocate(    niimm3d_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(    nidep3d_c) )    allocate(    nidep3d_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(     wsub3d_c) )    allocate(     wsub3d_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(       ccn1_c) )    allocate(       ccn1_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(       ccn2_c) )    allocate(       ccn2_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(       ccn3_c) )    allocate(       ccn3_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(       ccn4_c) )    allocate(       ccn4_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(       ccn5_c) )    allocate(       ccn5_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(       ccn6_c) )    allocate(       ccn6_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(    nsource_c) )    allocate(    nsource_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(     rndst1_c) )    allocate(     rndst1_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(     rndst2_c) )    allocate(     rndst2_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(     rndst3_c) )    allocate(     rndst3_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(     rndst4_c) )    allocate(     rndst4_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(     nacon1_c) )    allocate(     nacon1_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(     nacon2_c) )    allocate(     nacon2_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(     nacon3_c) )    allocate(     nacon3_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated(     nacon4_c) )    allocate(     nacon4_c(ims:ime, kms:kme, jms:jme            ) )

    if(.not.allocated(      ddvel_c) )    allocate(      ddvel_c(ims:ime,          jms:jme, num_chem  ) )
 end subroutine allocate_drydep

!=================================================================================================================
 subroutine deallocate_drydep
!=================================================================================================================
! mixactivate
    if(allocated(      f_ice_c) )    deallocate(      f_ice_c)
    if(allocated(     f_rain_c) )    deallocate(     f_rain_c)
    if(allocated(    tke_pbl_c) )    deallocate(    tke_pbl_c)
    if(allocated(   relhum3d_c) )    deallocate(   relhum3d_c)
    if(allocated(     nihf3d_c) )    deallocate(     nihf3d_c)
    if(allocated(    niimm3d_c) )    deallocate(    niimm3d_c)
    if(allocated(    nidep3d_c) )    deallocate(    nidep3d_c)
    if(allocated(     wsub3d_c) )    deallocate(     wsub3d_c)
    if(allocated(       ccn1_c) )    deallocate(       ccn1_c)
    if(allocated(       ccn2_c) )    deallocate(       ccn2_c)
    if(allocated(       ccn3_c) )    deallocate(       ccn3_c)
    if(allocated(       ccn4_c) )    deallocate(       ccn4_c)
    if(allocated(       ccn5_c) )    deallocate(       ccn5_c)
    if(allocated(       ccn6_c) )    deallocate(       ccn6_c)
    if(allocated(    nsource_c) )    deallocate(    nsource_c)
    if(allocated(     rndst1_c) )    deallocate(     rndst1_c)
    if(allocated(     rndst2_c) )    deallocate(     rndst2_c)
    if(allocated(     rndst3_c) )    deallocate(     rndst3_c)
    if(allocated(     rndst4_c) )    deallocate(     rndst4_c)
    if(allocated(     nacon1_c) )    deallocate(     nacon1_c)
    if(allocated(     nacon2_c) )    deallocate(     nacon2_c)
    if(allocated(     nacon3_c) )    deallocate(     nacon3_c)
    if(allocated(     nacon4_c) )    deallocate(     nacon4_c)

    if(allocated(      ddvel_c))     deallocate(      ddvel_c)
 end subroutine deallocate_drydep

!=================================================================================================================
 subroutine drydep_from_MPAS(state,diag,diag_physics,its,ite)    !!may need src_input
!=================================================================================================================
    !input arguments:
     type(mpas_pool_type),intent(in):: state
     type(mpas_pool_type),intent(in):: diag
     type(mpas_pool_type),intent(in):: diag_physics
     integer,intent(in):: its,ite
    !local variables:
     integer:: i,k,j,h
    !local pointers:
     real(kind=RKIND),dimension(:,:)  ,pointer :: tke_pbl
     real(kind=RKIND),dimension(:,:)  ,pointer ::  relhum
     real(kind=RKIND),dimension(:,:)  ,pointer ::  nihf3d
     real(kind=RKIND),dimension(:,:)  ,pointer :: niimm3d
     real(kind=RKIND),dimension(:,:)  ,pointer :: nidep3d
     real(kind=RKIND),dimension(:,:)  ,pointer ::  wsub3d
     real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn1
     real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn2
     real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn3
     real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn4
     real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn5
     real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn6
     real(kind=RKIND),dimension(:,:)  ,pointer :: nsource
     real(kind=RKIND),dimension(:,:)  ,pointer ::  rndst1
     real(kind=RKIND),dimension(:,:)  ,pointer ::  rndst2
     real(kind=RKIND),dimension(:,:)  ,pointer ::  rndst3
     real(kind=RKIND),dimension(:,:)  ,pointer ::  rndst4
     real(kind=RKIND),dimension(:,:)  ,pointer ::  nacon1
     real(kind=RKIND),dimension(:,:)  ,pointer ::  nacon2
     real(kind=RKIND),dimension(:,:)  ,pointer ::  nacon3
     real(kind=RKIND),dimension(:,:)  ,pointer ::  nacon4
!-----------------------------------------------------------------------------------------------------------------
     call mpas_pool_get_array(diag        ,'relhum'     ,relhum        )
     call mpas_pool_get_array(diag        ,'nihf3d'     ,nihf3d        )
     call mpas_pool_get_array(diag        ,'niimm3d'    ,niimm3d       )
     call mpas_pool_get_array(diag        ,'nidep3d'    ,nidep3d       )
     call mpas_pool_get_array(diag        ,'wsub3d'     ,wsub3d        )
     call mpas_pool_get_array(diag        ,'ccn1'       ,ccn1          )
     call mpas_pool_get_array(diag        ,'ccn2'       ,ccn2          )
     call mpas_pool_get_array(diag        ,'ccn3'       ,ccn3          )
     call mpas_pool_get_array(diag        ,'ccn4'       ,ccn4          )
     call mpas_pool_get_array(diag        ,'ccn5'       ,ccn5          )
     call mpas_pool_get_array(diag        ,'ccn6'       ,ccn6          )
     call mpas_pool_get_array(diag        ,'nsource'    ,nsource       )
     call mpas_pool_get_array(diag        ,'rndst1'     ,rndst1        )
     call mpas_pool_get_array(diag        ,'rndst2'     ,rndst2        )
     call mpas_pool_get_array(diag        ,'rndst3'     ,rndst3        )
     call mpas_pool_get_array(diag        ,'rndst4'     ,rndst4        )
     call mpas_pool_get_array(diag        ,'nacon1'     ,nacon1        )
     call mpas_pool_get_array(diag        ,'nacon2'     ,nacon2        )
     call mpas_pool_get_array(diag        ,'nacon3'     ,nacon3        )
     call mpas_pool_get_array(diag        ,'nacon4'     ,nacon4        )

     do h = 1,num_chem
     do j = jts,jte
     do i = its,ite
         ddvel_c(i,j,h)                  = 0._RKIND
     enddo
     enddo
     enddo

     do j = jts,jte
     do k = kts,kte
     do i = its,ite
          tke_pbl_c(i,k,j)    = 0._RKIND ! Don't need tke_pbl
         relhum3d_c(i,k,j)    =  relhum(k,i)
           nihf3d_c(i,k,j)    =  nihf3d(k,i)
          niimm3d_c(i,k,j)    = niimm3d(k,i)
          nidep3d_c(i,k,j)    = nidep3d(k,i)
           wsub3d_c(i,k,j)    =  wsub3d(k,i)
             ccn1_c(i,k,j)    =    ccn1(k,i)
             ccn2_c(i,k,j)    =    ccn2(k,i)
             ccn3_c(i,k,j)    =    ccn3(k,i)
             ccn4_c(i,k,j)    =    ccn4(k,i)
             ccn5_c(i,k,j)    =    ccn5(k,i)
             ccn6_c(i,k,j)    =    ccn6(k,i)
          nsource_c(i,k,j)    = nsource(k,i)
           rndst1_c(i,k,j)    =  rndst1(k,i)
           rndst2_c(i,k,j)    =  rndst2(k,i)
           rndst3_c(i,k,j)    =  rndst3(k,i)
           rndst4_c(i,k,j)    =  rndst4(k,i)
           nacon1_c(i,k,j)    =  nacon1(k,i)
           nacon2_c(i,k,j)    =  nacon2(k,i)
           nacon3_c(i,k,j)    =  nacon3(k,i)
           nacon4_c(i,k,j)    =  nacon4(k,i)
     enddo
     enddo
     enddo

!initialization:
 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    f_ice_c(i,k,j)  = 0.0_RKIND
    f_rain_c(i,k,j) = 0.0_RKIND
 enddo
 enddo
 enddo


 end subroutine drydep_from_MPAS

!=================================================================================================================
 subroutine drydep_to_MPAS(         state,diag,diag_physics,tend_chemistry,its,ite)
!=================================================================================================================

!  logical :: debug_on
!inout arguments:
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(inout):: diag
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_chemistry

 integer,intent(in):: its,ite
!local variables:
 integer:: i,k,j,h

!local pointers:
 real(kind=RKIND),dimension(:,:,:),pointer :: tend_mixdry
 real(kind=RKIND),dimension(:,:,:),pointer :: chem
 real(kind=RKIND),dimension(:,:,:),pointer :: scalars
 real(kind=RKIND),dimension(:,:)  ,pointer :: ddvel
 real(kind=RKIND),dimension(:,:)  ,pointer :: cldfrac
!real(kind=RKIND),dimension(:,:)  ,pointer :: cldfrac_old
 real(kind=RKIND),dimension(:,:)  ,pointer :: tke_pbl
 real(kind=RKIND),dimension(:,:)  ,pointer ::  relhum
 real(kind=RKIND),dimension(:,:)  ,pointer ::  nihf3d
 real(kind=RKIND),dimension(:,:)  ,pointer :: niimm3d
 real(kind=RKIND),dimension(:,:)  ,pointer :: nidep3d
 real(kind=RKIND),dimension(:,:)  ,pointer ::  wsub3d
 real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn1
 real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn2
 real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn3
 real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn4
 real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn5
 real(kind=RKIND),dimension(:,:)  ,pointer ::    ccn6
 real(kind=RKIND),dimension(:,:)  ,pointer :: nsource
 real(kind=RKIND),dimension(:,:)  ,pointer ::  rndst1
 real(kind=RKIND),dimension(:,:)  ,pointer ::  rndst2
 real(kind=RKIND),dimension(:,:)  ,pointer ::  rndst3
 real(kind=RKIND),dimension(:,:)  ,pointer ::  rndst4
 real(kind=RKIND),dimension(:,:)  ,pointer ::  nacon1
 real(kind=RKIND),dimension(:,:)  ,pointer ::  nacon2
 real(kind=RKIND),dimension(:,:)  ,pointer ::  nacon3
 real(kind=RKIND),dimension(:,:)  ,pointer ::  nacon4

!-------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_array(diag          ,'relhum'         ,relhum        )
 call mpas_pool_get_array(diag          ,'nihf3d'         ,nihf3d        )
 call mpas_pool_get_array(diag          ,'niimm3d'        ,niimm3d       )
 call mpas_pool_get_array(diag          ,'nidep3d'        ,nidep3d       )
 call mpas_pool_get_array(diag          ,'wsub3d'         ,wsub3d        )
 call mpas_pool_get_array(diag          ,'ccn1'           ,ccn1          )
 call mpas_pool_get_array(diag          ,'ccn2'           ,ccn2          )
 call mpas_pool_get_array(diag          ,'ccn3'           ,ccn3          )
 call mpas_pool_get_array(diag          ,'ccn4'           ,ccn4          )
 call mpas_pool_get_array(diag          ,'ccn5'           ,ccn5          )
 call mpas_pool_get_array(diag          ,'ccn6'           ,ccn6          )
 call mpas_pool_get_array(diag          ,'nsource'        ,nsource       )
 call mpas_pool_get_array(diag          ,'rndst1'         ,rndst1        )
 call mpas_pool_get_array(diag          ,'rndst2'         ,rndst2        )
 call mpas_pool_get_array(diag          ,'rndst3'         ,rndst3        )
 call mpas_pool_get_array(diag          ,'rndst4'         ,rndst4        )
 call mpas_pool_get_array(diag          ,'nacon1'         ,nacon1        )
 call mpas_pool_get_array(diag          ,'nacon2'         ,nacon2        )
 call mpas_pool_get_array(diag          ,'nacon3'         ,nacon3        )
 call mpas_pool_get_array(diag          ,'nacon4'         ,nacon4        )
 call mpas_pool_get_array(tend_chemistry,'tend_mixdry'    ,tend_mixdry   )
 call mpas_pool_get_array(tend_chemistry,'ddvel'          ,ddvel         )
 call mpas_pool_get_array(state         ,'chem'           ,chem          )
 call mpas_pool_get_array(state         ,'scalars'        ,scalars       )

 do h = 1,num_chem 
 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    if( .not. chem_c(i,k,j,h).gt.0.0_RKIND ) chem_c(i,k,j,h) = 0.0_RKIND
     tend_mixdry(h,k,i) = (chem_c(i,k,j,h) - chem(h,k,i))/dt_dyn
     chem(h,k,i)        = chem_c(i,k,j,h)
 enddo
 enddo
 enddo
 enddo

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
     scalars(p_qndrop,k,i)        = scalars_c(i,k,j,p_qndrop)
     scalars(p_qniceaer,k,i)      = scalars_c(i,k,j,p_qniceaer)
     scalars(p_qnicemix,k,i)      = scalars_c(i,k,j,p_qnicemix)
 enddo
 enddo
 enddo


 do j = jts,jte
 do k = kts,kte
 do i = its,ite
      relhum(k,i)       = relhum3d_c(i,k,j)
      nihf3d(k,i)       =   nihf3d_c(i,k,j)
     niimm3d(k,i)       =  niimm3d_c(i,k,j)
     nidep3d(k,i)       =  nidep3d_c(i,k,j)
      wsub3d(k,i)       =   wsub3d_c(i,k,j)
        ccn1(k,i)       =     ccn1_c(i,k,j)
        ccn2(k,i)       =     ccn2_c(i,k,j)
        ccn3(k,i)       =     ccn3_c(i,k,j)
        ccn4(k,i)       =     ccn4_c(i,k,j)
        ccn5(k,i)       =     ccn5_c(i,k,j)
        ccn6(k,i)       =     ccn6_c(i,k,j)
     nsource(k,i)       =  nsource_c(i,k,j)
      rndst1(k,i)       =   rndst1_c(i,k,j)
      rndst2(k,i)       =   rndst2_c(i,k,j)
      rndst3(k,i)       =   rndst3_c(i,k,j)
      rndst4(k,i)       =   rndst4_c(i,k,j)
      nacon1(k,i)       =   nacon1_c(i,k,j)
      nacon2(k,i)       =   nacon2_c(i,k,j)
      nacon3(k,i)       =   nacon3_c(i,k,j)
      nacon4(k,i)       =   nacon4_c(i,k,j)
 enddo
 enddo
 enddo

 do h = 1,num_chem
 do j = jts,jte
 do i = its,ite
     ddvel(h,i)        = ddvel_c(i,j,h)
 enddo
 enddo
 enddo

 end subroutine drydep_to_MPAS

!=================================================================================================================
 subroutine driver_drydep(itimestep,configs,mesh,state,sfc_input,diag,diag_physics,tend_chemistry,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input

 integer,intent(in):: its,ite
 integer,intent(in):: itimestep

!inout arguments:

 type(mpas_pool_type),intent(inout):: diag
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_chemistry
 type(mpas_pool_type),intent(inout):: state
 
!local variables:

 logical,pointer :: debug_on
 real, dimension(ims:ime,kms:kme,jms:jme,num_chem ) :: chem_old
 integer :: i,k,j,h,nv,ll,l
 integer :: itf,ktf,jtf
 integer :: k_a, k_c, kmax, m_mam
 integer :: idrydep_onoff,aer_mech_id
 integer :: iphase, itype
 real,parameter  :: m2cm = 100.
 real,dimension(ims:ime,config_flags%kdepvel,jms:jme,num_chem ) :: dep_vel 
 real,dimension(kms:kme) :: dryrho_1d
 real,dimension(ims:ime,jms:jme) :: frac_removed
 real,dimension(ims:ime,kms:kme,jms:jme) :: dryrho_phy
!turbulent transport
 real :: pblst(kts:kte),ekmfull(kts:kte+1),zzfull(kts:kte+1),zz(kts:kte)

! Used to pass parameters for subroutines.
 integer :: id,julday,numgas,ktau
 integer :: current_month,current_secs
 real    :: dtstep
 real, dimension(its:ite,        jts:jte          ) ::                    &
                                                      aer_res,            &
                                                      aer_res_def,        &
                                                      aer_res_zcen
 real, dimension(ims:ime,kms:kme,jms:jme,num_moist) ::                    &
                                                      moist  
 real, dimension(ims:ime,kms:kme,jms:jme,num_scalars) ::                  &
                                                      scalar  
 real, dimension(ims:ime,kms:kme,jms:jme,num_chem ) ::                    &
                                                      chem           
 real, dimension(ims:ime,        jms:jme,num_chem ) ::                    &
                                                      ddvel                     
 real, dimension(ims:ime,kms:kme,jms:jme          ) ::                    &
                                                      t_phy,              &
                                                      p_phy,              &
                                                      dz8w,               &
                                                      z,                  &
                                                      w,                  &
                                                      t8w,                &
                                                      p8w,                &
                                                      z_at_w,             &
                                                      alt,                &
                                                      exch_h,             &
                                                      h2oai,              &
                                                      h2oaj,              &
                                                      cldfra,             &
                                                      f_ice_phy,          &
                                                     f_rain_phy,          &
                                                      cldfra_old,         &
                                                      tke_pbl,            &
                                                      relhum3d,           &
                                                        nihf3d,           &
                                                       niimm3d,           &
                                                       nidep3d,           &
                                                        wsub3d,           &
                                                          ccn1,           &
                                                          ccn2,           &
                                                          ccn3,           &
                                                          ccn4,           &
                                                          ccn5,           &
                                                          ccn6,           &
                                                       nsource,           &
                                                        rndst1,           &
                                                        rndst2,           &
                                                        rndst3,           &
                                                        rndst4,           &
                                                        nacon1,           &
                                                        nacon2,           &
                                                        nacon3,           &
                                                        nacon4,           &
                                                       rho_phy
 integer,dimension(ims:ime,      jms:jme          ) ::                    &
                                                      ivgtyp
 real,   dimension(ims:ime,      jms:jme          ) ::                    &
                                                      tsk,                &
                                                      gsw,                &
                                                      vegfra,             &
                                                      pbl,                &
                                                      rmol,               &
                                                      ust,                &
                                                      hfx,                &
                                                      xland,              &
                                                      xlat,               &
                                                      xlong,              &
                                                      dx,                 &
                                                      raincv,             &
                                                      snowh,              &
                                                      znt

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_print_debug_info' ,debug_on )
 call drydep_from_MPAS(state,diag,diag_physics,its,ite)
 if(debug_on) call chemistry_message('drydep from MPAS success')

!-----------------------------------------------------------------------------------------------------------------
! Passing parameters for subroutines.
 id              = itimestep
 ktau            = itimestep
 dtstep          = dt_dyn
 dx              = dx_c
 current_month   = curr_month
 curr_secs       = curr_secs
 julday          = int(curr_julday)
 t_phy           = t_phy_c
 moist           = moist_c
 scalar          = scalars_c
 p8w             = pres2_hyd_c
 t8w             = t2_c
 raincv          = raincv_c
 p_phy           = pres_hyd_c
 chem            = chem_c
 chem_old        = chem
 rho_phy         = rho_phy_c
 alt             = alt_c
 dz8w            = dz_c
 ddvel           = ddvel_c
 ivgtyp          = ivgtyp_c
 xland           = xland_c
 tsk             = tsk_c
 gsw             = gsw_c
 vegfra          = vegfra_c
 pbl             = hpbl_c
 rmol            = rmol_c
 ust             = ust_c
 znt             = znt_c
 xlat            = xlat_c
 xlong           = xlong_c
 z               = z_c
 w               = w_c
 z_at_w          = z_c
 snowh           = snowh_c
 numgas          = num_gas
 hfx             = hfx_c
 aer_res_def     = 0
 aer_res_zcen    = 0
 aer_res         = 0
 exch_h          = exch_h_c
 cldfra          = cldfra_c
 cldfra_old      = cldfra_old_c
 f_ice_phy       = f_ice_c
 f_rain_phy      = f_rain_c

 tke_pbl         = tke_pbl_c
 relhum3d        = relhum3d_c
 nihf3d          = nihf3d_c
 niimm3d         = niimm3d_c
 nidep3d         = nidep3d_c
 wsub3d          = wsub3d_c
 ccn1            = ccn1_c
 ccn2            = ccn2_c
 ccn3            = ccn3_c
 ccn4            = ccn4_c
 ccn5            = ccn5_c
 ccn6            = ccn6_c
 nsource         = nsource_c
 rndst1          = rndst1_c
 rndst2          = rndst2_c
 rndst3          = rndst3_c
 rndst4          = rndst4_c
 nacon1          = nacon1_c
 nacon2          = nacon2_c
 nacon3          = nacon3_c
 nacon4          = nacon4_c

 jtf=min0(jte,jde-1)
 ktf=min0(kte,kde-1)
 itf=min0(ite,ide-1)

!-----------------------------------------------------------------------------------------------------------------
 
!
! compute dry deposition velocities = ddvel
!
! 28-jun-2005 rce - initialize ddvel=0; call aerosol drydep routine
!           only when gas_drydep_opt == WESELY
!       the wesely_driver routine computes aer_res, and currently
!	    you cannot compute aerosol drydep without it !!
! 08-jul-2005 rce - pass idrydep_onoff to mixactivate routines
!

 ddvel(its:ite,jts:jte,num_chem) = 0.0_RKIND

 idrydep_onoff = 0

!jwfeng: You can use wesely to calculate aer_res, if you call wesely_driver
!        in the future, you don't need calculate_aer_res to get aer_res.

 call wesely_driver(id,ktau,dtstep,                          &
      config_flags,current_month,                                &
      gmt,julday,t_phy,moist,p8w,t8w,raincv,                     &
      p_phy,chem,rho_phy,dz8w,ddvel,aer_res_def,aer_res_zcen,    &
      ivgtyp,tsk,gsw,vegfra,pbl,rmol,ust,znt,xlat,xlong,z,z_at_w,&
      snowh,numgas,                                              &
      ids,ide, jds,jde, kds,kde,                                 &
      ims,ime, jms,jme, kms,kme,                                 &
      its,ite, jts,jte, kts,kte                                  )

 drydep_select: SELECT CASE(config_flags%gas_drydep_opt)

     CASE ( WESELY )

     IF(config_flags%chem_opt /= CHEM_TRACER    .and.                   &
          config_flags%chem_opt /= CHEM_TRACE2    .and.                   &
          config_flags%chem_opt /= CO2_TRACER     .and.                   & 
          config_flags%chem_opt /= GHG_TRACER     .and.                   &
          config_flags%chem_opt /= CHEM_VASH      .and.                   &
          config_flags%chem_opt /= CHEM_VOLC_4BIN .and.                   &
          config_flags%chem_opt /= DUST           .and.                   &
          config_flags%chem_opt /= GOCART_SIMPLE  .and.                   &
          config_flags%chem_opt /= GOCARTRACM_KPP ) THEN     

! this does aerosol species (dust,seas, bc,oc) for gocart only

          IF(config_flags%chem_opt == MOZCART_KPP ) then
            call gocart_drydep_driver( dtstep,                            &
                  config_flags, numgas,                                   &
                  t_phy, moist, p8w, t8w, rmol,aer_res_def,               &
                  p_phy, chem, rho_phy, dz8w, ddvel, xland, hfx,          &
                  ivgtyp, tsk, vegfra, pbl, ust, znt, xlat, xlong,        &
                  ids,ide, jds,jde, kds,kde,                              &
                  ims,ime, jms,jme, kms,kme,                              &
                  its,ite, jts,jte, kts,kte )
          ENDIF
     ELSEIF ( config_flags%chem_opt == GOCART_SIMPLE ) then
!
! this does aerosol species (dust,seas, bc,oc) for gocart only
!
          call gocart_drydep_driver(dtstep,                                &
               config_flags,numgas,                                       &
               t_phy,moist,p8w,t8w,rmol,aer_res_def,                      &
               p_phy,chem,rho_phy,dz8w,ddvel,xland,hfx,                   &
               ivgtyp,tsk,vegfra,pbl,ust,znt,xlat,xlong,                  &
               ids,ide, jds,jde, kds,kde,                                 &
               ims,ime, jms,jme, kms,kme,                                 &
               its,ite, jts,jte, kts,kte                                  )

     ELSEIF ( config_flags%chem_opt == DUST ) then
!
! this does aerosol species (dust) for gocart only
!
          call gocart_drydep_driver(dtstep,                                &
               config_flags,numgas,                                        &
               t_phy,moist,p8w,t8w,rmol,aer_res_def,                       &
               p_phy,chem,rho_phy,dz8w,ddvel,xland,hfx,                    &
               ivgtyp,tsk,vegfra,pbl,ust,znt,xlat,xlong,                   &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )

     ELSEIF ( config_flags%chem_opt == GOCARTRACM_KPP) then
!
! this does aerosol species (dust,seas, bc,oc) for gocart only
!
          call gocart_drydep_driver(dtstep,                                &
               config_flags,numgas,                                       &
               t_phy,moist,p8w,t8w,rmol,aer_res_def,                      &
               p_phy,chem,rho_phy,dz8w,ddvel,xland,hfx,                   &
               ivgtyp,tsk,vegfra,pbl,ust,znt,xlat,xlong,                  &
               ids,ide, jds,jde, kds,kde,                                 &
               ims,ime, jms,jme, kms,kme,                                 &
               its,ite, jts,jte, kts,kte                                  )

     ELSE
          !Set dry deposition velocity to zero when using the
          !chemistry tracer mode.
          ddvel(its:ite,jte:jte,1:num_chem) = 0.
     END IF

     aer_res(its:ite,jts:jte) = aer_res_def(its:ite,jts:jte)
     if(debug_on) call mpas_log_write("--- Min/Max of aer_res = $r/$r",realArgs=(/maxval(aer_res(its:ite,jts:jte)),minval(aer_res(its:ite,jts:jte))/))
!    end if

     idrydep_onoff = 1
     aer_mech_id = 3

     if ((config_flags%aer_drydep_opt <= 0) .or. (aer_mech_id <= 0)) then
        if(debug_on) call chemistry_message('AEROSOL DRY DEP VELOCITIES  = 0.0')
     else if (config_flags%aer_drydep_opt <= 99) then
                                                    
!    END SELECT adrydep_select

     else ! (config_flags%aer_drydep_opt > 99)
        if(debug_on) call chemistry_message('DOING DRY DEP VELOCITIES THRU AER_DRYDEP_DRIVER')
        call aer_drydep_driver(debug_on,                                 &
             id, ktau, dtstep, config_flags, aer_mech_id,             &
             gmt, julday,                                             &
             t_phy, rho_phy, p_phy,                                   &
             alt, p8w, t8w, dz8w, z, z_at_w,                          &
             ust, aer_res, ivgtyp, vegfra, pbl, rmol, znt,            &
             moist, chem, ddvel,                                      &
             h2oai, h2oaj, numgas,                                    &
             ids,ide, jds,jde, kds,kde,                               &
             ims,ime, jms,jme, kms,kme,                               &
             its,ite, jts,jte, kts,kte                                )

     endif

     if (config_flags%aer_drydep_opt > 0) then
        if ((aer_mech_id > 0) .and. (aer_mech_id <= 4)) then
             ddvel(its:ite,jts:jte,numgas+1:num_chem) = min(0.50,ddvel(its:ite,jts:jte,numgas+1:num_chem) )
        endif
     endif

    CASE DEFAULT 
                                                     
 END SELECT drydep_select                              

 ll = num_chem    ! not confirmed now

 dep_vel(:,:,:,:) = 0.
 do l=1,ll
 do j=jts,jte
 do k=1,config_flags%kdepvel
 do i=its,ite
   dep_vel(i,k,j,l) = ddvel(i,j,l)
 enddo
 enddo
 enddo
 enddo


!   This will be called later from subgrd_transport_driver.F !!!!!!!!
!
!
 if (num_vert_mix == 0) then
    do j=jts,jte
    do i=its,ite
    pblst=0.
!
!====================================================================================
!Start with vertical turbulent mixing    -czhao
!====================================================================================
!
    do k=kts,kte+1
       zzfull(k)=z_at_w(i,k,j)-z_at_w(i,kts,j)
    enddo
    do k=kts,kte
       ekmfull(k)=max(1.e-6,exch_h(i,k,j))
    enddo
    ekmfull(kts)=0.
    ekmfull(kte+1)=0.

    do k=kts,kte
       zz(k)=(z_at_w(i,k+1,j)+z_at_w(i,k,j))/2.0 -z_at_w(i,kts,j)
    enddo

    do nv=1,num_chem-0
       do k=kts,kte
          pblst(k)=max(epsilc,chem(i,k,j,nv))
          dryrho_1d(k) = 1./alt(i,k,j)
       enddo

       if(.not.is_aerosol(nv)) call vertmx(dtstep,pblst,ekmfull,dryrho_1d, zzfull,zz,ddvel(i,j,nv),kts,kte)


       do k=kts,kte-1
          chem(i,k,j,nv)=max(epsilc,pblst(k))
       enddo
    enddo
    enddo !j loop
    enddo !i loop

 endif ! num_vert_mix = 0
!
!  vertical mixing and activation of aerosol
!
 where( alt(its:ite,kts:kte,jts:jte) /= 0. )  !get dry density to conserve mass in mixactivate, wig, 24-apr-2006
    dryrho_phy(its:ite,kts:kte,jts:jte) = 1./alt(its:ite,kts:kte,jts:jte)
 elsewhere
    dryrho_phy(its:ite,kts:kte,jts:jte) = 0.
 end where


!----------------------------------------
!Vertical turbulent mixing with dry deposition and aerosol activation for aerosols only  -czhao
!----------------------------------------
 if (cw_phase > 0) call mosaic_mixactivate (                        &
      id, ktau, dtstep, config_flags, idrydep_onoff,   &
      p_phy,dryrho_phy, t_phy, w, cldfra, cldfra_old,  &
      ddvel, z, dz8w, p8w, t8w, exch_h,                &
      moist(ims,kms,jms,P_QV), moist(ims,kms,jms,P_QC),&
      moist(ims,kms,jms,P_QI), moist(ims,kms,jms,P_QS),&
      f_ice_phy,f_rain_phy, scalar(ims,kms,jms,P_QNDROP),&
      scalar(ims,kms,jms,P_QNICEAER),&
      f_qv,f_qc, f_qi,f_qs, chem,   &
      ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,      &
      nihf3d,niimm3d,nidep3d,scalar(ims,kms,jms,p_qnicemix),  &
      rndst1,rndst2,rndst3,rndst4,nacon1,nacon2,nacon3,nacon4,  &
      tke_pbl,relhum3d,wsub3d,   &
      ids,ide, jds,jde, kds,kde,                        &
      ims,ime, jms,jme, kms,kme,                        &
      its,ite, jts,jte, kts,kte                         )

 do l=1,num_chem
 do j=jts,jte
 do k=kts,kte
 do i=its,ite
        chem(i,k,j,l) = max(epsilc,chem(i,k,j,l))
 enddo
 enddo
 enddo
 enddo

 if(debug_on) call chemistry_message('call mosaic settling routine')
 
 call mosaic_settling_driver(debug_on,dtstep,config_flags,t_phy,moist,  & 
      chem,rho_phy,dz8w,p8w,p_phy,dx,g, &
      ids,ide, jds,jde, kds,kde,                                        &
      ims,ime, jms,jme, kms,kme,                                        &
      its,ite, jts,jte, kts,kte                                         )

     scalars_c =     scalar
        chem_c =       chem
       ddvel_c =      ddvel
     tke_pbl_c =    tke_pbl
    relhum3d_c =   relhum3d
      nihf3d_c =     nihf3d
     niimm3d_c =    niimm3d
     nidep3d_c =    nidep3d
      wsub3d_c =     wsub3d
        ccn1_c =       ccn1
        ccn2_c =       ccn2
        ccn3_c =       ccn3
        ccn4_c =       ccn4
        ccn5_c =       ccn5
        ccn6_c =       ccn6
     nsource_c =    nsource
      rndst1_c =     rndst1
      rndst2_c =     rndst2
      rndst3_c =     rndst3
      rndst4_c =     rndst4
      nacon1_c =     nacon1
      nacon2_c =     nacon2
      nacon3_c =     nacon3
      nacon4_c =     nacon4

 if(debug_on) call chemistry_message('drydep driver success')
 call drydep_to_MPAS(         state,diag,diag_physics,tend_chemistry,its,ite)
 if(debug_on) call chemistry_message('drydep to MPAS success')

 end subroutine driver_drydep

!=================================================================================================================
END MODULE mpas_atmchem_driver_drydep
!=================================================================================================================

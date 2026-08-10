!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to compute optical arguments for chem(:,:,:,:). 
! * Author: Jiawang Feng from AEMOL USTC
! Last update: jwfeng 04/2020
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

module mpas_atmchem_driver_optical
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines
 use mpas_atmchem_vars
 use mpas_atmchem_configs
 use mpas_atmchem_utilities
 use mpas_atmphys_vars

 use module_optical_averaging

 implicit none
 private
 public:: allocate_optical,   &
          deallocate_optical, &
          driver_optical

CONTAINS

!=================================================================================================================
 subroutine allocate_optical
!=================================================================================================================
    if(.not.allocated( relhum3d_c) )      allocate(      relhum3d_c(ims:ime, kms:kme, jms:jme            ) )
    if(.not.allocated( tauaersw_c) )      allocate(      tauaersw_c(ims:ime, kms:kme, jms:jme, 1:nswbands) )
    if(.not.allocated( absaersw_c) )      allocate(      absaersw_c(ims:ime, kms:kme, jms:jme, 1:nswbands) )
    if(.not.allocated(   gaersw_c) )      allocate(        gaersw_c(ims:ime, kms:kme, jms:jme, 1:nswbands) )
    if(.not.allocated(   waersw_c) )      allocate(        waersw_c(ims:ime, kms:kme, jms:jme, 1:nswbands) )
    if(.not.allocated( extaersw_c) )      allocate(      extaersw_c(ims:ime, kms:kme, jms:jme, 1:nswbands) )
    if(.not.allocated( bscoefsw_c) )      allocate(      bscoefsw_c(ims:ime, kms:kme, jms:jme, 1:nswbands) )
    if(.not.allocated( tauaerlw_c) )      allocate(      tauaerlw_c(ims:ime, kms:kme, jms:jme, 1:nlwbands) )
    if(.not.allocated( tauaerlw_c) )      allocate(      tauaerlw_c(ims:ime, kms:kme, jms:jme, 1:nlwbands) )
    if(.not.allocated( extaerlw_c) )      allocate(      extaerlw_c(ims:ime, kms:kme, jms:jme, 1:nlwbands) )
    if(.not.allocated(    l2aer_c) )      allocate(         l2aer_c(ims:ime, kms:kme, jms:jme, 1:nswbands) )
    if(.not.allocated(    l3aer_c) )      allocate(         l3aer_c(ims:ime, kms:kme, jms:jme, 1:nswbands) )
    if(.not.allocated(    l4aer_c) )      allocate(         l4aer_c(ims:ime, kms:kme, jms:jme, 1:nswbands) )
    if(.not.allocated(    l5aer_c) )      allocate(         l5aer_c(ims:ime, kms:kme, jms:jme, 1:nswbands) )
    if(.not.allocated(    l6aer_c) )      allocate(         l6aer_c(ims:ime, kms:kme, jms:jme, 1:nswbands) )
    if(.not.allocated(    l7aer_c) )      allocate(         l7aer_c(ims:ime, kms:kme, jms:jme, 1:nswbands) )
 end subroutine allocate_optical

!=================================================================================================================
 subroutine deallocate_optical
!=================================================================================================================
    if(allocated( relhum3d_c) )      deallocate(      relhum3d_c )
    if(allocated( tauaersw_c) )      deallocate(      tauaersw_c )
    if(allocated( absaersw_c) )      deallocate(      absaersw_c )
    if(allocated(   gaersw_c) )      deallocate(        gaersw_c )
    if(allocated(   waersw_c) )      deallocate(        waersw_c )
    if(allocated( extaersw_c) )      deallocate(      extaersw_c )
    if(allocated( bscoefsw_c) )      deallocate(      bscoefsw_c )
    if(allocated( tauaerlw_c) )      deallocate(      tauaerlw_c )
    if(allocated( tauaerlw_c) )      deallocate(      tauaerlw_c )
    if(allocated( extaerlw_c) )      deallocate(      extaerlw_c )
    if(allocated(    l2aer_c) )      deallocate(         l2aer_c )
    if(allocated(    l3aer_c) )      deallocate(         l3aer_c )
    if(allocated(    l4aer_c) )      deallocate(         l4aer_c )
    if(allocated(    l5aer_c) )      deallocate(         l5aer_c )
    if(allocated(    l6aer_c) )      deallocate(         l6aer_c )
    if(allocated(    l7aer_c) )      deallocate(         l7aer_c )
 end subroutine deallocate_optical

!=================================================================================================================
 subroutine optical_from_MPAS(diag, its, ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: diag
 integer,intent(in):: its,ite



!local variables:
 integer:: i,k,j,h
 real(kind=RKIND),dimension(:,:)  ,pointer ::  relhum
!-----------------------------------------------------------------------------------------------------------------

    call mpas_pool_get_array(diag        ,'relhum'     ,relhum        )
    do j = jts,jte
    do k = kts,kte
    do i = its,ite
        relhum3d_c(i,k,j)    =  relhum(k,i)
    end do
    end do
    end do

         tauaersw_c(its:ite,kts:kte,jts:jte,1:nswbands) = 0.0 
         absaersw_c(its:ite,kts:kte,jts:jte,1:nswbands) = 0.0 
           gaersw_c(its:ite,kts:kte,jts:jte,1:nswbands) = 0.0 
           waersw_c(its:ite,kts:kte,jts:jte,1:nswbands) = 0.0 
         extaersw_c(its:ite,kts:kte,jts:jte,1:nswbands) = 0.0
         bscoefsw_c(its:ite,kts:kte,jts:jte,1:nswbands) = 0.0
         extaerlw_c(its:ite,kts:kte,jts:jte,1:nlwbands) = 0.0 
         tauaerlw_c(its:ite,kts:kte,jts:jte,1:nlwbands) = 0.0
            l2aer_c(its:ite,kts:kte,jts:jte,1:nswbands) = 0.0
            l3aer_c(its:ite,kts:kte,jts:jte,1:nswbands) = 0.0
            l4aer_c(its:ite,kts:kte,jts:jte,1:nswbands) = 0.0
            l5aer_c(its:ite,kts:kte,jts:jte,1:nswbands) = 0.0
            l6aer_c(its:ite,kts:kte,jts:jte,1:nswbands) = 0.0
            l7aer_c(its:ite,kts:kte,jts:jte,1:nswbands) = 0.0

 end subroutine optical_from_MPAS

!=================================================================================================================
 subroutine optical_to_MPAS(debug_on, diag_chemistry, its, ite)
!=================================================================================================================

 logical:: debug_on
!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_chemistry
 integer,intent(in):: its,ite

!local variables:
 integer:: i,k,j,h
 character(len=StrKIND):: msg
!local pointers:
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaersw, absaersw, extaersw, gaersw, waersw, bscoefsw
 real(kind=RKIND),dimension(:,:,:),pointer :: tauaerlw, extaerlw
 real(kind=RKIND),dimension(:,:,:),pointer :: l2aer, l3aer, l4aer, l5aer, l6aer, l7aer
 real(kind=RKIND),dimension(:,:  ),pointer :: aodsw, aodlw
!-----------------------------------------------------------------------------------------------------------------
 call mpas_pool_get_array(diag_chemistry,'tauaersw'     , tauaersw      )
 call mpas_pool_get_array(diag_chemistry,'absaersw'     , absaersw      )
 call mpas_pool_get_array(diag_chemistry,'gaersw'       , gaersw        )
 call mpas_pool_get_array(diag_chemistry,'waersw'       , waersw        )
 call mpas_pool_get_array(diag_chemistry,'aodsw'        , aodsw         )
 call mpas_pool_get_array(diag_chemistry,'aodlw'        , aodlw         )
 call mpas_pool_get_array(diag_chemistry,'extaersw'     , extaersw      )
 call mpas_pool_get_array(diag_chemistry,'bscoefsw'     , bscoefsw      )
 call mpas_pool_get_array(diag_chemistry,'tauaerlw'     , tauaerlw      )
 call mpas_pool_get_array(diag_chemistry,'extaerlw'     , extaerlw      )
 call mpas_pool_get_array(diag_chemistry,'l2aer'        , l2aer         )
 call mpas_pool_get_array(diag_chemistry,'l3aer'        , l3aer         )
 call mpas_pool_get_array(diag_chemistry,'l4aer'        , l4aer         )
 call mpas_pool_get_array(diag_chemistry,'l5aer'        , l5aer         )
 call mpas_pool_get_array(diag_chemistry,'l6aer'        , l6aer         )
 call mpas_pool_get_array(diag_chemistry,'l7aer'        , l7aer         )

 aodsw(1:nswbands,its:ite) = 0.0

 do h = 1,nswbands
 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    tauaersw(h,k,i)      = tauaersw_c(i,k,j,h)
    absaersw(h,k,i)      = absaersw_c(i,k,j,h)
    gaersw(h,k,i)        = gaersw_c(i,k,j,h)
    waersw(h,k,i)        = waersw_c(i,k,j,h) 
    aodsw(h,i)           = aodsw(h,i) + tauaersw_c(i,k,j,h)
    extaersw(h,k,i)      = extaersw_c(i,k,j,h)
    bscoefsw(h,k,i)      = bscoefsw_c(i,k,j,h)
    l2aer(h,k,i)         = l2aer_c(i,k,j,h)
    l3aer(h,k,i)         = l3aer_c(i,k,j,h)
    l4aer(h,k,i)         = l4aer_c(i,k,j,h)
    l5aer(h,k,i)         = l5aer_c(i,k,j,h)
    l6aer(h,k,i)         = l6aer_c(i,k,j,h)
    l7aer(h,k,i)         = l7aer_c(i,k,j,h)
 enddo
 enddo
 enddo
 enddo

 if (debug_on) then
    write(msg,*) 'max/min of tauaersw_c      =',maxval(tauaersw_c(its:ite,kts:kte,jts:jte,1:nswbands)),minval(tauaersw_c(its:ite,kts:kte,jts:jte,1:nswbands))
    call chemistry_message(msg)
 end if

 aodlw(1:nlwbands,its:ite) = 0.0
 do h = 1, nlwbands
 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    tauaerlw(h,k,i) = tauaerlw_c(i,k,j,h)
    aodlw(h,i)      = aodlw(h,i) + tauaerlw_c(i,k,j,h)
    extaerlw(h,k,i) = extaerlw_c(i,k,j,h)
 enddo
 enddo
 enddo
 enddo

 end subroutine optical_to_MPAS

subroutine driver_optical(diag, diag_chemistry, configs, its, ite)
    use  module_data_mosaic_therm, only: nbin_a
!------------------------------------------------------------------------
 IMPLICIT NONE
 type(mpas_pool_type),intent(in   ):: diag
 type(mpas_pool_type),intent(in   ):: configs
 type(mpas_pool_type),intent(inout):: diag_chemistry
 integer,intent(in):: its,ite

!         
! local variables
!
 logical :: haveaer = .true.
 integer :: option_method, option_mie
 integer :: i,j,k
 character(len=StrKIND),pointer :: aer_op_opt
 logical,pointer :: debug_on

 call mpas_pool_get_config(configs,'config_aer_op_opt',aer_op_opt)
 call mpas_pool_get_config(configs,'config_print_debug_info' ,debug_on )

 if(debug_on) call chemistry_message('optical averaging')

 aer_op_opt_select: SELECT CASE(trim(aer_op_opt))
 CASE ('volume_approx')
     option_method=1
     option_mie=1
 CASE ('maxwell_approx')
     option_method=2
     option_mie=1
 CASE ('volume_exact')
     option_method=1
     option_mie=2
 CASE ('maxwell_exact')
     option_method=2
     option_mie=2
 CASE ('shell_exact')
     option_method=3
     option_mie=2
 CASE DEFAULT
     call chemistry_message('WARNING: Invalid aer_op_opt. Defaulting to volume_approx.')
     option_method=1
     option_mie=1
 END SELECT aer_op_opt_select

 call optical_from_MPAS(diag, its, ite)

 if(debug_on) call chemistry_message('optical from MPAS success')
 if(debug_on) call chemistry_message('optical driver: call optical averaging')

 !BSINGH(PNNL)- The following variables had an undefined behavior at the boundaries.

 call optical_averaging(  debug_on = debug_on,                        &
         id = 1, curr_secs = real(curr_secs,kind=R8KIND), dtstep = dt_dyn,              &
         config_flags = config_flags, nbin_o = nbin_a,                &
         haveaer = haveaer, option_method = option_method, option_mie = option_mie, &
         chem = chem_c, dz8w = dz_c, alt = alt_c, relhum = relhum3d_c,&
      !  h2oai,h2oaj,h2oai_cz,h2oaj_cz,h2oac_cz,                      &
         tauaersw = tauaersw_c, absaersw = absaersw_c,                &
         extaersw = extaersw_c,   gaersw = gaersw_c,                  &
           waersw = waersw_c,   bscoefsw = bscoefsw_c,                &
         l2aer = l2aer_c, l3aer = l3aer_c, l4aer = l4aer_c,           &
         l5aer = l5aer_c, l6aer = l6aer_c, l7aer = l7aer_c,           &
         tauaerlw = tauaerlw_c, extaerlw = extaerlw_c,                &
         ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,        &
         ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,        &
         its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte         &
         )
 if(debug_on) call chemistry_message('optical averaging success')

 call optical_to_MPAS(debug_on, diag_chemistry, its, ite)
 if(debug_on) call chemistry_message('optical to MPAS success')

end subroutine driver_optical

END MODULE mpas_atmchem_driver_optical

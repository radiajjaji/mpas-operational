!************************************************************************
!
! Module to check if chemistry namelists are set up correctly.
! * Author: Jiawang Feng and Chun Zhao / AEMOL USTC
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
!   mpas_atmphys_control.F by Laura D. Fowler
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

module mpas_atmchem_control
 use mpas_dmpar
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmchem_utilities

 implicit none
 private
 public:: chemistry_namelist_check

 contains


!=================================================================================================================
 subroutine chemistry_namelist_check(configs)
!=================================================================================================================

!===========================================================================
!IN THIS SUBROUTINE, IT INCLUDES FOLLOWING PARTS:
!  + Check if chemistry namelists are set correctly. 
!                                         ---jwfeng   12/2019
!NOTE:
!    If you add a new chemistry namelist config in Registry.xml,
!    you'd better edit its "possible_values" in the same time.
!===========================================================================


!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: config_chem_opt, &
                                  config_dust_opt, &
                                  config_convection_scheme, &
                                  config_ctrans_opt

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine chemistry_namelist_check:'

 call mpas_pool_get_config(configs,'config_chem_opt'    ,config_chem_opt    )
 call mpas_pool_get_config(configs,'config_dust_opt'    ,config_dust_opt    )
 call mpas_pool_get_config(configs,'config_ctrans_opt'  ,config_ctrans_opt  )
 call mpas_pool_get_config(configs,'config_convection_scheme'  ,config_convection_scheme  )

!aerosol scheme:
 if(.not. (config_chem_opt .eq. 'no_chemistry'         .or. &
           config_chem_opt .eq. 'bin10aer_dust'        .or. &
           config_chem_opt .eq. 'bin9aer_dust'         .or. &
           config_chem_opt .eq. 'bin8aer_dust'         .or. &
           config_chem_opt .eq. 'bin5aer_dust'         .or. &
           config_chem_opt .eq. 'bin4aer_dust'         .or. &
           config_chem_opt .eq. 'bin10aer_dust_aq'     .or. &
           config_chem_opt .eq. 'bin9aer_dust_aq'      .or. &
           config_chem_opt .eq. 'bin8aer_dust_aq'      .or. &
           config_chem_opt .eq. 'bin5aer_dust_aq'      .or. &
           config_chem_opt .eq. 'bin4aer_dust_aq'))    then

    write(mpas_err_message,'(A,A10)') 'illegal value for config_chem_opt:', &
          trim(config_chem_opt)
    call chemistry_error_fatal(mpas_err_message)
 endif

!dust scheme:
 if(.not. (config_dust_opt .eq. 'gocart'     .or. &
           config_dust_opt .eq. 'mosaic'))    then

    write(mpas_err_message,'(A,A10)') 'illegal value for config_dust_opt:', &
          trim(config_dust_opt)
    call chemistry_error_fatal(mpas_err_message)
 endif

!ctrans scheme:
 if(.not. (config_ctrans_opt .eq. 'no_ctrans' .or. &
           config_ctrans_opt .eq. 'kfeta'     .or. &
           config_ctrans_opt .eq. 'general')) then

    write(mpas_err_message,'(A,A10)') 'illegal value for config_ctrans_opt:', &
          trim(config_ctrans_opt)
    call chemistry_error_fatal(mpas_err_message)
 endif
!Check if ctrans_opt matches convection_scheme
 if(trim(config_ctrans_opt) .eq. 'kfeta' .and. trim(config_convection_scheme) .ne. 'cu_kain_fritsch') then
    mpas_err_message = 'config_ctrans_opt: kfeta only works when config_convection_scheme is cu_kain_fritsch'
    call chemistry_error_fatal(mpas_err_message)
 endif

 end subroutine chemistry_namelist_check

!=================================================================================================================
 end module mpas_atmchem_control
!=================================================================================================================


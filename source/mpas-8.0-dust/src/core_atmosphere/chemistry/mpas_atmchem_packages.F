!************************************************************************
!
! Module to Set Up Definitions of All Chemistry Packages Through Configs.
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
!   mpas_atmphys_packages.F by Laura D. Fowler
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

module mpas_atmchem_packages
 use mpas_kind_types
 use mpas_derived_types,only : mpas_pool_type,mpas_io_context_type,MPAS_LOG_ERR
 use mpas_pool_routines,only : mpas_pool_get_config,mpas_pool_get_package
!! use mpas_io_units,only : stderrUnit
 use mpas_log,only : mpas_log_write


 implicit none
 private
 public:: atmchem_setup_packages

 contains


!=================================================================================================================
 function atmchem_setup_packages(configs,packages,iocontext) result(ierr)
!=================================================================================================================

!===========================================================================
!IN THIS FUNCTION, IT INCLUDES FOLLOWING PARTS:
!  + Definitions of all chemistry packages. 
!                                         ---jwfeng   12/2019
!NOTE:
!    If you want to add a package,you need to add a "package" in
!    Registry.xml and add a new config_chem_opt relating to your
!    new package.Don't forget to update package under chem vars.
!===========================================================================

!inout arguments:
type (mpas_pool_type), intent(inout) :: configs
type (mpas_pool_type), intent(inout) :: packages
type (mpas_io_context_type), intent(inout) :: iocontext

!local variables:
character(len=StrKIND),pointer:: config_chem_opt
logical,pointer:: bin10aer_dust_aq, bin9aer_dust_aq, bin8aer_dust_aq, bin5aer_dust_aq, bin4aer_dust_aq
logical,pointer:: bin10aer_dust, bin9aer_dust, bin8aer_dust, bin5aer_dust, bin4aer_dust
logical,pointer:: no_chemistry

integer :: ierr

!-----------------------------------------------------------------------------------------------------------------
!Debug Information,Not Necessary.
call mpas_log_write( ' ')
call mpas_log_write( '--- enter subroutine atmchem_setup_packages:')

ierr = 0

call mpas_log_write( '----- Setting up package variables -----')
call mpas_log_write( ' ')

!--- initialization of all packages for dust aerosol schemes:

call mpas_pool_get_config(configs,'config_chem_opt',config_chem_opt)

nullify(no_chemistry)
call mpas_pool_get_package(packages,'no_chemistryActive',no_chemistry)

nullify(bin10aer_dust)
call mpas_pool_get_package(packages,'bin10aer_dustActive',bin10aer_dust)

nullify(bin9aer_dust)
call mpas_pool_get_package(packages,'bin9aer_dustActive',bin9aer_dust)

nullify(bin8aer_dust)
call mpas_pool_get_package(packages,'bin8aer_dustActive',bin8aer_dust)

nullify(bin5aer_dust)
call mpas_pool_get_package(packages,'bin5aer_dustActive',bin5aer_dust)

nullify(bin4aer_dust)
call mpas_pool_get_package(packages,'bin4aer_dustActive',bin4aer_dust)

nullify(bin10aer_dust_aq)
call mpas_pool_get_package(packages,'bin10aer_dust_aqActive',bin10aer_dust_aq)

nullify(bin9aer_dust_aq)
call mpas_pool_get_package(packages,'bin9aer_dust_aqActive',bin9aer_dust_aq)

nullify(bin8aer_dust_aq)
call mpas_pool_get_package(packages,'bin8aer_dust_aqActive',bin8aer_dust_aq)

nullify(bin5aer_dust_aq)
call mpas_pool_get_package(packages,'bin5aer_dust_aqActive',bin5aer_dust_aq)

nullify(bin4aer_dust_aq)
call mpas_pool_get_package(packages,'bin4aer_dust_aqActive',bin4aer_dust_aq)

if(.not.associated(no_chemistry)     .or. &
   .not.associated(bin10aer_dust)    .or. &
   .not.associated(bin9aer_dust)     .or. &
   .not.associated(bin8aer_dust)     .or. &
   .not.associated(bin5aer_dust)     .or. &
   .not.associated(bin4aer_dust)     .or. &
   .not.associated(bin10aer_dust_aq) .or. &
   .not.associated(bin9aer_dust_aq)  .or. &
   .not.associated(bin8aer_dust_aq)  .or. &
   .not.associated(bin5aer_dust_aq)  .or. &
   .not.associated(bin4aer_dust_aq)) then
   call mpas_log_write('=================================================',messageType=MPAS_LOG_ERR)
   call mpas_log_write('* Error while setting up packages for aerosol schemes options(config_chem_opt). *',messageType=MPAS_LOG_ERR)
   call mpas_log_write('=================================================',messageType=MPAS_LOG_ERR)
   ierr = 1
   return
endif

no_chemistry     = .false.
bin10aer_dust    = .false.
bin9aer_dust     = .false.
bin8aer_dust     = .false.
bin5aer_dust     = .false.
bin4aer_dust     = .false.
bin10aer_dust_aq = .false.
bin9aer_dust_aq  = .false.
bin8aer_dust_aq  = .false.
bin5aer_dust_aq  = .false.
bin4aer_dust_aq  = .false.


CHEM_OPT:select case (config_chem_opt)
    case ('no_chemistry')
        no_chemistry    = .true.
    case ('bin10aer_dust')
        bin10aer_dust   = .true.
    case ('bin9aer_dust')
        bin9aer_dust    = .true.
    case ('bin8aer_dust')
        bin8aer_dust    = .true.
    case ('bin5aer_dust')
        bin5aer_dust    = .true.
    case ('bin4aer_dust')
        bin4aer_dust    = .true.
    case ('bin10aer_dust_aq')
        bin10aer_dust_aq   = .true.
    case ('bin9aer_dust_aq')
        bin9aer_dust_aq    = .true.
    case ('bin8aer_dust_aq')
        bin8aer_dust_aq    = .true.
    case ('bin5aer_dust_aq')
        bin5aer_dust_aq    = .true.
    case ('bin4aer_dust_aq')
        bin4aer_dust_aq    = .true.
    case default
     call mpas_log_write('=================================================',messageType=MPAS_LOG_ERR)

     call mpas_log_write('* Error while setting up packages for aerosol schemes options(config_chem_opt). *',messageType=MPAS_LOG_ERR)
     call mpas_log_write('=================================================',messageType=MPAS_LOG_ERR)
end select CHEM_OPT

call mpas_log_write( '    no_chemistry       = $l', logicArgs=(/no_chemistry     /) )
call mpas_log_write( '    bin10aer_dust      = $l', logicArgs=(/bin10aer_dust    /) )
call mpas_log_write( '    bin9aer_dust       = $l', logicArgs=(/bin9aer_dust     /) )
call mpas_log_write( '    bin8aer_dust       = $l', logicArgs=(/bin8aer_dust     /) )
call mpas_log_write( '    bin5aer_dust       = $l', logicArgs=(/bin5aer_dust     /) )
call mpas_log_write( '    bin4aer_dust       = $l', logicArgs=(/bin4aer_dust     /) )
call mpas_log_write( '    bin10aer_dust_aq   = $l', logicArgs=(/bin10aer_dust_aq /) )   
call mpas_log_write( '    bin9aer_dust_aq    = $l', logicArgs=(/bin9aer_dust_aq  /) )  
call mpas_log_write( '    bin8aer_dust_aq    = $l', logicArgs=(/bin8aer_dust_aq  /) )  
call mpas_log_write( '    bin5aer_dust_aq    = $l', logicArgs=(/bin5aer_dust_aq  /) )  
call mpas_log_write( '    bin4aer_dust_aq    = $l', logicArgs=(/bin4aer_dust_aq  /) )  

end function atmchem_setup_packages

!=================================================================================================================
 end module mpas_atmchem_packages
!=================================================================================================================

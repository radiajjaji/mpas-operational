!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to define config variables for chemistry configs. 
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

MODULE mpas_atmchem_configs
    use mpas_kind_types
    implicit none
    !config_chem_opt : default = 'bin8aer_dust_aq'
    character(len=StrKIND),parameter :: CHEM_TRACER                 = 'Not_Included'  ,&
                                        CHEM_TRACE2                 = 'Not_Included'  ,&
                                        CO2_TRACER                  = 'Not_Included'  ,&
                                        GHG_TRACER                  = 'Not_Included'  ,&
                                        CHEM_VASH                   = 'Not_Included'  ,& 
                                        CHEM_VOLC_4BIN              = 'Not_Included'  ,&
                                        DUST                        = 'Not_Included'  ,&
                                        GOCART_SIMPLE               = 'Not_Included'  ,& 
                                        GOCARTRACM_KPP              = 'Not_Included'  ,&
                                        MOZART_KPP                  = 'Not_Included'  ,&
                                        MOZCART_KPP                 = 'Not_Included'  ,&
                                        MOZART_MOSAIC_4BIN_VBS0_KPP = 'Not_Included'  ,&
                                        CBMZ_MOSAIC_DMS_8BIN_AQ     = 'bin8aer_dust_aq'

    integer,parameter :: RRTMG_SWSCHEME = 4    ,&
                         CAMUWPBLSCHEME = 9
    !config_gas_drydep_opt : default = 1
    integer,parameter :: WESELY           = 1

    !config_aer_drydep_opt : default = 101

    !config_diagnostic_chem : default = 1
    integer,parameter :: DEPVEL1          = 1
    integer,parameter :: DEPVEL2          = 2

    type MPAS_chem_config
        character(len=StrKIND) :: chem_opt
        character(len=StrKIND) :: dust_opt
        integer :: chem_in_opt
        integer :: gas_drydep_opt
        integer :: aer_drydep_opt
        integer :: diagnostic_chem
        integer :: kdepvel
        integer :: settling_option

        ! in mixactivate  some of these is in physics 
        !!! All of these need to be added to Registry.xml
        integer :: ra_sw_physics
        integer :: bl_pbl_physics
        integer :: wsubscheme
        integer :: ice_para
        integer :: mixcloud
    end type MPAS_chem_config

END MODULE mpas_atmchem_configs

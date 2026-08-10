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


module mpas_atmchem_binaer_init
use mpas_kind_types
use mpas_pool_routines

implicit none
public :: init_data_binaer_asect
save

contains
!-----------------------------------------------------------------------
! *** note - eventually is_aerosol will be a subr argument
subroutine init_data_binaer_asect( state, is_aerosol )
!   subroutine init_data_binaer_asect( id,n2o5_hetchem, is_aerosol )
!    subroutine init_data_mosaic_asect( )

    use module_data_binaer_asect
    use module_data_binaer_other, only:  lunerr, lunout,   &
        aboxtest_testmode, aboxtest_units_convert,   &
        aboxtest_rh_method, aboxtest_map_method,   &
        aboxtest_gases_fixed, aboxtest_min_temp,   &
        aboxtest_min_relhum, aboxtest_max_relhum
    use module_data_mosaic_therm, only:  nbin_a, nbin_a_maxd, n2o5_flag
    ! use module_mosaic_csuesat, only:  init_csuesat
    ! use module_mosaic_movesect, only:  move_sections, test_move_sections
    ! use module_peg_util, only:  peg_error_fatal
    use mpas_atmchem_utilities, only: chemistry_error_fatal, chemistry_message

    use mpas_atmchem_vars, only:  param_first_scalar, num_chem , &
        p_num_a01,   p_num_a02,   p_num_a03,   p_num_a04,   p_num_a05   , &
        p_num_a06,   p_num_a07,   p_num_a08,   p_num_a09,   p_num_a10   , &
        p_num_cw01,  p_num_cw02,  p_num_cw03,  p_num_cw04,  p_num_cw05  , &
        p_num_cw06,  p_num_cw07,  p_num_cw08,  p_num_cw09,  p_num_cw10

    implicit none

! *** note - eventually is_aerosol will be a subr argument
    type(mpas_pool_type),intent(in) :: state
    logical, intent(out) :: is_aerosol(num_chem)

!   local variables
    integer idum, itype, l, ldum, n, nhi, nsize_aer_dum
    real dum
    real, parameter :: pi = 3.14159265
    character*200 msg


!
!   set indexs to MPAS chem-array species
!
    call init_data_binaer_index(state)

!
!   set some "pegasus" control variables
!
    msectional = 20
    maerocoag = -2
    maerchem = 1
    maeroptical = 1
    maerchem_boxtest_output = -1

!
!   set ntype_aer = 1
!
    ntype_aer = 1

!
!   set number of aerosol bins using the wrf-chem sulfate pointers
!
    nsize_aer(:) = 0
    itype=1
! jwfeng from AEMOL USTC
! Assume that all schemes have num_a0*.p_num_a0* will be -1 if num_a0* is not included in.
    if (p_num_a01 .ge. param_first_scalar) nsize_aer(itype) = 1
    if (p_num_a02 .ge. param_first_scalar) nsize_aer(itype) = 2
    if (p_num_a03 .ge. param_first_scalar) nsize_aer(itype) = 3
    if (p_num_a04 .ge. param_first_scalar) nsize_aer(itype) = 4
    if (p_num_a05 .ge. param_first_scalar) nsize_aer(itype) = 5
    if (p_num_a06 .ge. param_first_scalar) nsize_aer(itype) = 6
    if (p_num_a07 .ge. param_first_scalar) nsize_aer(itype) = 7
    if (p_num_a08 .ge. param_first_scalar) nsize_aer(itype) = 8
    if (p_num_a09 .ge. param_first_scalar) nsize_aer(itype) = 9
    if (p_num_a10 .ge. param_first_scalar) nsize_aer(itype) = 10

    if (nsize_aer(itype) .le. 0) then
        call chemistry_error_fatal(&
            'init_data_binaer_asect - nsize_aer = 0')
    else if (nsize_aer(itype) .gt. maxd_asize) then
        call chemistry_error_fatal(&
            'init_data_binaer_asect - nsize_aer > maxd_asize' )
    end if

!
!   set nbin_a to total number of aerosol bins (for all types)
!
    nbin_a = 0
    do itype = 1, ntype_aer
        nbin_a = nbin_a + nsize_aer(itype)
    end do
    write(msg,'(a,(1x,I4))') '--- nbin_a = ' , nbin_a
    call chemistry_message(msg)
    if (nbin_a .gt. nbin_a_maxd) then
        call chemistry_error_fatal(&
            'init_data_binaer_asect - nbin_a > nbin_a_maxd' )
    end if

!
!   set nphase_aer (number of active aerosol species phases),
!   the xx_phase, and maerosolincw
!
    nphase_aer = 0
    maerosolincw = 0
    if (nsize_aer(1) .gt. 0) then
        nphase_aer = 1
        ai_phase = 1

        ! jwfeng from AEMOL USTC:
        ! Assume that all schemes have num_a0*.p_num_a0* will be -1 if num_a0* is not included in.
        if (p_num_cw01 .ge. param_first_scalar) then
            nphase_aer = 2
            cw_phase = 2
            maerosolincw = 1
        end if
    end if


! #if defined ( aboxtest_box_testing_active )
! ! *** start of "box testing" code section ***
! !     these code lines should be inactive when running wrf-chem
! !
! !   set some variables to "box test" values
!      call aboxtest_get_extra_args( 10,   &
!          msectional, maerosolincw, maerocoag,   &
!          maerchem, maeroptical, maerchem_boxtest_output,   &
!          lunerr, lunout, idum, dum, dum )
!      call aboxtest_get_extra_args( 11,   &
!          aboxtest_testmode, aboxtest_units_convert,   &
!          aboxtest_rh_method, aboxtest_map_method,   &
!          aboxtest_gases_fixed, nsize_aer_dum,   &
!          idum, idum, idum, dum, dum )
 
!      itype = 1
!      if (nsize_aer_dum > 0) nsize_aer(itype) = nsize_aer_dum
 
!      aboxtest_min_temp = 0.0
!      aboxtest_min_relhum = 0.0
!      aboxtest_max_relhum = 1.0
! ! ***  end  of "box testing" code section ***
! #endif

! copy over N2O5 control flag (DL 9/6/2013)
! n2o5_flag = n2o5_hetchem

!
!   set master aerosol chemical types
!
    !ntot_mastercomp_aer = 96  !czhao add dust or tr1r1-4
    ntot_mastercomp_aer = 3 

        l = 1 
        mastercompindx_seas_aer = l    
        name_mastercomp_aer( l ) = 'seas'
        dens_mastercomp_aer( l ) =  dens_seas_aer
        mw_mastercomp_aer(   l ) =    mw_seas_aer
        hygro_mastercomp_aer(l ) = hygro_seas_aer

        l = 2 
        mastercompindx_oc_aer    = l    
        name_mastercomp_aer( l ) = 'oc'
        dens_mastercomp_aer( l ) =  dens_oc_aer
        mw_mastercomp_aer(   l ) =    mw_oc_aer
        hygro_mastercomp_aer(l ) = hygro_oc_aer

        l = 3 
        mastercompindx_dust_aer = l    
        name_mastercomp_aer( l ) = 'dust'
        dens_mastercomp_aer( l ) =  dens_dust_aer
        mw_mastercomp_aer(   l ) =    mw_dust_aer
        hygro_mastercomp_aer(l ) = hygro_dust_aer
!
!   set section size arrays
! 
    do itype = 1, ntype_aer 
        nhi = nsize_aer(itype) 
        dlo_sect(1,itype) = 3.90625e-6 
        if (nhi.eq.4 .or. nhi.eq.8) then
            dhi_sect(nhi,itype) = 10.0e-4 
        else if (nhi.eq.5 .or. nhi.eq.10) then
            dhi_sect(nhi,itype) = 40.0e-4 
        else if (nhi.eq.9) then
            dhi_sect(nhi,itype) = 20.0e-4 
        else
            write( msg,*) 'Wrong nsize_aer:', nhi
            call chemistry_error_fatal(msg)
        end if
        dum = log( dhi_sect(nhi,itype)/dlo_sect(1,itype) ) / nhi !replaced alog by log by Manish Shrivastava on 11/28/2011. alog denoted natural log in fortran 77. log(x) is natural log in fortran 90
        do n = 2, nhi
                dlo_sect(n,itype) = dlo_sect(1,itype) * exp( (n-1)*dum )
                dhi_sect(n-1,itype) = dlo_sect(n,itype)
        end do
        do n = 1, nhi
            dcen_sect(n,itype) = sqrt( dlo_sect(n,itype)*dhi_sect(n,itype) )
            volumlo_sect(n,itype) = (pi/6.) * (dlo_sect(n,itype)**3)
            volumhi_sect(n,itype) = (pi/6.) * (dhi_sect(n,itype)**3)
            volumcen_sect(n,itype) = (pi/6.) * (dcen_sect(n,itype)**3)
            sigmag_aer(n,itype) = (dhi_sect(n,itype)/dlo_sect(n,itype))**0.289
        end do
        write(msg,'(a,10(1x,f8.3))') '--- dlo_sect  (um):' , dlo_sect(:,itype) *1.0e4
        call chemistry_message(msg)
        write(msg,'(a,10(1x,f8.3))') '--- dhi_sect  (um):' , dhi_sect(:,itype) *1.0e4
        call chemistry_message(msg)
        write(msg,'(a,10(1x,f8.3))') '--- dcen_sect (um):' , dcen_sect(:,itype)*1.0e4
        call chemistry_message(msg)
    end do

!
!   set pointers to MPAS chem-array species
!
    call init_data_binaer_ptr( is_aerosol )

!
!   csuesat initialization
!
!    call init_csuesat

!
!   move_sect initialization (and testing)
!
!    subr move_sections( iflag, iclm, jclm, k, m )
!    call move_sections(    -1,    1,    1, 1, 1 )
!
!    call test_move_sections( 1,   1,    1, 1, 1 )
    

end subroutine init_data_binaer_asect

subroutine init_data_binaer_index(state)
    use mpas_kind_types
    use mpas_pool_routines
    use mpas_atmchem_vars, only: &
        p_qv,        p_qc,        p_qi,        p_qs            , &
        p_qndrop,    p_qniceaer,  p_qnicemix                   , &
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


    implicit none
    type(mpas_pool_type),intent(in) :: state

!jwfeng: 12/31/2019 
!call mpas subroutine to pass index_* to p_* 
!physics
    call mpas_pool_get_dimension(state, 'index_qv'      , p_qv           )
    call mpas_pool_get_dimension(state, 'index_qc'      , p_qc           )
    call mpas_pool_get_dimension(state, 'index_qi'      , p_qi           )
    call mpas_pool_get_dimension(state, 'index_qs'      , p_qs           )
    call mpas_pool_get_dimension(state, 'index_qndrop'  , p_qndrop       )
    call mpas_pool_get_dimension(state, 'index_qniceaer', p_qniceaer     )
    call mpas_pool_get_dimension(state, 'index_qnicemix', p_qnicemix     )

!chemistry
    call mpas_pool_get_dimension(state, 'index_seas_a01', p_seas_a01     )
    call mpas_pool_get_dimension(state, 'index_seas_a02', p_seas_a02     )
    call mpas_pool_get_dimension(state, 'index_seas_a03', p_seas_a03     )
    call mpas_pool_get_dimension(state, 'index_seas_a04', p_seas_a04     )
    call mpas_pool_get_dimension(state, 'index_seas_a05', p_seas_a05     )
    call mpas_pool_get_dimension(state, 'index_seas_a06', p_seas_a06     )
    call mpas_pool_get_dimension(state, 'index_seas_a07', p_seas_a07     )
    call mpas_pool_get_dimension(state, 'index_seas_a08', p_seas_a08     )

    call mpas_pool_get_dimension(state, 'index_oc_a01',   p_oc_a01       )
    call mpas_pool_get_dimension(state, 'index_oc_a02',   p_oc_a02       )
    call mpas_pool_get_dimension(state, 'index_oc_a03',   p_oc_a03       )
    call mpas_pool_get_dimension(state, 'index_oc_a04',   p_oc_a04       )
    call mpas_pool_get_dimension(state, 'index_oc_a05',   p_oc_a05       )
    call mpas_pool_get_dimension(state, 'index_oc_a06',   p_oc_a06       )
    call mpas_pool_get_dimension(state, 'index_oc_a07',   p_oc_a07       )
    call mpas_pool_get_dimension(state, 'index_oc_a08',   p_oc_a08       )

    call mpas_pool_get_dimension(state, 'index_dust_a01', p_dust_a01     )
    call mpas_pool_get_dimension(state, 'index_dust_a02', p_dust_a02     )
    call mpas_pool_get_dimension(state, 'index_dust_a03', p_dust_a03     )
    call mpas_pool_get_dimension(state, 'index_dust_a04', p_dust_a04     )
    call mpas_pool_get_dimension(state, 'index_dust_a05', p_dust_a05     )
    call mpas_pool_get_dimension(state, 'index_dust_a06', p_dust_a06     )
    call mpas_pool_get_dimension(state, 'index_dust_a07', p_dust_a07     )
    call mpas_pool_get_dimension(state, 'index_dust_a08', p_dust_a08     )
    call mpas_pool_get_dimension(state, 'index_dust_a09', p_dust_a09     )
    call mpas_pool_get_dimension(state, 'index_dust_a10', p_dust_a10     )

    call mpas_pool_get_dimension(state, 'index_hysw_a01', p_hysw_a01     )
    call mpas_pool_get_dimension(state, 'index_hysw_a02', p_hysw_a02     )
    call mpas_pool_get_dimension(state, 'index_hysw_a03', p_hysw_a03     )
    call mpas_pool_get_dimension(state, 'index_hysw_a04', p_hysw_a04     )
    call mpas_pool_get_dimension(state, 'index_hysw_a05', p_hysw_a05     )
    call mpas_pool_get_dimension(state, 'index_hysw_a06', p_hysw_a06     )
    call mpas_pool_get_dimension(state, 'index_hysw_a07', p_hysw_a07     )
    call mpas_pool_get_dimension(state, 'index_hysw_a08', p_hysw_a08     )
    call mpas_pool_get_dimension(state, 'index_hysw_a09', p_hysw_a09     )
    call mpas_pool_get_dimension(state, 'index_hysw_a10', p_hysw_a10     )

    call mpas_pool_get_dimension(state, 'index_water_a01', p_water_a01   )
    call mpas_pool_get_dimension(state, 'index_water_a02', p_water_a02   )
    call mpas_pool_get_dimension(state, 'index_water_a03', p_water_a03   )
    call mpas_pool_get_dimension(state, 'index_water_a04', p_water_a04   )
    call mpas_pool_get_dimension(state, 'index_water_a05', p_water_a05   )
    call mpas_pool_get_dimension(state, 'index_water_a06', p_water_a06   )
    call mpas_pool_get_dimension(state, 'index_water_a07', p_water_a07   )
    call mpas_pool_get_dimension(state, 'index_water_a08', p_water_a08   )
    call mpas_pool_get_dimension(state, 'index_water_a09', p_water_a09   )
    call mpas_pool_get_dimension(state, 'index_water_a10', p_water_a10   )

    call mpas_pool_get_dimension(state, 'index_num_a01', p_num_a01     )
    call mpas_pool_get_dimension(state, 'index_num_a02', p_num_a02     )
    call mpas_pool_get_dimension(state, 'index_num_a03', p_num_a03     )
    call mpas_pool_get_dimension(state, 'index_num_a04', p_num_a04     )
    call mpas_pool_get_dimension(state, 'index_num_a05', p_num_a05     )
    call mpas_pool_get_dimension(state, 'index_num_a06', p_num_a06     )
    call mpas_pool_get_dimension(state, 'index_num_a07', p_num_a07     )
    call mpas_pool_get_dimension(state, 'index_num_a08', p_num_a08     )
    call mpas_pool_get_dimension(state, 'index_num_a09', p_num_a09     )
    call mpas_pool_get_dimension(state, 'index_num_a10', p_num_a10     )
 
    call mpas_pool_get_dimension(state, 'index_seas_cw01', p_seas_cw01     )
    call mpas_pool_get_dimension(state, 'index_seas_cw02', p_seas_cw02     )
    call mpas_pool_get_dimension(state, 'index_seas_cw03', p_seas_cw03     )
    call mpas_pool_get_dimension(state, 'index_seas_cw04', p_seas_cw04     )
    call mpas_pool_get_dimension(state, 'index_seas_cw05', p_seas_cw05     )
    call mpas_pool_get_dimension(state, 'index_seas_cw06', p_seas_cw06     )
    call mpas_pool_get_dimension(state, 'index_seas_cw07', p_seas_cw07     )
    call mpas_pool_get_dimension(state, 'index_seas_cw08', p_seas_cw08     )

    call mpas_pool_get_dimension(state, 'index_oc_cw01',   p_oc_cw01       )
    call mpas_pool_get_dimension(state, 'index_oc_cw02',   p_oc_cw02       )
    call mpas_pool_get_dimension(state, 'index_oc_cw03',   p_oc_cw03       )
    call mpas_pool_get_dimension(state, 'index_oc_cw04',   p_oc_cw04       )
    call mpas_pool_get_dimension(state, 'index_oc_cw05',   p_oc_cw05       )
    call mpas_pool_get_dimension(state, 'index_oc_cw06',   p_oc_cw06       )
    call mpas_pool_get_dimension(state, 'index_oc_cw07',   p_oc_cw07       )
    call mpas_pool_get_dimension(state, 'index_oc_cw08',   p_oc_cw08       )

    call mpas_pool_get_dimension(state, 'index_dust_cw01', p_dust_cw01     )
    call mpas_pool_get_dimension(state, 'index_dust_cw02', p_dust_cw02     )
    call mpas_pool_get_dimension(state, 'index_dust_cw03', p_dust_cw03     )
    call mpas_pool_get_dimension(state, 'index_dust_cw04', p_dust_cw04     )
    call mpas_pool_get_dimension(state, 'index_dust_cw05', p_dust_cw05     )
    call mpas_pool_get_dimension(state, 'index_dust_cw06', p_dust_cw06     )
    call mpas_pool_get_dimension(state, 'index_dust_cw07', p_dust_cw07     )
    call mpas_pool_get_dimension(state, 'index_dust_cw08', p_dust_cw08     )
    call mpas_pool_get_dimension(state, 'index_dust_cw09', p_dust_cw09     )
    call mpas_pool_get_dimension(state, 'index_dust_cw10', p_dust_cw10     )

    call mpas_pool_get_dimension(state, 'index_num_cw01', p_num_cw01     )
    call mpas_pool_get_dimension(state, 'index_num_cw02', p_num_cw02     )
    call mpas_pool_get_dimension(state, 'index_num_cw03', p_num_cw03     )
    call mpas_pool_get_dimension(state, 'index_num_cw04', p_num_cw04     )
    call mpas_pool_get_dimension(state, 'index_num_cw05', p_num_cw05     )
    call mpas_pool_get_dimension(state, 'index_num_cw06', p_num_cw06     )
    call mpas_pool_get_dimension(state, 'index_num_cw07', p_num_cw07     )
    call mpas_pool_get_dimension(state, 'index_num_cw08', p_num_cw08     )
    call mpas_pool_get_dimension(state, 'index_num_cw09', p_num_cw09     )
    call mpas_pool_get_dimension(state, 'index_num_cw10', p_num_cw10     )

end subroutine init_data_binaer_index


!-----------------------------------------------------------------------
    subroutine init_data_binaer_ptr( is_aerosol )

    ! use module_configure
    use mpas_atmchem_vars, only:  param_first_scalar,num_chem  , &
        p_qv,        p_qc,        p_qi,        p_qs            , &
        p_qndrop,    p_qniceaer,  p_qnicemix                   , &
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


    use module_data_binaer_asect
    use module_data_binaer_other, only:   &
        kh2so4, khno3, khcl, knh3,kmsa, ko3, kh2o, ktemp,   &
        kso2, kh2o2, khcho, koh, kho2,   &
        kno3, kno, kno2, khono, kpan,   &
                lmaxd, l2maxd, ltot, ltot2, lunout, lunerr, &
                name,kpcg1_b_c,kpcg2_b_c,kpcg3_b_c,kpcg4_b_c, &
                kpcg5_b_c,kpcg6_b_c,kpcg7_b_c,kpcg8_b_c, &
                kpcg9_b_c,kpcg1_b_o,kpcg2_b_o,kpcg3_b_o, &
                kpcg4_b_o,kpcg5_b_o,kpcg6_b_o,kpcg7_b_o, &
                kpcg8_b_o,kpcg9_b_o,kopcg1_b_c,kopcg2_b_c,&
                kopcg3_b_c, kopcg4_b_c,kopcg5_b_c,kopcg6_b_c,&
                kopcg7_b_c,kopcg8_b_c,kopcg1_b_o,kopcg2_b_o,&
                kopcg3_b_o,kopcg4_b_o,kopcg5_b_o,kopcg6_b_o,&
                kopcg7_b_o,kopcg8_b_o,&
                kpcg1_f_c,kpcg2_f_c,kpcg3_f_c,kpcg4_f_c, &
                kpcg5_f_c,kpcg6_f_c,kpcg7_f_c,kpcg8_f_c, &
                kpcg9_f_c,kpcg1_f_o,kpcg2_f_o,kpcg3_f_o, &
                kpcg4_f_o,kpcg5_f_o,kpcg6_f_o,kpcg7_f_o, &
                kpcg8_f_o,kpcg9_f_o,kopcg1_f_c,kopcg2_f_c,&
                kopcg3_f_c, kopcg4_f_c,kopcg5_f_c,kopcg6_f_c,&
                kopcg7_f_c,kopcg8_f_c,kopcg1_f_o,kopcg2_f_o,&
                kopcg3_f_o,kopcg4_f_o,kopcg5_f_o,kopcg6_f_o,&
                kopcg7_f_o,kopcg8_f_o,                      &
                ksmpa,ksmpbb,                      &
                kant1_c,kant2_c,kant3_c,kant4_c,kant1_o,kant2_o, &
                kant3_o,kant4_o,                            &
                kbiog1_c,kbiog2_c,kbiog3_c,kbiog4_c,kbiog1_o,kbiog2_o, &
                kbiog3_o,kbiog4_o, &
                kn2o5, kclno2                            

    ! use module_peg_util, only:  peg_error_fatal, peg_message
    use mpas_atmchem_utilities, only: chemistry_error_fatal, chemistry_message
    use module_mosaic_wetscav, only: initwet
        ! use module_scalar_tables, only:  chem_dname_table

    implicit none

!   subr arguments
        logical, intent(out) :: is_aerosol(num_chem)
!   local variables
    integer l, ll, n, p1st
    integer iaddto_ncomp, iaddto_ncomp_plustracer
    integer l_mastercomp, lptr_dum
    integer mcindx_dum
    integer isize, itype, iphase
    integer nphasetxt, nsizetxt, nspectxt, ntypetxt
    integer ncomp_dum(maxd_asize,maxd_aphase)
    integer ncomp_plustracer_dum(maxd_asize,maxd_aphase)

    integer y_seas, y_oc, y_dust, y_hysw, y_water, y_num
    integer y_cw_dust, y_cw_num

    character*200 msg
    character*20 phasetxt, sizetxt, spectxt, typetxt

    p1st = param_first_scalar
!
!   set up pointers to aerosol species in the wrf-chem "chem" array
!   note:  lptr=1 points to the first chem species which is "unused"
!
    itype=1
    lptr_seas_aer(:,itype,:)      = -1
    lptr_oc_aer(:,itype,:)        = -1
    lptr_dust_aer(:,itype,:)      = -1
    lptr_so4_aer(:,itype, :)      = -1
    lptr_msa_aer(:,itype, :)      = -1
    lptr_no3_aer(:,itype, :)      = -1
    lptr_cl_aer(:,itype, :)       = -1
    lptr_seas_aer(:,itype, :)     = -1
    lptr_co3_aer(:,itype, :)      = -1
    lptr_nh4_aer(:,itype, :)      = -1
    lptr_na_aer(:,itype, :)       = -1
    lptr_ca_aer(:,itype, :)       = -1
    lptr_oin_aer(:,itype, :)      = -1
    lptr_dust_aer(:,itype, :)     = -1
    lptr_tr1r1_aer(:,itype, :)    = -1
    lptr_tr1r2_aer(:,itype, :)    = -1
    lptr_tr1r3_aer(:,itype, :)    = -1
    lptr_tr1r4_aer(:,itype, :)    = -1
    lptr_oc_aer(:,itype, :)       = -1
    lptr_bc_aer(:,itype, :)       = -1
    lptr_pcg1_b_c_aer(:,itype, :) = -1
    lptr_pcg2_b_c_aer(:,itype, :) = -1
    lptr_pcg3_b_c_aer(:,itype, :) = -1
    lptr_pcg4_b_c_aer(:,itype, :) = -1
    lptr_pcg5_b_c_aer(:,itype, :) = -1
    lptr_pcg6_b_c_aer(:,itype, :) = -1
    lptr_pcg7_b_c_aer(:,itype, :) = -1
    lptr_pcg8_b_c_aer(:,itype, :) = -1
    lptr_pcg9_b_c_aer(:,itype, :) = -1
    lptr_pcg1_b_o_aer(:,itype, :) = -1
    lptr_pcg2_b_o_aer(:,itype, :) = -1
    lptr_pcg3_b_o_aer(:,itype, :) = -1
    lptr_pcg4_b_o_aer(:,itype, :) = -1
    lptr_pcg5_b_o_aer(:,itype, :) = -1
    lptr_pcg6_b_o_aer(:,itype, :) = -1
    lptr_pcg7_b_o_aer(:,itype, :) = -1
    lptr_pcg8_b_o_aer(:,itype, :) = -1
    lptr_pcg9_b_o_aer(:,itype, :) = -1
    lptr_opcg1_b_c_aer(:,itype, :)= -1
    lptr_opcg2_b_c_aer(:,itype, :)= -1
    lptr_opcg3_b_c_aer(:,itype, :)= -1
    lptr_opcg4_b_c_aer(:,itype, :)= -1
    lptr_opcg5_b_c_aer(:,itype, :)= -1
    lptr_opcg6_b_c_aer(:,itype, :)= -1
    lptr_opcg7_b_c_aer(:,itype, :)= -1
    lptr_opcg8_b_c_aer(:,itype, :)= -1
    lptr_opcg1_b_o_aer(:,itype, :)= -1
    lptr_opcg2_b_o_aer(:,itype, :)= -1
    lptr_opcg3_b_o_aer(:,itype, :)= -1
    lptr_opcg4_b_o_aer(:,itype, :)= -1
    lptr_opcg5_b_o_aer(:,itype, :)= -1
    lptr_opcg6_b_o_aer(:,itype, :)= -1
    lptr_opcg7_b_o_aer(:,itype, :)= -1
    lptr_opcg8_b_o_aer(:,itype, :)= -1
    lptr_pcg1_f_c_aer(:,itype, :) = -1
    lptr_pcg2_f_c_aer(:,itype, :) = -1
    lptr_pcg3_f_c_aer(:,itype, :) = -1
    lptr_pcg4_f_c_aer(:,itype, :) = -1
    lptr_pcg5_f_c_aer(:,itype, :) = -1
    lptr_pcg6_f_c_aer(:,itype, :) = -1
    lptr_pcg7_f_c_aer(:,itype, :) = -1
    lptr_pcg8_f_c_aer(:,itype, :) = -1
    lptr_pcg9_f_c_aer(:,itype, :) = -1
    lptr_pcg1_f_o_aer(:,itype, :) = -1
    lptr_pcg2_f_o_aer(:,itype, :) = -1
    lptr_pcg3_f_o_aer(:,itype, :) = -1
    lptr_pcg4_f_o_aer(:,itype, :) = -1
    lptr_pcg5_f_o_aer(:,itype, :) = -1
    lptr_pcg6_f_o_aer(:,itype, :) = -1
    lptr_pcg7_f_o_aer(:,itype, :) = -1
    lptr_pcg8_f_o_aer(:,itype, :) = -1
    lptr_pcg9_f_o_aer(:,itype, :) = -1
    lptr_opcg1_f_c_aer(:,itype, :)= -1
    lptr_opcg2_f_c_aer(:,itype, :)= -1
    lptr_opcg3_f_c_aer(:,itype, :)= -1
    lptr_opcg4_f_c_aer(:,itype, :)= -1
    lptr_opcg5_f_c_aer(:,itype, :)= -1
    lptr_opcg6_f_c_aer(:,itype, :)= -1
    lptr_opcg7_f_c_aer(:,itype, :)= -1
    lptr_opcg8_f_c_aer(:,itype, :)= -1
    lptr_opcg1_f_o_aer(:,itype, :)= -1
    lptr_opcg2_f_o_aer(:,itype, :)= -1
    lptr_opcg3_f_o_aer(:,itype, :)= -1
    lptr_opcg4_f_o_aer(:,itype, :)= -1
    lptr_opcg5_f_o_aer(:,itype, :)= -1
    lptr_opcg6_f_o_aer(:,itype, :)= -1
    lptr_opcg7_f_o_aer(:,itype, :)= -1
    lptr_opcg8_f_o_aer(:,itype, :)= -1
    lptr_smpa_aer(:,itype, :)     = -1
    lptr_smpbb_aer(:,itype, :)    = -1
    lptr_ant1_c_aer(:,itype, :)   = -1
    lptr_ant2_c_aer(:,itype, :)   = -1
    lptr_ant3_c_aer(:,itype, :)   = -1
    lptr_ant4_c_aer(:,itype, :)   = -1
    lptr_ant1_o_aer(:,itype, :)   = -1
    lptr_ant2_o_aer(:,itype, :)   = -1
    lptr_ant3_o_aer(:,itype, :)   = -1
    lptr_ant4_o_aer(:,itype, :)   = -1
    lptr_biog1_c_aer(:,itype, :)  = -1
    lptr_biog2_c_aer(:,itype, :)  = -1
    lptr_biog3_c_aer(:,itype, :)  = -1
    lptr_biog4_c_aer(:,itype, :)  = -1
    lptr_biog1_o_aer(:,itype, :)  = -1
    lptr_biog2_o_aer(:,itype, :)  = -1
    lptr_biog3_o_aer(:,itype, :)  = -1
    lptr_biog4_o_aer(:,itype, :)  = -1


    waterptr_aer(:,itype)         = -1
    hyswptr_aer(:,itype)          = -1
    numptr_aer(:,itype,:)         = -1
 

    if (nsize_aer(itype) .ge. 1) then
        if (p_seas_a01  .ge. p1st) lptr_seas_aer(01,itype,ai_phase) = p_seas_a01
        if (p_oc_a01    .ge. p1st) lptr_oc_aer(01,itype,ai_phase)   = p_oc_a01
        if (p_dust_a01  .ge. p1st) lptr_dust_aer(01,itype,ai_phase) = p_dust_a01
        if (p_hysw_a01  .ge. p1st) hyswptr_aer(01,itype)            = p_hysw_a01
        if (p_water_a01 .ge. p1st) waterptr_aer(01,itype)           = p_water_a01
        if (p_num_a01   .ge. p1st) numptr_aer(01,itype,ai_phase)    = p_num_a01
    end if

    if (nsize_aer(itype) .ge. 2) then
        if (p_seas_a02  .ge. p1st) lptr_seas_aer(02,itype,ai_phase) = p_seas_a02
        if (p_oc_a02    .ge. p1st) lptr_oc_aer(02,itype,ai_phase)   = p_oc_a02
        if (p_dust_a02  .ge. p1st) lptr_dust_aer(02,itype,ai_phase) = p_dust_a02
        if (p_hysw_a02  .ge. p1st) hyswptr_aer(02,itype)            = p_hysw_a02
        if (p_water_a02 .ge. p1st) waterptr_aer(02,itype)           = p_water_a02
        if (p_num_a02   .ge. p1st) numptr_aer(02,itype,ai_phase)    = p_num_a02
    end if

    if (nsize_aer(itype) .ge. 3) then
        if (p_seas_a03  .ge. p1st) lptr_seas_aer(03,itype,ai_phase) = p_seas_a03
        if (p_oc_a03    .ge. p1st) lptr_oc_aer(03,itype,ai_phase)   = p_oc_a03
        if (p_dust_a03  .ge. p1st) lptr_dust_aer(03,itype,ai_phase) = p_dust_a03
        if (p_hysw_a03  .ge. p1st) hyswptr_aer(03,itype)            = p_hysw_a03
        if (p_water_a03 .ge. p1st) waterptr_aer(03,itype)           = p_water_a03
        if (p_num_a03   .ge. p1st) numptr_aer(03,itype,ai_phase)    = p_num_a03
    end if

    if (nsize_aer(itype) .ge. 4) then
        if (p_seas_a04  .ge. p1st) lptr_seas_aer(04,itype,ai_phase) = p_seas_a04
        if (p_oc_a04    .ge. p1st) lptr_oc_aer(04,itype,ai_phase)   = p_oc_a04
        if (p_dust_a04  .ge. p1st) lptr_dust_aer(04,itype,ai_phase) = p_dust_a04
        if (p_hysw_a04  .ge. p1st) hyswptr_aer(04,itype)            = p_hysw_a04
        if (p_water_a04 .ge. p1st) waterptr_aer(04,itype)           = p_water_a04
        if (p_num_a04   .ge. p1st) numptr_aer(04,itype,ai_phase)    = p_num_a04
    end if

!jwfeng:
!I don't know why wrfchem has no "if" begore hysw,water and num pointer in bin 5,6,7,8
!So I just copy that.        
    if (nsize_aer(itype) .ge. 5) then
        if (p_seas_a05  .ge. p1st) lptr_seas_aer(05,itype,ai_phase) = p_seas_a05
        if (p_oc_a05    .ge. p1st) lptr_oc_aer(05,itype,ai_phase)   = p_oc_a05
        if (p_dust_a05  .ge. p1st) lptr_dust_aer(05,itype,ai_phase) = p_dust_a05
                                   hyswptr_aer(05,itype)            = p_hysw_a05
                                   waterptr_aer(05,itype)           = p_water_a05
                                   numptr_aer(05,itype,ai_phase)    = p_num_a05
    end if

    if (nsize_aer(itype) .ge. 6) then
        if (p_seas_a06  .ge. p1st) lptr_seas_aer(06,itype,ai_phase) = p_seas_a06
        if (p_oc_a06    .ge. p1st) lptr_oc_aer(06,itype,ai_phase)   = p_oc_a06
        if (p_dust_a06  .ge. p1st) lptr_dust_aer(06,itype,ai_phase) = p_dust_a06
                                   hyswptr_aer(06,itype)            = p_hysw_a06
                                   waterptr_aer(06,itype)           = p_water_a06
                                   numptr_aer(06,itype,ai_phase)    = p_num_a06
    end if

    if (nsize_aer(itype) .ge. 7) then
        if (p_seas_a07  .ge. p1st) lptr_seas_aer(07,itype,ai_phase) = p_seas_a07
        if (p_oc_a07    .ge. p1st) lptr_oc_aer(07,itype,ai_phase)   = p_oc_a07
        if (p_dust_a07  .ge. p1st) lptr_dust_aer(07,itype,ai_phase) = p_dust_a07
                                   hyswptr_aer(07,itype)            = p_hysw_a07
                                   waterptr_aer(07,itype)           = p_water_a07
                                   numptr_aer(07,itype,ai_phase)    = p_num_a07
    end if

    if (nsize_aer(itype) .ge. 8) then
        if (p_seas_a08  .ge. p1st) lptr_seas_aer(08,itype,ai_phase) = p_seas_a08
        if (p_oc_a08    .ge. p1st) lptr_oc_aer(08,itype,ai_phase)   = p_oc_a08
        if (p_dust_a08  .ge. p1st) lptr_dust_aer(08,itype,ai_phase) = p_dust_a08
                                   hyswptr_aer(08,itype)            = p_hysw_a08
                                   waterptr_aer(08,itype)           = p_water_a08
                                   numptr_aer(08,itype,ai_phase)    = p_num_a08
    end if

    if (nsize_aer(itype) .ge. 9) then
        if (p_dust_a09  .ge. p1st) lptr_dust_aer(09,itype,ai_phase) = p_dust_a09
                                   hyswptr_aer(09,itype)            = p_hysw_a09
                                   waterptr_aer(09,itype)           = p_water_a09
                                   numptr_aer(09,itype,ai_phase)    = p_num_a09
    end if

    if (nsize_aer(itype) .ge. 10) then
        if (p_dust_a10  .ge. p1st) lptr_dust_aer(10,itype,ai_phase) = p_dust_a10
                                   hyswptr_aer(10,itype)            = p_hysw_a10
                                   waterptr_aer(10,itype)           = p_water_a10
                                   numptr_aer(10,itype,ai_phase)    = p_num_a10
    end if

if (cw_phase .gt. 0) then
!   this code is "active" only when cw species are in the registry
    if (nsize_aer(itype) .ge. 1) then
        if(p_seas_cw01.ge.p1st) lptr_seas_aer(01,itype,cw_phase) = p_seas_cw01
        if(p_oc_cw01  .ge.p1st) lptr_oc_aer(01,itype,cw_phase)   = p_oc_cw01
        if(p_dust_cw01.ge.p1st) lptr_dust_aer(01,itype,cw_phase) = p_dust_cw01
        numptr_aer(01,itype,cw_phase)        = p_num_cw01
    end if

    if (nsize_aer(itype) .ge. 2) then
        if(p_seas_cw02.ge.p1st) lptr_seas_aer(02,itype,cw_phase) = p_seas_cw02
        if(p_oc_cw02  .ge.p1st) lptr_oc_aer(02,itype,cw_phase)   = p_oc_cw02
        if(p_dust_cw02.ge.p1st) lptr_dust_aer(02,itype,cw_phase) = p_dust_cw02
        numptr_aer(02,itype,cw_phase)        = p_num_cw02
    end if

    if (nsize_aer(itype) .ge. 3) then
        if(p_seas_cw03.ge.p1st) lptr_seas_aer(03,itype,cw_phase) = p_seas_cw03
        if(p_oc_cw03  .ge.p1st) lptr_oc_aer(03,itype,cw_phase)   = p_oc_cw03
        if(p_dust_cw03.ge.p1st) lptr_dust_aer(03,itype,cw_phase) = p_dust_cw03
        numptr_aer(03,itype,cw_phase)        = p_num_cw03
    end if

    if (nsize_aer(itype) .ge. 4) then
        if(p_seas_cw04.ge.p1st) lptr_seas_aer(04,itype,cw_phase) = p_seas_cw04
        if(p_oc_cw04.ge.p1st)   lptr_oc_aer(04,itype,cw_phase)   = p_oc_cw04
        if(p_dust_cw04.ge.p1st) lptr_dust_aer(04,itype,cw_phase) = p_dust_cw04
        numptr_aer(04,itype,cw_phase)        = p_num_cw04
    end if

    if (nsize_aer(itype) .ge. 5) then
        if(p_seas_cw05.ge.p1st) lptr_seas_aer(05,itype,cw_phase) = p_seas_cw05
        if(p_oc_cw05.ge.p1st)   lptr_oc_aer(05,itype,cw_phase)   = p_oc_cw05
        if(p_dust_cw05.ge.p1st) lptr_dust_aer(05,itype,cw_phase) = p_dust_cw05
        numptr_aer(05,itype,cw_phase)        = p_num_cw05
    end if

    if (nsize_aer(itype) .ge. 6) then
        if(p_seas_cw06.ge.p1st) lptr_seas_aer(06,itype,cw_phase) = p_seas_cw06
        if(p_oc_cw06  .ge.p1st) lptr_oc_aer(06,itype,cw_phase)   = p_oc_cw06
        if(p_dust_cw06.ge.p1st) lptr_dust_aer(06,itype,cw_phase) = p_dust_cw06
        numptr_aer(06,itype,cw_phase)        = p_num_cw06
    end if

    if (nsize_aer(itype) .ge. 7) then
        if(p_seas_cw07.ge.p1st) lptr_seas_aer(07,itype,cw_phase) = p_seas_cw07
        if(p_oc_cw07  .ge.p1st) lptr_oc_aer(07,itype,cw_phase)   = p_oc_cw07
        if(p_dust_cw07.ge.p1st) lptr_dust_aer(07,itype,cw_phase) = p_dust_cw07
        numptr_aer(07,itype,cw_phase)        = p_num_cw07
    end if

    if (nsize_aer(itype) .ge. 8) then
        if(p_seas_cw08.ge.p1st) lptr_seas_aer(08,itype,cw_phase) = p_seas_cw08
        if(p_oc_cw08.ge.p1st)   lptr_oc_aer(08,itype,cw_phase)   = p_oc_cw08
        if(p_dust_cw08.ge.p1st) lptr_dust_aer(08,itype,cw_phase) = p_dust_cw08
        numptr_aer(08,itype,cw_phase)        = p_num_cw08
    end if

    if (nsize_aer(itype) .ge. 9) then
        if(p_dust_cw09.ge.p1st) lptr_dust_aer(09,itype,cw_phase) = p_dust_cw09
        numptr_aer(09,itype,cw_phase)        = p_num_cw09
    end if

    if (nsize_aer(itype) .ge. 10) then
        if(p_dust_cw10.ge.p1st) lptr_dust_aer(10,itype,cw_phase) = p_dust_cw10
        numptr_aer(10,itype,cw_phase)        = p_num_cw10
    end if

end if

!
!   define the massptr_aer and mastercompptr_aer pointers
!   and the name() species names
!

!   first initialize
    do l = 1, l2maxd
        write( name(l), '(a,i4.4,15x)' ) 'r', l
    end do
    massptr_aer(:,:,:,:) = -999888777
    mastercompptr_aer(:,:) = -999888777

do itype = 1, ntype_aer
    if (itype .eq. 1) then
        typetxt = ' '
        ntypetxt = 1
        if (ntype_aer .gt. 1) then
        typetxt = '_t1'
        ntypetxt = 3
        end if
    else if (itype .le. 9) then
        write(typetxt,'(a,i1)') '_t', itype
        ntypetxt = 3
    else if (itype .le. 99) then
        write(typetxt,'(a,i2)') '_t', itype
        ntypetxt = 4
    else
        typetxt = '_t?'
        ntypetxt = 3
    end if

    ncomp_dum(:,:) = 0
    ncomp_plustracer_dum(:,:) = 0

    do isize = 1, nsize_aer(itype)
        n =isize

        if (isize .le. 9) then
            write(sizetxt,'(i1)') isize
            nsizetxt = 1
        else if (isize .le. 99) then
            write(sizetxt,'(i2)') isize
            nsizetxt = 2
        else if (isize .le. 999) then
            write(sizetxt,'(i3)') isize
            nsizetxt = 3
        else
            sizetxt = 's?'
            nsizetxt = 2
        end if


        do iphase = 1, nphase_aer

            if (iphase .eq. ai_phase) then
                phasetxt = 'a'
                nphasetxt = 1
            else if (iphase .eq. cw_phase) then
                phasetxt = 'cw'
                nphasetxt = 2
            else 
                phasetxt = 'p?'
                nphasetxt = 2
            end if

            do l_mastercomp = -2, ntot_mastercomp_aer

                iaddto_ncomp = 1
                iaddto_ncomp_plustracer = 1

                if (l_mastercomp .eq. -2) then
                    iaddto_ncomp = 0
                    iaddto_ncomp_plustracer = 0
                    lptr_dum = numptr_aer(n,itype,iphase)
                    mcindx_dum = -2
                    spectxt = 'num_'
                    nspectxt = 4
                
                else if (l_mastercomp .eq. -1) then
                    if (iphase .ne. ai_phase) cycle
                    iaddto_ncomp = 0
                    iaddto_ncomp_plustracer = 0
                    lptr_dum = waterptr_aer(n,itype)
                    mcindx_dum = -1
                    spectxt = 'water_'
                    nspectxt = 6
                else if (l_mastercomp .eq. 0) then
                    if (iphase .ne. ai_phase) cycle
                    iaddto_ncomp = 0
                    iaddto_ncomp_plustracer = 0
                    lptr_dum = hyswptr_aer(n,itype)
                    mcindx_dum = 0
                    spectxt = 'hysw_'
                    nspectxt = 5

                else if (l_mastercomp .eq. mastercompindx_seas_aer) then
                    lptr_dum = lptr_seas_aer(n,itype,iphase)
                    mcindx_dum = mastercompindx_seas_aer
                    spectxt = 'seas_'
                    nspectxt = 5

                else if (l_mastercomp .eq. mastercompindx_oc_aer) then
                    lptr_dum = lptr_oc_aer(n,itype,iphase)
                    mcindx_dum = mastercompindx_oc_aer
                    spectxt = 'oc_'
                    nspectxt = 3

                else if (l_mastercomp .eq. mastercompindx_dust_aer) then
                    lptr_dum = lptr_dust_aer(n,itype,iphase)
                    mcindx_dum = mastercompindx_dust_aer
                    spectxt = 'dust_'
                    nspectxt = 5
                else
                    ! goto 2500
                    cycle
                end if


                if (lptr_dum .gt. lmaxd) then
                    ! rce 2005-mar-14 - added check for lptr_dum > lmaxd
                    write( msg, '(a,3(1x,i4))' ) 'itype, isize, iphase =',   &
                    itype, isize, iphase
                    call chemistry_message( msg )
                    write( msg, '(a,3(1x,i4))' ) 'l_mastercomp, lptr_dum, lmaxd =',   &
                    l_mastercomp, lptr_dum, lmaxd
                    call chemistry_message( msg )
                    msg = '*** subr init_data_binaer_ptr error - lptr_dum > lmaxd'
                    call chemistry_error_fatal( msg )
                
                else if (lptr_dum .ge. p1st) then
                
                    ncomp_dum(isize,iphase) = ncomp_dum(isize,iphase) + iaddto_ncomp
                    ncomp_plustracer_dum(isize,iphase) =   &
                    ncomp_plustracer_dum(isize,iphase) + iaddto_ncomp_plustracer
                
                    name(lptr_dum) =   &
                    spectxt(1:nspectxt) // phasetxt(1:nphasetxt) //   &
                    sizetxt(1:nsizetxt) //  typetxt(1:ntypetxt)
                
                    if (l_mastercomp .eq. -2) then
                        !(numptr_aer is already set)
                        mprognum_aer(n,itype,iphase) = 1
                    
                    else if (l_mastercomp .eq. -1) then
                        !(waterptr_aer is already set)
                        continue
                    
                    else if (l_mastercomp .eq. 0) then
                        !(hyswptr_aer is already set)
                        continue
                    
                    else if (l_mastercomp .gt. 0) then
                        ll = ncomp_plustracer_dum(isize,iphase)
                        ! Note massptr_aer holds the species type l1, size bin n and is used at many places in movesect and binaer drydep code
                        massptr_aer(ll,n,itype,iphase) = lptr_dum
                        mastercompptr_aer(ll,itype) = mcindx_dum
                    
                        name_aer(ll,itype) = name_mastercomp_aer(mcindx_dum)
                        dens_aer(ll,itype) = dens_mastercomp_aer(mcindx_dum)
                        mw_aer(ll,itype) = mw_mastercomp_aer(mcindx_dum)
                        hygro_aer(ll,itype) = hygro_mastercomp_aer(mcindx_dum)
                    
                    end if
                end if ! if for lptr_dum
            end do    ! l_mastercomp = -1, ntot_mastercomp_aer 2500
        end do    ! iphase = 1, nphase_aer 2600
    end do    ! isize = 1, nsize_aer(itype) 2700


!   now set ncomp_aer and ncomp_plustracer_aer, 
!   *** and check that the values computed for each size and phase all match
    ncomp_aer(itype) = ncomp_dum(1,ai_phase)
    ncomp_plustracer_aer(itype) = ncomp_plustracer_dum(1,ai_phase)

    do iphase = 1, nphase_aer
        do isize = 1, nsize_aer(itype)
            if (ncomp_aer(itype) .ne. ncomp_dum(isize,iphase)) then
                msg =  '*** subr init_data_binaer_ptr - ' //   &
                'ncomp_aer .ne. ncomp_dum'
            call chemistry_message( msg )
            write(msg,9350) 'isize, itype, iphase =', isize, itype, iphase
            call chemistry_message( msg )
            write(msg,9350) 'ncomp_aer, ncomp_dum =',   &
                ncomp_aer(itype), ncomp_dum(isize,iphase)
            call chemistry_error_fatal( msg )
            end if
            if (ncomp_plustracer_aer(itype) .ne.   &
                ncomp_plustracer_dum(isize,iphase)) then
                msg = '*** subr init_data_binaer_ptr - ' //   &
                'ncomp_plustracer_aer .ne. ncomp_plustracer_dum'
            call chemistry_message( msg )
            write(msg,9350) 'isize, itype, iphase =', isize, itype, iphase
            call chemistry_message( msg )
            write(msg,9350)   &
                'ncomp_plustracer_aer, ncomp_plustracer_dum =',   &
                ncomp_plustracer_aer(itype),   &
                ncomp_plustracer_dum(isize,iphase)
            call chemistry_error_fatal( msg )
            end if
        end do
    end do


end do    ! itype = 1, ntype_aer 2800


9320    format( a, i1, i1, a, 8x )

!
!   output wrfch pointers
!

!jwfeng: 03/24/2020 
!modify to a nicer output format
9350    format( a, 32(1x,i4) )
    ! msg = ' '
    ! call chemistry_message( msg )
    ! msg = 'output from subr init_data_binaer_ptr'
    ! call chemistry_message( msg )
    write(msg,9350) '    nphase_aer = ', nphase_aer
    call chemistry_message( msg )

    do iphase=1,nphase_aer
    write(msg,'(a,i2)') '    Information of phase', iphase
    call chemistry_message( msg )
    msg = '    -------------------------------------------------------------------'
    call chemistry_message( msg )
    write(msg,9350) '            ntype_aer = ', ntype_aer
    call chemistry_message( msg )
    write(msg,9350) '            ncomp_aer = ', ncomp_aer
    call chemistry_message( msg )
      
    do itype=1,ntype_aer

    write(msg,9350) '                itype = ', itype
    call chemistry_message( msg )

    write(msg,9350) '            nsize_aer = ', nsize_aer(itype)
    call chemistry_message( msg )

    write(msg,9350) '        lptr_seas_aer = ',   &
        (lptr_seas_aer(n,itype,iphase), n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '          lptr_oc_aer = ',   &
        (lptr_oc_aer(n,itype,iphase),   n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '        lptr_dust_aer = ',   &
        (lptr_dust_aer(n,itype,iphase), n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '          hyswptr_aer = ',   &
        (hyswptr_aer(n,itype), n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '         waterptr_aer = ',   &
        (waterptr_aer(n,itype), n=1,nsize_aer(itype))
    call chemistry_message( msg )

    write(msg,9350) '           numptr_aer = ',   &
        (numptr_aer(n,itype,iphase), n=1,nsize_aer(itype))
    call chemistry_message( msg )

9352    format( 5a )

    do ll = 1, ncomp_plustracer_aer(itype)
    write(msg,9350) '          massptr_aer = ',   &
        (massptr_aer(ll,n,itype,iphase), n=1,nsize_aer(itype)),ll
        call chemistry_message( msg )


    end do
    end do ! type
    msg = '    -------------------------------------------------------------------'
    call chemistry_message( msg )
    end do ! phase

!
!   check aerosol species pointers for "validity"
!
    do iphase=1,nphase_aer
    do itype=1,ntype_aer

    y_seas = 0
    y_oc = 0
    y_dust = 0
    y_hysw = 0
    y_water = 0
    y_num = 0

    do n = 1, nsize_aer(itype)
        if (lptr_seas_aer(n,itype,iphase) .ge. p1st) y_seas = y_seas + 1
        if (lptr_oc_aer(n,itype,iphase)   .ge. p1st) y_oc   = y_oc + 1
        if (lptr_dust_aer(n,itype,iphase) .ge. p1st) y_dust = y_dust + 1
        if (hyswptr_aer(n,itype)          .ge. p1st) y_hysw = y_hysw + 1
        if (waterptr_aer(n,itype)         .ge. p1st) y_water = y_water + 1
        if (numptr_aer(n,itype,iphase)    .ge. p1st) y_num = y_num + 1
    end do

!   these must be defined for all aerosol bins

    if (y_dust .ne. nsize_aer(itype)) then
        msg = '*** subr init_data_binaer_ptr - ptr error for dust'
        call chemistry_message( msg )
        write(msg,9350) 'phase, type=', iphase,itype
        call chemistry_error_fatal( msg )
    else if (y_hysw .ne. nsize_aer(itype)) then
        msg = '*** subr init_data_binaer_ptr - ptr error for hysw'
        call chemistry_message( msg )
        write(msg,9350) 'phase, type=', iphase,itype
        call chemistry_error_fatal( msg )
    else if (y_water .ne. nsize_aer(itype)) then
        msg = '*** subr init_data_binaer_ptr - ptr error for water'
        call chemistry_message( msg )
        write(msg,9350) 'phase, type=', iphase,itype
        call chemistry_error_fatal( msg )
    else if (y_num .ne. nsize_aer(itype)) then
        msg = '*** subr init_data_binaer_ptr - ptr error for num'
        call chemistry_message( msg )
        write(msg,9350) 'phase, type=', iphase,itype
        call chemistry_error_fatal( msg )
    end if


!   these must be defined for all aerosol bins
!       or else undefined for all aerosol bins
!    if ((y_no3 .ne. 0) .and. (y_no3 .ne. nsize_aer(itype))) then
!        msg = '*** subr init_data_binaer_ptr - ptr error for no3'
!        call chemistry_message( msg )
!        write(msg,9350) 'phase, type=', iphase,itype
!        call chemistry_error_fatal( msg )
!   else if ((y_dust .ne. 0) .and. (y_dust .ne. nsize_aer(itype))) then
!   if ((y_dust .ne. 0) .and. (y_dust .ne. nsize_aer(itype))) then
!       msg = '*** subr init_data_binaer_ptr - ptr error for dust'
!       call chemistry_message( msg )
!       write(msg,9350) 'phase, type=', iphase,itype
!       call chemistry_error_fatal( msg )
!   end if

    enddo ! type
    enddo ! phase
!
!   calc ltot, ltot2, kh2o, ktemp
!
    is_aerosol(:) = .false.
    ltot = 0
    do iphase=1,nphase_aer
        do itype=1,ntype_aer
        do n = 1, nsize_aer(itype)
            do ll = 1, ncomp_plustracer_aer(itype)
               ltot = max( ltot, massptr_aer(ll,n,itype,iphase) )
               is_aerosol(massptr_aer(ll,n,itype,iphase))=.true.
            end do
            ltot = max( ltot, hyswptr_aer(n,itype) )
            ltot = max( ltot, waterptr_aer(n,itype) )
            ltot = max( ltot, numptr_aer(n,itype,iphase) )
            l = hyswptr_aer(n,itype)
            if (l .ge. p1st) is_aerosol(l)=.true.
            l = waterptr_aer(n,itype)
            if (l .ge. p1st) is_aerosol(l)=.true.
            l = numptr_aer(n,itype,iphase)
            if (l .ge. p1st) is_aerosol(l)=.true.
        end do
        end do
    end do


    kh2o = ltot + 1
    ktemp = ltot + 2
    ltot2 = ktemp

    write( msg, '(a,4(1x,i4))' ) 'ltot, ltot2, lmaxd, l2maxd =',   &
        ltot, ltot2, lmaxd, l2maxd
    call chemistry_message( msg )
    if ((ltot .gt. lmaxd) .or. (ltot2 .gt. l2maxd)) then
        msg = '*** subr init_data_binaer_ptr - ltot/ltot2 too big'
        call chemistry_error_fatal( msg )
    end if

    name(ktemp)  = 'temp'
    name(kh2o)   = 'h2o'

    call initwet(   &
    ntype_aer, nsize_aer, ncomp_aer,   &
    massptr_aer, dens_aer, numptr_aer,           &
    maxd_acomp, maxd_asize,maxd_atype, maxd_aphase,   &
    dcen_sect, sigmag_aer, &
    waterptr_aer, dens_water_aer, &
    scavimptblvol, scavimptblnum, nimptblgrow_mind,   &
    nimptblgrow_maxd, dlndg_nimptblgrow)

    return
end subroutine init_data_binaer_ptr

end module mpas_atmchem_binaer_init

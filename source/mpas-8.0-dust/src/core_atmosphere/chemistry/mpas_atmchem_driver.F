!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! General driver for all chemistry procecces.
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
!     The code structure is similar to the physics drivers
!     in MPAS-v5.2:src/core_atmosphere/physics/mpas_atmphys_driver.F
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

module mpas_atmchem_driver
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines
 use mpas_log, only : mpas_log_write

 use mpas_atmchem_driver_optical
 use mpas_atmchem_driver_emission
 use mpas_atmchem_driver_drydep
 use mpas_atmchem_driver_convection
 use mpas_atmchem_driver_wetscav
 use mpas_atmchem_interface
 use mpas_atmchem_vars
 use mpas_atmchem_binaer_init
 use mpas_atmchem_utilities
 use mpas_timer

 implicit none
 private
 public:: chemistry_driver,driver_chem_diag


!MPAS top chemistry driver.

! subroutine chemistry_driver is the top chemistry driver from which separate drivers.

 contains

!=================================================================================================================
 subroutine driver_chem_diag(domain)
!=================================================================================================================
 use module_data_binaer_asect

!inout arguments:
 type(domain_type),intent(inout):: domain

!local arguments:
 type(mpas_pool_type),pointer :: mesh
 type(mpas_pool_type),pointer :: state
 type(mpas_pool_type),pointer :: diag_chemistry
 type(block_type),pointer:: block

 integer:: thread
 integer,pointer:: nThreads, nVertLevels
 integer,dimension(:),pointer:: cellSolveThreadStart, cellSolveThreadEnd
 integer :: its,ite

 integer:: i,k,n
 integer:: l_dust,l_seas
 integer:: iphase,itype
 real(kind=RKIND),dimension(:,:),pointer :: pm10_dust, pm10_seas
 real(kind=RKIND),dimension(:,:,:),pointer :: chem

 itype = 1
 iphase = ai_phase

block => domain % blocklist
do while(associated(block))
    call mpas_pool_get_subpool(block%structs,'mesh'          ,mesh          )
    call mpas_pool_get_subpool(block%structs,'state'         ,state         )
    call mpas_pool_get_subpool(block%structs,'diag_chemistry',diag_chemistry)

    call mpas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)
    call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

    call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
    call mpas_pool_get_array(diag_chemistry,'pm10_dust',pm10_dust)
    call mpas_pool_get_array(diag_chemistry,'pm10_seas',pm10_seas)
    call mpas_pool_get_array(state,'chem',chem)

!jwfeng: Please CHECK if I have used OMP right
    do thread=1,nThreads
        its = cellSolveThreadStart(thread)
        ite = cellSolveThreadEnd(thread)
        do k = 1, nVertLevels
        do i = its,ite
           pm10_dust(k,i) = 0.0_RKIND
           pm10_seas(k,i) = 0.0_RKIND
        enddo
        enddo

        do k = 1, nVertLevels
        do i = its,ite
        do n = 1, nsize_aer(itype)
           l_dust = lptr_dust_aer(n,itype,iphase)
           if (l_dust > param_first_scalar .and. dhi_sect(n,itype) < 10.001e-4 ) then
               pm10_dust(k,i) = chem(l_dust,k,i) + pm10_dust(k,i)
           end if

           l_seas = lptr_seas_aer(n,itype,iphase)
           if (l_seas > param_first_scalar .and. dhi_sect(n,itype) < 10.001e-4 ) then
               pm10_seas(k,i) = chem(l_seas,k,i) + pm10_seas(k,i)
           end if
        enddo
        enddo
        enddo
    end do

    block => block % next
end do

 end subroutine driver_chem_diag

!=================================================================================================================
 subroutine chemistry_driver(domain,itimestep,xtime_d,xtime_s)
!=================================================================================================================

!input arguments:
 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: xtime_s,xtime_d

!inout arguments:
 type(domain_type),intent(inout):: domain

!local pointers:
 type(mpas_pool_type),pointer::  configs,        &
                                 mesh,           &
                                 state,          &
                                 sfc_input,      &
                                 diag,           &
                                 diag_physics,   &
                                 diag_chemistry, &
								 tend_chemistry

!local variables:
 type(block_type),pointer:: block
 type(dm_info),pointer:: dminfo

 integer:: time_lev
 integer:: thread

 integer,pointer:: nThreads
 integer,dimension(:),pointer:: cellSolveThreadStart, cellSolveThreadEnd
 logical,pointer:: config_drydep_on, config_wetscav_on, debug_on
 real,pointer:: config_photdt, config_dt
 logical :: do_photstep
 character(len=StrKIND),pointer:: config_dust_opt, config_aer_op_opt, config_ctrans_opt
!=================================================================================================================

call mpas_timer_start('chemistry driver')
call chemistry_message('----------------------------------------------------------------------')
call chemistry_message('----------------------Starting chemistry driver-----------------------')
call chemistry_message('----------------------------------------------------------------------')

block => domain % blocklist
do while(associated(block))

    call mpas_pool_get_subpool(block%structs,'mesh'          ,mesh          )
    call mpas_pool_get_subpool(block%structs,'state'         ,state         )
    call mpas_pool_get_subpool(block%structs,'diag'          ,diag          )
    call mpas_pool_get_subpool(block%structs,'diag_physics'  ,diag_physics  )
    call mpas_pool_get_subpool(block%structs,'diag_chemistry',diag_chemistry)
    call mpas_pool_get_subpool(block%structs,'sfc_input'     ,sfc_input     )
    call mpas_pool_get_subpool(block%structs,'tend_chemistry',tend_chemistry)

    call mpas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)

    call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(block % dimensions, 'cellSolveThreadEnd', cellSolveThreadEnd)

!	jwfeng from USTC AEMOL:  02/2020
!	Preparation for chemistry_driver

!	jwfeng:  03/2020
    call mpas_pool_get_config(block%configs,'config_print_debug_info' ,debug_on )
    call mpas_pool_get_config(block%configs,'config_aer_op_opt' ,config_aer_op_opt )
    call mpas_pool_get_config(block%configs,'config_drydep_on'  ,config_drydep_on  )
    call mpas_pool_get_config(block%configs,'config_ctrans_opt' ,config_ctrans_opt )
    call mpas_pool_get_config(block%configs,'config_wetscav_on' ,config_wetscav_on )
    call mpas_pool_get_config(block%configs,'config_dust_opt'   ,config_dust_opt   )
    call mpas_pool_get_config(block%configs,'config_photdt'     ,config_photdt     )
    call mpas_pool_get_config(block%configs,'config_dt'         ,config_dt         )

!-----------------------------------------------------------------------------------------!
! Allocate arrays shared by all chemistry parameterizations:
    call allocate_forall_chemistry
    if(debug_on) call chemistry_message('forall chemistry allocate success')
    time_lev = 1
!$OMP PARALLEL DO 
    do thread=1,nThreads
       call chemistry_from_MPAS(mesh,state,time_lev,sfc_input,diag,diag_physics,        &
                                cellSolveThreadStart(thread), cellSolveThreadEnd(thread))
    end do
!$OMP END PARALLEL DO 
!-----------------------------------------------------------------------------------------!
! Optical start 
    do_photstep = .false.
    if((xtime_d .lt. 1.e-12) .and. (xtime_s .lt. 1.e-12) ) then  !First step always do photstep
       do_photstep = .true.
    else
       if ( (config_photdt<=0) .or. &
            ( &
              ( xtime_s+real(config_dt,8)+0.01 ) >= &
              ( int( xtime_s/real(config_photdt,8)+1,8 )*real(config_photdt,8) ) &
            ) &
          ) then
          do_photstep = .true.
       endif
    endif
    if(do_photstep) then
        call mpas_log_write('--- time to run the dust optical driver config_aer_op_opt   = '//trim(config_aer_op_opt))
        call mpas_timer_start('optical driver')

        call allocate_optical
        if(debug_on) call chemistry_message('allocate optical success')
        do thread=1,nThreads
           call driver_optical(diag, diag_chemistry, block%configs, &
                               cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
        end do
        if(debug_on) call chemistry_message('optical driver success')
        call deallocate_optical
        if(debug_on) call chemistry_message('deallocate optical success')
! Optical end 
        call mpas_timer_stop('optical driver')
    else
        call mpas_log_write('--- time to run the dust optical driver config_aer_op_opt   = Dont do photstep at this timestep')
    endif
!-----------------------------------------------------------------------------------------!

!-----------------------------------------------------------------------------------------!
! Emisssion start 
    call mpas_log_write('--- time to run the dust emission driver config_dust_opt    = '//trim(config_dust_opt))
    call mpas_timer_start('emission driver')

dminfo => domain % dminfo

    call allocate_emission
    if(debug_on) call chemistry_message('allocate emission success')
!$OMP PARALLEL DO 
    do thread=1,nThreads
       call driver_emission(itimestep,dminfo,block%configs,mesh,state,sfc_input         , &
                            diag,diag_physics,diag_chemistry,tend_chemistry      , &
                            cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
    end do
!$OMP END PARALLEL DO 
    if(debug_on) call chemistry_message('emission driver success')
    call deallocate_emission
    if(debug_on) call chemistry_message('deallocate emission success')
! Emisssion end 
    call mpas_timer_stop('emission driver')
!-----------------------------------------------------------------------------------------!

!-----------------------------------------------------------------------------------------!
!   call mpas_timer_start('dry deposition')
if (config_drydep_on) then
    call mpas_log_write('--- time to run the dry deposition driver config_drydep_on  = $l',logicArgs = (/config_drydep_on/))
! Dry deposition start 
    call allocate_drydep
    if(debug_on) call chemistry_message('allocate drydep success')
!$OMP PARALLEL DO 
    do thread=1,nThreads
       call driver_drydep(itimestep,block%configs,mesh,state,                    &
                          sfc_input,diag,diag_physics,tend_chemistry,            &
                          cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
    end do
!$OMP END PARALLEL DO 
    if(debug_on) call chemistry_message('drydep driver success')
    call deallocate_drydep
    if(debug_on) call chemistry_message('deallocate drydep success')
! Dry deposition end 
else
    call mpas_log_write('--- time to run the dry deposition driver config_drydep_on  = $l',logicArgs = (/config_drydep_on/))
end if
!   call mpas_timer_stop('dry deposition')

!-----------------------------------------------------------------------------------------!
    call mpas_timer_start('chem convection')
       if(debug_on) call chemistry_message('05')
if (trim(config_ctrans_opt) .ne. 'no_ctrans') then
    call mpas_log_write('--- time to run the chem ctrans driver config_ctrans_opt    = '//trim(config_ctrans_opt))
! Chem convection start
    call allocate_convection
    if(debug_on) call chemistry_message('allocate ctrans success')
!$OMP PARALLEL DO 
    do thread=1,nThreads
       call driver_convection( &
            state, diag_physics, time_lev, &
            block%configs, diag_chemistry, tend_chemistry, &
            cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
    end do
!$OMP END PARALLEL DO 
    if(debug_on) call chemistry_message('ctrans driver success')
    call deallocate_convection
    if(debug_on) call chemistry_message('deallocate ctrans success')
! Chem convection end
else
    call mpas_log_write('--- time to run the chem ctrans driver config_ctrans_opt    = '//trim(config_ctrans_opt))
end if
    call mpas_timer_stop('chem convection')

!-----------------------------------------------------------------------------------------!
       if(debug_on) call chemistry_message('06')
    call mpas_timer_start('wet scavenging')
if (config_wetscav_on) then
    call mpas_log_write('--- time to run the wet scavenging driver config_wetscav_on = $l',logicArgs = (/config_wetscav_on/))
! wet scavenging start 
    call allocate_wetscav
    if(debug_on) call chemistry_message('allocate wetscav success')
!$OMP PARALLEL DO 
    do thread=1,nThreads
       call driver_wetscav(state, diag_physics, time_lev, block%configs, tend_chemistry, &
                          cellSolveThreadStart(thread),cellSolveThreadEnd(thread))
    if(debug_on) call chemistry_message('wetscav driver success')
    end do
!$OMP END PARALLEL DO 
    call deallocate_wetscav
    if(debug_on) call chemistry_message('deallocate wetscav success')
! wet scavenging end 
else
    call mpas_log_write('--- time to run the wet scavenging driver config_wetscav_on = $l',logicArgs = (/config_wetscav_on/))
end if
    call mpas_timer_stop('wet scavenging')



!-----------------------------------------------------------------------------------------!

! Deallocate arrays shared by all chemistry parameterizations:
!-----------------------------------------------------------------------------------------!
    call deallocate_forall_chemistry
    if(debug_on) call chemistry_message('forall chemistry deallocate success')

    block => block % next
end do 

call mpas_timer_stop('chemistry driver')

 end subroutine chemistry_driver

!=================================================================================================================
 end module mpas_atmchem_driver
!=================================================================================================================

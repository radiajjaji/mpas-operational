!************************************************************************
!
! Module to initialize varaiables for gocart dust scheme.
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
!   WRF-Chem-v3.5:chem/module_gocart_dust_init.F
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************
#define DM_BCAST_CHAR(A)     call mpas_dmpar_bcast_char(dminfo,A)
#define DM_BCAST_INTEGER(A)  call mpas_dmpar_bcast_int(dminfo,A)
#define DM_BCAST_INTEGERS(A) call mpas_dmpar_bcast_ints(dminfo,size(A),A)
#define DM_BCAST_REAL(A)     call mpas_dmpar_bcast_real(dminfo,A)
#define DM_BCAST_REALS(A)    call mpas_dmpar_bcast_reals(dminfo,size(A),A)

Module mpas_atmchem_gocart_dust_init

    use module_data_gocart_dust
    use mpas_kind_types
    use mpas_derived_types
    use mpas_pool_routines
    use mpas_atmchem_utilities
    use mpas_dmpar


    implicit none
    private
    public:: gocart_dust_init_forMPAS

    logical,public,save :: has_done_gocart_init = .false.         !Greenwich mean time hour of model start (hr)

    contains

    subroutine gocart_dust_init_forMPAS(dminfo, configs)
        type(dm_info), intent(in):: dminfo
        integer :: istat,lc,slcats,iindex,lumatch=1
        character(len=StrKIND),pointer:: mminsl
        character(4) :: sltype
        type(mpas_pool_type),intent(in):: configs
        character*128:: mess , message
        integer,parameter:: open_ok = 0

        call mpas_pool_get_config(configs,'input_soil_data',mminsl)

! read in soil properties from soilparm.tbl:

        if(dminfo % my_proc_id == IO_NODE ) then


           open(16,file='SOILPARM.TBL',form='FORMATTED',status='OLD',iostat=istat)
           if(istat /= open_ok) &
              call chemistry_error_fatal('mpas_atmchem_gocart_dust_init.F: gocart_dust_init_forMPAS: ' // &
                                       'failure opening SOILPARM.TBL' )
       
           write(mess,*) '   input soil texture classification = ', trim (mminsl)
           call chemistry_message(mess)
       
           lumatch=0
           read(16,*)
           read(16,2000,end=2003) sltype
           2000 format(a4)
           read(16,*)slcats,iindex
           if(sltype.eq.mminsl)then
              write(mess,*) '   soil texture classification = ', trim ( sltype ) , ' found', &
                         slcats,' categories'
              call chemistry_message ( mess )
              lumatch=1
           endif
       
!prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008:

           if(size(bb) < slcats .or. &
              size(drysmc) < slcats .or. &
              size(f11   ) < slcats .or. &
              size(maxsmc) < slcats .or. &
              size(refsmc) < slcats .or. &
              size(satpsi) < slcats .or. &
              size(satdk ) < slcats .or. &
              size(satdw ) < slcats .or. &
              size(wltsmc) < slcats .or. &
              size(qtz   ) < slcats) then
              call wrf_error_fatal('table sizes too small for value of slcats in mpas_atmchem_gocart_dust_init.F')
           endif
       
           if(sltype.eq.mminsl) then
              do lc = 1, slcats
                 read(16,*) iindex,bb(lc),drysmc(lc),f11(lc),maxsmc(lc),refsmc(lc),satpsi(lc), &
                             satdk(lc),satdw(lc),wltsmc(lc),qtz(lc)
              enddo
           endif
       
           2003 continue
           close(16)


           if(lumatch.eq.0)then
              call chemistry_message( 'soil texture in input file does not ' )
              call chemistry_message( 'match soilparm table'                 )
              call chemistry_error_fatal( 'inconsistent or missing soilparm file' )
           endif
        endif
       
!distribute data to all processors:

        DM_BCAST_INTEGER(lumatch)
        DM_BCAST_CHAR(sltype)
        DM_BCAST_CHAR(mminsl)
        DM_BCAST_INTEGER(slcats)
        DM_BCAST_INTEGER(iindex)
        DM_BCAST_REALS(bb)
        DM_BCAST_REALS(drysmc)
        DM_BCAST_REALS(f11)
        DM_BCAST_REALS(maxsmc)
        DM_BCAST_REALS(refsmc)
        DM_BCAST_REALS(satpsi)
        DM_BCAST_REALS(satdk)
        DM_BCAST_REALS(satdw)
        DM_BCAST_REALS(wltsmc)
        DM_BCAST_REALS(qtz)

       
        call mpas_log_write(' LUMATCH=$i',intArgs=(/lumatch/))
        call mpas_log_write(' SLTYPE ='//trim(sltype))
        call mpas_log_write(' MMINSL ='//mminsl)
        call mpas_log_write(' SLCATS =$i',intArgs=(/slcats/))
        call mpas_log_write(' IINDEX =$i',intArgs=(/iindex/))
        call mpas_log_write('')

        call mpas_log_write('    end read SOILPARM.TBL')
       
    end subroutine

END Module mpas_atmchem_gocart_dust_init

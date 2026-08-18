










      SUBROUTINE INITPOST
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    INITPOST    INITIALIZE POST FOR RUN
!   PRGRMMR: RUSS TREADON    ORG: W/NP2      DATE: 93-11-10
!     
! ABSTRACT:  THIS ROUTINE INITIALIZES CONSTANTS AND
!   VARIABLES AT THE START OF AN ETA MODEL OR POST 
!   PROCESSOR RUN.
!
!   THIS ROUTINE ASSUMES THAT INTEGERS AND REALS ARE THE SAME SIZE
!   .     
!     
! PROGRAM HISTORY LOG:
!   93-11-10  RUSS TREADON - ADDED DOCBLOC
!   98-05-29  BLACK - CONVERSION OF POST CODE FROM 1-D TO 2-D
!   99-01 20  TUCCILLO - MPI VERSION
!   01-10-25  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!   02-06-19  MIKE BALDWIN - WRF VERSION
!   02-08-15  H CHUANG - UNIT CORRECTION AND GENERALIZE PROJECTION OPTIONS
!     
! USAGE:    CALL INIT
!   INPUT ARGUMENT LIST:
!     NONE     
!
!   OUTPUT ARGUMENT LIST: 
!     NONE
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       NONE
!     LIBRARY:
!       COMMON   - CTLBLK
!                  LOOKUP
!                  SOILDEPTH
!
!    
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : CRAY C-90
!$$$  
      use vrbls4d, only: dust, smoke
      use vrbls3d, only: t, u, uh, v, vh, wh, q, pmid, t, omga, pint, alpint,plume,&
              qqr, qqs, qqg, qqnw, qqni,qqnr, cwm, qqw, qqi, qqr, qqs, extcof55,&
              f_ice, f_rain, f_rimef, q2, zint, zmid, ttnd, cfr, cfr_raw, REF_10CM,&
              qqnwfa,qqnifa,taod5503d,aextc55,DUEM,DUSD,DUDP,DUWT,TOT_DUST
      use vrbls2d, only: tmax, qrmax, htop, hbot, cuppt, fis, cfrach, cfracl,      &
              sr, cfrach, cfracm, wspd10max, w_up_max, w_dn_max, w_mean, refd_max, &
              up_heli_max, up_heli_max16, grpl_max, up_heli, up_heli16,            &
              up_heli_min,up_heli_min16,up_heli_max02,up_heli_min02,               &
              up_heli_max03,up_heli_min03,rel_vort_max,rel_vort_max01,             &
              wspd10umax,wspd10vmax,refdm10c_max,hail_max2d,hail_maxk1,ltg1_max,   &
              ltg2_max, ltg3_max, nci_ltg, nca_ltg, nci_wq, nca_wq, nci_refd,      &
              u10, v10, th10, q10, tshltr, mrshltr,dusmass,                        &
              nca_refd, qv2m, qshltr, smstav, smstot, ssroff, bgroff, sfcevp,      &
              sfcexc, vegfrc, acsnow, cmc, sst, thz0, qz0, uz0, vz0, qs, qvg,      &
              z0, ustar, akhs, akms, radot, ths, acsnom, cuprec, ancprc, acprec,   &
              rainc_bucket, pcp_bucket, cprate, prec, snownc, snow_bucket,         &
              graup_bucket, swddni, swddif, mean_frp, acgraup, acfrain,            &
              graupelnc, albedo, rswin, rswout, swdnbc, swddnic,                   &
              swddifc, swupbc, swupt, czen, czmean, rlwin, lwdnbc, lwupbc,         &
              rainnc_bucket, taod5502d, aerasy2d, aerssa2d, lwp, iwp,              &
              sigt4, rlwtoa, rswinc, aswin, aswout, alwin, alwout, alwtoa, aswtoa, &
              tg, soiltb, twbs, qwbs,grnflx, sfcshx, sfclhx, subshx, snopcx,       &
              sfcuvx, potevp, ncfrcv, ncfrst, sno, si, pctsno, snonc, tsnow,       &
              ivgtyp, isltyp, islope, pblh, pblhgust, f,                           &
              QVl1,REFC_10CM,REF1KM_10CM,REF4KM_10CM,                              &
              SWRADmean,U10mean,V10mean,SPDUV10mean,SWNORMmean,SNFDEN,SNDEPAC,     &
              hbotd,hbots,htopd,htops,tot_edust, epsr,VIS_DUST
      use soil, only: smc, sh2o, stc, sldpth, sllevel, erod
      use masks, only: lmv, lmh, vtm, sice, gdlat, gdlon, sm, dx, dy, htm
      use ctlblk_mod, only: jsta_2l, jend_2u, filename, datahandle, datestr,       &
              ihrst, imin, idat, sdat, ifhr, ifmin, imp_physics, jsta, jend,       &
              spval,gdsdegr, modelname, pt, icu_physics, jsta_m, jend_m, nsoil,    &
              isf_surface_physics, nsoil, ardlw, ardsw, asrfc, me, mpi_comm_comp,  &
              nphs, smflag, spl, lsm, dt, prec_acc_dt, dtq2, tsrfc, trdlw,         &
              trdsw, theat, tclod, tprec, nprec, alsl, im, jm, lm, grib
      use params_mod, only: capa, g, rd, d608, tfrz, ad05, cft0, stbol,            &
              p1000, pi, rtd, lheat, dtr, erad
      use lookup_mod, only: thl, plq, ptbl, ttbl, rdq, rdth, rdp, rdthe, pl,       &
              qs0, sqs, sthe, the0, ttblq, rdpq, rdtheq, stheq, the0q
      use gridspec_mod, only: gridtype, dxval, latstart, latlast, lonstart,        &
              lonlast, dyval, cenlat, cenlon, maptype, truelat1, truelat2,         &
              standlon, psmapf
      use wrf_io_flags_mod, only:
!- - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - -
      implicit none
!
!     INCLUDE/SET PARAMETERS.
!     
      INCLUDE "mpif.h"
!
! This version of INITPOST shows how to initialize, open, read from, and
! close a NetCDF dataset. In order to change it to read an internal (binary)
! dataset, do a global replacement of _ncd_ with _int_. 

      character(len=31) :: VarName
      integer :: Status
      character startdate*19,SysDepInfo*80
! 
!     NOTE: SOME INTEGER VARIABLES ARE READ INTO DUMMY ( A REAL ). THIS IS OK
!     AS LONG AS REALS AND INTEGERS ARE THE SAME SIZE.
!
!     ALSO, EXTRACT IS CALLED WITH DUMMY ( A REAL ) EVEN WHEN THE NUMBERS ARE
!     INTEGERS - THIS IS OK AS LONG AS INTEGERS AND REALS ARE THE SAME SIZE.

      INTEGER IDATE(8),JDATE(8)
!
!     DECLARE VARIABLES.
!     
      REAL RINC(5)
      REAL DUMMY ( IM, JM )
      REAL DUMMY2 ( IM, JM ), MSFT(IM,JM)
      INTEGER IDUMMY ( IM, JM )
      REAL,allocatable ::  DUM3D ( :, :, : )
        real, allocatable::  pvapor(:,:)
        real, allocatable::  pvapor_orig(:,:)      
      REAL, ALLOCATABLE :: THV(:,:,:)
 
      integer js,je,jev,iyear,imn,iday,itmp,ioutcount,istatus,          &
              ii,jj,ll,i,j,l,nrdlw,nrdsw,n,igdout,irtn,idyvald,         &
              idxvald,NSRFC , lflip, k, k1
      real DZ,TSPH,TMP,QMEAN,PVAPORNEW,DUMCST,TLMH,RHO,ZSF,ZPBLTOP
      real t2,th2,x2m,p2m,tsk, fact, temp
      real LAT

      integer jdn, numr, ic, jc, ierr
      integer, external :: iw3jdn
      real sun_zenith,sun_azimuth, ptop_low, ptop_mid, ptop_high
      real watericetotal, cloud_def_p, radius
      real totcount, cloudcount
      real delta_theta4gust

      real(kind=8), parameter :: alpha_d = 0.6d0   ! m2/kg
      real(kind=8), parameter :: H_dust  = 1500.d0 ! m
      real(kind=8), parameter :: eps     = 1.0d-8

      real(kind=8) :: dp, beta_ext
      integer :: lbot
      real(8) :: dustlev, dust10max

!
!
!***********************************************************************
!     START INIT HERE.
!
      ALLOCATE ( THV(IM,JSTA_2L:JEND_2U,LM) )
      ALLOCATE ( DUM3D ( IM+1, JM+1, LM+1 ) )
      WRITE(6,*)'INITPOST:  ENTER INITPOST'
!     
      gridtype='A'
      hbotd=0
      hbots=0
      htopd=0
      htops=0
      ttnd=0
      qs=0
!     STEP 1.  READ MODEL OUTPUT FILE
!
! LMH always = LM for sigma-type vert coord
! LMV always = LM for sigma-type vert coord

      do j = jsta_2l, jend_2u
        do i = 1, im
            LMV ( i, j ) = lm
            LMH ( i, j ) = lm
        end do
      end do

! HTM VTM all 1 for sigma-type vert coord

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            HTM ( i, j, l ) = 1.0
            VTM ( i, j, l ) = 1.0
        end do
       end do
      end do
!
      call ext_ncd_ioinit(SysDepInfo,Status)
      call ext_ncd_open_for_read( trim(fileName), 0, 0, " ",        &  
           DataHandle, Status)
      if ( Status /= 0 ) then
         print*,'error opening ',fileName, ' Status = ', Status ; stop
      endif
      print *,'DateStr before calling ext_ncd_get_next_time=',DateStr

!  The end j row is going to be jend_2u for all variables except for V.
      JS=JSTA_2L
      JE=JEND_2U
      IF (JEND_2U.EQ.JM) THEN
       JEV=JEND_2U+1
      ELSE
       JEV=JEND_2U
      ENDIF
!
! Getting start time
      call ext_ncd_get_dom_ti_char(DataHandle,'SIMULATION_START_DATE',  &
        startdate, status )
      jdate=0
      idate=0
      read(startdate,15)iyear,imn,iday,ihrst,imin      
 15   format(i4,1x,i2,1x,i2,1x,i2,1x,i2)
      print*,'processing yr mo day hr min='                             &
        ,idat(3),idat(1),idat(2),idat(4),idat(5)
      idate(1)=iyear
      idate(2)=imn
      idate(3)=iday
      idate(5)=ihrst
      idate(6)=imin
      SDAT(1)=imn
      SDAT(2)=iday
      SDAT(3)=iyear
      jdate(1)=idat(3)
      jdate(2)=idat(1)
      jdate(3)=idat(2)
      jdate(5)=idat(4)
      jdate(6)=idat(5)
      CALL W3DIFDAT(JDATE,IDATE,0,RINC)
      ifhr=nint(rinc(2)+rinc(1)*24.)
      ifmin=nint(rinc(3))

!  OK, since all of the variables are dimensioned/allocated to be
!  the same size, this means we have to be careful int getVariable
!  to not try to get too much data.  For example, 
!  DUM3D is dimensioned IM+1,JM+1,LM+1 but there might actually
!  only be im,jm,lm points of data available for a particular variable.  

      call ext_ncd_get_dom_ti_integer(DataHandle,'MP_PHYSICS'          &
       ,itmp,1,ioutcount,istatus)
      imp_physics=itmp

! Initializes constants for Ferrier microphysics       

      if(imp_physics==5 .or. imp_physics==85 .or. imp_physics==95)then
       CALL MICROINIT(imp_physics)
      end if
      
      call ext_ncd_get_dom_ti_integer(DataHandle,'CU_PHYSICS'          &
       ,itmp,1,ioutcount,istatus)
        icu_physics=itmp
	
      call ext_ncd_get_dom_ti_integer(DataHandle,'SF_SURFACE_PHYSICS'          &
       ,itmp,1,ioutcount,istatus)
        isf_surface_physics=itmp
	
! get 3-D variables
      ii=im/2
      jj=(jsta+jend)/2
      ll=lm

! Read Temperature (theta)

      VarName='T'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,       &
        IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            t ( i, j, l ) = dum3d ( i, j, l ) 
        end do
       end do
      end do

! Read Zonal Wind

      do l=1,lm
      end do
      VarName='U'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,         &
        IM+1,1,JM+1,LM+1,IM+1,JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im+1
            u ( i, j, l ) = dum3d ( i, j, l )
        end do
       end do
!  fill up UH which is U at P-points including 2 row halo
       do j = jsta_2l, jend_2u
        do i = 1, im
            UH (I,J,L) = (dum3d(I,J,L)+dum3d(I+1,J,L))*0.5
        end do
       end do
      end do

! Read Meridional Wind

      VarName='V'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,         &
        IM+1,1,JM+1,LM+1,IM, JS,JEV,LM)
      do l = 1, lm
       do j = jsta_2l, jev
        do i = 1, im
            v ( i, j, l ) = dum3d ( i, j, l )
        end do
       end do
!  fill up VH which is V at P-points including 2 row halo
       do j = jsta_2l, jend_2u
        do i = 1, im
          VH(I,J,L) = (dum3d(I,J,L)+dum3d(I,J+1,L))*0.5
        end do
       end do
      end do

! Read Vertical velocity

      VarName='W'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,         &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM+1)
!     do l = 1, lm+1
!      do j = jsta_2l, jend_2u
!       do i = 1, im
!           w ( i, j, l ) = dum3d ( i, j, l )
!       end do
!      end do
!     end do
!  fill up WH which is W at P-points including 2 row halo
      DO L=1,LM
        DO I=1,IM
         DO J=JSTA_2L,JEND_2U 
          WH(I,J,L) = (DUM3D(I,J,L)+DUM3D(I,J,L+1))*0.5
         ENDDO
        ENDDO
      ENDDO

! Read qv

      VarName='QVAPOR'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,         &
        IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            q ( i, j, l ) = dum3d ( i, j, l )/(1.0+dum3d ( i, j, l ))
        end do
       end do
      end do

! 3d pressure 

      VarName='P_HYD'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,         &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            PMID(I,J,L)=DUM3D(I,J,L)
            thv ( i, j, l ) = T(I,J,L)*(Q(I,J,L)*0.608+1.)
            T ( i, j, l ) = T(I,J,L)*(PMID(I,J,L)*1.E-5)**CAPA
             if(abs(t( i, j, l )).gt.1.0e-3)                              &
              omga(I,J,L) = -WH(I,J,L)*pmid(i,j,l)*G/                     &
                              (RD*t(i,j,l)*(1.+D608*q(i,j,l)))
        end do
       end do
      end do
!tgs - 6 June 2012 - added check on monotonic PMID
      do l = 2, lm-1
        ll=lm-l+1
        do j = jsta_2l, jend_2u
          do i = 1, im
            if((PMID(I,J,ll-1) - PMID(I,J,ll)).ge.0.) then

              PMID(I,J,LL)=0.5*(PMID(I,J,LL+1)+PMID(I,J,LL-1))

            endif
          end do
        end do
      end do
!
! check the lowest level too: set l=1, then extroplate the lowest level

      fact=6.0/11.0
      ll=lm
      do j = jsta_2l, jend_2u
        do i = 1, im
          if((PMID(I,J,ll-1) - PMID(I,J,ll)).ge.0.) then

            PMID(I,J,LL)=PMID(I,J,ll-1) +                   &
                         fact*(PMID(I,J,LL-1)-PMID(I,J,LL-2))

          endif
        end do
      end do

! end check

      DO L=2,LM
         DO I=1,IM
            DO J=JSTA_2L,JEND_2U
              PINT(I,J,L)=(PMID(I,J,L-1)+PMID(I,J,L))*0.5 
              ALPINT(I,J,L)=ALOG(PINT(I,J,L))
            ENDDO
         ENDDO
      ENDDO

!---  Compute max temperature in the column up to level 20 
!---    to be used later in precip type computation

       do j = jsta_2l, jend_2u
        do i = 1, im
           tmax(i,j)=0.
        end do
       end do

      do l = 2,20 
       lflip = lm - l + 1
       do j = jsta_2l, jend_2u
        do i = 1, im
           tmax(i,j)=max(tmax(i,j),t(i,j,lflip))
        end do
       end do
      end do


! Brad comment out the output of individual species for Ferriers scheme within
! ARW in Registry file

      qqw=0.
      qqr=0.
      qqs=0.
      qqi=0.
      qqg=0. 
      qqni=0. 
      qqnr=0. 
      cwm=0.

! extinction coef for aerosol
      extcof55=0.
      
      if(imp_physics.ne.5 .and. imp_physics.ne.0)then

! Read qc

      VarName='QCLOUD'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,         &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM)

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
! partition cloud water and ice for WSM3 
	    if(imp_physics.eq.3)then 
             if(t(i,j,l) .ge. TFRZ)then  
              qqw ( i, j, l ) = dum3d ( i, j, l )
	     else
	      qqi  ( i, j, l ) = dum3d ( i, j, l )
	     end if
            else ! bug fix provided by J CASE
             qqw ( i, j, l ) = dum3d ( i, j, l )
	    end if  	     
        end do
       end do
      end do
      end if 
      


      if(imp_physics.ne.1 .and. imp_physics.ne.3                          &
        .and. imp_physics.ne.5 .and. imp_physics.ne.0)then

! Read qi

      VarName='QICE'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,         &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            qqi ( i, j, l ) = dum3d ( i, j, l )
        end do
       end do
      end do
      end if
      

      if(imp_physics.ne.5 .and. imp_physics.ne.0)then

! Read qr

      VarName='QRAIN'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,         &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
! partition rain and snow for WSM3 	
          if(imp_physics .eq. 3)then
	    if(t(i,j,l) .ge. TFRZ)then  
             qqr ( i, j, l ) = dum3d ( i, j, l )
	    else
	     qqs ( i, j, l ) = dum3d ( i, j, l )
	    end if
           else ! bug fix provided by J CASE
            qqr ( i, j, l ) = dum3d ( i, j, l )  
	   end if 
            dummy(i,j)=dum3d(i,j,l)
        end do
       end do
      end do
       do j = jsta_2l, jend_2u
        do i = 1, im
           qrmax(i,j)=0.
        end do
       end do

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
           qrmax(i,j)=max(qrmax(i,j),qqr(i,j,l))
        end do
       end do
      end do

      end if 
     

      if(imp_physics.ne.1 .and. imp_physics.ne.3 .and.                    &
         imp_physics.ne.5 .and. imp_physics.ne.0)then

! Read qs

      VarName='QSNOW'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,         &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            qqs ( i, j, l ) = dum3d ( i, j, l )
            dummy(i,j)=dum3d(i,j,l)
        end do
       end do
      end do
      end if
      

      if(imp_physics.eq.2 .or. imp_physics.eq.6 .or.                      &
         imp_physics.eq.8 .or. imp_physics.eq.9 .or. imp_physics.eq.28)then

! Read qg

      VarName='QGRAUP'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,         &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            qqg ( i, j, l ) = dum3d ( i, j, l )
        end do
       end do
      end do
      end if  

      if(imp_physics.eq.8 .or. imp_physics.eq.9 .or.imp_physics.eq.28)then

! read qni

      VarName='QNICE'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,        &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            qqni ( i, j, l ) = dum3d ( i, j, l )
        end do
       end do
      end do

! Read qnr

      VarName='QNRAIN'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,        &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            qqnr ( i, j, l ) = dum3d ( i, j, l )
        end do
       end do
      end do
      end if

! For aerosol aware microphyscis
      if(imp_physics.eq.28) then

! Read qnc

      VarName='QNCLOUD'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D, &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            qqnw ( i, j, l ) = dum3d ( i, j, l )
        end do
       end do
      end do

! read qnw

      VarName='QNWFA'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D, &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            qqnwfa ( i, j, l ) = dum3d ( i, j, l )
        end do
       end do
      end do

! Read qnifa

      VarName='QNIFA'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D, &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            qqnifa ( i, j, l ) = dum3d ( i, j, l )
        end do
       end do
      end do
      end if

      if(imp_physics.ne.5)then     
!HC SUM UP ALL CONDENSATE FOR CWM       
       do l = 1, lm
        do j = jsta_2l, jend_2u
         do i = 1, im
          IF(QQR(I,J,L).LT.SPVAL)THEN
           CWM(I,J,L)=QQR(I,J,L)
          END IF
          IF(QQI(I,J,L).LT.SPVAL)THEN
           CWM(I,J,L)=CWM(I,J,L)+QQI(I,J,L)
          END IF
          IF(QQW(I,J,L).LT.SPVAL)THEN
           CWM(I,J,L)=CWM(I,J,L)+QQW(I,J,L)
          END IF
          IF(QQS(I,J,L).LT.SPVAL)THEN
           CWM(I,J,L)=CWM(I,J,L)+QQS(I,J,L)
          END IF
          IF(QQG(I,J,L).LT.SPVAL)THEN
           CWM(I,J,L)=CWM(I,J,L)+QQG(I,J,L)
          END IF 
         end do
        end do
       end do
      else
      
! Read cwm

       VarName='CWM'
       call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,        &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM)
       do l = 1, lm
        do j = jsta_2l, jend_2u
         do i = 1, im
            CWM ( i, j, l ) = dum3d ( i, j, l )
         end do
        end do
       end do 

       VarName='F_ICE_PHY'
       call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,        &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM)
       do l = 1, lm
        do j = jsta_2l, jend_2u
         do i = 1, im
            F_ICE( i, j, l ) = dum3d ( i, j, l )
         end do
        end do
       end do 
       
       VarName='F_RAIN_PHY'
       call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,        &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM)
       do l = 1, lm
        do j = jsta_2l, jend_2u
         do i = 1, im
            F_RAIN( i, j, l ) = dum3d ( i, j, l )
         end do
        end do
       end do 

       VarName='F_RIMEF_PHY'
       call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,       &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM)
       do l = 1, lm
        do j = jsta_2l, jend_2u
         do i = 1, im
            F_RIMEF( i, j, l ) = dum3d ( i, j, l )
         end do
        end do
       end do 

      end if
      
       VarName='HTOP'
      IF(ICU_PHYSICS .EQ. 3 .or. ICU_PHYSICS .EQ. 5) VarName='CUTOP'
       call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,       &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            HTOP ( i, j ) = float(LM)-dummy(i,j)+1.0
        end do
       end do
       VarName='HBOT'
      IF(ICU_PHYSICS .EQ. 3 .or. ICU_PHYSICS .EQ. 5) VarName='CUBOT'
       call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,       &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            HBOT ( i, j ) = float(LM)-dummy(i,j)+1.0
        end do
       end do 
       
       VarName='CUPPT'
       call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,       &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            CUPPT ( i, j ) = dummy ( i, j )
        end do
       end do
       
      
      IF(MODELNAME == 'RAPR')THEN

      call getVariable(fileName,DateStr,DataHandle,'TKE_PBL',DUM3D,      &
        IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      ELSE

      call getVariable(fileName,DateStr,DataHandle,'TKE',DUM3D,          &
        IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      ENDIF

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            q2 ( i, j, l ) = dum3d ( i, j, l )
        end do
       end do
      end do

! Height

      VarName='SOILHGT'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            FIS ( i, j ) = dummy ( i, j ) * G
        end do
       end do

      VarName='GHT'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,         &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            ZINT(I,J,L)=DUM3D(I,J,L)
        enddo
       enddo
      enddo

!  get sfc pressure

      VarName='PSFC'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)

      DO J=Jsta,jend
      DO I=1,IM
        if((PINT(I,J,lm) - DUMMY(I,J)).ge.0.) then 
           DUMMY(I,J)=PMID(I,J,LM)*1.001
        endif
        PINT(I,J,LM+1)=DUMMY(I,J)
        ALPINT(I,J,LM+1)=ALOG(PINT(I,J,LM+1))
        PINT(I,J,1)=PMID(I,J,1)
        ALPINT(I,J,1)=ALOG(PINT(I,J,1))
      ENDDO
      ENDDO

      do j = js, je
       do i = 1, im
           ZINT(I,J,LM+1)=FIS(I,J)/G
           DUMMY(I,J)=FIS(I,J)
       end do
      end do
      DO L=LM,1,-1
       do j = js, je
        do i = 1, im
         DUMMY2(I,J)=HTM(I,J,L)*T(I,J,L)*(Q(I,J,L)*D608+1.0)*RD*          &
                   (ALPINT(I,J,L+1)-ALPINT(I,J,L))+DUMMY(I,J)
         DUM3D(I,J,L)=ZINT(I,J,L)-DUMMY2(I,J)/g
         ZINT(I,J,L)=DUMMY2(I,J)/G
         DUMMY(I,J)=DUMMY2(I,J)
        ENDDO
       ENDDO
      END DO
      DO L=LM,1,-1
       do j = js, je
        do i = 1, im
         if(i.eq.im/2.and.j.eq.(jsta+jend)/2) then
                      print*,'DIFF heights model-unipost= ',         &
         i,j,l,DUM3D(I,J,L)
         endif
        ENDDO
       ENDDO
      END DO


      IF(MODELNAME == 'RAPR')THEN

       DO L=1,LM-1
        DO I=1,IM
         DO J=JS,JE
          FACT=(ALOG(PMID(I,J,L))-ALPINT(I,J,L))/                   & 
               max(1.e-6,(ALPINT(I,J,L+1)-ALPINT(I,J,L)))  
          ZMID(I,J,L)=ZINT(I,J,L)+(ZINT(I,J,L+1)-ZINT(I,J,L))*FACT
          dummy(i,j)=ZMID(I,J,L)
         if((ALPINT(I,J,L+1)-ALPINT(I,J,L)) .lt. 1.e-6) print*,       &
                 'P(K+1) and P(K) are too close, i,j,L,',             &
                 'ALPINT(I,J,L+1),ALPINT(I,J,L),ZMID = ',             &
                  i,j,l,ALPINT(I,J,L+1),ALPINT(I,J,L),ZMID(I,J,L)
         ENDDO
        ENDDO
       print*,'max/min ZMID= ',l,maxval(dummy),minval(dummy)
       ENDDO

        DO I=1,IM
         DO J=JS,JE
          DO L=1,LM
           ZINT(I,J,L+1) =AMIN1(ZINT(I,J,L)-2.,ZINT(I,J,L+1))
           ZMID(I,J,L)=(ZINT(I,J,L+1)+ZINT(I,J,L))*0.5  ! ave of z
          ENDDO
          dummy(i,j)=ZMID(I,J,LM)
         ENDDO
        ENDDO
       print*,'max/min ZMID= ',lm,maxval(ZMID(1:im,js:je,lm)),       &
                                  minval(ZMID(1:im,js:je,lm))
     ELSE

       DO L=1,LM
        DO I=1,IM
         DO J=JS,JE
          ZMID(I,J,L)=(ZINT(I,J,L+1)+ZINT(I,J,L))*0.5  ! ave of z
         ENDDO
        ENDDO
       ENDDO

     ENDIF  ! IF(MODELNAME == 'RAPR')THEN

!
! E. James - 8 Dec 2017: this is for HRRR-smoke; it needs to be after ZINT
! is defined.
!
      if(imp_physics.eq.28) then
      VarName='AOD3D_SMOKE'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D, &
        IM+1,1,JM+1,LM+1,IM, JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            taod5503d ( i, j, l ) = dum3d ( i, j, l )
            dz = ZINT( i, j, l ) - ZINT( i, j, l+1 )
            aextc55 ( i, j, l ) = taod5503d ( i, j, l ) / dz
        end do
       end do
      end do
      end if

! get 3-d soil variables
      VarName='SMOIS'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,        &   
        IM+1,1,JM+1,LM+1,IM,JS,JE,NSOIL)
      do l = 1, nsoil
       do j = jsta_2l, jend_2u
        do i = 1, im
            smc ( i, j, l ) = dum3d ( i, j, nsoil-l+1)
            sh2o ( i, j, l ) = dum3d ( i, j, nsoil-l+1)
        end do
       end do
      end do
      
      
      VarName='SEAICE'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &  
        IM,1,JM,1,IM,JS,JE,1)
     
      do j = jsta_2l, jend_2u
        do i = 1, im
            SICE( i, j ) = dummy ( i, j )
        end do
       end do


      VarName='TSLB'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,        & 
        IM+1,1,JM+1,LM+1,IM,JS,JE,NSOIL)
      do l = 1, nsoil
       do j = jsta_2l, jend_2u
        do i = 1, im
            stc ( i, j, l ) = dum3d ( i, j, nsoil-l+1)
        end do
       end do
      end do

      VarName='EROD'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,        & 
        IM+1,1,JM+1,LM+1,IM,JS,JE,3)
      do l = 1, 3
       do j = jsta_2l, jend_2u
        do i = 1, im
            erod ( i, j, l ) = dum3d ( i, j, 4-l)
        end do
       end do
      end do

! bitmask out high, middle, and low cloud cover

       do j = jsta_2l, jend_2u
        do i = 1, im
            CFRACH ( i, j ) = SPVAL/100.
	    CFRACL ( i, j ) = SPVAL/100.
	    CFRACM ( i, j ) = SPVAL/100.
        end do
       end do

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            CFR( i, j, l ) = SPVAL
        end do
       end do
      end do
      
      VarName='SR'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &  
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SR( i, j ) = dummy ( i, j )
        end do
       end do
       
! WRF EM outputs 3D cloud cover now

      VarName='CLDFRA'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,        & 
         IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l=1,lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            CFR ( i, j, l ) = dum3d ( i, j, l )
            CFR_RAW ( i, j, l ) = dum3d ( i, j, l )
        end do
       end do
      end do 

      call ext_ncd_get_dom_ti_real(DataHandle,'DX',tmp,                &      
          1,ioutcount,istatus)
      dxval=nint(tmp)
      write(6,*) 'dxval= ', dxval
       IF(MODELNAME .EQ. 'NCAR' .OR. MODELNAME == 'RAPR')THEN
        if(imp_physics.ne.5 .and. imp_physics.ne.0)then
! Compute 3-D cloud fraction not provided from ARW
        Cloud_def_p = 0.0000001
        do k = 1,lm
          do j = jsta_2l, jend_2u
             do i = 1, im
               dummy(i,j)=QQW(i,j,k)
               dummy2(i,j)=QQI(i,j,k)
             enddo
          enddo

          call AllGETHERV(dummy)
          call AllGETHERV(dummy2)

          do j = jsta_2l, jend_2u
          do i = 1, im
            radius = 30000.
            numr = nint(radius/dxval)
              totcount = 0.
              cloudcount=0.
              cfr(i,j,k) = 0.
              do ic = max(1,I-numr),min(im,I+numr)
              do jc = max(1,J-numr),min(jm,J+numr)
                totcount = totcount+1.
                watericetotal = dummy(ic,jc) + dummy2(ic,jc)
                if ( watericetotal .gt. cloud_def_p) &
                     cloudcount=cloudcount+1.
              enddo
              enddo
              cfr(i,j,k) = min(1.,cloudcount/totcount)
          enddo
          enddo

        enddo 
        do k=1,lm
          print *,'min/max CFR, k= ',minval(CFR(:,:,k)),maxval(CFR(:,:,k)),k
        enddo
!LOW, MID and HIGH cloud fractions
        PTOP_LOW  = 64200.
        PTOP_MID  = 35000.
        PTOP_HIGH = 15000.
!LOW
        do j = jsta_2l, jend_2u
          do i = 1, im
             CFRACL(I,J)=0.
             CFRACM(I,J)=0.
             CFRACH(I,J)=0.

           do k = 1,lm
              if (PMID(I,J,K) .ge. PTOP_LOW) then
!LOW
                CFRACL(I,J)=max(CFRACL(I,J),cfr(i,j,k))
              elseif (PMID(I,J,K) .lt. PTOP_LOW .and. PMID(I,J,K) .ge. PTOP_MID) then
!MID
                CFRACM(I,J)=max(CFRACM(I,J),cfr(i,j,k))
              elseif (PMID(I,J,K) .lt. PTOP_MID .and. PMID(I,J,K) .ge. PTOP_HIGH) then
!HIGH
                CFRACH(I,J)=max(CFRACH(I,J),cfr(i,j,k))
              endif
           enddo

          enddo
        enddo

        print *,' MIN/MAX CFRACL ',minval(CFRACL),maxval(CFRACL)
        print *,' MIN/MAX CFRACM ',minval(CFRACM),maxval(CFRACM)
        print *,' MIN/MAX CFRACH ',minval(CFRACH),maxval(CFRACH)
        ENDIF  ! Not Ferrier or null mp_physics
      ENDIF   ! NCAR or RAPR

! GRAUP
! E. James - 8 Dec 2017: SMOKE from WRF-CHEM
!
     plume = 0.
     VarName='NCLPLM'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,        & 
         IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l = 1,lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            plume ( i, j, l) = dum3d ( i, j, l )
        end do
       end do
      end do

!-------------------------------------------------------------
! DUST
!-------------------------------------------------------------

      DUST     = 0.0d0
      DUSMASS  = 0.0d0
      VIS_DUST = 0.0d0

!-------------------------------------------------------------
! Read MPAS dust bins
!-------------------------------------------------------------
      VarName='DUST_1'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,        &
     &     IM,1,JM,LM,IM,JS,JE,LM)

      do l = 1, lm
        do j = jsta_2l, jend_2u
          do i = 1, im
            DUST(i,j,l,1) = DUM3D(i,j,l)
          end do
        end do
      end do

      VarName='DUST_2'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,        &
     &     IM+1,1,JM+1,LM+1,IM,JS,JE,LM)

      do l = 1, lm
        do j = jsta_2l, jend_2u
          do i = 1, im
            DUST(i,j,l,2) = DUM3D(i,j,l)
          end do
        end do
      end do

      VarName='DUST_3'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,        &
     &     IM+1,1,JM+1,LM+1,IM,JS,JE,LM)

      do l = 1, lm
        do j = jsta_2l, jend_2u
          do i = 1, im
            DUST(i,j,l,3) = DUM3D(i,j,l)
          end do
        end do
      end do

      VarName='DUST_4'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,        &
     &     IM+1,1,JM+1,LM+1,IM,JS,JE,LM)

      do l = 1, lm
        do j = jsta_2l, jend_2u
          do i = 1, im
            DUST(i,j,l,4) = DUM3D(i,j,l)
          end do
        end do
      end do

      print*,'Read DUST_1..DUST_4 at domain center :',                  &
     &        DUST(im/2,jm/2,lm,1), DUST(im/2,jm/2,lm,2),              &
     &        DUST(im/2,jm/2,lm,3), DUST(im/2,jm/2,lm,4)

!-------------------------------------------------------------
! Use the 10 lowest levels near the surface
!-------------------------------------------------------------
      lbot = max(1, lm-9)
      
      DUSMASS  = 0.0d0
      VIS_DUST = 0.0d0
      
      do j = jsta_2l, jend_2u
        do i = 1, im
      
          dust10max = 0.0d0
      
          do l = lbot, lm
      
            dustlev = max(dble(DUST(i,j,l,1)), 0.0d0) +                 &
           &          max(dble(DUST(i,j,l,2)), 0.0d0) +                 &
           &          max(dble(DUST(i,j,l,3)), 0.0d0) +                 &
           &          max(dble(DUST(i,j,l,4)), 0.0d0)
      
            dust10max = max(dust10max, dustlev)
      
          end do
      
          DUSMASS(i,j) = dust10max / 3.2d0
      
          VIS_DUST(i,j) = 10000.0d0 /                                   &
               &              (1.0d0 + (DUSMASS(i,j)/1310.0d0)**0.95d0)
          
          VIS_DUST(i,j) = min(10000.0d0, max(VIS_DUST(i,j), 200.0d0))
      
        end do
      end do

      print*,'MPAS DUSMASS/VIS_DUST at domain center :',                &
     &        DUSMASS(im/2,jm/2), VIS_DUST(im/2,jm/2)

      VarName='TOT_EDUST'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            TOT_EDUST ( i, j ) = dummy ( i, j )
        end do
       end do
      DUDP = 0.
      VarName='DRYDEP_1'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            DUDP ( i, j, 1 ) = dummy ( i, j )
        end do
       end do
      VarName='DRYDEP_2'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            DUDP ( i, j, 2 ) = dummy ( i, j )
        end do
       end do
      VarName='DRYDEP_3'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            DUDP ( i, j, 3 ) = dummy ( i, j )
        end do
       end do
      VarName='DRYDEP_4'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            DUDP ( i, j, 4 ) = dummy ( i, j )
        end do
       end do
      VarName='DRYDEP_5'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            DUDP ( i, j, 5 ) = dummy ( i, j )
        end do
       end do
      DUSD=0.
      VarName='GRASET_1'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            DUSD ( i, j, 1 ) = dummy ( i, j )
        end do
       end do
      VarName='GRASET_2'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            DUSD ( i, j, 2 ) = dummy ( i, j )
        end do
       end do
 
      VarName='GRASET_3'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            DUSD ( i, j, 3 ) = dummy ( i, j )
        end do
       end do
      VarName='GRASET_4'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            DUSD ( i, j, 4 ) = dummy ( i, j )
        end do
       end do
 
      VarName='GRASET_5'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            DUSD ( i, j, 5 ) = dummy ( i, j )
        end do
       end do
     VarName='PM2_5_DRY'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,&
         IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l=1,lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            DUST ( i, j, l, 6) = dum3d ( i, j, l )
        end do
       end do
      end do
     VarName='PM10'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,&
         IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l=1,lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            DUST ( i, j, l, 7) = dum3d ( i, j, l )
        end do
       end do
      end do

     VarName='so2'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,&
         IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l=1,lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            DUST ( i, j, l, 8) = dum3d ( i, j, l )
        end do
       end do
      end do


! Soil layer/depth - extract thickness of soil layers from wrf output

      SLDPTH=0.0 !Assign bogus value

      ! RUC LSM - use depths of center of soil layer
      IF(iSF_SURFACE_PHYSICS==3)then ! RUC LSM
       call getVariable(fileName,DateStr,DataHandle,'ZS',SLLEVEL,  &
        NSOIL,1,1,1,NSOIL,1,1,1)
       print*,'SLLEVEL= ',(SLLEVEL(N),N=1,NSOIL)
      ELSE
       call getVariable(fileName,DateStr,DataHandle,'DZS',SLDPTH,  &
        NSOIL,1,1,1,NSOIL,1,1,1)
! Radi ... I will force the layers to be conform with AFAD operational MPAS
       SLDPTH(:) = 0.0
       SLDPTH(1) = 0.1
       SLDPTH(2) = 0.3
       SLDPTH(3) = 0.6
       SLDPTH(4) = 1.0
       print*,'SLDPTH= ',(SLDPTH(N),N=1,NSOIL)
      END IF

! SRD
! get 2-d variables

      VarName='WSPD10MAX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            WSPD10MAX ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='WSPD10UMAX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            WSPD10UMAX ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='WSPD10VMAX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            WSPD10VMAX ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='W_UP_MAX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            W_UP_MAX ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='W_DN_MAX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            W_DN_MAX ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='W_MEAN'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            W_MEAN ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='REFD_MAX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            REFD_MAX ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='REFDM10C_MAX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            REFDM10C_MAX ( i, j ) = dummy ( i, j )
        end do
       end do


      VarName='UP_HELI_MAX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            UP_HELI_MAX ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='UP_HELI_MAX16'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            UP_HELI_MAX16 ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='UP_HELI_MIN'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            UP_HELI_MIN ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='UP_HELI_MIN16'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            UP_HELI_MIN16 ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='UP_HELI_MAX02'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            UP_HELI_MAX02 ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='UP_HELI_MIN02'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            UP_HELI_MIN02 ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='UP_HELI_MAX03'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            UP_HELI_MAX03 ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='UP_HELI_MIN03'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            UP_HELI_MIN03 ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='REL_VORT_MAX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            REL_VORT_MAX ( i, j ) = dummy ( i, j )
        end do
       end do
        
      VarName='REL_VORT_MAX01'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            REL_VORT_MAX01 ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='GRPL_MAX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            GRPL_MAX ( i, j ) = dummy ( i, j )
        end do
       end do


      VarName='HAIL_MAXK1'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            HAIL_MAXK1 ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='HAIL_MAX2D'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            HAIL_MAX2D ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='UH'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            UP_HELI ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='UH16'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            UP_HELI16 ( i, j ) = dummy ( i, j )
        end do
       end do 

      VarName='LTG1_MAX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            LTG1_MAX ( i, j ) = dummy ( i, j )
        end do
       end do 

      VarName='LTG2_MAX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            LTG2_MAX ( i, j ) = dummy ( i, j )
        end do
       end do
    
      VarName='LTG3_MAX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            LTG3_MAX ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='NCI_LTG'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            NCI_LTG ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='NCA_LTG'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            NCA_LTG ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='NCI_WQ'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            NCI_WQ ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='NCA_WQ'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            NCA_WQ ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='NCI_REFD'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            NCI_REFD ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='NCA_REFD'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            NCA_REFD ( i, j ) = dummy ( i, j )
        end do
       end do
! 
! SRD
!
! CRA REFLECTIVITY VARIABLES
      VarName='REFL_10CM'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,&
         IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
       do l=1,lm
        do j = jsta_2l, jend_2u
         do i = 1, im
            REF_10CM ( i, j, l) = dum3d ( i, j, l )
         end do
        end do
       end do
      deallocate(DUM3D)

      VarName='COMPOSITE_REFL_10CM'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            REFC_10CM ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='REFL_10CM_1KM'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            REF1KM_10CM ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='REFL_10CM_4KM'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            REF4KM_10CM ( i, j ) = dummy ( i, j )
        end do
       end do
! CRA
! get 2-d variables

      VarName='U10'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &  
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            U10 ( i, j ) = dummy( i, j )
        end do
       end do
      VarName='V10'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            V10 ( i, j ) = dummy( i, j )
        end do
       end do

! RAP/HRRR time-averaged wind
      VarName='U10MEAN'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY, &       
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            U10mean ( i, j ) = dummy( i, j )
        end do
       end do
!
      VarName='V10MEAN'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY, &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            V10mean ( i, j ) = dummy( i, j )
        end do
       end do
!
      VarName='SPDUV10MEAN'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY, &       
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SPDUV10mean ( i, j ) = dummy( i, j )
        end do
       end do
!

       do j = jsta_2l, jend_2u
        do i = 1, im
            TH10 ( i, j ) = SPVAL
	    Q10 ( i, j ) = SPVAL
        end do
       end do

! get 2-m theta 
      VarName='TH2'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &  
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            TSHLTR ( i, j ) = dummy ( i, j )
        end do
       end do
      VarName='Q2'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &  
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
          MRSHLTR ( i, j ) = dummy (i, j )  ! Shelter Mixing ratio
          IF(MODELNAME == 'RAPR')THEN
            QV2M ( i, j )   =  dummy ( i, j ) ! 2-m mix. ratio
            QSHLTR ( i, j ) = dummy ( i, j )/(1.0+dummy ( i, j )) ! 2-m spec. hum.
            QVl1 ( i, j )   =  q ( i, j, lm ) ! spec. humidity at lev. 1
          ELSE
            QV2M ( i, j ) =  dummy ( i, j )
            QSHLTR ( i, j ) = dummy ( i, j )/(1.0+dummy ( i, j ))
          ENDIF
        end do
       end do

      IF(MODELNAME == 'RAPR')THEN
        VarName='MAVAIL'
      ELSE
        VarName='SMSTAV'
      END IF

      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SMSTAV ( i, j ) = dummy ( i, j )
        end do
       end do
       
      VarName='SFROFF' 
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        & 
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SSROFF ( i, j ) = dummy ( i, j )
        end do
       end do
      VarName='UDROFF'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        & 
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            BGROFF ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='VEGFRA'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            VEGFRC ( i, j ) = dummy ( i, j )/100.
        end do
       end do
      VarName='ACSNOW'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        & 
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            ACSNOW ( i, j ) = dummy ( i, j )
        end do
       end do
      VarName='ACSNOM'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        & 
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            ACSNOM ( i, j ) = dummy ( i, j )
        end do
       end do
      VarName='CANWAT'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        & 
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            CMC ( i, j ) = dummy ( i, j )
        end do
       end do
      VarName='SST'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SST ( i, j ) = dummy ( i, j )
        end do
       end do
      VarName='THZ0'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        & 
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            THZ0 ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='Z0'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            Z0( i, j ) = dummy ( i, j )
        end do
       end do
      VarName='UST'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            USTAR( i, j ) = dummy ( i, j )
        end do
       end do
      VarName='TSK'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &  
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
           RADOT ( i, j ) = DUMMY(i,j)**4.0/STBOL
           THS ( i, j ) = dummy ( i, j )                                 &
                  *(P1000/PINT(I,J,NINT(LMH(I,J))+1))**CAPA
        end do
       end do

      VarName='EMISS'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,     &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            EPSR( i, j ) = dummy ( i, j )
        end do
       end do

      IF(MODELNAME .EQ. 'RAPR')THEN
        call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
          IM,1,JM,1,IM,JS,JE,1)
         do j = jsta_2l, jend_2u
          do i = 1, im
              RADOT ( i, j ) = RADOT (i, j) * dummy ( i, j )
          end do
         end do
      END IF

      VarName='RAINC'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &  
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            CUPREC ( i, j ) = dummy ( i, j ) * 0.001
        end do
       end do
      VarName='RAINNC'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            ANCPRC ( i, j ) = dummy ( i, j )* 0.001
        end do
       end do

       do j = jsta_2l, jend_2u
        do i = 1, im
            ACPREC(I,J)=ANCPRC(I,J)+CUPREC(I,J)
        end do
       end do  

!-- RAINC_bucket is "ACCUMULATED CUMULUS PRECIPITATION OVER BUCKET_DT PERIODS OF TIME"

      VarName='PREC_ACC_C'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            rainc_bucket ( i, j ) = dummy ( i, j )
        end do
       end do

!-- RAINNC_bucket  is "ACCUMULATED GRID SCALE  PRECIPITATION OVER BUCKET_DT PERIODS OF TIME"

      VarName='PREC_ACC_NC'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            rainnc_bucket ( i, j ) = dummy ( i, j )
        end do
       end do

       do j = jsta_2l, jend_2u
        do i = 1, im
            PCP_BUCKET(I,J)=rainc_bucket(I,J)+rainnc_bucket(I,J)
        end do
       end do

      VarName='RAINCV'
      DUMMY=0.0
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
!-- CPRATE is in [m] per time step
            CPRATE ( i, j ) = dummy ( i, j )* 0.001
        end do
       end do
     

      VarName='RAINNCV'
      DUMMY2=0.0
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,       &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
!-- PREC is in [m] per time step
            prec ( i, j ) = (dummy ( i, j )+dummy2(i,j))* 0.001
        end do
       end do

      VarName='SNOWNCV'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
!-- SNOW is in [m] per time sep
            snownc ( i, j ) = dummy ( i, j ) * 0.001
        end do
       end do

!-- SNOW_bucket  is "ACCUMULATED GRID SCALE SNOW OVER BUCKET_DT PERIODS OF TIME"

        write(6,*) 'getting SNOW_ACC_NC, [mm] '
!      VarName='SNOW_BUCKET'
      VarName='SNOW_ACC_NC'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            snow_bucket ( i, j ) = dummy ( i, j )
        end do
       end do

!-- GRAUP_bucket  is "ACCUMULATED GRID SCALE GRAUPEL OVER BUCKET_DT PERIODS OF TIME"

        write(6,*) 'getting GRAUP_ACC_NC, [mm] '
      VarName='GRAUP_ACC_NC'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,       &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            graup_bucket ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='ACGRAUP'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &   
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            acgraup ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='ACFRAIN'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &   
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            acfrain ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='GRAUPELNCV'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,       &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
!-- GRAUPEL in in [m] per time step
            graupelnc ( i, j ) = dummy ( i, j ) * 0.001
        end do
       end do


      VarName='ALBEDO'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            ALBEDO ( i, j ) = dummy ( i, j )
        end do
       end do
!
!  GSD trunk as 2014.4.14
!      VarName='GSW'
!      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,&
!        IM,1,JM,1,IM,JS,JE,1)
!       do j = jsta_2l, jend_2u
!        do i = 1, im
!            RSWIN ( i, j ) = dummy ( i, j )
! HCHUANG: GSW is actually net downward shortwave in ncar wrf
!             RSWIN ( i, j ) = dummy ( i, j )/(1.0-albedo(i,j))
!             RSWOUT ( i, j ) = RSWIN ( i, j ) - dummy ( i, j )
!        end do
!       end do

      VarName='SWDOWN'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
! HCHUANG: SWDOWN is actually net downward shortwave in ncar wrf
             RSWIN ( i, j ) = dummy ( i, j )
             RSWOUT ( i, j ) = RSWIN ( i, j ) * ALBEDO ( i, j )
        end do
       end do

      VarName='SWDDNI'
! Shortwave surface downward direct normal irradiance
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY, &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
             SWDDNI ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='SWDDIF'
! Shortwave surface downward diffuse horizontal irradiance
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY, &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
             SWDDIF ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='SWDNBC'
! Shortwave surface downward clear-sky shortwave irradiance
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY, &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
             SWDNBC ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='SWDDNIC'
! Clear-sky shortwave surface downward direct normal irradiance
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY, &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
             SWDDNIC ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='SWDDIFC'
! Clear-sky shortwave surface downward diffuse horizontal irradiance
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY, &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
             SWDDIFC ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='SWUPBC'
! Clear-sky surface upwelling shortwave flux
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY, &
        IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
       do i = 1, im
            SWUPBC ( i, j ) = dummy ( i, j )
       end do
      end do

      VarName='SWUPT'
! Upward shortwave flux at top of atmosphere
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY, &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
             SWUPT ( i, j ) = dummy ( i, j )
        end do
       end do

!
! E. James - 8 Dec 2017: Fire Radiative Power from HRRR-smoke
!
      VarName='MEAN_FRP'
! Mean fire radiative power
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY, &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
             MEAN_FRP ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='TAOD5502D'
! Total aerosol optical depth
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,     &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            TAOD5502D ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='AERASY2D'
! Aerosol asymmetry parameter
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,     &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            AERASY2D ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='AERSSA2D'
!  Aerosol single-scattering albedo
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,     &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            AERSSA2D ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='LWP'
!  Liquid water path
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,     &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            LWP ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='IWP'
!  Ice water path
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,     &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            IWP ( i, j ) = dummy ( i, j )
        end do
       end do

! time_averaged SWDOWN
      VarName='SWRADMEAN'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY, &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
! averaged incoming solar radiation
             SWRADmean ( i, j ) = dummy ( i, j )
        end do
       end do

! time_averaged SWNORM
      VarName='SWNORMMEAN'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY, &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
! averaged incoming solar radiation
             SWNORMmean ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='GLW'
! Downwelling longwave flux at surface
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &  
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            RLWIN ( i, j ) = dummy ( i, j )
        end do
       end do
! ncar wrf does not output sigt4 so make sig4=sigma*tlmh**4
       do j = jsta_2l, jend_2u
        do i = 1, im
             TLMH=T(I,J,NINT(LMH(I,J)))
             SIGT4 ( i, j ) =  5.67E-8*TLMH*TLMH*TLMH*TLMH
        end do
       end do

      VarName='LWDNBC'
! Clear-sky downwelling longwave flux at surface
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY, &
        IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
       do i = 1, im
            LWDNBC ( i, j ) = dummy ( i, j )
       end do
      end do

      VarName='LWUPBC'
! Clear-sky upwelling longwave flux at surface
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY, &
        IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
       do i = 1, im
            LWUPBC ( i, j ) = dummy ( i, j )
       end do
      end do

! Top of the atmosphere outgoing LW radiation
      VarName='OLR'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            RLWTOA ( i, j ) = dummy ( i, j )
        end do
       end do


! NCAR WRF does not output accumulated fluxes so set the bitmap of these fluxes to 0
      do j = jsta_2l, jend_2u
        do i = 1, im
!	   RLWTOA(I,J)=SPVAL
	   RSWINC(I,J)=SPVAL
           ASWIN(I,J)=SPVAL  
	   ASWOUT(I,J)=SPVAL
	   ALWIN(I,J)=SPVAL
	   ALWOUT(I,J)=SPVAL
	   ALWTOA(I,J)=SPVAL
	   ASWTOA(I,J)=SPVAL
	   ARDLW=1.0
	   ARDSW=1.0
	   NRDLW=1
	   NRDSW=1
        end do
       end do

      VarName='TMN'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        & 
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            TG ( i, j ) = dummy ( i, j )
            SOILTB ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='HFX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        & 
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            TWBS(I,J)= dummy ( i, j )
!            SFCSHX ( i, j ) = dummy ( i, j )
!            ASRFC=1.0
        end do
       end do

! latent heat flux
      IF(iSF_SURFACE_PHYSICS.NE.3) then
      VarName='LH'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            QWBS(I,J) = dummy ( i, j )
!            SFCLHX ( i, j ) = dummy ( i, j )
        end do
       end do
      else
       VarName='QFX'
       call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
        do j = jsta_2l, jend_2u
         do i = 1, im
             QWBS(I,J) = dummy ( i, j ) * LHEAT
         end do
        end do
      ENDIF

! ground heat fluxes       
      VarName='GRDFLX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &  
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            GRNFLX(I,J) = dummy ( i, j )
        end do
       end do       

! NCAR WRF does not output accumulated fluxes so bitmask out these fields
      do j = jsta_2l, jend_2u
        do i = 1, im
           SFCSHX(I,J)=SPVAL  
	   SFCLHX(I,J)=SPVAL
	   SUBSHX(I,J)=SPVAL
	   SNOPCX(I,J)=SPVAL
	   SFCUVX(I,J)=SPVAL
	   POTEVP(I,J)=SPVAL
	   NCFRCV(I,J)=SPVAL
	   NCFRST(I,J)=SPVAL
	   ASRFC=1.0
	   NSRFC=1
        end do
       end do

!      VarName='WEASD'
! Snow water equivalent
      VarName='SNOW'  ! WRF V2 replace WEASD with SNOW
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &  
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
           if( dummy ( i, j ) <= 5000.0 .and. dummy ( i, j ) >=0.0) then
              SNO ( i, j ) = dummy ( i, j )
           elseif( dummy ( i, j ) > 5000.0) then 
              SNO ( i, j ) = 5000.0
           elseif( dummy ( i, j ) < 0.0 ) then 
              SNO ( i, j ) = 0.0
           else 
              SNO ( i, j ) = 0.0
           endif
        end do
       end do
! Snow depth
      VarName='SNOWH'   
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &  
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
           if( dummy ( i, j ) <= 50.0 .and. dummy ( i, j ) >=0.0) then
              SI ( i, j ) = dummy ( i, j ) * 1000.
           elseif( dummy ( i, j ) > 50.0) then 
              SI ( i, j ) = 50.0 * 1000.
           elseif( dummy ( i, j ) < 0.0 ) then 
              SI ( i, j ) = 0.0
           else 
              SI ( i, j ) = 0.0
           endif
        end do
       end do

! snow cover
      VarName='SNOWC'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        & 
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            PCTSNO ( i, j ) = dummy ( i, j )
        end do
       end do 

! Accumulated grid-scale snow and ice precipitation
      VarName='SNOWNC'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SNONC  ( i, j ) = dummy ( i, j )
        end do
       end do

! snowfall density
      VarName='RHOSNF'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,  &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SNFDEN ( i, j ) = max(0.,dummy ( i, j ))
        end do
       end do

! snowfall accumulation
      VarName='SNOWFALLAC'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY, &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SNDEPAC ( i, j ) = dummy ( i, j )
        end do
       end do

! snow temperature at the interface of 2 snow layers
      VarName='SOILT1'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &  
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
           TSNOW ( i, j ) = dummy ( i, j )
        end do
       end do

! GET VEGETATION TYPE

      call getIVariableN(fileName,DateStr,DataHandle,'IVGTYP',IDUMMY,    &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            IVGTYP ( i, j ) = idummy ( i, j ) 
        end do
       end do
       
      VarName='ISLTYP'
      call getIVariableN(fileName,DateStr,DataHandle,VarName,IDUMMY,     & 
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            ISLTYP ( i, j ) = idummy ( i, j ) 
        end do
       end do

!      VarName='ISLOPE'
!      call getIVariableN(fileName,DateStr,DataHandle,VarName,IDUMMY,     &
!        IM,1,JM,1,IM,JS,JE,1)
!       do j = jsta_2l, jend_2u
!        do i = 1, im
!            ISLOPE( i, j ) = idummy ( i, j )
!        end do
!       end do
       

! XLAND 1 land 2 sea
      VarName='XLAND'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SM ( i, j ) = dummy ( i, j ) - 1.0
        end do
       end do

! PBL depth
      VarName='PBLH'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,     &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            PBLH ( i, j ) = dummy ( i, j )
        end do
       end do
      IF(MODELNAME .EQ. 'RAPR')THEN
! PBL depth from GSD
       delta_theta4gust=0.5
       do j = jsta_2l, jend_2u
        do i = 1, im
!!   Is there any mixed layer at all?
          if (thv(i,j,lm-1) .lt. (thv(i,j,lm) + delta_theta4gust)) then
            ZSF=ZINT(I,J,NINT(LMH(I,J))+1)
!!   Calculate k1 level as first above PBL top
            do 34 k=3,LM
              k1 = k
!! - give theta-v at the sfc a 0.5K boost in
!!         the PBL height definition
              if (thv(i,j,lm-k+1).gt.(thv(i,j,lm) + delta_theta4gust))  &
                      go to 341
 34         continue
 341        continue
              zpbltop = zmid(i,j,lm-k1+1) +                             &
                   ((thv(i,j,lm)+delta_theta4gust)-thv(i,j,lm-k1+1))    &
                 * (zmid(i,j,lm-k1+2)-zmid(i,j,lm-k1+1))                &
                 / (thv(i,j,lm-k1+2) - thv(i,j,lm-k1+1))
            PBLHGUST ( i, j ) = zpbltop - zsf
          else
            PBLHGUST ( i, j ) = 0.
          endif
        end do
       end do
      ENDIF

      VarName='XLAT_C'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            GDLAT ( i, j ) = dummy ( i, j )
            f(i,j) = 1.454441e-4*sin(gdlat(i,j)*DTR)
        end do
       end do
! pos north
      VarName='XLONG_C'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,         &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            GDLON ( i, j ) = dummy ( i, j )
        end do
       end do
! pos east
       call collect_loc(gdlat,dummy)
       if(me.eq.0)then
        latstart=nint(dummy(1,1)*gdsdegr)
        latlast=nint(dummy(im,jm)*gdsdegr)
       end if
       call mpi_bcast(latstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(latlast,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       write(6,*) 'laststart,latlast A calling bcast= ',latstart,latlast


       call collect_loc(gdlon,dummy)
       if(me.eq.0)then
        if (grib=="grib2" .and. dummy(1,1) <0.) dummy(1,1)=dummy(1,1)+360.
        if (grib=="grib2" .and. dummy(im,jm) <0.) dummy(im,jm)=dummy(im,jm)+360.
        lonstart=dummy(1,1)*gdsdegr
        lonlast=dummy(im,jm)*gdsdegr
       end if
       call mpi_bcast(lonstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(lonlast,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       write(6,*) 'lonstart,lonlast A calling bcast= ',lonstart,lonlast
!
! obtain map scale factor
!      VarName='msft'
      VarName='MAPFAC_M'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,        &  
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            MSFT ( i, j ) = dummy ( i, j ) 
        end do
       end do

! physics calling frequency
      VarName='STEPBL'
      call getIVariableN(fileName,DateStr,DataHandle,VarName,NPHS,       &
        1,1,1,1,1,1,1,1)
     

!        ncdump -h

! ncar wrf does not output zenith angle so make czen=czmean so that
! RSWIN can be output normally in SURFCE
      IF(MODELNAME .NE. 'RAPR')THEN
       do j = jsta_2l, jend_2u
        do i = 1, im
             CZEN ( i, j ) = 1.0
             CZMEAN ( i, j ) = CZEN ( i, j )
        end do
       end do
      ELSE

        jdn=iw3jdn(idat(3),idat(1),idat(2))
        do j=jsta,jend
         do i=1,im
             call zensun(jdn,float(idat(4)),gdlat(i,j),gdlon(i,j)     &
               ,pi,sun_zenith,sun_azimuth)
             temp=sun_zenith/rtd
             czen(i,j)=cos(temp)
             CZMEAN ( i, j ) = CZEN ( i, j )
         end do
        end do
      ENDIF


!!
!! 
!!
        write(6,*) 'filename in INITPOST=', filename,' is'

!	status=nf_open(filename,NF_NOWRITE,ncid)
!	        write(6,*) 'returned ncid= ', ncid
!        status=nf_get_att_real(ncid,varid,'DX',tmp)
!	dxval=int(tmp)
!        status=nf_get_att_real(ncid,varid,'DY',tmp)
!	dyval=int(tmp)
!        status=nf_get_att_real(ncid,varid,'CEN_LAT',tmp)
!	cenlat=int(1000.*tmp)
!        status=nf_get_att_real(ncid,varid,'CEN_LON',tmp)
!	cenlon=int(1000.*tmp)
!        status=nf_get_att_real(ncid,varid,'TRUELAT1',tmp)
!	truelat1=int(1000.*tmp)
!        status=nf_get_att_real(ncid,varid,'TRUELAT2',tmp)
!	truelat2=int(1000.*tmp)
!        status=nf_get_att_real(ncid,varid,'MAP_PROJ',tmp)
!        maptype=int(tmp)
!	status=nf_close(ncid)

!	dxval=30000.
! 	dyval=30000.
!
!        write(6,*) 'dxval= ', dxval
!        write(6,*) 'dyval= ', dyval
!        write(6,*) 'cenlat= ', cenlat
!        write(6,*) 'cenlon= ', cenlon
!        write(6,*) 'truelat1= ', truelat1
!        write(6,*) 'truelat2= ', truelat2
!        write(6,*) 'maptype is ', maptype
!
!tgs        call ext_ncd_get_dom_ti_real(DataHandle,'DX',tmp,                &      
!          1,ioutcount,istatus)
!        dxval=nint(tmp)
!        write(6,*) 'dxval= ', dxval
        call ext_ncd_get_dom_ti_real(DataHandle,'DY',tmp,                &
          1,ioutcount,istatus)
        dyval=nint(tmp)
        write(6,*) 'dyval= ', dyval
        call ext_ncd_get_dom_ti_real(DataHandle,'CEN_LAT',tmp,           &
          1,ioutcount,istatus)
        cenlat=nint(gdsdegr*tmp)
        write(6,*) 'cenlat= ', cenlat
        call ext_ncd_get_dom_ti_real(DataHandle,'CEN_LON',tmp,           &
          1,ioutcount,istatus)
        if(tmp < 0) tmp=360.0 + tmp
        cenlon=nint(gdsdegr*tmp)
        write(6,*) 'cenlon= ', cenlon
        call ext_ncd_get_dom_ti_integer(DataHandle,'MAP_PROJ',itmp,      &
          1,ioutcount,istatus)
        maptype=itmp
        write(6,*) 'maptype is ', maptype
        if(maptype.ne.6)then
        call ext_ncd_get_dom_ti_real(DataHandle,'TRUELAT1',tmp,          &
          1,ioutcount,istatus)
        truelat1=nint(gdsdegr*tmp)
        write(6,*) 'truelat1= ', truelat1
        if(maptype.ne.2)then  !PS projection excluded
          call ext_ncd_get_dom_ti_real(DataHandle,'TRUELAT2',tmp,          &
            1,ioutcount,istatus)
          truelat2=nint(gdsdegr*tmp)
          write(6,*) 'truelat2= ', truelat2
        endif
        endif
	call ext_ncd_get_dom_ti_real(DataHandle,'STAND_LON',tmp,         &
          1,ioutcount,istatus)
        if(tmp < 0) tmp=360.0 + tmp
        STANDLON=nint(gdsdegr*tmp)
        write(6,*) 'STANDLON= ', STANDLON

!MEB not sure how to get these 
       do j = jsta_2l, jend_2u
        do i = 1, im
            DX ( i, j ) = dxval/MSFT(I,J)  
            DY ( i, j ) = dyval/MSFT(I,J)  
        end do
       end do

! Convert DXVAL and DYVAL for ARW rotated latlon from meters to radian
        if(maptype==6)then
         dxval=(DXVAL * 360.)/(ERAD*2.*pi)*gdsdegr
         dyval=(DYVAL * 360.)/(ERAD*2.*pi)*gdsdegr
         print*,'dx and dy for arw rotated latlon= ', &
         dxval,dyval
        end if

!tgs Define smoothing flag for isobaric output 
              IF(MODELNAME == 'RAPR')THEN
                SMFLAG=.TRUE.
              ELSE
                SMFLAG=.TRUE.
              ENDIF

! generate look up table for lifted parcel calculations

      THL=210.
      PLQ=70000.

      CALL TABLE(PTBL,TTBL,PT,                                           &  
                RDQ,RDTH,RDP,RDTHE,PL,THL,QS0,SQS,STHE,THE0)

      CALL TABLEQ(TTBLQ,RDPQ,RDTHEQ,PLQ,THL,STHEQ,THE0Q)

!     
!     
      IF(ME.EQ.0)THEN
        WRITE(6,*)'  SPL (POSTED PRESSURE LEVELS) BELOW: '
        WRITE(6,51) (SPL(L),L=1,LSM)
   50   FORMAT(14(F4.1,1X))
   51   FORMAT(8(F8.1,1X))
      ENDIF
!     
!     COMPUTE DERIVED TIME STEPPING CONSTANTS.
!
!need to get DT
      call ext_ncd_get_dom_ti_real(DataHandle,'DT',tmp,1,ioutcount,istatus)
      DT=abs(tmp)

!need to get period of time for precipitation buckets
      call ext_ncd_get_dom_ti_real(DataHandle,'PREC_ACC_DT',tmp,1,ioutcount,istatus)
      prec_acc_dt=abs(tmp)
       
!      DT = 120. !MEB need to get DT
      NPHS = 1  !CHUANG SET IT TO 1 BECAUSE ALL THE INST PRECIP ARE ACCUMULATED 1 TIME STEP
      DTQ2 = DT * NPHS  !MEB need to get physics DT
      TSPH = 3600./DT   !MEB need to get DT
! Randomly specify accumulation period because WRF EM does not
! output accumulation fluxes yet and accumulated fluxes are bit
! masked out

      TSRFC=1.0
      TRDLW=1.0
      TRDSW=1.0
      THEAT=1.0
      TCLOD=1.0

      TPREC=float(NPREC)/TSPH
      IF(NPREC.EQ.0)TPREC=float(ifhr)  !in case buket does not get emptied

!tgs      TPREC=float(ifhr)  ! WRF EM does not empty precip buket at all

!      TSRFC=float(NSRFC)/TSPH
!      TRDLW=float(NRDLW)/TSPH
!      TRDSW=float(NRDSW)/TSPH
!      THEAT=float(NHEAT)/TSPH
!      TCLOD=float(NCLOD)/TSPH
!      TPREC=float(NPREC)/TSPH
!MEB need to get DT

!how am i going to get this information?
!      NPREC  = INT(TPREC *TSPH+D50)
!      NHEAT  = INT(THEAT *TSPH+D50)
!      NCLOD  = INT(TCLOD *TSPH+D50)
!      NRDSW  = INT(TRDSW *TSPH+D50)
!      NRDLW  = INT(TRDLW *TSPH+D50)
!      NSRFC  = INT(TSRFC *TSPH+D50)
!how am i going to get this information?
!     
!     IF(ME.EQ.0)THEN
!       WRITE(6,*)' '
!       WRITE(6,*)'DERIVED TIME STEPPING CONSTANTS'
!       WRITE(6,*)' NPREC,NHEAT,NSRFC :  ',NPREC,NHEAT,NSRFC
!       WRITE(6,*)' NCLOD,NRDSW,NRDLW :  ',NCLOD,NRDSW,NRDLW
!     ENDIF
!

!      VarName='RAINCV'
!      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
!     &  IM,1,JM,1,IM,JS,JE,1)
!       do j = jsta_2l, jend_2u
!        do i = 1, im
!            CUPPT ( i, j ) = dummy ( i, j )* 0.001*(TRDLW*3600.)	    
!        end do
!       end do


      
!     COMPUTE DERIVED MAP OUTPUT CONSTANTS.
      DO L = 1,LSM
         ALSL(L) = ALOG(SPL(L))
      END DO
! close up shop
       call ext_ncd_ioclose ( DataHandle, Status )
!
!HC WRITE IGDS OUT FOR WEIGHTMAKER TO READ IN AS KGDSIN
        if(me.eq.0)then
        igdout=110
!        open(igdout,file='griddef.out',form='unformatted'
!     +  ,status='unknown')
        if(maptype .eq. 1)THEN  ! Lambert conformal
          WRITE(igdout)3
          WRITE(6,*)'igd(1)=',3
          WRITE(igdout)im
          WRITE(igdout)jm
          WRITE(igdout)LATSTART
          WRITE(igdout)LONSTART
          WRITE(igdout)8
!          WRITE(igdout)CENLON
          WRITE(igdout)STANDLON
          WRITE(igdout)DXVAL
          WRITE(igdout)DYVAL
          WRITE(igdout)0
          WRITE(igdout)64
          WRITE(igdout)TRUELAT2
          WRITE(igdout)TRUELAT1
          WRITE(igdout)255
        ELSE IF(MAPTYPE .EQ. 2)THEN  !Polar stereographic
          WRITE(igdout)5
          WRITE(igdout)im
          WRITE(igdout)jm
          WRITE(igdout)LATSTART
          WRITE(igdout)LONSTART
          WRITE(igdout)8
          WRITE(igdout)CENLON
          WRITE(igdout)DXVAL
          WRITE(igdout)DYVAL
          WRITE(igdout)0
          WRITE(igdout)64
          WRITE(igdout)TRUELAT2  !Assume projection at +-90
          WRITE(igdout)TRUELAT1
          WRITE(igdout)255
        !  Note: The calculation of the map scale factor at the standard
        !        lat/lon and the PSMAPF
        ! Get map factor at 60 degrees (N or S) for PS projection, which will
        ! be needed to correctly define the DX and DY values in the GRIB GDS
          if (TRUELAT1 .LT. 0.) THEN
            LAT = -60.
          else
            LAT = 60.
          end if

          CALL MSFPS (LAT,TRUELAT1*0.001,PSMAPF)

        ELSE IF(MAPTYPE .EQ. 3)THEN  !Mercator
          WRITE(igdout)1
          WRITE(igdout)im
          WRITE(igdout)jm
          WRITE(igdout)LATSTART
          WRITE(igdout)LONSTART
          WRITE(igdout)8
          WRITE(igdout)latlast
          WRITE(igdout)lonlast
          WRITE(igdout)TRUELAT1
          WRITE(igdout)0
          WRITE(igdout)64
          WRITE(igdout)DXVAL
          WRITE(igdout)DYVAL
          WRITE(igdout)255
        ELSE IF(MAPTYPE.EQ.6 )THEN  ! ARW rotated lat/lon grid
          WRITE(igdout)205
          WRITE(igdout)im
          WRITE(igdout)jm
          WRITE(igdout)LATSTART
          WRITE(igdout)LONSTART+360.
          WRITE(igdout)136
          WRITE(igdout)CENLAT
          WRITE(igdout)CENLON
          WRITE(igdout)DXVAL
          WRITE(igdout)DYVAL
          WRITE(igdout)64
          WRITE(igdout)LATLAST
          WRITE(igdout)LONLAST
          WRITE(igdout)0

         END IF

! following for hurricane wrf post
     
          open(10,file='copygb_hwrf.txt',form='formatted',status='unknown')
          idxvald = abs(LONLAST-LONSTART)/(im-2)
          idyvald = abs(LATLAST-LATSTART)/(jm-2)
          write(10,1010) IM-1,JM-1,LATSTART,LONSTART,LATLAST,LONLAST,   &
                         idxvald,idyvald

1010      format('255 0 ',2(I4,x),I8,x,I9,x,'136 ',I8,x,I9,x,            &
                 2(I8,x),'0')
          close (10)
        end if
!     
      DEALLOCATE (THV)

      RETURN
      END

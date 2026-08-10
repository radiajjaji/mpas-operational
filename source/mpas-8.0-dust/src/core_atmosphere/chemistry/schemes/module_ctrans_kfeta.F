MODULE module_ctrans_kfeta

!************************************************************************
! This computer software was prepared by Battelle Memorial Institute,
! hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with
! the Department of Energy (DOE). NEITHER THE GOVERNMENT NOR THE
! CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
! LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to Compute Convective Transport following KFeta scheme 
! * Author: Chun Zhao 
! Last update: July 2011
!
! Contact:
! Chun Zhao, PhD
! Scientist
! Pacific Northwest National Laboratory
! P.O. Box 999, MSIN K9-30
! Richland, WA, 99352
! Phone: (509) 371-6372
! Email: Chun.Zhao@pnnl.gov
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
! * Zhao, C., Y. Wang, Y. Choi, T. Zeng: 
!   Summertime impact of convective transport and lightning NOx production
!   over North America: modeling dependence on meteorological simulations,
!   ACP, 9, 4315-4327, 2009. 
!
! Additional information:
! *  www.pnl.gov/atmospheric/research/wrf-chem
!
!************************************************************************
!===========================================================================
!CUMULUS CONVECTION PARAMETERIZATION USING NEW-KFETA SCHEME 
!
! New Kain-Fritsch Scheme: The scheme is based on Kain and Fritsh(1990) and 
! Kain and Fritsch(1993), but has been modified based on testing within the Eta
! model. As with the original KF scheme, it utilizes a simple cloud model with 
! moist updrafts and downdrafts, including the effects of detrainment, 
! entrainment, and relatively simple microphysics. It differs from the original
! KF scheme in the following ways:
!   1.A minimum entrainment rate is imposed to suppress widespread convection
!   in marginally unstable, relatively dry environments.
!   2.Shallow (non-precipitating) convection is allowed for any updraft that 
!   does not reach minimum cloud depth for precipitating clouds; this minimum 
!   depth varies as a function of cloud-base temperature.
!   3.The entrainment rate is allowed to vary as a function of low-level
!   convergence.
!   4.Downdraft changes:
!     ---Source layer is the entire 150-200 mb deep layer just above cloud base.
!     ---Mass flux is specified as a fraction of updraft mass flux at cloud
!     base. Fraction is a function of source layer RH rather than wind shear 
!     or other parameters, i.e., old precipitation efficiency relationship not
!     used.
!     ---Detrainment is specified to occur in updraft source layer and below.
!
!                                                   ----czhao  03/2006 
!IN THIS SUBROUTINE, IT INCLUDES FOLLOWING PARTS:
!  + CUMULUS SCAVENGING FOR SOLUBLE GASES  
!  + NEW_KFETA SCHEME
!                                         ---czhao    03/2006
!=============================================================================

USE module_mosaic_wetscav, only:e_ice, compute_l2g
USE mpas_atmchem_vars, only:is_aerosol
USE mpas_atmchem_utilities

CONTAINS

!-------------------------------------------------------------
   subroutine kfeta_cup_transport(dt,dx,rho_phy,alt,chem, &
              u,v,t_phy,dz8w,p_phy,z,                     &
              umf3d,uer3d,udr3d,dmf3d,der3d,ddr3d,wu3d,   &
              ldb2d,let2d,ksl2d,lfs2d,klc2d,              &
              cutop,cubot,nca2d,ishall2d,check2d,         &
              num_chem,chem_opt,scav_opt,                 &
              ids,ide, jds,jde, kds,kde,                  &
              ims,ime, jms,jme, kms,kme,                  &
              its,ite, jts,jte, kts,kte                   )

!-------------------------------------------------------------
   implicit none
!-------------------------------------------------------------
! umf3d      kg  s-1   cloud base mass flux (updraft)
! dmf3d      kg  s-1   downdraft flux(negative value)
! uer3d      kg  s-1   entrainment rate (updraft)
! udr3d      kg  s-1   detrainment rate (updraft)
! der3d      kg  s-1   entrainment rate (downdraft,negative)
! ddr3d      kg  s-1   detrainment rate (downdraft)
! wu3d       m   s-1   updraft velocity
! ldb2d      unitless     downdraft end level
! lfs2d      unitless     downdraft originating level
! ksl2d      unitless     updraft source top layer
! klc2d      unitless     updraft source bot layer
! let2d      unitless     topest cloud detrainment layer
! cutop      unitless     cloud top level
! cubot      unitless     cloud bot level
! ishall2d   unitless     switch for deep or shallow convection
! check2d    unitless     check if deep convection is really occur
! nca2d      unitless     convection relaxation counter 
!-------------------------------------------------------------

   integer,      intent(in   ) ::                               &
                                  ids,ide, jds,jde, kds,kde,    &
                                  ims,ime, jms,jme, kms,kme,    &
                                  its,ite, jts,jte, kts,kte,    &
                                  num_chem

   character(len=512), intent(in) :: chem_opt
   logical, intent(in) :: scav_opt

   real,  dimension( ims:ime , kms:kme , jms:jme )         ,    &
          intent(in   ) ::                                      &
                         u,    &
                         v,    &
                     t_phy,    &
                     z,        &       !heights of model levels in m
                     p_phy,    &
                      dz8w,    &       !layer thickness in m
                      alt,     &       !inversed dry air density in m3/kg
                   rho_phy             !moisture air density in kg/m3
   real, intent(in   ) :: dt
   real, dimension( ims:ime , jms:jme ), intent(in   ) :: dx
   real, dimension( ims:ime , kms:kme , jms:jme, num_chem ),    &
         intent(inout) ::               chem

   real, dimension (ims:ime,kms:kme,jms:jme), intent(in) ::    &
     umf3d,uer3d,udr3d,dmf3d,der3d,ddr3d,wu3d
   integer, dimension(ims:ime,jms:jme),intent(in) ::    &
     ldb2d,let2d,ksl2d,lfs2d,klc2d,cutop,cubot,ishall2d,check2d

   real, dimension(ims:ime,jms:jme),intent(in) :: nca2d

!  local variable

   character(len=512) :: msg

   real dtin
   real,dimension(kms:kme)  ::  airdens_gas   ! in mole/m3
   real,dimension(kms:kme)  ::  airdens_aer   ! in kg/m3
   real,dimension(kms:kme)  ::  airdens   ! in kg/m3
   real,dimension(kms:kme)  ::  air_ratio   ! dry_air to moisture_air ratio 
!  real,dimension(kms:kme)  ::     &
!    fluxup0,fluxdown0,uentr0,udetr0,dentr0,ddetr0
   real,dimension(kms:kme)  ::     &
     fluxup,fluxdown,uentr,udetr,dentr,ddetr,difmass,submass
   real,dimension(kms:kme)  ::    tracup,tracdw,tracsub1,tracsub2
   real   tracsum
   real uairmass,uenmass,udemass,denmass,ddemass,tracermass1,tracermass2, &
    tracermass,volume,density,airmass,airmass1,airmass2,fluxtra,fluxrat
   real,dimension(kms:kme)  :: xs0  !original mixing ratio
   real,dimension(kms:kme)  :: xs1  !new mixing ratio after transport
   real,dimension(kms:kme)  :: txs  !temporary mixing ratio
   real, parameter :: avogadro = 6.0225e23  ! molecules/mole
   real, parameter :: airmw = 28.98  ! g/mole
   integer i, j, k,  nv 
   integer ztop              ! cloud top layer
   integer zb                ! cloud bot layer
   integer zlc               ! updraft originating layer
   integer zldb              ! downdraft end layer   
   integer zlfs              ! downdraft originating layer
   integer zlet              ! topest cloud detrainment layer
   integer zksl              ! the updraft source top layer 
   integer ishall            ! switch for deep or shallow convection
   real   sumass1,sumass2,totracer  ! check mass balance 
   real,dimension(kms:kme,1:2)  :: scavf  !updraft and downdraft scavenging coefficient
   real :: dx2
   real,dimension(kms:kme)  :: dzq,tmp,pmid
   real,dimension(kms:kme)  :: vud,vdd 
   real,dimension(kms:kme,num_chem)  :: sol_f 
   real,dimension(num_chem) ::  sol_f0 
   real    :: epsilc

   epsilc=1.e-30
   dtin=dt

   do 1000 j=jts,jte
   do 1000 i=its,ite
   !obtain the soluble information for scavenging 
   dx2=dx(i,j)*dx(i,j)
   do k=kts,kte
      sol_f(k,1:num_chem)=1.0
      sol_f0(1:num_chem)=1.0
      soluble_select: select case(trim(chem_opt))
      case ( "bin4aer_dust_aq" , "bin4aer_dust"  , &
             "bin5aer_dust_aq" , "bin5aer_dust"  , &
             "bin8aer_dust_aq" , "bin8aer_dust"  , &
             "bin9aer_dust_aq" , "bin9aer_dust"  , &
             "bin10aer_dust_aq", "bin10aer_dust" )
          call mosaic_soluble(sol_f0,num_chem,chem(i,k,j,1:num_chem))
          sol_f(k,1:num_chem)=sol_f0
      case default
      end select soluble_select
   enddo !k

   do 900 nv=1,num_chem

   do k=kts,kte
   if(.not.is_aerosol(nv))then !gases in ppmv
       xs0(k)=1.e-6*max(epsilc,chem(i,k,j,nv))  ! ppmv to mole/mole
       xs1(k)=xs0(k)
       txs(k)=xs1(k)
   else !aerosol in ug/kg dry-air 
       xs0(k)=max(epsilc,chem(i,k,j,nv))
       xs1(k)=xs0(k)
       txs(k)=xs1(k)
   endif
   airdens_gas(k)=rho_phy(i,k,j)*1000.0/airmw !kg/m3 to mole/m3
   airdens_aer(k)=1.0/alt(i,k,j) !dry air density kg/m3
   air_ratio(k)=(1.0/alt(i,k,j))/rho_phy(i,k,j) !dry air to moist air ratio
   dzq(k)=dz8w(i,k,j)
   enddo

         !----------------------------------------------------------------   
         !before convection,check tracer mass balance
   sumass1=0.
   do k=kts,kte
   if(.not.is_aerosol(nv))then !gases in mole/mole 
       sumass1=sumass1+dzq(k)*airdens_gas(k)*dx2*txs(k)
       if (sumass1 .gt. 1.e20) then  ! larger than 1.e10 mole
           write(msg,*) 'tracer mass > 1.e20 mole',nv,i,j,k
           call chemistry_message( msg )
           write(msg,*) sumass1
           call chemistry_message( msg )
       endif 
   else   ! aerosol in ug/kg-dryair 
       sumass1=sumass1+dzq(k)*airdens_aer(k)*dx2*txs(k)
       if (sumass1 .gt. 1.e30) then  ! larger than 1.e30 ug 
           write(msg,*) 'tracer mass > 1.e30 ug',nv,i,j,k
           print *, 'tracer mass > 1.e30 ug',nv,i,j,k
           write(msg,*) sumass1
           print *, sumass1
       endif
   endif
   enddo

   !----------------------------------------------------------------   
   !0 pre-initialize 
   !----------------------------------------------------------------   
   do k=kts,kte
   if(.not.is_aerosol(nv))then !for gas, in mole/m3
       airdens(k)=airdens_gas(k)
   else  !for aerosol, in kg/m3 dry-air 
       airdens(k)=airdens_aer(k)
   endif
   tmp(k)=t_phy(i,k,j)
   pmid(k)=p_phy(i,k,j)
   enddo

   do k=kts,kte
   scavf(k,1:2)=0.0
   tracup(k)=0.
   tracdw(k)=0.
   tracsub1(k)=0.
   tracsub2(k)=0.
   submass(k)=0.
   enddo

   ztop=cutop(i,j)
   zb=cubot(i,j)
   zlc=klc2d(i,j)
   zlfs =lfs2d(i,j)
   zldb =ldb2d(i,j)
   zlet =let2d(i,j)
   zksl= ksl2d(i,j)
   ishall=ishall2d(i,j)

   do k=kts,kte
   if(.not.is_aerosol(nv))then !gases in mole/mole 
       fluxup(k)=umf3d(i,k,j)*1.e3/airmw !kg/s to mole/s
       fluxdown(k)=-1.0*dmf3d(i,k,j)*1.e3/airmw !kg/s to mole/s
              !if (fluxdown(k).lt.0) fluxdown(k)=0.0
       uentr(k)=uer3d(i,k,j)*1.e3/airmw !kg/s to mole/s
       udetr(k)=udr3d(i,k,j)*1.e3/airmw !kg/s to mole/s
       dentr(k)=-1.0*der3d(i,k,j)*1.e3/airmw !kg/s to mole/s
       ddetr(k)=ddr3d(i,k,j)*1.e3/airmw !kg/s to mole/s
   else 
       fluxup(k)=umf3d(i,k,j)*air_ratio(k) !kg/s air to kg/s dry-air
       fluxdown(k)=-1.0*dmf3d(i,k,j)*air_ratio(k) !kg/s air to kg/s dry-air
              !if (fluxdown(k).lt.0) fluxdown(k)=0.0
       uentr(k)=uer3d(i,k,j)*air_ratio(k) !kg/s air to kg/s dry-air
       udetr(k)=udr3d(i,k,j)*air_ratio(k) !kg/s air to kg/s dry-air
       dentr(k)=-1.0*der3d(i,k,j)*air_ratio(k) !kg/s air to kg/s dry-air
       ddetr(k)=ddr3d(i,k,j)*air_ratio(k) !kg/s air to kg/s dry-air
   endif
   if (fluxup(k).gt.0) then 
       vud(k)=min(50.0,max(wu3d(i,k,j),0.0))
       vdd(k)=min(50.0,max(vud(k)*fluxdown(k)/fluxup(k),0.0))
   else 
       vud(k)=0.0
       vdd(k)=0.0
   endif
   enddo

             !make convective output consistent
   if (maxval(abs(fluxup(kts:kte)))==0 ) then
       ztop=0 
       zb=ztop
   endif
   if (maxval(abs(fluxdown(kts:kte)))==0) then
       zlfs=0
       zldb=zlfs
   endif

!------------------------------------------------------------------
!...do not allow any cloud from this layer if:
!
!...            1.) if there is no cape, or 
!...            2.) cloud top is at model level just above lcl, or
!...            3.) cloud top is within updraft source layer, or
!...            4.) cloud-top detrainment layer begins within 
!...                updraft source layer.
!------------------------------------------------------------------

   if((zb.lt.ztop).and.(check2d(i,j).ne.0))then !allowed convection block

   !================================================================   
   !1 start deep convective transport 
   !================================================================   
   !note:
   !ishall: switch for deep convection or shallow convection
   !for deep and shallow convection, the difference is just for downdraft
   !flux. the calculation of updraft is same. shallow convection has no downdraft.
   !----------------------------------------------------------------   

       if (ishall==0) then  
     !----------------------------------------------------------------   
     !1.0 calculate the scavenging coefficient for deep convective transport 
     !----------------------------------------------------------------   
           if (scav_opt) then 

               call cuscav_kf(i,j,nv,dzq,dx2,airdens,vud,vdd,tmp,pmid, &
                    num_chem,chem(i,kms:kme,j,1:num_chem),scavf,  &
                    ids,ide, jds,jde, kds,kde,                       &
                    ims,ime, jms,jme, kms,kme,                       &
                    its,ite, jts,jte, kts,kte                        ) 
          !apply the internal-mixing soluble fraction for scavenging aerosol 
               if(is_aerosol(nv))then
                   do k=kts,kte
                       scavf(k,1:2)=scavf(k,1:2)*sol_f(k,nv)
                   enddo
               endif
          !set downdraft scavenging zero for now
               do k=kts,kte
                   scavf(k,2)=0.0
               enddo
           endif

     !----------------------------------------------------------------   
     !1.1 compute updraft effect in deep convection                  
     !----------------------------------------------------------------   
           do k=kts,kte
               difmass(k) = 0.0
               submass(k) = 0.
           enddo

           do k=zlc+1,ztop-1
               if (fluxup(k) .eq.0 ) then
                   write(msg,*) 'wrong fluxup',i,j,k,ztop,zlc
                   print *, 'wrong fluxup',i,j,k,ztop,zlc
                   stop
               endif
          !    fluxup ajust for mass consistence
               uenmass=uentr(k)*dtin ! mole of entrain air
               udemass=udetr(k)*dtin ! mole of detrain air

               if(uenmass.ne.udemass) then 
                   if((fluxup(k)-fluxup(k-1))/(uenmass-udemass).gt.1.e-5.and.  &
                       (fluxup(k)-fluxup(k-1))/(uenmass-udemass).lt.0 ) then
                       write(msg,*) 'updraft ajustment at',i,j,k
                       print *, 'updraft ajustment at',i,j,k
                       write(msg,*) fluxup(k),fluxup(k+1),uenmass,udemass
                       print *, fluxup(k),fluxup(k+1),uenmass,udemass
                       fluxup(k)=fluxup(k-1)+uenmass-udemass
                   endif
               endif
           enddo !k

     !----------------------------------------------------------------   
     !1.1.1 calculate tracer tendence in lc layer
     !----------------------------------------------------------------   
         !gases in mole/mole or aerosol in ug/kg-dry air
           fluxtra=fluxup(zlc)*txs(zlc)*dtin  ! mole or ug of tracer
           uairmass=fluxup(zlc)*dtin  ! mole or kg of air mass
           volume=dzq(zlc)*dx2     !m3
           density=airdens(k)  !mole or kg /m3
           tracup(zlc)=-fluxtra
           airmass=volume*density-uairmass !mole or kg
           difmass(zlc)=-uairmass  !save for subsidence
           fluxrat=txs(zlc)*(1.0-scavf(zlc,1)) ! mole/mole or ug/kg dry-air

     !----------------------------------------------------------------   
     !1.1.2 calculate other under cloud top layers
     !----------------------------------------------------------------   
           do k=zlc+1,ztop-1
               uairmass=fluxup(k-1)*dtin ! mole or kg of air mass from lower level
               uenmass=uentr(k)*dtin ! mole or kg of entrain air
               udemass=udetr(k)*dtin ! mole or kg of detrain air
               volume=dzq(k)*dx2 ! m3
               density=airdens(k) ! mole or kg /m3
               airmass=volume*density-uenmass+udemass  !mole or kg
               difmass(k)=-uenmass+udemass  !save for subsidence
               fluxtra=uairmass*fluxrat+uenmass*txs(k)*(1.-scavf(k,1))-udemass*fluxrat !mole gas or ug aerosol
               tracup(k)=udemass*fluxrat-uenmass*txs(k) !mole gas or ug aerosol
               fluxrat=fluxtra*(1.-scavf(k,1))/(fluxup(k)*dtin) !mole/mole or ug/kg dry-air
           enddo

     !----------------------------------------------------------------   
     !1.1.3 calculate the top layer of cloud, no entrainment
     !----------------------------------------------------------------   
           uairmass=fluxup(ztop-1)*dtin  !mole or kg airmass from lower level
           volume=dzq(ztop)*dx2
           density=airdens(ztop)
           airmass=volume*density+uairmass ! mole or kg
           difmass(ztop)=uairmass !mole or kg dry-air
           submass(ztop)=difmass(ztop) !mole or kg
           tracup(ztop)=uairmass*fluxrat !mole gas or ug aerosol

     !----------------------------------------------------------------   
     !1.1.4 check mass balance for updraft
     !----------------------------------------------------------------   
           tracsum=0.0
           do k=zlc,ztop
           tracsum=tracsum+tracup(k)
           enddo

           if( (tracsum/sumass1) .gt. 1.e-5) then
               write(msg,*) 'tracer mass is not balance during updraft!!!'
               print *, 'tracer mass is not balance during updraft!!!'
               write(msg,*) tracsum,i,j,nv
               print *, tracsum,i,j,nv
               write(msg,*) tracsum/sumass1
               print *, tracsum/sumass1
               stop
           endif

!          do k=kts,kte
!           volume=dzq(k)*dx2
!           density=airdens(k)
!           airmass=volume*density+difmass(k) ! mole or kg air
!           tracermass=volume*density*txs(k)+tracup(k) ! mole gas or ug aerosol
!          enddo

     !------------------------------------------------------------
     !1.1.5 calculate tracer concentration due to updraft in deep convection
     !------------------------------------------------------------
           do k=zlc,ztop
               !gases in mole/mole and aeroso in ug/kg dry-air 
               volume=dzq(k)*dx2
               density=airdens(k)
               airmass=volume*density+difmass(k) ! mole or kg air
               tracermass=volume*density*txs(k)+tracup(k) ! mole gas or ug aerosol
               xs1(k)=tracermass/airmass ! mole/mole or ug/kg
               txs(k)=xs1(k)
           enddo

     !----------------------------------------------------------------   
     !1.1.6 check mass balance for column 
     !----------------------------------------------------------------   
           sumass2=0.
           do k=kts,kte
               !gas in mole/mole and aerosol in ug/kg
               if (.not.is_aerosol(nv)) then 
                   sumass2=sumass2+xs1(k)*(dzq(k)*airdens(k)*dx2+difmass(k))
                   if (sumass2 .gt. 1.e20) then  ! larger than 1.e10 mole
                       write(msg,*) 'tracer mass2 > 1.e20 mole',nv,i,j,k
                       print *, 'tracer mass2 > 1.e20 mole',nv,i,j,k
                       write(msg,*) sumass2
                       print *, sumass2
                   endif
               else   ! aerosol in ug/kg-dryair 
                   sumass2=sumass2+xs1(k)*(dzq(k)*airdens(k)*dx2+difmass(k))
                   if (sumass2 .gt. 1.e30) then  ! larger than 1.e30 ug 
                       write(msg,*) 'tracer mass2 > 1.e30 ug',nv,i,j,k
                       call chemistry_message(msg)
                       write(msg,*) sumass2
                       call chemistry_message(msg)
                   endif
               endif
           enddo

           if ( ((sumass2-sumass1)/sumass1) .gt. 1.e-5) then
               write(msg,*) 'tracer mass is not balance after updraft!!!!!'
               call chemistry_message(msg)
               write(msg,*) sumass1,sumass2
               call chemistry_message(msg)
               write(msg,*) (sumass2-sumass1)/sumass1
               call chemistry_message(msg)

               totracer=0.
               do k=kts,kte
                  totracer=totracer+tracup(k)
                  write(msg,*)nv,is_aerosol(nv),k,tracup(k)
                  call chemistry_message(msg)
               enddo
               write(msg,*) 'total tracer tendence',i,j,ishall,totracer
               call chemistry_message(msg)
               stop
           endif

     !----------------------------------------------------------------   
     !1.2 calculate effect of compensating due to updraft in deep convection 
     !----------------------------------------------------------------   

     !-----------------------------------------------   
     !1.2.1 calculate the top layer of cloud
     !-----------------------------------------------   
           tracsub1(ztop)=submass(ztop)*txs(ztop)
           xs1(ztop)=xs1(ztop)
           txs(ztop)= xs1(ztop)

     !-----------------------------------------------   
     !1.2.2 calclulate iner layers between cloud top and original level
     !-----------------------------------------------   
           do k=ztop-1,zlc+1,-1
            !a. mixing with upper layer downward flux
           volume=dzq(k)*dx2  ! m3
           density=airdens(k) ! mole or kg/m3
           airmass=volume*density+difmass(k)+submass(k+1) ! mole or kg
            !tracermass=(volume*density+difmass(k))*txs(k)+tracsub1(k+1)*(1.-scavf(k,1))
           tracermass=(volume*density+difmass(k))*txs(k)+tracsub1(k+1) !mole or ug
           xs1(k)=tracermass/airmass !mole/mole or ug/kg
           txs(k)=xs1(k)
            !b. consider the downward flux in the layer
           airmass=volume*density
           submass(k)=difmass(k)+submass(k+1)
            !check submass 
           if (submass(k).lt.0) then
               write(msg,*) 'error subside mass flux'
               call chemistry_message(msg)
               write(msg,*) i,j,k,submass(k),difmass(k),submass(k+1)
               call chemistry_message(msg)
               stop
           endif
           tracsub1(k) =submass(k)*txs(k)
           tracermass=tracermass-tracsub1(k)
           xs1(k)=tracermass/airmass
           txs(k)=xs1(k)
           enddo

     !-----------------------------------
     !1.2.3 calculate original updraft level
     !-----------------------------------
           volume=dzq(zlc)*dx2 !m3
           density=airdens(zlc) !mole/m3 or kg/m3
           airmass=volume*density !mole or kg
            !tracermass=(volume*density+difmass(zlc))*txs(zlc)+tracsub1(zlc+1)*(1.-scavf(k,1))
           tracermass=(volume*density+difmass(zlc))*txs(zlc)+tracsub1(zlc+1) !mole or ug
           xs1(zlc)=tracermass/airmass
           txs(zlc)=xs1(zlc)

     !--------------------------------------
     !1.2.4 check mass balance for column 
     !--------------------------------------
           sumass2=0.
           do k=kts,kte
           sumass2=sumass2+xs1(k)*(dzq(k)*airdens(k)*dx2)
           enddo

           if ( ((sumass2-sumass1)/sumass1) .gt. 1.e-5) then
               write(msg,*) 'tracer mass is not balance after subsidence!!!!!'
               call chemistry_message(msg)
               write(msg,*) sumass1, sumass2,(sumass2-sumass1)/sumass1
               call chemistry_message(msg)

               totracer=0.
               do k=kts,kte
               totracer=totracer+tracsub1(k)
               write(msg,*) k,tracsub1(k)
               call chemistry_message(msg)
               enddo
               write(msg,*) 'total tracer tendence',i,j,ishall,nv,is_aerosol(nv),totracer
               call chemistry_message(msg)
               stop
           endif

     !----------------------------------------------------------------   
     !1.3 compute downdraft effect in deep convection 
     !----------------------------------------------------------------   

           if(zlfs.gt.zldb ) then  ! allowed downdraft for deep convec.   

     !----------------------------------------------------------------   
     !1.3.0 initializing 
     !----------------------------------------------------------------   
               do k=kts,kte
               difmass(k) = 0.0
               submass(k) = 0.0
               enddo

               do k=zlfs-1,zldb+1,-1
               if (fluxdown(k).eq.0 ) then
                   write(msg,*) 'wrong fluxdown',i,j,k,zlfs,zldb
                   print *, 'wrong fluxdown',i,j,k,zlfs,zldb
                   stop
               endif
           !fluxdown ajust for mass consistence
               denmass=dentr(k)*dtin ! mole of entrain air
               ddemass=ddetr(k)*dtin ! mole of detrain air
               if(denmass.ne.ddemass) then
                   if((fluxdown(k-1)-fluxdown(k))/(denmass-ddemass).gt.1.e-5.and.  &
                       (fluxdown(k-1)-fluxdown(k))/(denmass-ddemass).lt.0 ) then
                       write(msg,*) 'downdraft ajustment at',i,j,k
                       print *, 'downdraft ajustment at',i,j,k
                       write(msg,*) fluxdown(k),fluxdown(k+1),denmass,ddemass
                       print *, fluxdown(k),fluxdown(k+1),denmass,ddemass
                       fluxdown(k-1)=fluxdown(k)+denmass-ddemass
                   endif
               endif
               enddo !k


     !----------------------------------------------------------------   
     !1.3.1 calculate the lfs layer tracer tendence 
     !----------------------------------------------------------------   
               uairmass=fluxdown(zlfs)*dtin  ! mole or kg of air mass
               volume=dzq(zlfs)*dx2 !m3
               density=airdens(zlfs) !mole or kg /m3
               airmass=volume*density-uairmass !mole or kg
               difmass(zlfs)=-uairmass !mole or kg, save for uplifting
               fluxrat=txs(zlfs)*(1.-scavf(zlfs,2))
               tracdw(zlfs)=-uairmass*txs(zlfs) !mole or ug of tracer

     !----------------------------------------------------------------   
     !1.3.2 calculate layer between lfs and ldb 
     !----------------------------------------------------------------   
               do k=zlfs-1,zldb+1,-1
                   uairmass=fluxdown(k+1)*dtin ! mole or kg airmass from upper level
                   denmass=dentr(k)*dtin ! mole or kg entrain air
                   ddemass=ddetr(k)*dtin ! mole or kg detrain air
                   volume=dzq(k)*dx2  ! m3
                   density=airdens(k)  !mole or kg /m3
                   airmass=volume*density-denmass+ddemass !mole or kg
                   difmass(k)=-denmass+ddemass !save for uplifting
                   fluxtra=uairmass*fluxrat+denmass*txs(k)*(1.-scavf(k,2))-ddemass*fluxrat !mole gas or ug aerosol
                   tracdw(k)=ddemass*fluxrat-denmass*txs(k) !mole gas or ug aerosol
                   fluxrat=(fluxtra/(fluxdown(k)*dtin))*(1.-scavf(k,2))
               enddo

     !----------------------------------------------------------------   
     !1.3.3 calculate ldb layer
     !----------------------------------------------------------------   
               uairmass=fluxdown(zldb+1)*dtin
               volume=dzq(zldb)*dx2
               density=airdens(zldb)
               airmass=volume*density+uairmass
               difmass(zldb)=uairmass
               submass(zldb)=difmass(zldb)
               tracdw(zldb)=uairmass*fluxrat 

     !----------------------------------------------------------------   
     !1.3.4 check mass balance for downdraft 
     !----------------------------------------------------------------   
               tracsum=0.0
               do k=zldb,zlfs
                   tracsum=tracsum+tracdw(k)
               enddo

               if( (tracsum/sumass1).gt.1.e-5) then
                   write(msg,*) 'tracer mass is not balance during downdraft!!!'
                   print *, 'tracer mass is not balance during downdraft!!!'
                   write(msg,*) tracsum,i,j,nv
                   print *, tracsum,i,j,nv
                   write(msg,*) tracsum/sumass1
                   print *, tracsum/sumass1
                   stop
               endif

     !----------------------------------------------------------------   
     !1.3.5 calculate tracer concentration due to downdraft in deep convection 
     !----------------------------------------------------------------   
               do k=zlfs,zldb,-1
                   volume=dzq(k)*dx2 !m3
                   density=airdens(k) !mole or kg /m3
                   airmass=volume*density+difmass(k) ! mole or kg
                   tracermass=volume*density*txs(k)+tracdw(k) ! mole or ug
                   xs1(k)=tracermass/airmass ! mole/mole or ug/kg
                   txs(k)=xs1(k)
               enddo

     !----------------------------------------------------------------   
     !1.3.6 check mass balance for column 
     !----------------------------------------------------------------   
               sumass2=0.
               do k=kts,kte
                   sumass2=sumass2+xs1(k)*(dzq(k)*airdens(k)*dx2+difmass(k))
               enddo

               if ( ((sumass2-sumass1)/sumass1) .gt. 1.e-5) then
                   write(msg,*) 'tracer mass is not balance after downdraft!!!!'
                   write(msg,*) sumass1, sumass2,(sumass2-sumass1)/sumass1

                   totracer=0.
                   do k=kts,kte
                       totracer=totracer+tracdw(k)
                   write(msg,*) k,tracdw(k)
                   enddo
                   write(msg,*) 'total tracer tendence',i,j,nv,is_aerosol(nv),ishall, totracer
                   stop
               endif

     !----------------------------------------------------------------   
     !1.4 calculate the effect of compensating due to downdraft in deep convection 
     !----------------------------------------------------------------   

     !----------------------------------------------------------------   
     !1.4.1 calculate the ldb layer 
     !----------------------------------------------------------------   
               tracsub2(zldb)=submass(zldb)*txs(zldb)
               xs1(zldb)=xs1(zldb)
               txs(zldb)=txs(zldb)

     !----------------------------------------------------------------   
     !1.4.2 calculate the iner layer between cloud top and lfs layer 
     !----------------------------------------------------------------   
               do k=zldb+1,zlfs-1
               !a. mixing with lower layer upward flux
                   volume=dzq(k)*dx2
                   density=airdens(k)
                   airmass=volume*density+difmass(k)+submass(k-1) ! mole or kg
                   !tracermass=(volume*density+difmass(k))*txs(k)+tracsub2(k-1)*(1.-scavf(k,2))
                   tracermass=(volume*density+difmass(k))*txs(k)+tracsub2(k-1) !mole or ug
                   xs1(k)=tracermass/airmass !mole/mole or ug/kg
                   txs(k)=xs1(k)
                   !b. consider the upward flux in the layer
                   airmass=volume*density
                   submass(k)=difmass(k)+submass(k-1)
                   ! check submass 
                   if (submass(k).lt.0) then
                       write(msg,*) 'error uplift mass flux'
                       print*,'error uplift mass flux'
                       write(msg,*) i,j,k,nv,submass(k),difmass(k),submass(k-1)
                       print*,i,j,k,nv,submass(k),difmass(k),submass(k-1)
                       stop
                   endif
                   tracsub2(k)=submass(k)*txs(k)
                   tracermass=tracermass-tracsub2(k)
                   xs1(k)=tracermass/airmass
                   txs(k)=xs1(k)
               enddo

     !----------------------------------------------------------------   
     !1.4.3 calculate the lfs layer
     !----------------------------------------------------------------   
               volume=dzq(zlfs)*dx2
               density=airdens(zlfs)
               airmass=volume*density ! mole or kg 
            !tracermass=(volume*density+difmass(zlfs))*txs(zlfs)+tracsub2(zlfs-1)*(1.-scavf(k,2))
               tracermass=(volume*density+difmass(zlfs))*txs(zlfs)+tracsub2(zlfs-1) !mole or ug
               xs1(zlfs)=tracermass/airmass
               txs(zlfs)=xs1(zlfs)

     !----------------------------------------------------------------   
     !1.4.4 check mass balance for column 
     !----------------------------------------------------------------   
               sumass2=0.
               do k=kts,kte
                   sumass2=sumass2+xs1(k)*(dzq(k)*airdens(k)*dx2)
               enddo

               if ( ((sumass2-sumass1)/sumass1) .gt. 1.e-5) then
                   write(msg,*) 'tracer mass is not balance after uplifting!!!!'
                   print *, 'tracer mass is not balance after uplifting!!!!'
                   write(msg,*) sumass1, sumass2, (sumass2-sumass1)/sumass1
                   print *, sumass1, sumass2, (sumass2-sumass1)/sumass1

                   totracer=0.
                   do k=kts,kte
                   totracer=totracer+tracsub2(k)
                   write(msg,*) k,tracsub2(k)
                   print*,k,tracsub2(k)
                   enddo
                   write(msg,*) 'total tracer tendence',i,j,ishall,nv,is_aerosol(nv),totracer
                   print *, 'total tracer tendence',i,j,ishall,nv,is_aerosol(nv),totracer
                   stop
               endif

           endif !lfs, downdraft

       endif ! ishall, for deep convection

   !================================================================   
   !2 start shallow convective transport 
   !================================================================   
       if ( ishall==1) then  ! for shallow convection

     !----------------------------------------------------------------   
     !2.1 calculate tracer tendence between lc layer and cloud top 
     !----------------------------------------------------------------   
     !note:the mixing between each two layers every time 
     !----------------------------------------------------------------   

           do k=zlc,ztop-1
               uairmass=fluxup(k)*dtin  ! mole or kg
               airmass1=dzq(k)*dx2*airdens(k) !mole or kg, current layer
               airmass2=dzq(k+1)*dx2*airdens(k+1) !mole or kg, upper layer
               tracermass1=airmass1*txs(k)       & ! current layer
                   -uairmass*txs(k)       & ! upward flux
                   +uairmass*txs(k+1)       ! downward flux
               tracermass2=airmass2*txs(k+1)     & ! upper layer
                   -uairmass*txs(k+1)     & ! downward flux
                   +uairmass*txs(k)         ! upward flux
               xs1(k)=tracermass1/airmass1   ! mole/mole or ug/kg
               xs1(k+1)=tracermass2/airmass2
               txs(k:k+1)=xs1(k:k+1)
           end do

     !----------------------------------------------------------------   
     !2.2 check tracer mass balance 
     !----------------------------------------------------------------   
           sumass2=0.
           do k=kts,kte
               sumass2=sumass2+xs1(k)*(dzq(k)*airdens(k)*dx2) ! mole or ug
           enddo

           if ( (sumass2-sumass1).gt.sumass1*1.e-5) then
               write(msg,*) 'tracer mass is not balance after shallow subside!'
               call chemistry_message(msg)
               write(msg,*) sumass1, sumass2
               call chemistry_message(msg)
               stop
           endif

       endif   !ishall, for shallow convection

   !================================================================   
   !3 check tracer mass balance due to convective transport 
   !================================================================   
       sumass2=0.
       do k=kts,kte
           sumass2=sumass2+xs1(k)*(dzq(k)*airdens(k)*dx2)
           !gas in mole/mole and aerosol in ug/kg
           if (.not.is_aerosol(nv)) then
               if (sumass2 .gt. 1.e20) then  ! larger than 1.e10 mole
                   write(msg,*) 'tracer mass2 > 1.e20 mole',nv,i,j,k,sumass2
                   call chemistry_message(msg)
               endif
           else   ! aerosol in ug/kg-dryair 
               if (sumass2 .gt. 1.e30) then  ! larger than 1.e30 ug 
                   write(msg,*) 'tracer mass2 > 1.e30 ug',nv,i,j,k,sumass2
                   call chemistry_message(msg)
               endif
           endif
       enddo


       if ( ((sumass2-sumass1)/sumass1) .gt. 1.e-5) then
           write(msg,*) 'tracer mass is not balance after convection!!!!!'
           call chemistry_message(msg)
           write(msg,*) sumass1, sumass2, (sumass2-sumass1)/sumass1
           call chemistry_message(msg)

           totracer=0.
           do k=zlc,ztop
               totracer=totracer+tracup(k)+tracsub1(k)+tracdw(k)+tracsub2(k)
               write(msg,*) k,tracup(k),tracsub1(k),tracdw(k),tracsub2(k)
               call chemistry_message(msg)
           enddo
           write(msg,*) 'total tracer tendence',i,j,ishall,nv,totracer
           call chemistry_message(msg)
           stop
       endif


   endif !check2d, for convective transport allowed

   do k=kts,kte
       if(.not.is_aerosol(nv))then !gases in ppmv
           chem(i,k,j,nv)=xs1(k)*1.e6  ! mole/mole to ppmv
       else !aerosol in ug/kg dry-air 
           chem(i,k,j,nv)=xs1(k) 
       endif
   enddo


   900   continue  !num_chem
   1000   continue  !i,j

 end subroutine kfeta_cup_transport 

!===========================================================================
 subroutine sorgam_soluble(sol_f0,num_chem,chem0)

 !calculate the soluble fraction of the internal-mixing aerosol based on hygroscopicity
 !this is not accurate, just for temporary application !!!! 

 use module_data_sorgam

 implicit none

   real, dimension(num_chem ),intent(inout)  :: sol_f0 
   integer, intent(in)  :: num_chem
   real, dimension(num_chem ),intent(in) :: chem0
   integer n,m,l,lmass
   real sumhygro,sumvol,hygro

          do n=1,ntype_aer
          do m=1,nsize_aer(n)
           sumhygro=0.
           sumvol=0. 
           do l=1,ncomp_aer(n)
             lmass=massptr_aer(l,m,n,ai_phase)
             sumhygro = sumhygro+hygro_aer(l,n)*chem0(lmass)/dens_aer(l,n)
             sumvol = sumvol+chem0(lmass)/dens_aer(l,n)
           enddo ! comp
           hygro=sumhygro/sumvol
           do l=1,ncomp_aer(n)
             lmass=massptr_aer(l,m,n,ai_phase)
             sol_f0(lmass)=min(max(hygro/1.0,0.0),1.0)
           enddo
          enddo
          enddo

 end subroutine sorgam_soluble

!===========================================================================
subroutine mosaic_soluble(sol_f0,num_chem,chem0)

!calculate the soluble fraction of the internal-mixing aerosol based on hygroscopicity
!this is not accurate, just for temporary application !!!! 

use module_data_binaer_asect

implicit none

real, dimension(num_chem ),intent(inout)  :: sol_f0
integer, intent(in)  :: num_chem
real, dimension(num_chem ),intent(in) :: chem0
integer n,m,l,lmass
real sumhygro,sumvol,hygro

do n=1,ntype_aer
    do m=1,nsize_aer(n)
        sumhygro=0.
        sumvol=0.
        do l=1,ncomp_aer(n)
            lmass=massptr_aer(l,m,n,ai_phase)
            sumhygro = sumhygro+hygro_aer(l,n)*chem0(lmass)/dens_aer(l,n)
            sumvol = sumvol+chem0(lmass)/dens_aer(l,n)
        enddo ! comp
        hygro=sumhygro/sumvol
        do l=1,ncomp_aer(n)
            lmass=massptr_aer(l,m,n,ai_phase)
            sol_f0(lmass)=min(max(hygro/1.0,0.0),1.0)
        enddo
    enddo
enddo

end subroutine mosaic_soluble

!===========================================================================
 subroutine cuscav_kf (ii,jj,n,dzq,dx2,airdens,vud,vdd,tmp,pmid,        &
                       num_chem,chem0,f,   &
                       ids,ide, jds,jde, kds,kde,                       &
                       ims,ime, jms,jme, kms,kme,                       &
                       its,ite, jts,jte, kts,kte                        )

!******************************************************************************
!  subroutine cuscav_kf computes f, the fraction of soluble tracer lost by 
!  scavenging in convective cloud updrafts. 
!
!  arguments as input:
!  (1 ) n    (integer) : tracer number
!
!  arguments as output:
!  (2 ) f    (real)  : fraction of tracer scavenged in cloud updraft [0-1]
!
   implicit none


   ! arguments
   integer,intent(in)  ::    &
                          ii,jj,n, &
                          ids,ide, jds,jde, kds,kde,    &
                          ims,ime, jms,jme, kms,kme,    &
                          its,ite, jts,jte, kts,kte
   real,intent(in),dimension(kms:kme)  :: dzq,tmp,pmid,airdens
   real,intent(in)  :: dx2
   real,intent(in),dimension(kms:kme)  :: vud,vdd
   integer, intent(in)  :: num_chem

   real,intent(inout),dimension(kms:kme,num_chem) :: chem0
   real,intent(inout),dimension(kms:kme,2) :: f

   ! local variables 
   integer            :: i, j, k,nn
   real               :: l2g, i2g, c_tot, f_l, f_i, kr, so2loss

   ! kc is the conversion rate from cloud condensate to precip [s^-1]
   real, parameter    :: kc        = 5.e-3

   real csr_conv(3) ! for three modes, aerosol activation rate
   real csr

   real, dimension (kms:kme) :: cldliq,cldice 
   real c_h2o

   csr_conv      = (/0.60, 0.99, 0.99/) ! see module_mosaic_wetscav.f

   !calculate the updraft velocity (m/s)
   !this is grid-average updraft velocity, 
   !now use real updraft velocity from convective scheme
!  do k = kts,kte
!    !fluxup: mole or kg/s  vud,vdd: m/s
!    vud(k) = fluxup(k)/dx2/airdens(k)
!    vdd(k) = fluxdown(k)/dx2/airdens(k)
!    vud(k) = 5.0 ! over ocean, or 10.0 over land 
!    vdd(k) = 5.0 
!  enddo

   !-------------------------------------------------------
   !calculate the convective cloud property
   ! cldliq, the cloud liquid water content [cm3 h2o/cm3 air],
   !         is a function of the local kelvin temperature:
   ! cldliq = 2e-6                    [     t >= 268 k    ]
   ! cldliq = 2e-6 * ((t - 248) / 20) [ 248 k < t < 268 k ]
   ! cldliq = 0                       [     t <= 248 k    ]
   !-------------------------------------------------------
   do k = kts,kte
       if ( tmp(k) >= 268.0 ) then
           cldliq(k) = 2e-6
       else if ( tmp(k) > 248.0 .and. tmp(k) < 268.0 ) then
           cldliq(k) = 2e-6 * ( ( tmp(k) - 248.0 ) / 20.0 )
       else
           cldliq(k) = 0.0
       endif
       !-------------------------------------------------------
       ! cldice, the cloud ice content [cm3 ice/cm3 air] is given by:
       ! cldice = 2e-6 - cldliq
       !-------------------------------------------------------
       cldice(k) = 2e-6 - cldliq(k)
   enddo

   !-------------------------------
   ! compute_f begins here!
   !-------------------------------
   if (is_aerosol(n)) then 
       csr=csr_conv(2)  ! now assume aerosol has the same activation rate for different size
                        ! in future, it should pass massptr_aer here to decide the size 
       call f_aerosol( kc,dzq,csr,f,vud,vdd,      &
           ids,ide, jds,jde, kds,kde,       &
           ims,ime, jms,jme, kms,kme,       &
           its,ite, jts,jte, kts,kte        )
   endif

 end subroutine cuscav_kf
!===========================================================================

 subroutine f_aerosol( kc,dzq,csr,f,vud,vdd,         &
                       ids,ide, jds,jde, kds,kde,       &
                       ims,ime, jms,jme, kms,kme,       &
                       its,ite, jts,jte, kts,kte        )
!
!******************************************************************************
!  subroutine f_aerosol returns the fraction of aerosol scavenged in updrafts
!
!  arguments as input:
!  (1 ) kc (real) : conversion rate from cloud condensate to precip [s^-1]
!
!  arguments as output:
!  (2 ) f  (real) : fraction of aerosol scavenged in updrafts [unitless]
!******************************************************************************
!
   implicit none


   ! arguments
   integer,intent(in)  ::    &
                          ids,ide, jds,jde, kds,kde,    &
                          ims,ime, jms,jme, kms,kme,    &
                          its,ite, jts,jte, kts,kte
   real, intent(in)  :: kc
   real, intent(in)  :: csr !aerosol activation rate 
   real, intent(in)  :: dzq(kms:kme)
   real, intent(in)  :: vud(kms:kme)
   real, intent(in)  :: vdd(kms:kme)
   real, intent(inout) :: f(kms:kme,2)

   ! local variables
   integer           :: i, j, k
   real              :: kr

   !=================================================================
   ! f_aerosol begins here!
   ! we set k = kc*csr, and compute f accordingly (cf jacob et al 2000 )    
   !=================================================================
   kr=kc*csr

   ! turn off scavenging in the first level by setting f = 0
   f(kts,1:2) = 0.0

   ! apply scavenging in levels 2 and higher
   do k = kts+1, kte 
     ! (eq. 2, jacob et al, 2000, with k = kc)
     if (vud(k).gt.0) then 
      f(k,1) = 1.0 - exp( -kr * dzq(k) / vud(k) )
     else 
      f(k,1) = 0.0
     endif
     if (vdd(k).gt.0) then 
      f(k,2) = 1.0 - exp( -kr * dzq(k) / vdd(k) )
     else  
      f(k,2) = 0.0
     endif
   enddo

 end subroutine f_aerosol

END MODULE module_ctrans_kfeta


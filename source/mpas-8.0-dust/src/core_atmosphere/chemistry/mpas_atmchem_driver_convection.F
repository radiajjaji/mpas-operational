!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC)
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to compute convection for chem.
! * Author: Jiawang Feng from AEMOL USTC
! Last update: jwfeng 03/2021
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

module mpas_atmchem_driver_convection

   use mpas_kind_types
   use mpas_pool_routines
   use mpas_atmchem_vars
   use mpas_atmchem_utilities
   use mpas_atmphys_vars
   use mpas_atmchem_constants, only: xlv, g, r_v
   use mpas_constants, only : cp
   !chemistry schemes:
   use module_ctrans_kfeta, only: kfeta_cup_transport
   use module_ctrans_general, only: general_cup_transport

   implicit none

   private

   public :: allocate_convection, &
             deallocate_convection, &
             driver_convection

contains

   subroutine allocate_convection

      ! For "cu_kain_fritsch" and "cu_kain_fritsch_sca"
      if(.not.allocated(umf_c        )) allocate(umf_c(ims:ime,kms:kme,jms:jme))
      if(.not.allocated(uer_c        )) allocate(uer_c(ims:ime,kms:kme,jms:jme))
      if(.not.allocated(udr_c        )) allocate(udr_c(ims:ime,kms:kme,jms:jme))
      if(.not.allocated(dmf_c        )) allocate(dmf_c(ims:ime,kms:kme,jms:jme))
      if(.not.allocated(der_c        )) allocate(der_c(ims:ime,kms:kme,jms:jme))
      if(.not.allocated(ddr_c        )) allocate(ddr_c(ims:ime,kms:kme,jms:jme))
      if(.not.allocated( wu_c        )) allocate( wu_c(ims:ime,kms:kme,jms:jme))
      if(.not.allocated(ldb_c        )) allocate(ldb_c(ims:ime,        jms:jme))
      if(.not.allocated(let_c        )) allocate(let_c(ims:ime,        jms:jme))
      if(.not.allocated(ksl_c        )) allocate(ksl_c(ims:ime,        jms:jme))
      if(.not.allocated(lfs_c        )) allocate(lfs_c(ims:ime,        jms:jme))
      if(.not.allocated(klc_c        )) allocate(klc_c(ims:ime,        jms:jme))
      if(.not.allocated(ishall_c     )) allocate(ishall_c(ims:ime,     jms:jme))
      if(.not.allocated(check_kfeta_c)) allocate(check_kfeta_c(ims:ime,jms:jme))
      if(.not.allocated(cutop_c      )) allocate(cutop_c(ims:ime,      jms:jme))
      if(.not.allocated(cubot_c      )) allocate(cubot_c(ims:ime,      jms:jme))
      if(.not.allocated(  nca_c      )) allocate(  nca_c(ims:ime,      jms:jme))

      ! For "cu_grell_freitas"
      if(.not.allocated(k22_shalc_c  )) allocate(k22_shalc_c(ims:ime,  jms:jme))
      if(.not.allocated(kbcon_shalc_c)) allocate(kbcon_shalc_c(ims:ime,jms:jme))
      if(.not.allocated(ktop_shalc_c )) allocate(ktop_shalc_c(ims:ime, jms:jme))
      if(.not.allocated(xmb_shalc_c  )) allocate(xmb_shalc_c(ims:ime,  jms:jme))

      ! For tendency
      if(.not.allocated(tend_ctrans_c  )) allocate(tend_ctrans_c(ims:ime,kms:kme,jms:jme,1:num_chem))
      if(.not.allocated(tend_wetconv_c )) allocate(tend_wetconv_c(ims:ime,jms:jme,1:num_chem))
      ! For diagnostics
      if(.not.allocated(aq_gas_ratio_c )) allocate(aq_gas_ratio_c(ims:ime,kms:kme,jms:jme))
      if(.not.allocated(ktop_deepc_c   )) allocate(ktop_deepc_c(ims:ime,jms:jme))
      if(.not.allocated(kbcon_deepc_c  )) allocate(kbcon_deepc_c(ims:ime,jms:jme))

   end subroutine allocate_convection

   subroutine deallocate_convection

      ! For "cu_kain_fritsch" and "cu_kain_fritsch_sca"
      if(allocated(umf_c        )) deallocate(umf_c         )
      if(allocated(uer_c        )) deallocate(uer_c         )
      if(allocated(udr_c        )) deallocate(udr_c         )
      if(allocated(dmf_c        )) deallocate(dmf_c         )
      if(allocated(der_c        )) deallocate(der_c         )
      if(allocated(ddr_c        )) deallocate(ddr_c         )
      if(allocated( wu_c        )) deallocate( wu_c         )
      if(allocated(ldb_c        )) deallocate(ldb_c         )
      if(allocated(let_c        )) deallocate(let_c         )
      if(allocated(ksl_c        )) deallocate(ksl_c         )
      if(allocated(lfs_c        )) deallocate(lfs_c         )
      if(allocated(klc_c        )) deallocate(klc_c         )
      if(allocated(ishall_c     )) deallocate(ishall_c      )
      if(allocated(check_kfeta_c)) deallocate(check_kfeta_c )
      if(allocated(cutop_c      )) deallocate(cutop_c       )
      if(allocated(cubot_c      )) deallocate(cubot_c       )
      if(allocated(  nca_c      )) deallocate(  nca_c       )

      ! For "cu_grell_freitas"
      if(allocated(k22_shalc_c  )) deallocate(k22_shalc_c  )
      if(allocated(kbcon_shalc_c)) deallocate(kbcon_shalc_c)
      if(allocated(ktop_shalc_c )) deallocate(ktop_shalc_c )
      if(allocated(xmb_shalc_c  )) deallocate(xmb_shalc_c  )

      ! For tendency
      if(allocated(tend_ctrans_c  )) deallocate(tend_ctrans_c  )
      if(allocated(tend_wetconv_c )) deallocate(tend_wetconv_c )
      ! For diagnostics
      if(allocated(aq_gas_ratio_c )) deallocate(aq_gas_ratio_c )
      if(allocated(ktop_deepc_c   )) deallocate(ktop_deepc_c   )
      if(allocated(kbcon_deepc_c  )) deallocate(kbcon_deepc_c  )

   end subroutine deallocate_convection

   subroutine convection_from_MPAS(convection_scheme, diag_physics, its, ite)

      type(mpas_pool_type), intent(in) :: diag_physics
      integer, intent(in) :: its, ite
      character(len = StrKIND) :: convection_scheme

      real(kind = RKIND), dimension(:, :), pointer :: umf, uer, udr
      real(kind = RKIND), dimension(:, :), pointer :: dmf, der, ddr, wu
      real(kind = RKIND), dimension(:), pointer :: cutop, cubot, nca
      integer, dimension(:), pointer :: k22_shalc, kbcon_shalc, ktop_shalc
      real(kind = RKIND), dimension(:), pointer :: xmb_shalc
      integer, dimension(:), pointer :: ldb, let, ksl, lfs, klc
      integer, dimension(:), pointer :: ishall, check_kfeta
      !local variables:
      integer :: i, k, j

      convection_select: select case( trim(convection_scheme) )
      case ("cu_kain_fritsch", "cu_kain_fritsch_sca")
         call mpas_pool_get_array(diag_physics, 'umf'        , umf         )
         call mpas_pool_get_array(diag_physics, 'uer'        , uer         )
         call mpas_pool_get_array(diag_physics, 'udr'        , udr         )
         call mpas_pool_get_array(diag_physics, 'dmf'        , dmf         )
         call mpas_pool_get_array(diag_physics, 'der'        , der         )
         call mpas_pool_get_array(diag_physics, 'ddr'        , ddr         )
         call mpas_pool_get_array(diag_physics, 'wu'         , wu          )
         call mpas_pool_get_array(diag_physics, 'ldb'        , ldb         )
         call mpas_pool_get_array(diag_physics, 'let'        , let         )
         call mpas_pool_get_array(diag_physics, 'ksl'        , ksl         )
         call mpas_pool_get_array(diag_physics, 'lfs'        , lfs         )
         call mpas_pool_get_array(diag_physics, 'klc'        , klc         )
         call mpas_pool_get_array(diag_physics, 'cutop'      , cutop       )
         call mpas_pool_get_array(diag_physics, 'cubot'      , cubot       )
         call mpas_pool_get_array(diag_physics, 'nca'        , nca         )
         call mpas_pool_get_array(diag_physics, 'ishall'     , ishall      )
         call mpas_pool_get_array(diag_physics, 'check_kfeta', check_kfeta )

         do j = jts,jte
            do i = its,ite
               ktop_deepc_c(i,j)  = 0
               kbcon_deepc_c(i,j) = 0
               ldb_c(i,j)   = ldb(i)
               let_c(i,j)   = let(i)
               ksl_c(i,j)   = ksl(i)
               lfs_c(i,j)   = lfs(i)
               klc_c(i,j)   = klc(i)
               nca_c(i,j)   = nca(i)
               cutop_c(i,j) = cutop(i)
               cubot_c(i,j) = cubot(i)
               ishall_c(i,j)= ishall(i)
               check_kfeta_c(i,j) = check_kfeta(i)
               do k = kts,kte
                  aq_gas_ratio_c(i,k,j) = 0.0
                  umf_c(i,k,j) = umf(k,i)
                  uer_c(i,k,j) = uer(k,i)
                  udr_c(i,k,j) = udr(k,i)
                  dmf_c(i,k,j) = dmf(k,i)
                  der_c(i,k,j) = der(k,i)
                  ddr_c(i,k,j) = ddr(k,i)
                  wu_c(i,k,j)  = wu(k,i)
               enddo
            enddo
         enddo
      case ("cu_grell_freitas")
         call mpas_pool_get_array(diag_physics, 'k22_shalc',   k22_shalc  )
         call mpas_pool_get_array(diag_physics, 'kbcon_shalc', kbcon_shalc)
         call mpas_pool_get_array(diag_physics, 'ktop_shalc',  ktop_shalc )
         call mpas_pool_get_array(diag_physics, 'xmb_shalc',   xmb_shalc  )
         do j = jts, jte
            do i = its,ite
               ktop_deepc_c(i,j)  = 0
               kbcon_deepc_c(i,j) = 0

               k22_shalc_c(i,j)   = k22_shalc(i)
               kbcon_shalc_c(i,j) = kbcon_shalc(i)
               ktop_shalc_c(i,j)  = ktop_shalc(i)
               xmb_shalc_c(i,j)   = xmb_shalc(i)
            enddo
         enddo
      case default
         call chemistry_message('There are no corresponding convection ' &
                                  //'scheme for: ['//trim(convection_scheme)//'] But going on' )
!         call chemistry_error_fatal('There are no corresponding convection ' &
!                                  //'scheme for: ['//trim(convection_scheme)//']' )
      end select convection_select

   end subroutine convection_from_MPAS

   subroutine convection_to_MPAS(state, time_lev, configs, diag_physics, diag_chemistry, tend_chemistry, its, ite)

      !in arguments:
      type(mpas_pool_type), intent(in   ) :: configs
      !inout arguments:
      type(mpas_pool_type), intent(inout) :: state
      type(mpas_pool_type), intent(inout) :: tend_chemistry
      type(mpas_pool_type), intent(inout) :: diag_physics
      type(mpas_pool_type), intent(inout) :: diag_chemistry
      integer, intent(in) :: its, ite, time_lev

      !local variables:
      integer :: i, k, j, h

      !local pointers:
      real(kind = RKIND), dimension(:, :, :), pointer :: chem
      real(kind = RKIND), dimension(:, :, :), pointer :: tend_ctrans
      real(kind = RKIND), dimension(:, :   ), pointer :: tend_wetconv
      real(kind = RKIND), dimension(:, :), pointer :: aq_gas_ratio
      integer, dimension(:), pointer :: ktop_deepc, kbcon_deepc
      character(len = StrKIND), pointer :: ctrans_opt

      call mpas_pool_get_array(state, 'chem', chem, time_lev)
      call mpas_pool_get_array(tend_chemistry, 'tend_ctrans', tend_ctrans)
      call mpas_pool_get_array(tend_chemistry, 'tend_wetconv',tend_wetconv)
      call mpas_pool_get_array(diag_chemistry, 'aq_gas_ratio', aq_gas_ratio)
      call mpas_pool_get_array(diag_physics, 'kbcon_deepc', kbcon_deepc)
      call mpas_pool_get_array(diag_physics, 'ktop_deepc', ktop_deepc)
      call mpas_pool_get_config(configs, 'config_ctrans_opt', ctrans_opt)

      convection_select: select case(trim(ctrans_opt))
      case ("kfeta")
         do h = 1, num_chem
         do j = jts, jte
         do k = kts, kte
         do i = its, ite
            if (.not. chem_c(i, k, j, h) .gt. 0.0_RKIND) chem_c(i, k, j, h) = 0.0_RKIND
            tend_ctrans(h, k, i) = (chem_c(i, k, j, h) - chem(h, k, i)) / dt_dyn
            chem(h, k, i) = chem_c(i, k, j, h)
         enddo
         enddo
         enddo
         enddo
      case ("general")
         do h = 1, num_chem
         do j = jts, jte
         do k = kts, kte
         do i = its, ite
            if (.not. chem_c(i, k, j, h) .gt. 0.0_RKIND) chem_c(i, k, j, h) = 0.0_RKIND
            tend_ctrans(h, k, i) = tend_ctrans_c(i, k, j, h)
            chem(h, k, i) = chem_c(i, k, j, h)
         enddo
         enddo
         enddo
         enddo
         do j = jts, jte
         do k = kts, kte
         do i = its, ite
            aq_gas_ratio(k, i) = aq_gas_ratio_c(i, k, j)
         enddo
         enddo
         enddo
         do j = jts, jte
         do i = its, ite
            ktop_deepc(i)  = ktop_deepc_c(i, j)
            kbcon_deepc(i) = kbcon_deepc_c(i, j)
         enddo
         enddo
         do h = 1, num_chem
         do j = jts, jte
         do i = its, ite
            tend_wetconv(h, i) = tend_wetconv_c(i, j, h)
         enddo
         enddo
         enddo

      case default
         call chemistry_error_fatal('Wrong ctrans_opt'//trim(ctrans_opt)//']' )
      end select convection_select


   end subroutine convection_to_MPAS

   subroutine driver_convection(state, diag_physics, time_lev, configs, diag_chemistry, tend_chemistry, its, ite)

      type(mpas_pool_type), intent(in) :: configs
      type(mpas_pool_type), intent(inout) :: state
      type(mpas_pool_type), intent(inout) :: tend_chemistry
      type(mpas_pool_type), intent(inout) :: diag_physics
      type(mpas_pool_type), intent(inout) :: diag_chemistry
      integer, intent(in) :: its, ite, time_lev

      !local vars
      logical, pointer :: debug_on, scav_on
      character(len = StrKIND), pointer :: convection_scheme, chem_opt, ctrans_opt

      call mpas_pool_get_config(configs, 'config_print_debug_info', debug_on)
      call mpas_pool_get_config(configs, 'config_chem_opt', chem_opt)
      call mpas_pool_get_config(configs, 'config_cu_scav_on', scav_on)
      call mpas_pool_get_config(configs, 'config_convection_scheme', convection_scheme)
      call mpas_pool_get_config(configs, 'config_ctrans_opt', ctrans_opt)
      !copy MPAS arrays to local arrays:
      call convection_from_MPAS(convection_scheme, diag_physics, its, ite)
      if (debug_on) call chemistry_message('convection from MPAS success')

      convection_select: select case(trim(ctrans_opt))
      case ("kfeta")
         if (debug_on) CALL chemistry_message('convection_DRIVER CALLING kfeta')
         call kfeta_cup_transport(                                     &
              dt=dt_dyn, dx=dx_c, rho_phy=rho_phy_c, alt=alt_c       , &
              chem=chem_c, u=u_phy_c, v=v_phy_c, t_phy=t_phy_c       , &
              dz8w=dz_c, p_phy=pres_hyd_c, z=z_c                     , &
              umf3d=umf_c, uer3d=uer_c, udr3d=udr_c, dmf3d=dmf_c     , &
              der3d=der_c, ddr3d=ddr_c,  wu3d= wu_c, ldb2d=ldb_c     , &
              let2d=let_c, ksl2d=ksl_c, lfs2d=lfs_c, klc2d=klc_c     , &
              cutop=int(cutop_c), cubot=int(cubot_c),nca2d=nca_c     , &
              ishall2d=ishall_c, check2d=check_kfeta_c               , &
              num_chem=num_chem, chem_opt=chem_opt, scav_opt=scav_on , &
              ids=ids, ide=ide, kds=kds, kde=kde, jds=jds, jde=jde   , &
              ims=ims, ime=ime, kms=kms, kme=kme, jms=jms, jme=jme   , &
              its=its, ite=ite, kts=kts, kte=kte, jts=jts, jte=jte     )
      case ("general")
         !jwfeng: *_shallow vars only present in cu_grell_freitas scheme
         !        if you want to use general ctrans module, you must set
         !        ishallow to 0 to avoid using *_shallow vars.
         if (trim(convection_scheme) .ne. "cu_grell_freitas") then
            call general_cup_transport(                                     &
                 dt=dt_dyn, rho_phy=rho_phy_c, raincv=raincv_c            , &
                 chem=chem_c, wetdep_2d = tend_wetconv_c                  , &
                 u=u_phy_c, v=v_phy_c, t_phy=t_phy_c                      , &
                 qv=qv_c, dz8w=dz_c, p_phy=pres_hyd_c, xlv=xlv            , &
                 cp=cp, g=g, r_v=r_v, z=z_c                               , &
                 tracert_3d = tend_ctrans_c, aq_gas_ratio = aq_gas_ratio_c, &
                 ktop_2d = ktop_deepc_c, kbcon_2d = kbcon_deepc_c         , &
                 ishallow=0, num_gas=num_gas, num_chem=num_chem           , &
                 chem_opt=chem_opt, cu_scav_on=scav_on, conv_tr_aqchem=0  , &
                 ids=ids, ide=ide, kds=kds, kde=kde, jds=jds, jde=jde     , &
                 ims=ims, ime=ime, kms=kms, kme=kme, jms=jms, jme=jme     , &
                 its=its, ite=ite, kts=kts, kte=kte, jts=jts, jte=jte       )
         else
            ! jwfeng: GF can calculate shallow convection through setting ishallow=1.
            !         But we don't do shallow for now to compared to KF
            call general_cup_transport(                                     &
                 dt=dt_dyn, rho_phy=rho_phy_c, raincv=raincv_c            , &
                 chem=chem_c, wetdep_2d = tend_wetconv_c                  , &
                 u=u_phy_c, v=v_phy_c, t_phy=t_phy_c                      , &
                 qv=qv_c, dz8w=dz_c, p_phy=pres_hyd_c, xlv=xlv            , &
                 cp=cp, g=g, r_v=r_v, z=z_c                               , &
                 tracert_3d = tend_ctrans_c, aq_gas_ratio = aq_gas_ratio_c, &
                 ktop_2d = ktop_deepc_c, kbcon_2d = kbcon_deepc_c         , &
                 k22_shallow=k22_shalc_c, kbcon_shallow=kbcon_shalc_c     , &
                 ktop_shallow=ktop_shalc_c, xmb_shallow=xmb_shalc_c       , &
                 ishallow=0, num_gas=num_gas, num_chem=num_chem           , &
                 chem_opt=chem_opt, cu_scav_on=scav_on, conv_tr_aqchem=0  , &
                 ids=ids, ide=ide, kds=kds, kde=kde, jds=jds, jde=jde     , &
                 ims=ims, ime=ime, kms=kms, kme=kme, jms=jms, jme=jme     , &
                 its=its, ite=ite, kts=kts, kte=kte, jts=jts, jte=jte       )
         end if
      case default
         call chemistry_error_fatal('Wrong ctrans_opt'//trim(ctrans_opt)//']' )
      end select convection_select

      if (debug_on) call chemistry_message('convection scheme success')
      !copy local arrays to MPAS grid:
      call convection_to_MPAS(state, time_lev, configs, diag_physics, diag_chemistry, tend_chemistry, its, ite)
      if (debug_on) call chemistry_message('convection to MPAS success')

   end subroutine driver_convection

end module mpas_atmchem_driver_convection

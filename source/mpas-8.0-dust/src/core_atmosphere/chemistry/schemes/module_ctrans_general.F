!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC)
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to Compute Convective Transport following KFeta scheme
! * Author: Jiawang Feng and Chun Zhao
! Last update: jwfeng 11/2021
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

module module_ctrans_general

   use mpas_atmchem_utilities, only: chemistry_message, chemistry_error_fatal
   use module_cu_gd, only: cup_dd_edt, cup_dd_he, cup_dd_moisture      , &
                           cup_dd_nms, cup_env, cup_env_clev, cup_kbcon, &
                           cup_ktop, cup_maximi, cup_minimi, cup_up_he , &
                           cup_up_moisture, cup_up_nms
   use mpas_atmchem_vars, only: p_co, p_so2, p_hno3, p_hno4, p_dms     , &
                                p_n2o5, p_nh3, p_h2o2, p_o3, p_ora1    , &
                                p_op1, p_paa, p_sulf, p_so4aj, p_nh4aj , &
                                p_no3aj, p_bc1, p_bc2, p_oc1, p_oc2    , &
                                p_seas_1, p_seas_2, p_seas_3, p_seas_4
   implicit none

   private

   real, parameter :: qcldwtr_cutoff = 1.0e-6 ! kg / kg
   real, parameter :: mwdry = 28.966  ! molecular mass of dry air (g/mol)
   real, parameter :: mwso4 = 96.00   ! molecular mass of so4-- (g/mol)
   real, parameter :: mwno3 = 62.0    ! molecular mass of no3- (g/mol)
   real, parameter :: mwnh4 = 18.0985 ! molecular mass of nh4+ (g/mol)

   public :: general_cup_transport

contains

   subroutine general_cup_transport(                  &
              dt, rho_phy, raincv, chem, wetdep_2d  , &
              u, v, t_phy, qv, dz8w, p_phy          , &
              xlv, cp, g, r_v, z                    , &
              tracert_3d, aq_gas_ratio              , &
              ktop_2d, kbcon_2d                     , &
              k22_shallow, kbcon_shallow            , &
              ktop_shallow, xmb_shallow, ishallow   , &
              num_gas, num_chem, chem_opt           , &
              cu_scav_on, conv_tr_aqchem            , &
              ids, ide, kds, kde, jds, jde          , &
              ims, ime, kms, kme, jms, jme          , &
              its, ite, kts, kte, jts, jte            )

      integer, intent(in) :: num_gas, ishallow, num_chem
      integer, intent(in) :: ids, ide, jds, jde, kds, kde
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: its, ite, jts, jte, kts, kte
      integer, intent(in) :: conv_tr_aqchem
      logical, intent(in) :: cu_scav_on
      character(len=*), intent(in) :: chem_opt

      real, intent(in) :: dt
      real, intent(in) :: xlv, r_v
      real, intent(in) :: cp, g
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: qv
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: u
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: v
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: t_phy
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: z
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: p_phy
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: dz8w
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: rho_phy
      real, dimension(ims:ime, jms:jme), intent(in) :: raincv
      real, dimension(ims:ime, kms:kme, jms:jme, num_chem), intent(inout) :: chem
      real, dimension(ims:ime, jms:jme, num_chem), intent(out) :: wetdep_2d
      real, dimension(ims:ime, jms:jme), intent(in), optional :: xmb_shallow

      integer, dimension(ims:ime, jms:jme), intent(in), optional :: k22_shallow
      integer, dimension(ims:ime, jms:jme), intent(in), optional :: kbcon_shallow
      integer, dimension(ims:ime, jms:jme), intent(in), optional :: ktop_shallow
      integer, dimension(ims:ime, jms:jme), intent(out), optional :: kbcon_2d
      integer, dimension(ims:ime, jms:jme), intent(out), optional :: ktop_2d
      real, dimension(ims:ime, kms:kme, jms:jme, num_chem), intent(out) :: tracert_3d
      real, dimension(ims:ime, kms:kme, jms:jme), intent (inout) :: aq_gas_ratio ! Output to diagnostics
      ! local vars
      real, dimension(its:ite, kts:kte) :: outt, outq, outqc
      real, dimension(its:ite) :: pret, ter11
      real, dimension(its:ite, num_chem) :: wetdep_1d
      ! basic environmental input includes moisture convergence (mconv)
      ! omega (omeg), windspeed (us, vs), and a flag (aaeq) to turn off
      ! convection for this call only and at that particular gridpoint
      real, dimension(its:ite, kts:kte) :: T, TN, q, qo, PO, P, US, VS, hstary
      real, dimension(its:ite, kts:kte, num_chem) :: tracer, tracert, tracert3
      real, dimension(its:ite) :: Z1, PSUR, AAEQ, xmb3
      integer, dimension(its:ite) :: ktop, kbcon, k23, kbcon3, ktop3

      integer :: nv, i, j, k, icldck, ipr, jpr, npr
      integer :: itf, jtf, ktf, iopt
      real :: tcrit, dp, dq, epsilc
      character(len=300) :: message
      real :: tracert_ratio


      ! jwfeng: *_shallow vars should present if ishallow equals 1
      if (ishallow .eq. 1) then
         if (.not. (present(xmb_shallow)   .and. &
                    present(k22_shallow)   .and. &
                    present(kbcon_shallow) .and. &
                    present(ktop_shallow)) ) then
            call chemistry_error_fatal("*_shallow vars should present if ishallow equals 1")
         end if
      end if

      epsilc = 1.e-30
      wetdep_1d(:, :) = 0.0
      wetdep_2d(:,:,:) = 0.0
      ktop_2d(:,:) = 0.0
      kbcon_2d(:,:) = 0.0
      tracert_3d(:,:,:,:) = 0.0
      ipr = 0
      jpr = 0
      npr = 1
      if (p_co .gt. 0) npr = p_co
      tcrit = 258.
      iopt = 0
      itf = min(ite, ide - 1)
      ktf = min(kte, kde - 1)
      jtf = min(jte, jde - 1)

      do j = jts, jtf
         do i = its, itf
            ktop(i) = 0
            kbcon(i)= 0
            psur(i) = p_phy(i, kts, j) * .01
            ter11(i) = z(i, kts, j)
            aaeq(i) = 0.
            ! rainrate is input for this transport/wet-deposition routine
            pret(i) = raincv(i, j) / dt
            if (pret(i) .le. 0.) aaeq(i) = 20.
         enddo
         if (ishallow .eq. 1) then
            do i = its, itf
               xmb3(i) = xmb_shallow(i, j)
               ktop3(i) = ktop_shallow(i, j)
               k23(i) = k22_shallow(i, j)
               kbcon3(i) = kbcon_shallow(i, j)
            enddo
         else
            do i = its, itf
               xmb3(i) = 0.
               ktop3(i) = 0
               k23(i) = 0
               kbcon3(i) = 0
            enddo
         endif
         do k = kts, ktf
            do i = its, itf
               po(i, k) = p_phy(i, k, j) * .01
               p(i, k) = po(i, k)
               us(i, k) = u(i, k, j)
               vs(i, k) = v(i, k, j)
               t(i, k) = t_phy(i, k, j)
               q(i, k) = qv(i, k, j)
               if (q(i, k) .lt. 1.e-08) q(i, k) = 1.e-08
            enddo
         enddo
         do nv=1,num_chem
            do K=kts,ktf
               do I=ITS,ITF
                  tracer(i,k,nv)=max(epsilc,chem(i,k,j,nv))
                  tracert(i,k,nv)=0.
                  tracert3(i,k,nv)=0.
               enddo
            enddo
         enddo
         do k=kts,ktf
            do i=its,itf
               ! cu_co_ten(i,k,j)=0.
               ! hstary(i,k)=hstar4(nv)*exp(dhr(nv)*(1./t(i,k)-1./298.))
               if(i.eq.ipr.and.j.eq.jpr)then
                  write(message, "(I3,x,4F9.3)") k, pret(i), tracer(i, k, npr), p(i, k), z(i, k, j)
                  call chemistry_message(trim(message))
               endif
            enddo
         enddo

         !---- CALL NON_RESOLVED CONVECTIVE TRANSPORT
         call cup_ct(ktop, kbcon, k23, kbcon3, ktop3, xmb3, &
              tracer, j, aaeq, t, q, ter11, pret, p, tracert, tracert3, &
              hstary, dt, psur, us, vs, tcrit, wetdep_1d, &
              xlv, r_v, cp, g, ipr, jpr, npr, num_chem, chem_opt, &
              cu_scav_on, conv_tr_aqchem, chem, aq_gas_ratio, &
              ishallow, num_gas, ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, &
              its, ite, jts, jte, kts, kte)

         do nv = 1, num_chem
            do i = its, itf
               if (pret(i) .le. 0.) then
                  do k = kts, ktf
                     tracert(i, k, nv) = 0.
                  enddo
               endif
            enddo
         enddo
         if (ishallow .eq. 1) then
            call neg_check_ct('shallow', pret, ktop3, epsilc, dt, &
                 tracer, tracert3, iopt, num_chem, &
                 its, ite, kts, kte, itf, ktf, ipr, jpr, npr, j)
            do nv = 1, num_chem
               do i = its, itf
                  do k = kts, ktf
                     chem(i, k, j, nv) = max(epsilc, chem(i, k, j, nv) &
                                       + tracert3(i, k, nv) * dt)
                     tracert_3d(i, k, j, nv) = tracert_3d(i, k, j, nv) + tracert3(i, k, nv)
                  enddo
               enddo
            enddo
         endif

         ! now deep convection
         call neg_check_ct('deep', pret, ktop, epsilc, dt, &
              tracer, tracert, iopt, num_chem, &
              its, ite, kts, kte, itf, ktf, ipr, jpr, npr, j)
         do nv = 1, num_chem
            do i = its, itf
               if (pret(i) .gt. 0.) then
                  do k = kts, ktf
                     chem(i, k, j, nv) = max(epsilc, chem(i, k, j, nv) &
                                       + tracert(i, k, nv) * dt)
                     ! if (nv .eq. npr) then
                     !    cu_co_ten(i, k, j) = tracert(i, k, npr) * dt
                     ! endif
                     tracert_3d(i, k, j, nv) = tracert_3d(i, k, j, nv) + tracert(i, k, nv)
                  enddo
               else
                  do k = kts, ktf
                     tracert(i, k, nv) = 0.
                     ! if (nv .eq. npr)cu_co_ten(i, k, j) = 0.
                  enddo
               endif
            enddo
         enddo

         ! jwfeng: Check if tracert is conserved vertically when cu_scav_on==.false.
         if (.not. cu_scav_on) then
            do nv = 1, num_chem
               do i = its, itf
                  tracert_ratio = dt*sum(tracert_3d(i,:,j,nv)*dz8w(i,:,j)*rho_phy(i,:,j))/ &
                                     sum(      chem(i,:,j,nv)*dz8w(i,:,j)*rho_phy(i,:,j))*100.0
                  if (abs(tracert_ratio) .ge. 1.0 .and. sum(chem(i,:,j,nv)) .ge. 1.0e-10) then
                     write(message, "(a,I5,x,a,I3,x,a,F9.4,a)") "iCell=", i, "iChem=", nv, "diff_ratio=", tracert_ratio, "%"
                     call chemistry_message("Tendencies of ctrans are not conserved at:")
                     call chemistry_message(trim(message))
                  end if
               enddo
            enddo
         endif
         wetdep_2d(its:ite, j, :) = wetdep_1d(its:ite, :)
         ktop_2d(its:ite, j) = ktop(its:ite)
         kbcon_2d(its:ite, j)= kbcon(its:ite)

      enddo ! end of loop: j = jts, jtf

   end subroutine general_cup_transport


   subroutine cup_ct(ktop,kbcon,k23,kbcon3,ktop3,xmb3,tracer,j,aaeq,t,q,z1,  &
              pre,p,tracert,tracert3,hstary,dtime,psur,us,vs,tcrit,    &
              wetdep,xl,rv,cp,g,ipr,jpr,npr,num_chem,chem_opt,         &
              cu_scav_on,conv_tr_aqchem, chem, aq_gas_ratio,           &
              ishallow,num_gas,ids,ide, jds,jde, kds,kde,              &
              ims,ime, jms,jme, kms,kme,                               &
              its,ite, jts,jte, kts,kte                                )

      integer, intent(in) :: num_chem, ids, ide, jds, jde, kds, kde
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: conv_tr_aqchem, ishallow
      integer, intent(in) :: its, ite, jts, jte, kts, kte
      integer, intent(in) :: ipr, jpr, npr, num_gas
      integer, intent(in) :: j
      logical, intent(in) :: cu_scav_on
      character(len=*), intent(in) :: chem_opt
      ! basic environmental input includes moisture convergence (mconv)
      ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
      ! convection for this call only and at that particular gridpoint
      real, intent (in) :: dtime, tcrit, xl, cp, rv, g
      real, dimension(its:ite, kts:kte), intent(in) :: T, P, US, VS, HSTARY
      real, dimension (its:ite), intent (in) :: Z1, PSUR, AAEQ, xmb3

      !tracert = output temp tendency (per s)
      !pre = input precip
      real, dimension (its:ite, kts:kte, num_chem), intent(inout) :: tracert
      real, dimension (its:ite, kts:kte, num_chem), intent(inout) :: tracer
      real, dimension (its:ite, kts:kte, num_chem), intent(inout) :: tracert3
      real, dimension (ims:ime, kms:kme, jms:jme, num_chem), intent(inout) :: chem
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: aq_gas_ratio ! Output to diagnostics
      real, dimension (its:ite, kts:kte), intent(inout) :: Q
      real, dimension (its:ite), intent(inout) :: pre
      integer, dimension (its:ite), intent(inout) :: ktop, kbcon, k23, kbcon3, ktop3

      ! Local Varibles
      real, dimension (its:ite, 1:3) :: edtc
      !************* the following are your basic environmental
      !              variables. They carry a "_cup" if they are
      !              on model cloud levels (staggered). They carry
      !              an "o"-ending (z becomes zo), if they are the forced
      !              variables. They are preceded by x (z becomes xz)
      !              to indicate modification by some typ of cloud
      !
      ! z          = heights of model levels
      ! q          = environmental mixing ratio
      ! qes        = environmental saturation mixing ratio
      ! t          = environmental temp
      ! p          = environmental pressure
      ! he         = environmental moist static energy
      ! hes        = environmental saturation moist static energy
      ! z_cup      = heights of model cloud levels
      ! q_cup      = environmental q on model cloud levels
      ! qes_cup    = saturation q on model cloud levels
      ! t_cup      = temperature (Kelvin) on model cloud levels
      ! p_cup      = environmental pressure
      ! he_cup     = moist static energy on model cloud levels
      ! hes_cup    = saturation moist static energy on model cloud levels
      ! gamma_cup  = gamma on model cloud levels
      ! hcd        = moist static energy in downdraft
      ! zd         = normalized downdraft mass flux
      ! dby        = buoancy term
      ! entr       = entrainment rate
      ! zd         = downdraft normalized mass flux
      ! entr       = entrainment rate
      ! hcd        = h in model cloud
      ! bu         = buoancy term
      ! zd         = normalized downdraft mass flux
      ! gamma_cup  = gamma on model cloud levels
      ! mentr_rate = entrainment rate
      ! qcd        = cloud q (including liquid water) after entrainment
      ! qrch       = saturation q in cloud
      ! pwd        = evaporate at that level
      ! pwev       = total normalized integrated evaoprate (I2)
      ! entr       = entrainment rate
      ! z1         = terrain elevation
      ! entr       = downdraft entrainment rate
      ! jmin       = downdraft originating level
      ! kdet       = level above ground where downdraft start detraining
      ! psur       = surface pressure
      ! z1         = terrain elevation
      ! zd         = downdraft normalized mass flux
      ! zu         = updraft normalized mass flux
      ! mbdt       = arbitrary numerical parameter
      ! dtime      = dt over which forcing is applied
      ! kbcon      = LFC of parcel from k22
      ! k22        = updraft originating level
      ! dby        = buoancy term
      ! ktop       = cloud top (output)
      ! xmb        = total base mass flux
      ! hc         = cloud moist static energy
      ! hkb        = moist static energy at originating level
      ! mentr_rate = entrainment rate

      real, dimension(its:ite, kts:kte) :: he, hes, qes, z, pwdper
      real, dimension(its:ite, kts:kte) :: qes_cup, q_cup, he_cup, hes_cup
      real, dimension(its:ite, kts:kte) :: z_cup, p_cup, gamma_cup, t_cup
      real, dimension(its:ite, kts:kte) :: dby, qc, qrcd, pwd, pw, hcd, qcd
      real, dimension(its:ite, kts:kte) :: dbyd, hc, qrc, zu, zd, clw_all, zd3
      real, dimension(its:ite, kts:kte) :: clw_all3, cd3, pw3, zu3
      ! cd  = detrainment function for updraft
      ! cdd = detrainment function for downdraft
      real, dimension(its:ite, kts:kte) :: cd, cdd, cdd3, scr1
      real, dimension(its:ite, kts:kte) :: DELLAH, DELLAQ, DELLAT, DELLAQC
      ! edt = epsilon
      real, dimension(its:ite) :: edt, HKB, QKB, edt3
      real, dimension(its:ite) :: XMB, PWAV, PWEV, BU
      real, dimension(its:ite) :: cap_max, cap_max_increment
      real, dimension(its:ite, kts:kte, num_chem) :: tr3_c, tr3_up, tr3_pw
      real, dimension(its:ite, num_chem) :: trkb3
      real, dimension(its:ite, kts:kte, num_chem) :: tr_c, tr_up, tr_dd
      real, dimension(its:ite, kts:kte, num_chem) :: tr_dd3, tre_cup
      real, dimension(its:ite, kts:kte, num_chem) :: tr_pw, tr_pwd
      real, dimension (its:ite, num_chem) :: trkb, wetdep
      integer, dimension(its:ite) :: kzdown, KDET, K22, KB
      integer, dimension(its:ite) :: JMIN, kstabi, kstabm
      integer, dimension(its:ite) :: ierr, KBMAX, ierr5
      integer :: ki, I, K, KK
      integer :: itf, jtf, ktf
      real :: day, dz, mbdt, entr_rate, radius
      real :: entrd_rate, mentr_rate, mentrd_rate
      real :: zcutdown, edtmax, edtmin, depth_min
      real :: zkbmax, z_detr, zktop
      real :: dh, cap_maxs, entr_rate3, mentr_rate3
      character(len=3000) message

      itf = MIN(ite, ide - 1)
      ktf = MIN(kte, kde - 1)
      jtf = MIN(jte, jde - 1)

      !sms$distribute end
      day = 86400.
      !--- specify entrainmentrate and detrainmentrate
      radius = 12000.
      !--- gross entrainment rate (these may be changed later on in the
      !--- program, depending what your detrainment is!!)
      entr_rate = .2 / radius
      entr_rate3 = .2 / 200.
      !--- entrainment of mass
      mentrd_rate = 0.
      mentr_rate = entr_rate
      mentr_rate3 = entr_rate3
      !--- initial detrainmentrates
      do k = kts, ktf
         do i = its, itf
            cd(i, k) = 0.1 * entr_rate
            cd3(i, k) = entr_rate3
            cdd(i, k) = 0.
            cdd3(i, k) = 0.
            clw_all(i, k) = 0.
         enddo
      enddo
      do ki = 1, num_chem
         do k = kts, ktf
            do i = its, itf
               tr_dd3(i, k, ki) = 0.
            enddo
         enddo
      enddo
      !--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
      !    base mass flux
      edtmax = .8
      edtmin = .2
      !--- minimum depth (m), clouds must have
      depth_min = 500.
      !--- maximum depth (mb) of capping
      !--- inversion (larger cap = no convection)
      cap_maxs = 175.
!sms$to_local(grid_dh: <1, mix :size>, <2, mjx :size>) begin
      do i = its, itf
         kbmax(i) = 1
         cap_max_increment(i) = 0.
         edt(i) = 0.
         edt3(i) = 0.
         kstabm(i) = ktf - 1
         IERR(i) = 0
         IERR5(i) = 0
         if (ktop3(i) .lt. 2)ierr5(i) = 20
         if (xmb3(i) .eq. 0.)ierr5(i) = 21
      enddo
      do i = its, itf
         cap_max(i) = cap_maxs
         do k = kts, kte
            zd3(i, k) = 0.
         enddo
      enddo
      !--- max height(m) above ground where updraft air can originate
      zkbmax = 4000.
      !--- height(m) above which no downdrafts are allowed to originate
      zcutdown = 3000.
      !--- depth(m) over which downdraft detrains all its mass
      z_detr = 1250.
      mbdt=dtime*4.E-03
      !--- calculate moist static energy, heights, qes
      call cup_env(z,qes,he,hes,t,q,p,z1, &
           psur,ierr,tcrit,0,xl,cp,   &
           itf,     jtf,     ktf,     &
           its,ite, jts,jte, kts,kte)
      !--- environmental values on cloud levels
      call cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup, &
           hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
           ierr,z1,xl,rv,cp,          &
           itf, jtf, ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_env_clev_tr(tracer,tre_cup,num_chem,ierr, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      do i = its, itf
         if (ierr(i) .eq. 0) then
            do k = kts, ktf - 2
               if (z_cup(i, k) .gt. zkbmax + z1(i)) then
                  kbmax(i) = k
                  exit ! jump out loop k
               endif
            enddo
            !--- level where detrainment for downdraft starts
            do k = kts, ktf
               if (z_cup(i, k) .gt. z_detr + z1(i)) then
                  kdet(i) = k
                  exit ! jump out loop k
               endif
            enddo
         endif
      enddo
      !--- determine level with highest moist static energy content - k22
      call cup_maximi(he_cup,3,kbmax,k22,ierr, &
           itf, jtf, ktf, &
           its,ite, jts,jte, kts,kte)
      do i = its, itf
         if (ierr(I) .eq. 0.) then
            if (K22(I) .ge. KBMAX(i)) ierr(i) = 2
         endif
      enddo
      ! done with basic stuff that is needed everywhere, from here on do not do
      ! deep convection where aaeq .ne.0
      do i = its, itf
         if (aaeq(i) .ne. 0.) then
            ierr(i) = 20
         endif
      enddo
      !--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
      call cup_kbcon(cap_max_increment,1,k22,kbcon,he_cup,hes_cup, &
           ierr,kbmax,p_cup,cap_max, &
           itf, jtf, ktf, &
           its,ite, jts,jte, kts,kte)
      !--- increase detrainment in stable layers
      CALL cup_minimi(HEs_cup,Kbcon,kstabm,kstabi,ierr,  &
           itf, jtf, ktf, &
           its,ite, jts,jte, kts,kte)
      do i = its, itf
         if (ierr(i) .eq. 0.) then
            if (kstabm(i) - 1 .gt. kstabi(i)) then
               do k = kstabi(i), kstabm(i) - 1
                  cd(i, k) = cd(i, k - 1) + 1.5 * entr_rate
                  if (cd(i, k) .gt. 10.0 * entr_rate) then
                     cd(i, k) = 10.0 * entr_rate
                  endif
               enddo
            endif
         endif
      enddo
      !--- calculate incloud moist static energy
      call cup_up_he(k22,hkb,z_cup,cd,mentr_rate,he_cup,hc, &
           kbcon,ierr,dby,he,hes_cup, &
           itf,     jtf,     ktf,     &
           its,ite, jts,jte, kts,kte)
      !--- DETERMINE CLOUD TOP - KTOP
      call cup_ktop(1,dby,kbcon,ktop,ierr, &
           itf, jtf, ktf, &
           its,ite, jts,jte, kts,kte)
      do i = its, itf
         kzdown(i) = 0
         if (ierr(i) .eq. 0) then
            zktop = (z_cup(i, ktop(i)) - z1(i)) * .6
            zktop = min(zktop + z1(i), zcutdown + z1(i))
            do k = kts, ktf
               if (z_cup(i, k) .gt. zktop) then
                  kzdown(i) = k
                  exit
               endif
            enddo
         endif
      end do
      !--- doWNDRAFT ORIGINATING LEVEL - JMIN
      call cup_minimi(HEs_cup,K22,kzdown,JMIN,ierr, &
           itf, jtf, ktf, &
           its,ite, jts,jte, kts,kte)
      do 100 i=its,ite
         if (ierr(I) .eq. 0.) THEN
            if (jmin(i) .le. 3) then
               ierr(i) = 9
               go to 100
            endif
            !--- check whether it would have buoyancy, if there where
            !--- no entrainment/detrainment
101 continue
            if (jmin(i) - 1 .lt. KDET(I)) kdet(i) = jmin(i) - 1
            if (jmin(i) .ge. Ktop(I) - 1) jmin(i) = ktop(i) - 2
            ki = jmin(i)
            hcd(i, ki) = hes_cup(i, ki)
            DZ = Z_cup(i, Ki + 1) - Z_cup(i, Ki)
            dh = dz * (HCD(i, Ki) - hes_cup(i, ki))
            dh = 0.
            do k = ki - 1, 1, -1
               hcd(i, k) = hes_cup(i, jmin(i))
               dz = z_cup(i, k + 1) - z_cup(i, k)
               dh = dh + dz * (hcd(i, k) - hes_cup(i, k))
               if (dh .gt. 0.) then
                  jmin(i) = jmin(i) - 1
                  if (jmin(i) .gt. 3) then
                     go to 101
                  else if (jmin(i) .le. 3) then
                     ierr(i) = 9
                     go to 100
                  endif
               endif
            enddo
            if (jmin(i) .le. 3) then
               ierr(i) = 4
            endif
         endif
100 continue ! end loop for i
      ! - Must have at least depth_min m between cloud convective base
      !     and cloud top.
      do i = its, itf
         if (ierr(i) .eq. 0.) then
            if (jmin(i) .le. 3) then
               ierr(i) = 9
            endif
            if (z_cup(i, ktop(i)) - z_cup(i, kbcon(i)) .lt. depth_min) then
               ierr(i) = 6
            endif
         endif
      enddo

      !--- normalized updraft mass flux profile
      call cup_up_nms(zu,z_cup,mentr_rate,cd,kbcon,ktop,ierr,k22, &
           itf,     jtf,     ktf,     &
           its,ite, jts,jte, kts,kte)
      !--- normalized downdraft mass flux profile,also work on
      !--- bottom detrainment in this routine
      call cup_dd_nms(zd,z_cup,cdd,mentrd_rate,jmin,ierr, &
           0,kdet,z1,                 &
           itf,     jtf,     ktf,     &
           its,ite, jts,jte, kts,kte)
      !--- downdraft moist static energy
      call cup_dd_he(hes_cup,zd,hcd,z_cup,cdd,mentrd_rate, &
           jmin,ierr,he,dbyd,he_cup,  &
           itf,     jtf,     ktf,     &
           its,ite, jts,jte, kts,kte)
      !--- calculate moisture properties of downdraft
      call cup_dd_moisture(zd,hcd,hes_cup,qcd,qes_cup, &
           pwd,q_cup,z_cup,cdd,mentrd_rate,jmin,ierr,gamma_cup, &
           pwev,bu,qrcd,q,he,t_cup,1,xl, &
           itf,     jtf,     ktf,        &
           its,ite, jts,jte, kts,kte)
      !--- calculate moisture properties of updraft
      call cup_up_moisture(ierr,z_cup,qc,qrc,pw,pwav, &
           kbcon,ktop,cd,dby,mentr_rate,clw_all,      &
           q,GAMMA_cup,zu,qes_cup,k22,q_cup,xl, &
           itf,     jtf,     ktf,               &
           its,ite, jts,jte, kts,kte)
      !--- DETERMINE doWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
      call cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
           pwev,edtmax,edtmin,3,edtc, &
           itf,     jtf,     ktf,     &
           its,ite, jts,jte, kts,kte)
      do i = its, itf
         if (ierr(i) .eq. 0) then
            edt(i) = edtc(i, 2)
         endif
      enddo
      ! massflux from precip and normalized cloud properties
      pwdper = 0.
      do i = its, itf
         if (ierr(i) .gt. 0) then
            pre(i) = 0.
         endif
         if (ierr(i) .eq. 0) then
            xmb(i) = pre(i) / (pwav(i) + edt(i) * pwev(i))
            !--- percent of that that is evaporated (pwd is negative)
            if (i .eq. ipr .and. j .eq. jpr) then
               write(message, "(a,x,3F9.3)") 'xmb, edt, pwav = ', xmb(i), edt(i), pwav(i)
               call chemistry_message(trim(message))
               write(message, "(a,x,F9.3)") 'k, pwdper(i, k), pw, pwd(i, k)', z1(i)
               call chemistry_message(trim(message))
            endif
            do k = 1, ktop(i)
               pwdper(i, k) = - edt(i) * pwd(i, k) / pwav(i)
               if (i .eq. ipr .and. j .eq. jpr) then
                  write(message, "(I3,x,3F9.3)")  k, pwdper(i, k), pw(i, k), pwd(i, k)
                  call chemistry_message(trim(message))
               endif
            enddo
         endif
      enddo
!!! NOW WE HAVE EVREYTHING TO CALCULATE TRACER TRANSPORT AND WET DEPOSITION !!!
      !--- calculate incloud tracer distribution
      call cup_up_tracer(ierr,tcrit,t,pre,z_cup,p,tracer,tre_cup,&
           tr_up,tr_pw,tr_c,hstary,pw,&
           clw_all,kbcon,ktop,cd,mentr_rate,zu,k22,&
           num_gas,chem_opt,cu_scav_on,conv_tr_aqchem, chem, aq_gas_ratio, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte, &
           ipr,jpr,j,npr,num_chem,'deep')
      ! for shallow convection, only 3 subroutines required
      if (ishallow .eq. 1) then
         call cup_up_nms(zu3,z_cup,mentr_rate3,cd3,kbcon3,ktop3,ierr5,k23, &
              itf,     jtf,     ktf,     &
              its,ite, jts,jte, kts,kte)
         call cup_up_tracer(ierr5,tcrit,t,pre,z_cup,p,tracer,&
              tre_cup,tr3_up,tr3_pw,tr3_c,hstary,pw3,&
              clw_all3,kbcon3,ktop3,cd3,mentr_rate3,zu3,k23,&
              num_gas,chem_opt,cu_scav_on,conv_tr_aqchem, chem, aq_gas_ratio, &
              ids,ide, jds,jde, kds,kde, &
              ims,ime, jms,jme, kms,kme, &
              its,ite, jts,jte, kts,kte, &
              ipr,jpr,j,npr,num_chem,'shallow')
         call cup_dellas_tr(ierr5,z_cup,p_cup,tr_dd3,edt3,zd3,cdd3,&
              tracer,tracert3,j,mentrd_rate,zu3,g,xmb3,            &
              cd3,tr3_up,ktop3,k23,kbcon3,mentr_rate3,jmin,tre_cup,kdet, &
              k23,ipr,jpr,npr,'shallow',num_chem,                  &
              ids,ide, jds,jde, kds,kde,                           &
              ims,ime, jms,jme, kms,kme,                           &
              its,ite, jts,jte, kts,kte                            )
      endif

      call cup_dd_tracer(ierr,z_cup,qrcd,tracer,tre_cup,tr_up,tr_dd, &
           tr_pw,tr_pwd,jmin,cdd,mentrd_rate,zd,pwdper,wetdep,xmb,k22, &
           num_gas,num_chem,ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)

      !--- calculate transport tendencies
      !--- 1. in bottom layer
      call cup_dellabot_tr(ipr,jpr,tre_cup,ierr,z_cup,p,tr_dd,edt, &
           zd,cdd,tracer,tracert,j,mentrd_rate,z,g,xmb, &
           num_chem,ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      !--- 2. everywhere else
      call cup_dellas_tr(ierr,z_cup,p_cup,tr_dd,edt,zd,cdd,      &
           tracer,tracert,j,mentrd_rate,zu,g,xmb,                &
           cd,tr_up,ktop,k22,kbcon,mentr_rate,jmin,tre_cup,kdet, &
           k22,ipr,jpr,npr,'deep',num_chem,                      &
           ids,ide, jds,jde, kds,kde,                            &
           ims,ime, jms,jme, kms,kme,                            &
           its,ite, jts,jte, kts,kte                             )
      ! may need more below for wet deposition......
      !
      !
      ! call cup_output_wd (   &
      !      ids,ide, jds,jde, kds,kde, &
      !      ims,ime, jms,jme, kms,kme, &
      !      its,ite, jts,jte, kts,kte)

   end subroutine cup_ct

   subroutine cup_dellabot_tr(ipr, jpr, tre_cup, ierr, z_cup, p_cup, &
                              tr_dd, edt, zd, cdd, tracer, tracert, j, &
                              mentrd_rate, z, g, xmb, &
                              num_chem, ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              its, ite, jts, jte, kts, kte)

      integer, intent (in) :: num_chem
      integer, intent (in) :: ids, ide, jds, jde, kds, kde
      integer, intent (in) :: ims, ime, jms, jme, kms, kme
      integer, intent (in) :: its, ite, jts, jte, kts, kte
      integer, intent (in) :: j, ipr, jpr
      ! ierr error value, maybe modified in this routine
      real, dimension (its:ite, kts:kte, 1:num_chem), intent (out) :: tracert
      real, dimension (its:ite, kts:kte, 1:num_chem), intent (in) :: tre_cup, tracer, tr_dd
      real, dimension (its:ite, kts:kte), intent (in) :: z_cup, p_cup, zd, cdd, z
      real, dimension (its:ite), intent (in) :: edt, xmb
      real, intent (in) :: g, mentrd_rate
      integer, dimension (its:ite), intent (inout) :: ierr
      ! local variables in this routine
      integer i
      real detdo, detdo1, detdo2, entdo, dp, dz, subin, totmas
      integer :: itf, ktf, nv, npr
      character(len=300) message
      npr=24
      itf=MIN(ite,ide-1)
      ktf=MIN(kte,kde-1)
      if (j .eq. jpr) call chemistry_message('in cup dellabot ')
      tracert = 0.
      do i = its, itf
         if (ierr(i) .ne. 0) cycle
         dz = z_cup(i, 2) - z_cup(i, 1)
         DP = 100. * (p_cup(i, 1) - P_cup(i, 2))
         detdo1 = edt(i) * zd(i, 2) * CDD(i, 1) * DZ
         detdo2 = edt(i) * zd(i, 1)
         entdo = edt(i) * zd(i, 2) * mentrd_rate * dz
         subin = - EDT(I) * zd(i, 2)
         detdo = detdo1 + detdo2 - entdo + subin
         do nv = 1, num_chem
            tracert(I, 1, nv) = (detdo1 * .5 * (tr_dd(i, 1, nv) + tr_dd(i, 2, nv)) &
                              + detdo2 * tr_dd(i, 1, nv) &
                              + subin * tre_cup(i, 2, nv) &
                              - entdo * tracer(i, 1, nv)) * g / dp * xmb(i)
         enddo
         if (j .eq. jpr .and. i .eq. ipr) then
            write(message, "(a,x,8F9.3)") 'in cup dellabot ', tracert(I, 1, npr) , &
            detdo1, detdo2, subin, entdo, tr_dd(i, 1, npr), tr_dd(i, 2, npr), tracer(i, 1, npr)
            call chemistry_message(trim(message))
         end if
      enddo
   end subroutine cup_dellabot_tr

   subroutine cup_dellas_tr(ierr, z_cup, p_cup, tr_dd, edt, zd, cdd, &
                           tracer, tracert, j, mentrd_rate, zu, g, xmb, &
                           cd, tr_up, ktop, k22, kbcon, mentr_rate, jmin, &
                           tre_cup, kdet, kpbl, &
                           ipr, jpr, npr, name, num_chem, &
                           ids, ide, jds, jde, kds, kde, &
                           ims, ime, jms, jme, kms, kme, &
                           its, ite, jts, jte, kts, kte)

      integer, intent (in) :: num_chem, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte
      integer, intent (in) :: j, ipr, jpr, npr
      ! ierr error value, maybe modified in this routine
      real, dimension (its:ite, kts:kte, 1:num_chem), intent (inout) :: tracert
      real, dimension (its:ite, kts:kte, 1:num_chem), intent (in) :: tr_up, tr_dd, tre_cup, tracer
      real, dimension (its:ite, kts:kte), intent (in) :: z_cup, p_cup, zd, cdd, cd, zu
      real, dimension (its:ite), intent (in) :: edt, xmb
      real, intent (in) :: g, mentrd_rate, mentr_rate
      integer, dimension (its:ite), intent (in) :: kbcon, ktop, k22, jmin, kdet, kpbl
      integer, dimension (its:ite), intent (inout) :: ierr
      character * (*), intent (in) :: name
      ! local variables in this routine

      integer i, k, nv
      real detdo1, detdo2, entdo, dp, dz, subin, detdo, entup, detup, subdown, entdoj, entupk, detupk, totmas
      integer :: itf, ktf, kstart
      character(len=300) message
      ! npr=24
      itf = MIN(ite, ide - 1)
      ktf = MIN(kte, kde - 1)
      kstart = kts + 1
      if (name .eq. 'shallow') kstart = kts
      i = ipr
      if (j .eq. jpr) then
         call chemistry_message('in dellas kpbl(i),k22(i),kbcon(i),ktop(i),jmin(i)')
         write(message, "(5I4)") kpbl(i),k22(i),kbcon(i),ktop(i),jmin(i)
         call chemistry_message(trim(message))
      endif
      do nv = 1, num_chem
         do k = kstart, kte
            do i = its, itf
               tracert(i, k, nv) = 0.
            enddo
         enddo
      enddo
      do k = kts + 1, ktf - 1
         do i = its, ite
            if (ierr(i) .ne. 0) cycle
            if (K .gt. KTOP(I)) cycle
            !--- specify detrainment of downdraft, has to be consistent
            !--- with zd calculations in soundd.
            DZ = Z_cup(I, K + 1) - Z_cup(I, K)
            detdo = edt(i) * CDD(i, K) * DZ * ZD(i, k + 1)
            entdo = edt(i) * mentrd_rate * dz * zd(i, k + 1)
            subin = zu(i, k + 1) - zd(i, k + 1) * edt(i)
            entup = 0.
            detup = 0.
            if (k .ge. kbcon(i) .and. k .lt. ktop(i)) then
               entup = mentr_rate * dz * zu(i, k)
               detup = CD(i, K + 1) * DZ * ZU(i, k)
            endif
            subdown = (zu(i, k) - zd(i, k) * edt(i))
            entdoj = 0.
            entupk = 0.
            detupk = 0.
!
            if (k .eq. jmin(i)) then
               entdoj = edt(i) * zd(i, k)
            endif

            if (k .eq. k22(i) - 1) then
               entupk = zu(i, kpbl(i))
            endif

            if (k .gt. kdet(i)) then
               detdo = 0.
            endif

            if (k .eq. ktop(i) - 0) then
               detupk = zu(i, ktop(i))
               subin = 0.
            endif
            if (k .lt. kbcon(i)) then
               detup = 0.
            endif
            !--- CHANGED DUE TO SUBSIDENCE AND ENTRAINMENT
            totmas = subin - subdown + detup - entup - entdo + &
                     detdo - entupk - entdoj + detupk
            if (j .eq. jpr .and. i .eq. ipr) then
               write(message, "(a,x,I3,x,3F9.3)") 'k, totmas, sui, sud = ', k, totmas, subin, subdown
               call chemistry_message(trim(message))
            end if
            if (abs(totmas) .gt. 1.e-6) then
               write(message, "(a,3I3,x,F9.3,a)") '*********************',i,j,k,totmas,name
               call chemistry_message(trim(message))
               write(message, "(4I5)") kpbl(i), k22(i), kbcon(i), ktop(i)
               call chemistry_message(trim(message))
               call chemistry_error_fatal( 'cup_dellas_tr: TOTMAS > CRITICAL VALUE')
            endif
            dp = 100. * (p_cup(i, k - 1) - p_cup(i, k))
            do nv = 1, num_chem
               ! tracert(i,k,nv)=(subin*tre_cup(i,k+1,nv) &
               !            -subdown*tre_cup(i,k,nv) &
               tracert(i,k,nv)=(subin*tracer(i,k+1,nv) &
                  -subdown*tracer(i,k,nv) &
                  +detup*.5*(tr_up(i,K+1,nv)+tr_up(i,K,nv)) &
                  +detdo*.5*(tr_dd(i,K+1,nv)+tr_dd(i,K,nv)) &
                  -entup*tracer(i,k,nv) &
                  -entdo*tracer(i,k,nv) &
                  -entupk*tre_cup(i,k22(i),nv) &
                  -entdoj*tre_cup(i,jmin(i),nv) &
                  +detupk*tr_up(i,ktop(i),nv) &
                  )*g/dp*xmb(i)
            enddo
            if(i.eq.ipr.and.j.eq.jpr)then
               write(message, "(I3,x,4F9.3)") k, tracert(i, k, npr), subin * tre_cup(i, k + 1, npr), subdown * tre_cup(i, k, npr), &
                                 detdo * .5 * (tr_dd(i, K + 1, npr) + tr_dd(i, K, npr))
               call chemistry_message(trim(message))
               write(message, "(I3,x,4F9.3)") k, detup * .5 * (tr_up(i, K + 1, npr) + tr_up(i, K, npr)), detupk * tr_up(i, ktop(i), npr), &
                                 entup * tracer(i, k, npr), entdo * tracer(i, k, npr)
               call chemistry_message(trim(message))
               write(message, "(I3,x,3F9.3)") k, entupk * tre_cup(i, k, npr), detupk, tr_up(i, ktop(i), npr)
               call chemistry_message(trim(message))
            endif
         enddo
      enddo
   end subroutine cup_dellas_tr

   subroutine cup_env_clev_tr(tracer,tre_cup,num_chem,ierr, &
                              ids,ide, jds,jde, kds,kde, &
                              ims,ime, jms,jme, kms,kme, &
                              its,ite, jts,jte, kts,kte)

      integer, intent(in) :: num_chem
      integer, intent(in) :: ids, ide, jds, jde, kds, kde
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: its, ite, jts, jte, kts, kte
      integer, dimension (its:ite), intent(in) :: ierr
      real, dimension (its:ite, kts:kte, 1:num_chem), intent(in) :: tracer
      real, dimension (its:ite, kts:kte, 1:num_chem), intent(out) :: tre_cup
      !  local variables in this routine
      integer :: i, k, nv, itf, ktf
      itf = MIN(ite, ide - 1)
      ktf = MIN(kte, kde - 1)
      do nv = 1, num_chem
         do k = kts + 1, ktf
            do i = its, ite
               if (ierr(i) .eq. 0) then
                  tre_cup(i, k, nv) = .5 * (tracer(i, k - 1, nv) + tracer(i, k, nv))
               endif
            enddo
         enddo
      enddo
      do nv = 1, num_chem
         do i = its, ite
            if (ierr(i) .eq. 0) then
               tre_cup(i, kts, nv) = tracer(i, kts, nv)
            endif
         enddo
      enddo

   end subroutine cup_env_clev_tr

   subroutine  cup_up_tracer(ierr, tcrit, t, pre, z_cup, p, tracer, tre_cup, tr_up, &
               tr_pw, tr_c, hstary, cupclw, clw_all, kbcon, ktop, cd, mentr_rate, zu, k22, &
               num_gas, chem_opt, cu_scav_on, conv_tr_aqchem, chem, aq_gas_ratio, &
               ids, ide, jds, jde, kds, kde, &
               ims, ime, jms, jme, kms, kme, &
               its, ite, jts, jte, kts, kte, ipr, jpr, j, npr, num_chem, name)
      ! USE module_configure
!       use module_state_description, only: RADM2SORG,RADM2SORG_AQ,RACMSORG_AQ,RACMSORG_KPP,   &
!                                           RADM2SORG_KPP,RACM_ESRLSORG_KPP,RACM_SOA_VBS_KPP,  &
!                                           RADM2SORG_AQCHEM,RACMSORG_AQCHEM_KPP,RACM_ESRLSORG_AQCHEM_KPP
      use module_ctrans_aqchem, only: aqchem
      use module_data_binaer_asect, only: nsize_aer, ncomp_aer, massptr_aer, dens_aer, &
          mastercompindx_seas_aer, mastercompindx_dust_aer, mastercompindx_oc_aer

      ! Aqeuous species pointers INCLUDE File
      !...........PARAMETERS and their descriptions:
      integer, parameter :: NGAS = 12     ! number of gas-phase species for AQCHEM
      integer, parameter :: NAER = 36     ! number of aerosol species for AQCHEM
      integer, parameter :: NLIQS = 41    ! number of liquid-phase species in AQCHEM
      !...pointers for the AQCHEM array GAS
      integer, parameter :: LSO2    =  1  ! Sulfur Dioxide
      integer, parameter :: LHNO3   =  2  ! Nitric Acid
      integer, parameter :: LN2O5   =  3  ! Dinitrogen Pentoxide
      integer, parameter :: LCO2    =  4  ! Carbon Dioxide
      integer, parameter :: LNH3    =  5  ! Ammonia
      integer, parameter :: LH2O2   =  6  ! Hydrogen Perioxide
      integer, parameter :: LO3     =  7  ! Ozone
      integer, parameter :: LFOA    =  8  ! Formic Acid
      integer, parameter :: LMHP    =  9  ! Methyl Hydrogen Peroxide
      integer, parameter :: LPAA    = 10  ! Peroxyacidic Acid
      integer, parameter :: LH2SO4  = 11  ! Sulfuric Acid
      integer, parameter :: LHCL    = 12  ! Hydrogen Chloride
      !...pointers for the aqchEM array AEROSOL
      integer, parameter :: LSO4AKN  =  1  ! Aitken-mode Sulfate
      integer, parameter :: LSO4ACC  =  2  ! Accumulation-mode Sulfate
      integer, parameter :: LSO4COR  =  3  ! Coarse-mode Sulfate
      integer, parameter :: LNH4AKN  =  4  ! Aitken-mode Ammonium
      integer, parameter :: LNH4ACC  =  5  ! Accumulation-mode Ammonium
      integer, parameter :: LNO3AKN  =  6  ! Aitken-mode Nitrate
      integer, parameter :: LNO3ACC  =  7  ! Accumulation-mode Nitrate
      integer, parameter :: LNO3COR  =  8  ! Coarse-mode Nitrate
      integer, parameter :: LORGAAKN =  9  ! Aitken-mode anthropogenic SOA
      integer, parameter :: LORGAACC = 10  ! Accumulation-mode anthropogenic SOA
      integer, parameter :: LORGPAKN = 11  ! Aitken-mode primary organic aerosol
      integer, parameter :: LORGPACC = 12  ! Accumulation-mode primary organic aerosol
      integer, parameter :: LORGBAKN = 13  ! Aitken-mode biogenic SOA
      integer, parameter :: LORGBACC = 14  ! Accumulation-mode biogenic SOA
      integer, parameter :: LECAKN   = 15  ! Aitken-mode elemental carbon
      integer, parameter :: LECACC   = 16  ! Accumulation-mode elemental carbon
      integer, parameter :: LPRIAKN  = 17  ! Aitken-mode primary aerosol
      integer, parameter :: LPRIACC  = 18  ! Accumulation-mode primary aerosol
      integer, parameter :: LPRICOR  = 19  ! Coarse-mode primary aerosol
      integer, parameter :: LNAAKN   = 20  ! Aitken-mode Sodium
      integer, parameter :: LNAACC   = 21  ! Accumulation-mode Sodium
      integer, parameter :: LNACOR   = 22  ! Coarse-mode Sodium
      integer, parameter :: LCLAKN   = 23  ! Aitken-mode Chloride ion
      integer, parameter :: LCLACC   = 24  ! Accumulation-mode Chloride ion
      integer, parameter :: LCLCOR   = 25  ! Coarse-mode Chloride ion
      integer, parameter :: LNUMAKN  = 26  ! Aitken-mode number
      integer, parameter :: LNUMACC  = 27  ! Accumulation-mode number
      integer, parameter :: LNUMCOR  = 28  ! Coarse-mode number
      integer, parameter :: LSRFAKN  = 29  ! Aitken-mode surface area
      integer, parameter :: LSRFACC  = 30  ! Accumulation-mode surface area
      integer, parameter :: LNACL    = 31  ! Sodium Chloride aerosol for AE3 only {depreciated in AE4}
      integer, parameter :: LCACO3   = 32  ! Calcium Carbonate aerosol (place holder)
      integer, parameter :: LMGCO3   = 33  ! Magnesium Carbonate aerosol (place holder)
      integer, parameter :: LA3FE    = 34  ! Iron aerosol (place holder)
      integer, parameter :: LB2MN    = 35  ! Manganese aerosol (place holder)
      integer, parameter :: LK       = 36  ! Potassium aerosol (Cl- tracked separately) (place holder)
      ! on input
      ! only local wrf dimensions are need as of now in this routine
      integer, intent(in) :: ids, ide, jds, jde, kds, kde
      integer, intent(in) :: num_gas, conv_tr_aqchem
      integer, intent(in) :: num_chem, ims, ime, jms, jme, kms, kme
      integer, intent(in) :: its,ite, jts,jte, kts,kte,ipr,jpr,j,npr
      logical, intent(in) :: cu_scav_on
      real, dimension (its:ite, kts:kte), intent (in) :: z_cup, cd, zu, p, hstary, t
      real, dimension (its:ite, kts:kte), intent (inout) :: cupclw, clw_all
      real, dimension (its:ite, kts:kte, 1:num_chem), intent (inout) :: tr_up, tr_c, tr_pw
      real, dimension (its:ite, kts:kte, 1:num_chem), intent (in) :: tre_cup, tracer
      real, dimension (its:ite), intent (in) :: pre
      character(len = *), intent(in) :: name, chem_opt
      ! entr = entrainment rate
      real, intent (in) :: mentr_rate, tcrit
      integer, dimension (its:ite), intent (in) :: kbcon, ktop, k22
      ! chem is used to calculate volume of aerosols
      real, dimension(ims:ime, kms:kme, jms:jme, num_chem), intent(in) :: chem
      real, dimension(ims:ime, kms:kme, jms:jme), intent (inout) :: aq_gas_ratio ! Output to diagnostics
      ! ierr error value, maybe modified in this routine

      integer, dimension (its:ite), intent (inout) :: ierr
      !  local variables in this routine
      real :: conc_equi, conc_mxr, partialp, taucld
      real :: trcc, trch, dh, c0, dz, radius, airm, dens
      integer :: iall, i, k, iwd, nv
      integer :: itf, ktf, iaer, igas
      integer :: iComp, iChem, iType, iSize, iPhase
      ! Gas / aqueous phase partitioning for wet scavenging / deposition of gas
      ! phase and aerosol species:
      real, dimension(1:ncomp_aer(1)) :: aq_gas_ratio_aer ! jwfeng
      ! I/O for AQCHEM:
      real precip ! Precipitation rate (mm/h)
      real, dimension (ngas) :: gas, gaswdep
      real, dimension (naer) :: aerosol, aerwdep
      real, dimension (nliqs) :: liquid
      real hpwdep
      real alfa0, alfa2, alfa3 ! Aerosol scavenging coefficients for Aitken mode
      ! jwfeng: Total aerosol mass inside grid box.
      real :: total_ratio, total_vol, temp_vol
      character(len=300) message

      itf = MIN(ite, ide - 1)
      ktf = MIN(kte, kde - 1)
      iall = 0
      c0 = .002
      c0 = .8 ! test by jwfeng
      aq_gas_ratio(:,:,:) = 0.0
      iwd = 0
!--- no precip for small clouds
!      if(mentr_rate.gt.0.)then
!        radius=.2/mentr_rate
!        if(radius.lt.900.)c0=0.
      if (name .eq. 'shallow') c0 = 0.
!        if(radius.lt.900.)iall=0
!      endif
      do nv = 1, num_chem
         do k = kts, ktf
            do i = its, itf
        ! Initialize species mass in rain water:
               tr_pw(i,k,nv)=0.
        ! Initialize species mass in updraft:
               if(ierr(i).eq.0)tr_up(i,k,nv)=tre_cup(i,k,nv)
        ! Initialize species mass in cloud + rain water:
               tr_c(i,k,nv)=0.
            enddo
         enddo
      enddo

      do nv=1,num_chem
         do i=its,itf
            if(ierr(i).eq.0.)then
               do k=k22(i),kbcon(i)-1
                  tr_up(i,k,nv)=tre_cup(i,k22(i),nv)
               enddo
            endif
         enddo
      enddo

      if (cu_scav_on) then
         iType = 1
         iPhase = 1 ! for now, just calc for iphase = ai_phase
         aq_gas_ratio_aer(1:ncomp_aer(iType)) = 0.0
         aq_gas_ratio_aer(mastercompindx_seas_aer) = 1.0
         aq_gas_ratio_aer(mastercompindx_dust_aer) = 0.1
         aq_gas_ratio_aer(mastercompindx_oc_aer)   = 0.1
         do k = kts, ktf
            do i = its, itf
               total_vol = 0.0
               total_ratio= 0.0
               do iSize = 1, nsize_aer(iType)
                  do iComp = 1, ncomp_aer(iType)
                     iChem = massptr_aer(iComp,iSize,iType,iPhase)
                     temp_vol = max(chem(i, k, j, iChem), 0.0) / dens_aer(iComp, iType)
                     total_vol = total_vol + temp_vol
                     total_ratio = total_ratio + temp_vol * aq_gas_ratio_aer(iComp)
                  end do
               end do
               aq_gas_ratio(i,k,j) = total_ratio / total_vol
               if (aq_gas_ratio(i,k,j)<0.0 .or. aq_gas_ratio(i,k,j)>1.0) then
                  write(message, "(a,x,f9.4)") "Wrong aq_gas_ratio = ", aq_gas_ratio(i,k,j)
                  call chemistry_error_fatal(message)
               end if
            end do
         end do
      end if


      do k = kts + 1, ktf - 1
         do i = its, itf
            !if (k == kts + 1) then
            !   write(message, "(a, I3, a, I3)") "KBCON/KTOP = ", kbcon(i), "/", ktop(i)
            !   call chemistry_message(trim(message))
            !end if
            if (ierr(i) .ne. 0) cycle
            if (K .lt. KBCON(I)) cycle
            if (K .gt. KTOP(I) + 1) cycle
            DZ = Z_cup(i, K) - Z_cup(i, K - 1)
            if (cupclw(i, k) .le. 0.) cupclw(i, k) = 0.
            if (clw_all(i, k) .le. 0.) clw_all(i, k) = 0.
      ! Steady state plume equation, for what could be in cloud before anything
      ! happens (kg/kg). tr_up would be the concentration if tracers were conserved.
      !

            do nv = 1, num_chem
               if (i .eq. ipr .and. j .eq. jpr .and. nv .eq. npr) then
                  write(message, "(I3,x,5F9.3)") k, tr_up(i, k - 1, nv), tr_up(i, k, nv), tr_pw(i, k - 1, nv), clw_all(i, k), cupclw(i, k)
                  call chemistry_message(trim(message))
               end if
               tr_up(i, k, nv) = (tr_up(i, k - 1, nv) * (1. - .5 * cd(i, k) * dz) + mentr_rate * &
               dz * tracer(i, k - 1, nv)) / (1. + mentr_rate * dz - .5 * cd(i, k) * dz)
               tr_up(i, k, nv) = max(1.e-16, tr_up(i, K, nv))
            enddo
            ! Aqueous chemistry
            ! if ((chem_opt .EQ. RADM2SORG .OR. chem_opt .EQ. RADM2SORG_AQ .OR. chem_opt .EQ. RACMSORG_AQ .OR. & 
            !    chem_opt .EQ. RACMSORG_KPP .OR. chem_opt .EQ. RADM2SORG_KPP .OR. chem_opt .EQ. RACM_ESRLSORG_KPP .OR. & 
            !    chem_opt .EQ. RACM_SOA_VBS_KPP .OR. chem_opt .EQ. RADM2SORG_AQCHEM .OR. chem_opt .EQ. RACMSORG_AQCHEM_KPP .OR. &
            !    chem_opt .EQ. RACM_ESRLSORG_AQCHEM_KPP) &
            !    .and. (conv_tr_aqchem == 1)) then
            if (conv_tr_aqchem == 1) then
               ! For MADE/SORGAM derived schemes with aqueous chemistry
               ! Air mass density
               dens = 0.1*p(i,k)/t(i,k)*mwdry/8.314472 ! kg/m3
               ! Column air number density:
               airm = 1000.0*dens*dz/mwdry ! mol/m2
               ! Wet scavenging initialization for AQCHEM
               GASWDEP = 0.0
               AERWDEP = 0.0
               HPWDEP  = 0.0
               ! We provide a precipitation rate and aerosol scavenging rates of zero,
               ! in order to prevent wet scavenging in AQCHEM (it is treated later):
               precip = 0.0 ! mm/hr
               alfa0 = 0.0
               alfa2 = 0.0
               alfa3 = 0.0
               ! Gas phase concentrations before aqueous phase chemistry
               ! (with units conversion ppm -> mol/mol)
               gas(:) = 0.0
               gas(lco2)   = 380.0e-6
               gas(lso2)   = tr_up(i,k,p_so2)*1.0e-6
               gas(lhno3)  = tr_up(i,k,p_hno3)*1.0e-6
               gas(ln2o5)  = tr_up(i,k,p_n2o5)*1.0e-6
               gas(lnh3)   = tr_up(i,k,p_nh3)*1.0e-6
               gas(lh2o2)  = tr_up(i,k,p_h2o2)*1.0e-6
               gas(lo3)    = tr_up(i,k,p_o3)*1.0e-6
               gas(lfoa)   = tr_up(i,k,p_ora1)*1.0e-6
               gas(lmhp)   = tr_up(i,k,p_op1)*1.0e-6
               gas(lpaa)   = tr_up(i,k,p_paa)*1.0e-6
               gas(lh2so4) = tr_up(i,k,p_sulf)*1.0e-6
               ! Aerosol mass concentrations before aqueous phase chemistry
               ! (with units conversion ug/kg -> mol/mol). Although AQCHEM
               ! accounts for much of the aerosol compounds in MADE, they are
               ! not treated at the moment by AQCHEM, as the mapping between
               ! the organic compound groups in MADE and AQCHEM is not obvious.
               aerosol(:) = 0.0
               ! We assume all accumulation mode particles are activated in cumulus clouds:
               aerosol(lso4acc) = tr_up(i,k,p_so4aj)*1.0e-9*mwdry/mwso4
               aerosol(lnh4acc) = tr_up(i,k,p_nh4aj)*1.0e-9*mwdry/mwnh4
               aerosol(lno3acc) = tr_up(i,k,p_no3aj)*1.0e-9*mwdry/mwno3
               ! Cloud lifetime:
               taucld = 1800.0
               if (clw_all(i, k) * dens .gt. qcldwtr_cutoff) then ! Cloud water > threshold
                  CALL AQCHEM(          &
                           t(i,k)     , &
                           p(i,k)*100., &
                           taucld     , &
                           precip     , &
                           clw_all(i,k)*dens, &
                           clw_all(i,k)*dens, &
                           airm , &
                           ALFA0, &
                           ALFA2, &
                           ALFA3, &
                           GAS  , &
                           AEROSOL, &
                           LIQUID , &
                           GASWDEP, &
                           AERWDEP, &
                           HPWDEP)
               endif
               ! Gas phase concentrations after aqueous phase chemistry
               ! (with units conversion mol/mol -> ppm)
               tr_up(i,k,p_so2)  =  gas(lso2)*1.0e6
               tr_up(i,k,p_hno3) =  gas(lhno3)*1.0e6
               tr_up(i,k,p_n2o5) =  gas(ln2o5)*1.0e6
               tr_up(i,k,p_nh3)  =  gas(lnh3)*1.0e6
               tr_up(i,k,p_h2o2) =  gas(lh2o2)*1.0e6
               tr_up(i,k,p_o3)   =  gas(lo3)*1.0e6
               tr_up(i,k,p_ora1) =  gas(lfoa)*1.0e6
               tr_up(i,k,p_op1)  =  gas(lmhp)*1.0e6
               tr_up(i,k,p_paa)  =  gas(lpaa)*1.0e6
               tr_up(i,k,p_sulf) =  gas(lh2so4)*1.0e6

               ! Aerosol mass concentrations
               ! (with units conversion mol/mol -> ug/kg)
               tr_up(i,k,p_so4aj) = aerosol(lso4acc)*1.0e9*mwso4/mwdry
               tr_up(i,k,p_nh4aj) = aerosol(lnh4acc)*1.0e9*mwnh4/mwdry
               tr_up(i,k,p_no3aj) = aerosol(lno3acc)*1.0e9*mwno3/mwdry
            endif
            ! wet scavenging option (turn off by setting
            if (cu_scav_on) then
               !if (k == 20) then
               !   write(message, "(a,I3,a,I3,a,I3,a,F9.3,a,F9.3)") "i/k/j = ", i, "/", k, "/", j, "total_vol = ", total_ratio, ", aq_gas_ratio = ", aq_gas_ratio(i,k,j)
               !   call chemistry_message(trim(message))
               !end if
               do nv = 1, num_chem
                  ! aq_gas_ratio = 0.0
                  ! ! Fraction of gas phase species that partions into the liquid phase:
                  ! if (nv.eq.p_so2)  aq_gas_ratio = 1.0
                  ! if (nv.eq.p_sulf) aq_gas_ratio = 1.0
                  ! if (nv.eq.p_nh3)  aq_gas_ratio = 1.0
                  ! if (nv.eq.p_hno3) aq_gas_ratio = 1.0
                  ! if (nv.gt.num_gas) aq_gas_ratio = 0.5
                  ! if (nv.eq.p_bc1 .or. nv.eq.p_oc1 .or. nv.eq.p_dms) aq_gas_ratio=0.
                  ! if (nv.eq.p_sulf .or. nv.eq.p_seas_1 .or. nv.eq.p_seas_2) aq_gas_ratio=1.
                  ! if (nv.eq.p_seas_3 .or. nv.eq.p_seas_4) aq_gas_ratio=1.
                  ! if (nv.eq.p_bc2 .or. nv.eq.p_oc2) aq_gas_ratio=0.8
                  if (aq_gas_ratio(i,k,j) > 0.0) then
                     tr_c(i,k,nv)  = aq_gas_ratio(i,k,j)*tr_up(i,k,nv) ! Amount of species cloud + rain water
                     trch          = tr_up(i,k,nv)-tr_c(i,k,nv) ! Amount of species remaining in gas phase
                     ! trcc          = (tr_up(i,k,nv)-trch)/(1.+c0*dz*zu(i,k)) ! Amount of species cloud + rain water
                     ! tr_pw(i,k,nv) = c0*dz*trcc*zu(i,k) ! Amount of species in rain water
                     tr_pw(i,k,nv) = c0*dz*tr_c(i,k,nv)*zu(i,k) ! Amount of species in rain water
                     trcc          = tr_c(i,k,nv)*(1 - c0) ! Amount of species cloud + rain water
                     tr_up(i,k,nv) = trcc+trch ! Total amount of species in updraft = conc in liq water (trcc) + conc in gas phase (trch)
                  endif
               enddo
            endif ! parameterization wetscavenging option
         enddo
      enddo
   end subroutine cup_up_tracer


   SUBROUTINE CUP_DD_TRACER(IERR, Z_CUP, QRCD, TRACER, TRE_CUP, TR_UP, TR_DD, &
         tr_pw,tr_pwd,jmin,cdd,entr,zd,pwdper,wetdep,xmb,k22, &
         num_gas,num_chem,ids,ide, jds,jde, kds,kde, &
         ims,ime, jms,jme, kms,kme, &
         its,ite, jts,jte, kts,kte)

      use mpas_atmchem_constants, only: mwdry


      ! only local wrf dimensions are need as of now in this routine

      integer, intent(in) :: num_gas, num_chem
      integer, intent(in) :: ids, ide, jds, jde, kds, kde
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: its, ite, jts, jte, kts, kte
      real, dimension(its:ite, kts:kte), intent(in) :: pwdper, zd, cdd, qrcd, z_cup
      real, dimension(its:ite), intent(in) :: xmb
      real, dimension(its:ite, kts:kte, 1:num_chem), intent(inout) :: tr_dd, tr_pwd, tr_up
      real, dimension(its:ite, kts:kte, 1:num_chem), intent(in) :: tre_cup, tracer, tr_pw
      real, dimension(its:ite, 1:num_chem), intent(out) :: wetdep

      ! entr = entrainment rate
      real, intent(in) :: entr
      integer, dimension (its:ite), intent(in) :: jmin
      ! ierr error value, maybe modified in this routine
      integer, dimension (its:ite), intent(inout) :: ierr, k22
      ! local variables in this routine
      integer :: iall, i, k, nv, ki, kj
      real :: dh, c0, dz, radius
      integer :: itf, ktf
      real :: evaporate, condensate

      itf = MIN(ite, ide - 1)
      ktf = MIN(kte, kde - 1)

      ! Initialize the tracer amount that evaporated from rain water:
      tr_pwd(its:ite, kts:kte, 1:num_chem) = 0.0

      ! Initialize wet deposition:
      wetdep(its:ite, 1:num_chem) = 0.0

      ! Calculate wet deposition with re-evaporation (based on wet scavenging in cup_up_tracer);
      ! assume no gas takeup by rain during fall for now

      do i=its,ite
         if(ierr(i).eq.0)then
            do nv=1,num_gas ! only gas phase species evaporate along with rain water
               do k=ktf,kts,-1 ! descending loop over all levels
                  ! start with initializing the condensate with the tracer amount in rain water on this level:
                  condensate = tr_pw(i,k,nv)
                  do kj=k,kts,-1 ! descending loop over this and lower levels
                     ! evaporated tracer amount at the current level:
                     evaporate = condensate*pwdper(i,kj)
                     ! accumulate the evaporated tracer amount:
                     tr_pwd(i,kj,nv) = tr_pwd(i,kj,nv) + evaporate
                     ! Remove the evaporated tracer amount from the condensate:
                     condensate = max(0.0,condensate - evaporate)
                  enddo
               enddo
               ! Calculate the wet deposition as the column integral over the
               ! tracer amount in rain water - evaporated tracer amount:
               do k=kts,ktf
                  wetdep(i,nv) = wetdep(i,nv) + tr_pw(i,k,nv) - tr_pwd(i,k,nv)
               enddo
               wetdep(i,nv) = max(0.0,wetdep(i,nv))
            enddo
         endif
      enddo
      ! In downdraft, do only transport of tracers
      do i=its,ite
         if(ierr(i).eq.0)then
            do nv=1,num_chem
               tr_dd(i,jmin(i),nv)=tre_cup(i,jmin(i),nv) ! Tracer amount in downdraft
            enddo
            do ki=jmin(i)-1,1,-1
               DZ=Z_cup(i,ki+1)-Z_cup(i,ki)
               do nv=1,num_chem
                  tr_dd(i,ki,nv)=(tr_dd(i,ki+1,nv)*(1.-.5*CDD(i,ki)*DZ) &
                     +entr*DZ*tracer(i,ki,nv) &
                     )/(1.+entr*DZ-.5*CDD(i,ki)*DZ)
               enddo
            enddo
         endif
      enddo
      ! Add evaporation from rain water:
      do i=its,ite
         IF(ierr(I).eq.0)then
            do nv=1,num_chem
               do k=kts,ktf
                  tr_dd(i,k,nv) = tr_dd(i,k,nv) + tr_pwd(i,k,nv)
               enddo
            enddo
         endif
      enddo
      ! To obtain wet deposition, multiply with base mass flux; units are given
      ! assuming that the base mass flux units [xmb] = kg(air)/m2/s:
      do nv=1,num_chem
         if (nv <= num_gas) then
            do i=its,itf
               if(ierr(I).eq.0)then
                  wetdep(i,nv)=wetdep(i,nv)*xmb(i)/mwdry ! mmol/m2/s for gas phase species
                  wetdep(i,nv) = max(0.0,wetdep(i,nv))
               endif
            enddo
         else
            do i=its,itf
               if(ierr(I).eq.0)then
                  wetdep(i,nv)=wetdep(i,nv)*xmb(i) ! ug/m2/s for aerosol mass, #/m2/s for aerosol number
                  wetdep(i,nv) = max(0.0,wetdep(i,nv))
               endif
            enddo
         endif
      enddo

   end subroutine cup_dd_tracer


   subroutine neg_check_ct(name, pret, ktop, epsilc, dt, q, outq, iopt, &
              num_chem, its, ite, kts, kte, itf, ktf, ipr, jpr, npr, j)

      integer, intent(in) :: iopt, num_chem, its, ite, kts, kte
      integer, intent(in) :: itf, ktf, ipr, jpr, npr, j

      real, dimension(its:ite, kts:kte, num_chem), intent(inout) :: q, outq
      real, dimension(its:ite), intent(in) :: pret
      integer, dimension(its:ite), intent(in) :: ktop
      real, intent(in) :: dt, epsilc
      real :: tracermin, tracermax, thresh, qmem, qmemf, qmem2, qtest, qmem1
      character*(*) name
      character(len=300) message
      integer nv, i, k
      ! check whether routine produces negative q's. This can happen, since
      ! tendencies are calculated based on forced q's. This should have no
      ! influence on conservation properties, it scales linear through all
      ! tendencies. Use iopt=0 to test for each tracer seperately, iopt=1
      ! for a more severe limitation...
      thresh = epsilc
      !thresh=1.e-30
      if (iopt .eq. 0) then
         do nv = 1, num_chem
            do i = its, itf
               if (pret(i) .le. 0 .and. name .eq. 'deep') cycle
               tracermin = q(i, kts, nv)
               tracermax = q(i, kts, nv)
               do k = kts + 1, kte - 1
                  tracermin = min(tracermin, q(i, k, nv))
                  tracermax = max(tracermax, q(i, k, nv))
               enddo
               tracermin = max(tracermin, thresh)
               qmemf = 1.
               ! first check for minimum restriction
               do k = kts, ktop(i)
               ! tracer tendency
                  qmem = outq(i, k, nv)
                  ! only necessary if there is a tendency
                  if (qmem .lt. 0.) then
                     qtest = q(i, k, nv) + outq(i, k, nv) * dt
                     if (qtest .lt. tracermin) then
                        ! qmem2 would be the maximum allowable tendency
                        qmem1 = outq(i, k, nv)
                        qmem2 = (tracermin - q(i, k, nv)) / dt
                        qmemf = min(qmemf, qmem2 / qmem1)
                        if (qmemf .gt. 1.) then
                           write(message, "(a,x,2F9.3)") 'something wrong in negct_1', qmem2, qmem1
                           call chemistry_message(trim(message))
                        endif
                        if (i .eq. ipr .and. j .eq. jpr .and. nv .eq. npr) then
                           write(message, "(I3,x,4F9.3)") k, qtest, qmem2, qmem1, qmemf
                           call chemistry_message(trim(message))
                        endif
                        qmemf = max(qmemf, 0.)
                     endif
                  endif
               enddo
               do k = kts, ktop(i)
                  outq(i, k, nv) = outq(i, k, nv) * qmemf
               enddo
               ! now check max
               qmemf = 1.
               do k = kts, ktop(i)
                  ! tracer tendency
                  qmem = outq(i, k, nv)
                  ! only necessary if there is a tendency
                  if (qmem .gt. 0.) then
                     qtest = q(i, k, nv) + outq(i, k, nv) * dt
                     if (qtest .gt. tracermax) then
                        ! qmem2 would be the maximum allowable tendency
                        qmem1 = outq(i, k, nv)
                        qmem2 = (tracermax - q(i, k, nv)) / dt
                        qmemf = min(qmemf, qmem2 / qmem1)
                        if (qmemf .gt. 1.) then
                           write(message, "(a,x,2F9.3)") 'something wrong in negct_2', qmem2, qmem1
                           call chemistry_error_fatal(trim(message))
                        endif
                        if (i .eq. ipr .and. j .eq. jpr .and. nv .eq. npr) then
                           write(message, "(I3,x,4F9.3)") k, qtest, qmem2, qmem1, qmemf
                           call chemistry_message(trim(message))
                        endif
                        qmemf = max(qmemf, 0.)
                     endif
                  endif
               enddo
               do k = kts, ktop(i)
                  outq(i, k, nv) = outq(i, k, nv) * qmemf
               enddo
            enddo
         enddo
      ! ELSE
      elseif (iopt .eq. 1) then
         do i = its, itf
            qmemf = 1.
            do k = kts, ktop(i)
               do nv = 1, num_chem
                  ! tracer tendency
                  qmem = outq(i, k, nv)
                  ! only necessary if tendency is larger than zero
                  if (qmem .lt. 0.) then
                     qtest = q(i, k, nv) + outq(i, k, nv) * dt
                     if (qtest .lt. thresh) then
                        ! qmem2 would be the maximum allowable tendency
                        qmem1 = outq(i, k, nv)
                        qmem2 = (thresh - q(i, k, nv)) / dt
                        qmemf = min(qmemf, qmem2 / qmem1)
                        qmemf = max(0., qmemf)
                     endif
                  endif
               enddo
            enddo
            do nv = 1, num_chem
               do k = kts, ktop(i)
                  outq(i, k, nv) = outq(i, k, nv) * qmemf
               enddo
            enddo
         enddo
      endif
   end subroutine neg_check_ct
!-------------------------------------------------------
END MODULE module_ctrans_general

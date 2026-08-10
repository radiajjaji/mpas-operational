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

module module_cu_gd

   implicit none

contains

   subroutine cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
              pwev,edtmax,edtmin,maxens2,edtc,               &
              itf,     jtf,     ktf,                         &
              its,ite, jts,jte, kts,kte                     )

      integer, intent(in) :: itf, jtf, ktf
      integer, intent(in) :: its, ite, jts, jte, kts, kte
      integer, intent(in) :: maxens2
       ! ierr error value, maybe modified in this routine
      real, dimension (its:ite, kts:kte), intent(in) :: us, vs, z, p
      real, dimension (its:ite, 1:maxens2), intent(out) :: edtc
      real, dimension (its:ite), intent(out) :: edt
      real, dimension (its:ite), intent(in) :: pwav, pwev
      real, intent(in) :: edtmax, edtmin
      integer, dimension (its:ite), intent(in) :: ktop, kbcon
      integer, dimension (its:ite), intent(inout) :: ierr

      !  local variables in this routine
      integer i, k, kk
      real    einc, pef, pefb, prezk, zkbc
      real, dimension (its:ite) :: vshear, sdp, vws
      !--- determine downdraft strength in terms of windshear
      ! */ calculate an average wind shear over the depth of the cloud
      do i = its, itf
         edt(i) = 0.
         vws(i) = 0.
         sdp(i) = 0.
         vshear(i) = 0.
      enddo
      do kk = kts, ktf - 1
         do i = its, itf
            if (ierr(i) .ne. 0) cycle
            if (kk .le. min0(ktop(i), ktf) .and. kk .ge. kbcon(i)) then
               vws(i) = vws(i) + &
                  (abs((us(i,kk+1)-us(i,kk))/(z(i,kk+1)-z(i,kk))) &
                  +   abs((vs(i,kk+1)-vs(i,kk))/(z(i,kk+1)-z(i,kk)))) * &
                  (p(i,kk) - p(i,kk+1))
               sdp(i) = sdp(i) + p(i,kk) - p(i,kk+1)
            endif
            if (kk .eq. ktf-1)vshear(i) = 1.e3 * vws(i) / sdp(i)
         enddo
      enddo
      do i = its, itf
         if (ierr(i) .eq. 0) then
            pef=(1.591-.639*vshear(i)+.0953*(vshear(i)**2) &
               -.00496*(vshear(i)**3))
            if (pef .gt. edtmax) pef = edtmax
            if (pef .lt. edtmin) pef = edtmin
!
!--- cloud base precip efficiency
!
            zkbc=z(i,kbcon(i))*3.281e-3
            prezk = .02
            if (zkbc .gt. 3.) then
               prezk=.96729352+zkbc*(-.70034167+zkbc*(.162179896+zkbc &
                  *(- 1.2569798e-2+zkbc*(4.2772e-4-zkbc*5.44e-6))))
            endif
            if (zkbc .gt. 25) then
               prezk = 2.4
            endif
            pefb = 1. / (1. + prezk)
            if (pefb .gt. edtmax) pefb = edtmax
            if (pefb .lt. edtmin) pefb = edtmin
            edt(i) = 1. - .5 * (pefb + pef)
            !--- edt here is 1-precipeff!
            !           einc=(1.-edt(i))/float(maxens2)
            !           einc=edt(i)/float(maxens2+1)
            !--- 20 percent
            einc = .2 * edt(i)
            do k = 1, maxens2
               edtc(i, k) = edt(i) + float(k - 2) * einc
            enddo
         endif
      enddo
      do i = its, itf
         if (ierr(i) .eq. 0) then
            do k = 1, maxens2
               edtc(i, k) = -edtc(i, k) * pwav(i) / pwev(i)
               if (edtc(i, k) .gt. edtmax)edtc(i, k) = edtmax
               if (edtc(i, k) .lt. edtmin)edtc(i, k) = edtmin
            enddo
         endif
      enddo

   end subroutine cup_dd_edt

   subroutine cup_dd_he(hes_cup,zd,hcd,z_cup,cdd,entr, &
         jmin,ierr,he,dby,he_cup,                 &
         itf,     jtf,     ktf,                   &
         its,ite, jts,jte, kts,kte                )

      integer, intent (in) :: itf, jtf, ktf
      integer, intent (in) :: its, ite, jts, jte, kts, kte
      ! hcd = downdraft moist static energy
      ! he = moist static energy on model levels
      ! he_cup = moist static energy on model cloud levels
      ! hes_cup = saturation moist static energy on model cloud levels
      ! dby = buoancy term
      ! cdd = detrainment function
      ! z_cup = heights of model cloud levels
      ! entr = entrainment rate
      ! zd = downdraft normalized mass flux
      !
      real, dimension (its:ite, kts:kte), intent (in) :: he, he_cup, hes_cup, z_cup, cdd, zd
      ! entr = entrainment rate
      real, intent (in) :: entr
      integer, dimension (its:ite), intent (in) :: jmin
      ! input and output
      ! ierr error value, maybe modified in this routine
      integer, dimension (its:ite), intent (inout) :: ierr
      real, dimension (its:ite, kts:kte), intent (out) :: hcd, dby
      !  local variables in this routine
      integer :: i, k, ki
      real :: dz

      do k = kts + 1, ktf
         do i = its, itf
            dby(i, k) = 0.
            if (ierr(i) .eq. 0) then
               hcd(i, k) = hes_cup(i, k)
            endif
         enddo
      enddo

      do i = its, itf
         if (ierr(i) .eq. 0) then
            k = jmin(i)
            hcd(i, k) = hes_cup(i, k)
            dby(i, k) = hcd(i, jmin(i)) - hes_cup(i, k)
!
            do ki = jmin(i) - 1, 1, - 1
               dz = z_cup(i, ki + 1) - z_cup(i, ki)
               hcd(i,ki)=(hcd(i,ki+1)*(1.-.5*cdd(i,ki)*dz) &
                  +entr*dz*he(i,ki) &
                  )/(1.+entr*dz-.5*cdd(i,ki)*dz)
               dby(i, ki) = hcd(i, ki) - hes_cup(i, ki)
            enddo
         endif
      enddo
   end subroutine cup_dd_he

   subroutine cup_dd_moisture(zd,hcd,hes_cup,qcd,qes_cup, &
              pwd,q_cup,z_cup,cdd,entr,jmin,ierr,         &
              gamma_cup,pwev,bu,qrcd,                     &
              q,he,t_cup,iloop,xl,                        &
              itf,     jtf,     ktf,                      &
              its,ite, jts,jte, kts,kte                   )

      integer, intent (in) :: itf, jtf, ktf
      integer, intent (in) :: its, ite, jts, jte, kts, kte
      ! cdd = detrainment function
      ! q = environmental q on model levels
      ! q_cup = environmental q on model cloud levels
      ! qes_cup = saturation q on model cloud levels
      ! hes_cup = saturation h on model cloud levels
      ! hcd = h in model cloud
      ! bu = buoancy term
      ! zd = normalized downdraft mass flux
      ! gamma_cup = gamma on model cloud levels
      ! mentr_rate = entrainment rate
      ! qcd = cloud q (including liquid water) after entrainment
      ! qrch = saturation q in cloud
      ! pwd = evaporate at that level
      ! pwev = total normalized integrated evaoprate (i2)
      ! entr = entrainment rate
      real, dimension(its:ite, kts:kte), intent(in) :: zd, t_cup, hes_cup
      real, dimension(its:ite, kts:kte), intent(in) :: hcd, qes_cup, q_cup
      real, dimension(its:ite, kts:kte), intent(in) :: z_cup, cdd, gamma_cup
      real, dimension(its:ite, kts:kte), intent(in) :: q, he
      real, intent (in) :: entr, xl
      integer, intent (in) :: iloop
      integer, dimension (its:ite), intent (in) :: jmin
      integer, dimension (its:ite), intent (inout) :: ierr
      real, dimension (its:ite, kts:kte), intent (out) :: qcd, qrcd, pwd
      real, dimension (its:ite), intent (out) :: pwev, bu
      !  local variables in this routine
      integer :: i, k, ki
      real :: dh, dz, dqeva

      do i = its, itf
         bu(i) = 0.
         pwev(i) = 0.
      enddo
      do k = kts, ktf
         do i = its, itf
            qcd(i, k) = 0.
            qrcd(i, k) = 0.
            pwd(i, k) = 0.
         enddo
      enddo

      do i = its, itf
         if (ierr(i) .eq. 0) then
            k = jmin(i)
            dz = z_cup(i, k + 1) - z_cup(i, k)
            qcd(i, k) = q_cup(i, k)
            ! qcd(i, k) = .5 * (qes_cup(i, k) + q_cup(i, k))
            qrcd(i, k) = qes_cup(i, k)
            pwd(i, jmin(i)) = min(0., qcd(i, k) - qrcd(i, k))
            pwev(i) = pwev(i) + pwd(i, jmin(i))
            qcd(i, k) = qes_cup(i, k)
            dh = hcd(i, k) - hes_cup(i, k)
            bu(i) = dz * dh
            do ki = jmin(i) - 1, 1, - 1
               dz = z_cup(i, ki + 1) - z_cup(i, ki)
               qcd(i,ki)=(qcd(i,ki+1)*(1.-.5*cdd(i,ki)*dz) &
                  +entr*dz*q(i,ki) &
                  )/(1.+entr*dz-.5*cdd(i,ki)*dz)
               !--- to be negatively buoyant, hcd should be smaller than hes!
               dh = hcd(i, ki) - hes_cup(i, ki)
               bu(i) = bu(i) + dz * dh
               qrcd(i,ki)=qes_cup(i,ki)+(1./xl)*(gamma_cup(i,ki) &
                  /(1.+gamma_cup(i,ki)))*dh
               dqeva = qcd(i, ki) - qrcd(i, ki)
               if (dqeva .gt. 0.) dqeva = 0.
               pwd(i, ki) = zd(i, ki) * dqeva
               qcd(i, ki) = qrcd(i, ki)
               pwev(i) = pwev(i) + pwd(i, ki)
!        if(iloop.eq.1.and.i.eq.102.and.j.eq.62)then
!         print *,'in cup_dd_moi ', hcd(i,ki),hes_cup(i,ki),dh,dqeva
!        endif
            enddo
            if (pwev(i) .eq. 0 .and. iloop .eq. 1) then
               ierr(i) = 7
            endif
            if (bu(i) .ge. 0 .and. iloop .eq. 1) then
               ierr(i) = 7
            endif
         endif
      enddo

   end subroutine cup_dd_moisture

   subroutine cup_dd_nms(zd,z_cup,cdd,entr,jmin,ierr, &
              itest,kdet,z1,                          &
              itf,     jtf,     ktf,                  &
              its,ite, jts,jte, kts,kte               )

      integer, intent (in) :: itf, jtf, ktf
      integer, intent (in) :: its, ite, jts, jte, kts, kte
      ! z_cup = height of cloud model level
      ! z1 = terrain elevation
      ! entr = downdraft entrainment rate
      ! jmin = downdraft originating level
      ! kdet = level above ground where downdraft start detraining
      ! itest = flag to whether to calculate cdd

      real, dimension (its:ite, kts:kte), intent (in) :: z_cup
      real, dimension (its:ite), intent (in) :: z1
      real, intent (in) :: entr
      integer, dimension (its:ite), intent (in) :: jmin, kdet
      integer, intent (in) :: itest
      ! ierr error value, maybe modified in this routine
      integer, dimension (its:ite), intent (inout) :: ierr
      ! zd is the normalized downdraft mass flux
      ! cdd is the downdraft detrainmen function

      real, dimension (its:ite, kts:kte), intent (out) :: zd
      real, dimension (its:ite, kts:kte), intent (inout) :: cdd
      !  local variables in this routine
      integer :: i, k, ki
      real :: a, perc, dz
      !--- perc is the percentage of mass left when hitting the ground
      perc = .03
      do k = kts, ktf
         do i = its, itf
            zd(i, k) = 0.
            if (itest .eq. 0)cdd(i, k) = 0.
         enddo
      enddo
      a = 1. - perc

      do i = its, itf
         if (ierr(i) .eq. 0) then
            zd(i, jmin(i)) = 1.
            !--- integrate downward, specify detrainment(cdd)!
            do ki = jmin(i) - 1, 1, - 1
               dz = z_cup(i, ki + 1) - z_cup(i, ki)
               if (ki .le. kdet(i) .and. itest .eq. 0) then
                  cdd(i,ki)=entr+(1.- (a*(z_cup(i,ki)-z1(i)) &
                     +perc*(z_cup(i,kdet(i))-z1(i)) ) &
                     /(a*(z_cup(i,ki+1)-z1(i)) &
                     +perc*(z_cup(i,kdet(i))-z1(i))))/dz
               endif
               zd(i,ki)=zd(i,ki+1)*(1.+(entr-cdd(i,ki))*dz)
            enddo
         endif
      enddo

   end subroutine cup_dd_nms

   subroutine cup_env(z,qes,he,hes,t,q,p,z1, &
         psur,ierr,tcrit,itest,xl,cp,        &
         itf,     jtf,     ktf,              &
         its,ite, jts,jte, kts,kte           )

      integer, intent (in) :: itf, jtf, ktf
      integer, intent (in) :: its, ite, jts, jte, kts, kte
      !
      ! ierr error value, maybe modified in this routine
      ! q = environmental mixing ratio
      ! qes = environmental saturation mixing ratio
      ! t = environmental temp
      ! tv = environmental virtual temp
      ! p = environmental pressure
      ! z = environmental heights
      ! he = environmental moist static energy
      ! hes = environmental saturation moist static energy
      ! psur = surface pressure
      ! z1 = terrain elevation
      !
      real, dimension (its:ite, kts:kte), intent (in) :: p, t
      real, dimension (its:ite, kts:kte), intent (out) :: he, hes, qes
      real, dimension (its:ite, kts:kte), intent (inout) :: z, q
      real, dimension (its:ite), intent (in) :: psur, z1
      real, intent (in) :: xl, cp
      integer, dimension (its:ite), intent (inout) :: ierr
      integer, intent (in) :: itest

      !  local variables in this routine
      integer :: i, k, iph
      real, dimension (1:2) :: ae,be,ht
      real, dimension (its:ite,kts:kte) :: tv
      real :: tcrit,e,tvbar

      ht(1) = xl / cp
      ht(2) = 2.834e6 / cp
      be(1) = .622 * ht(1) / .286
      ae(1) = be(1) / 273. + alog(610.71)
      be(2) = .622 * ht(2) / .286
      ae(2) = be(2) / 273. + alog(610.71)
      do k = kts, ktf
         do i = its, itf
            if (ierr(i) .eq. 0) then
               !csgb - iph is for phase, dependent on tcrit (water or ice)
               iph = 1
               if (t(i, k) .le. tcrit) iph = 2
               e = exp(ae(iph) - be(iph) / t(i, k))
               qes(i, k) = .622 * e / (100. * p(i, k) - e)
               if (qes(i, k) .le. 1.e-08)qes(i, k) = 1.e-08
               if (q(i, k) .gt. qes(i, k))q(i, k) = qes(i, k)
               tv(i, k) = t(i, k) + .608 * q(i, k) * t(i, k)
            endif
         enddo
      enddo
      !--- z's are calculated with changed h's and q's and t's
      !--- if itest=2
      if (itest .ne. 2) then
         do i = its, itf
            if (ierr(i) .eq. 0) then
               z(i, 1) = max(0., z1(i)) - (alog(p(i, 1)) - &
                  alog(psur(i))) * 287. * tv(i, 1) / 9.81
            endif
         enddo

         ! --- calculate heights
         do k = kts + 1, ktf
            do i = its, itf
               if (ierr(i) .eq. 0) then
                  tvbar = .5 * tv(i, k) + .5 * tv(i, k - 1)
                  z(i, k) = z(i, k - 1) - (alog(p(i, k)) - &
                     alog(p(i, k - 1))) * 287. * tvbar / 9.81
               endif
            enddo
         enddo
      else
         do k = kts, ktf
            do i = its, itf
               if (ierr(i) .eq. 0) then
                  z(i,k)=(he(i,k)-1004.*t(i,k)-2.5e6*q(i,k))/9.81
                  z(i, k) = max(1.e-3, z(i, k))
               endif
            enddo
         enddo
      endif
      !--- calculate moist static energy - he
      !    saturated moist static energy - hes
      do k = kts, ktf
         do i = its, itf
            if (ierr(i) .eq. 0) then
               if(itest .eq. 0) he(i,k)=9.81*z(i,k)+1004.*t(i,k)+2.5e06*q(i,k)
               hes(i,k)=9.81*z(i,k)+1004.*t(i,k)+2.5e06*qes(i,k)
               if (he(i, k) .ge. hes(i, k))he(i, k) = hes(i, k)
            endif
         enddo
      enddo

   end subroutine cup_env

   subroutine cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,   &
              he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
              ierr,z1,xl,rv,cp,                                &
              itf,     jtf,     ktf,                           &
              its,ite, jts,jte, kts,kte                        )

      integer, intent (in) :: itf, jtf, ktf
      integer, intent (in) :: its, ite, jts, jte, kts, kte
      !
      ! ierr error value, maybe modified in this routine
      ! q = environmental mixing ratio
      ! q_cup = environmental mixing ratio on cloud levels
      ! qes = environmental saturation mixing ratio
      ! qes_cup = environmental saturation mixing ratio on cloud levels
      ! t = environmental temp
      ! t_cup = environmental temp on cloud levels
      ! p = environmental pressure
      ! p_cup = environmental pressure on cloud levels
      ! z = environmental heights
      ! z_cup = environmental heights on cloud levels
      ! he = environmental moist static energy
      ! he_cup = environmental moist static energy on cloud levels
      ! hes = environmental saturation moist static energy
      ! hes_cup = environmental saturation moist static energy on cloud levels
      ! gamma_cup = gamma on cloud levels
      ! psur = surface pressure
      ! z1 = terrain elevation
      !
      real, dimension(its:ite, kts:kte), intent(in) :: qes, q, he
      real, dimension(its:ite, kts:kte), intent(in) :: hes, z, p, t
      real, dimension(its:ite, kts:kte), intent(out) :: qes_cup, q_cup
      real, dimension(its:ite, kts:kte), intent(out) :: he_cup, hes_cup
      real, dimension(its:ite, kts:kte), intent(out) :: z_cup, p_cup
      real, dimension(its:ite, kts:kte), intent(out) :: gamma_cup, t_cup
      real, dimension(its:ite), intent(in) :: psur, z1
      real, intent (in) :: xl, rv, cp
      integer, dimension (its:ite), intent (inout) :: ierr
      !  local variables in this routine
      integer :: i, k

      do k = kts + 1, ktf
         do i = its, itf
            if (ierr(i) .eq. 0) then
               qes_cup(i,k)=.5*(qes(i,k-1)+qes(i,k))
               q_cup(i,k)=.5*(q(i,k-1)+q(i,k))
               hes_cup(i,k)=.5*(hes(i,k-1)+hes(i,k))
               he_cup(i,k)=.5*(he(i,k-1)+he(i,k))
               if (he_cup(i, k) .gt. hes_cup(i, k))he_cup(i, k) = hes_cup(i, k)
               z_cup(i,k)=.5*(z(i,k-1)+z(i,k))
               p_cup(i,k)=.5*(p(i,k-1)+p(i,k))
               t_cup(i,k)=.5*(t(i,k-1)+t(i,k))
               gamma_cup(i,k)=(xl/cp)*(xl/(rv*t_cup(i,k) &
                  *t_cup(i,k)))*qes_cup(i,k)
            endif
         enddo
      enddo
      do i=its,itf
         if(ierr(i).eq.0)then
            qes_cup(i,1)=qes(i,1)
            q_cup(i,1)=q(i,1)
            hes_cup(i,1)=hes(i,1)
            he_cup(i,1)=he(i,1)
            z_cup(i,1)=.5*(z(i,1)+z1(i))
            p_cup(i,1)=.5*(p(i,1)+psur(i))
            t_cup(i,1)=t(i,1)
            gamma_cup(i,1)=xl/cp*(xl/(rv*t_cup(i,1) &
               *t_cup(i,1)))*qes_cup(i,1)
         endif
      enddo

   end subroutine cup_env_clev

   subroutine cup_kbcon(cap_inc,iloop,k22,kbcon,he_cup,hes_cup, &
                        ierr, kbmax, p_cup, cap_max, &
                        itf, jtf, ktf, &
                        its, ite, jts, jte, kts, kte)


      integer, intent (in) :: itf, jtf, ktf
      integer, intent (in) :: its, ite, jts, jte, kts, kte
      real, dimension (its:ite, kts:kte), intent (in) :: he_cup, hes_cup, p_cup
      real, dimension (its:ite), intent (in) :: cap_max, cap_inc
      integer, dimension (its:ite), intent (in) :: kbmax
      integer, dimension (its:ite), intent (inout) :: kbcon, k22, ierr
      integer, intent (in) :: iloop

      !local variables
      integer :: i
      real :: pbcdif, plus
      ! --- determine the level of convective cloud base  - kbcon
      do i = its, itf
         kbcon(i)=1
         if(ierr(i).ne.0)cycle
         kbcon(i)=k22(i)
         go to 32
31 continue
         kbcon(i) = kbcon(i) + 1
         if (kbcon(i) .gt. kbmax(i) + 2) then
            if (iloop .lt. 4)ierr(i) = 3
!        if (iloop .lt. 4)ierr(i) = 997
            cycle
         endif
32 continue
         if (he_cup(i, k22(i)) .lt. hes_cup(i, kbcon(i))) go to 31

         ! cloud base pressure and max moist static energy pressure
         ! i.e., the depth (in mb) of the layer of negative buoyancy
         if (kbcon(i) - k22(i) .eq. 1) cycle
         pbcdif = - p_cup(i, kbcon(i)) + p_cup(i, k22(i))
         plus = max(25., cap_max(i) - float(iloop - 1) * cap_inc(i))
         if (iloop .eq. 4) plus = cap_max(i)
         if (pbcdif .gt. plus) then
            k22(i) = k22(i) + 1
            kbcon(i) = k22(i)
            go to 32
         endif
      enddo

   end subroutine cup_kbcon

   subroutine cup_ktop(ilo,dby,kbcon,ktop,ierr, &
         itf, jtf, ktf,                    &
         its,ite, jts,jte, kts,kte         )

      integer, intent (in) :: itf, jtf, ktf
      integer, intent (in) :: its, ite, jts, jte, kts, kte
      ! dby = buoancy term
      ! ktop = cloud top (output)
      ! ilo = flag
      ! ierr error value, maybe modified in this routine
      real, dimension (its:ite, kts:kte), intent (inout) :: dby
      integer, dimension (its:ite), intent (in) :: kbcon
      integer, intent (in) :: ilo
      integer, dimension (its:ite), intent (out) :: ktop
      integer, dimension (its:ite), intent (inout) :: ierr
      !  local variables
      integer :: i, k

      do 42 i = its, itf
         ktop(i) = 1
         if (ierr(i) .eq. 0) then
            do 40 k = kbcon(i) + 1, ktf - 1
               if (dby(i, k) .le. 0.) then
                  ktop(i) = k - 1
                  go to 41
               endif
            40 continue

               if (ilo .eq. 1) ierr(i) = 5
               !if (ilo .eq. 2)ierr(i) = 998
               go to 42
41 continue
            do k = ktop(i) + 1, ktf
               dby(i, k) = 0.
            enddo
         endif
      42 continue

   end subroutine cup_ktop

   subroutine cup_maximi(array,ks,ke,maxx,ierr, &
              itf,     jtf,     ktf,            &
              its,ite, jts,jte, kts,kte         )

      integer, intent (in) :: itf, jtf, ktf
      integer, intent (in) :: its, ite, jts, jte, kts, kte
      ! array input array
      ! x output array with return values
      ! kt output array of levels
      ! ks, kend  check - range
      real, dimension (its:ite, kts:kte), intent (in) :: array
      integer, dimension (its:ite), intent (in) :: ierr, ke
      integer, intent (in) :: ks
      integer, dimension (its:ite), intent (out) :: maxx
      real, dimension (its:ite) :: x
      real :: xar
      integer :: i, k

      do i = its, itf
         maxx(i) = ks
         if (ierr(i) .eq. 0) then
            x(i) = array(i, ks)

            do k = ks, ke(i)
               xar = array(i, k)
               if (xar .ge. x(i)) then
                  x(i) = xar
                  maxx(i) = k
               endif
            enddo
         endif
      enddo

   end subroutine cup_maximi

   subroutine cup_minimi(array,ks,kend,kt,ierr, &
              itf,     jtf,     ktf,            &
              its,ite, jts,jte, kts,kte         )

      integer, intent (in) :: itf, jtf, ktf
      integer, intent (in) :: its, ite, jts, jte, kts, kte
      ! array input array
      ! x output array with return values
      ! kt output array of levels
      ! ks, kend  check - range
      real, dimension (its:ite, kts:kte), intent (in) :: array
      integer, dimension (its:ite), intent (in) :: ierr, ks, kend
      integer, dimension (its:ite), intent (out) :: kt
      real, dimension (its:ite) :: x
      integer :: i, k, kstop

      do i = its, itf
         kt(i) = ks(i)
         if (ierr(i) .eq. 0) then
            x(i) = array(i, ks(i))
            kstop = max(ks(i) + 1, kend(i))

            do k = ks(i) + 1, kstop
               if (array(i, k) .lt. x(i)) then
                  x(i) = array(i, k)
                  kt(i) = k
               endif
            enddo
         endif
      enddo

   end subroutine cup_minimi

   subroutine cup_up_he(k22,hkb,z_cup,cd,entr,he_cup,hc, &
              kbcon,ierr,dby,he,hes_cup,                 &
              itf,     jtf,     ktf,                     &
              its,ite, jts,jte, kts,kte                  )


      integer, intent (in) :: itf, jtf, ktf
      integer, intent (in) :: its, ite, jts, jte, kts, kte
      ! hc = cloud moist static energy
      ! hkb = moist static energy at originating level
      ! he = moist static energy on model levels
      ! he_cup = moist static energy on model cloud levels
      ! hes_cup = saturation moist static energy on model cloud levels
      ! dby = buoancy term
      ! cd = detrainment function
      ! z_cup = heights of model cloud levels
      ! entr = entrainment rate
      real, dimension (its:ite, kts:kte), intent (in) :: he, he_cup
      real, dimension (its:ite, kts:kte), intent (in) :: hes_cup, z_cup, cd
  ! entr = entrainment rate
      real, intent (in) :: entr
      integer, dimension (its:ite), intent (in) :: kbcon, k22
      integer, dimension (its:ite), intent (inout) :: ierr
      real, dimension (its:ite, kts:kte), intent (out) :: hc, dby
      real, dimension (its:ite), intent (out) :: hkb

      ! local variables
      integer :: i, k
      real :: dz

      !--- moist static energy inside cloud
      do i = its, itf
         if (ierr(i) .eq. 0.) then
            hkb(i) = he_cup(i, k22(i))
            do k = 1, k22(i)
               hc(i, k) = he_cup(i, k)
               dby(i, k) = 0.
            enddo
            do k = k22(i), kbcon(i) - 1
               hc(i, k) = hkb(i)
               dby(i, k) = 0.
            enddo
            k = kbcon(i)
            hc(i, k) = hkb(i)
            dby(i, kbcon(i)) = hkb(i) - hes_cup(i, k)
         endif
      enddo
      do k = kts + 1, ktf
         do i = its, itf
            if (k .gt. kbcon(i) .and. ierr(i) .eq. 0.) then
               dz = z_cup(i, k) - z_cup(i, k - 1)
               hc(i,k)=(hc(i,k-1)*(1.-.5*cd(i,k)*dz)+entr* &
                  dz*he(i,k-1))/(1.+entr*dz-.5*cd(i,k)*dz)
               dby(i, k) = hc(i, k) - hes_cup(i, k)
            endif
         enddo
      enddo

   end subroutine cup_up_he

   subroutine cup_up_moisture(ierr,z_cup,qc,qrc,pw,pwav,    &
              kbcon,ktop,cd,dby,mentr_rate,clw_all,         &
              q,gamma_cup,zu,qes_cup,k22,qe_cup,xl,         &
              itf,     jtf,     ktf,                        &
              its,ite, jts,jte, kts,kte                     )

      integer, intent (in) :: itf, jtf, ktf
      integer, intent (in) :: its, ite, jts, jte, kts, kte
      ! cd = detrainment function
      ! q = environmental q on model levels
      ! qe_cup = environmental q on model cloud levels
      ! qes_cup = saturation q on model cloud levels
      ! dby = buoancy term
      ! cd = detrainment function
      ! zu = normalized updraft mass flux
      ! gamma_cup = gamma on model cloud levels
      ! mentr_rate = entrainment rate
      real, dimension (its:ite, kts:kte), intent (in) :: q, zu
      real, dimension (its:ite, kts:kte), intent (in) ::gamma_cup, qe_cup
      real, dimension (its:ite, kts:kte), intent (in) ::dby, qes_cup
      real, dimension (its:ite, kts:kte), intent (in) ::z_cup, cd
      ! entr = entrainment rate
      real, intent (in) :: mentr_rate, xl
      integer, dimension (its:ite), intent (in) :: kbcon, ktop, k22
      integer, dimension (its:ite), intent (inout) :: ierr
      ! qc = cloud q (including liquid water) after entrainment
      ! qrch = saturation q in cloud
      ! qrc = liquid water content in cloud after rainout
      ! pw = condensate that will fall out at that level
      ! pwav = totan normalized integrated condensate (i1)
      ! c0 = conversion rate (cloud to rain)
      real, dimension (its:ite, kts:kte), intent (out) :: qc, qrc
      real, dimension (its:ite, kts:kte), intent (out) :: pw, clw_all
      real, dimension (its:ite), intent (out) :: pwav

      ! local variables
      integer :: iall, i, k
      real :: dh, qrch, c0, dz, radius

      iall = 0
      c0 = .002
      !--- no precip for small clouds
      if (mentr_rate .gt. 0.) then
         radius = .2 / mentr_rate
         if (radius .lt. 900.) c0 = 0.
!         if (radius .lt. 900.) iall = 0
      endif
      do i = its, itf
         pwav(i) = 0.
      enddo
      do k = kts, ktf
         do i = its, itf
            pw(i, k) = 0.
            if (ierr(i) .eq. 0) qc(i, k) = qes_cup(i, k)
            clw_all(i, k) = 0.
            qrc(i, k) = 0.
         enddo
      enddo
      do i = its, itf
         if (ierr(i) .eq. 0.) then
            do k = k22(i), kbcon(i) - 1
               qc(i, k) = qe_cup(i, k22(i))
            enddo
         endif
      enddo

      do k = kts + 1, ktf
         do i = its, itf
            if (ierr(i) .ne. 0) cycle
            if (k .lt. kbcon(i)) cycle
            if (k .gt. ktop(i)) cycle
            dz = z_cup(i, k) - z_cup(i, k - 1)
!------    1. steady state plume equation, for what could
!------       be in cloud without condensation
            qc(i,k)=(qc(i,k-1)*(1.-.5*cd(i,k)*dz)+mentr_rate* &
               dz*q(i,k-1))/(1.+mentr_rate*dz-.5*cd(i,k)*dz)
!
!--- saturation  in cloud, this is what is allowed to be in it
!
            qrch=qes_cup(i,k)+(1./xl)*(gamma_cup(i,k) &
               /(1.+gamma_cup(i,k)))*dby(i,k)
!
!------- liquid water content in cloud after rainout
!
            clw_all(i, k) = qc(i, k) - qrch
            qrc(i,k)=(qc(i,k)-qrch)/(1.+c0*dz*zu(i,k))
            if (qrc(i, k) .lt. 0.) then
               qrc(i, k) = 0.
            endif
!
!-------   3.condensation
!
            pw(i, k) = c0 * dz * qrc(i, k) * zu(i, k)
            if (iall .eq. 1) then
               qrc(i, k) = 0.
               pw(i, k) = (qc(i, k) - qrch) * zu(i, k)
               if (pw(i, k) .lt. 0.)pw(i, k) = 0.
            endif
!
!----- set next level
!
            qc(i, k) = qrc(i, k) + qrch
!
!--- integrated normalized ondensate
!
            pwav(i) = pwav(i) + pw(i, k)
         enddo
      enddo

   end subroutine cup_up_moisture

   subroutine cup_up_nms(zu, z_cup, entr, cd, kbcon, ktop, ierr, k22, &
              itf, jtf, ktf, &
              its, ite, jts, jte, kts, kte)

      IMPLICIT NONE
      ! on input
      integer, intent (in) :: itf, jtf, ktf
      integer, intent (in) :: its, ite, jts, jte, kts, kte
      ! cd = detrainment function
      real, dimension (its:ite, kts:kte), intent (in) :: z_cup, cd
      ! entr = entrainment rate
      real, intent (in) :: entr
      integer, dimension (its:ite), intent (in) :: kbcon, ktop, k22
      ! input and output
      ! ierr error value, maybe modified in this routine
      integer, dimension (its:ite), intent (inout) :: ierr
      ! zu is the normalized mass flux
      real, dimension (its:ite, kts:kte), intent (out) :: zu
      ! local variables in this routine

      integer :: i, k
      real :: dz
      ! initialize for this go around
      do k = kts, ktf
          do i = its, itf
              zu(i, k) = 0.
          enddo
      enddo
      ! do normalized mass budget
      do i = its, itf
          if (ierr(i) .eq. 0) then
              do k = k22(i), kbcon(i)
                  zu(i, k) = 1.
              end do
              do k = kbcon(i) + 1, ktop(i)
                  dz = z_cup(i, k) - z_cup(i, k - 1)
                  zu(i, k) = zu(i, k - 1) * (1. + (entr - cd(i, k)) * dz)
              end do
          end if
      end do

   end subroutine cup_up_nms

!-------------------------------------------------------
end module module_cu_gd

!************************************************************************
!
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

 module mpas_atmchem_constants
 use mpas_kind_types
 use mpas_constants, only : pii, cp, gravity, R_d => rgas, rvord

 implicit none
 public
 save

 
!defines the constants needed for chemistry from WRF.
!Jiawang Feng.
!2019-09-26.
!=================================================================================================================

!real(kind=RKIND),parameter :: gravity = 9.80616
 real(kind=RKIND),parameter :: epsilc  = 1.e-16
 real(kind=RKIND),parameter :: r_v     = 461.6 !gas constant for water vapor  [J/kg/K]
!real(kind=RKIND),parameter :: r_d     = 287.0 !gas constant for dry air      [J/kg/K]
 real(kind=RKIND),parameter :: g       = 9.80616
 real(kind=RKIND),parameter :: pi      = pii

 real(kind=RKIND),parameter :: rhowater= 1000.
 real(kind=RKIND),parameter :: rvovrd  = rvord
 real(kind=RKIND),parameter :: mwdry   = 28.966
 real(kind=RKIND),parameter :: ep_2    = r_d/r_v
 real(kind=RKIND),parameter :: xlv     = 2.5e6
 real(kind=RKIND),parameter :: svp1    = 0.6112
 real(kind=RKIND),parameter :: svp2    = 17.67
 real(kind=RKIND),parameter :: svp3    = 29.65
 real(kind=RKIND),parameter :: karman  = 0.4             !Von Karman constant

 contains

!=================================================================================================================
 subroutine atmchem_constants_init()
!=================================================================================================================
!dummy subroutine that does not do anything.

 end subroutine atmchem_constants_init

!=================================================================================================================
 end module mpas_atmchem_constants
!=================================================================================================================

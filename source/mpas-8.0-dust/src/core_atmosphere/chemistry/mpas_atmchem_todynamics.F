!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to pass all chemistry tendencies from chemistry to dynamic modules.
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
!     MPAS-v5.2:src/core_atmosphere/physics/mpas_atmphys_todynamics.F
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

module mpas_atmchem_todynamics
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines
 use mpas_dmpar
 use mpas_atmchem_vars
 use mpas_atmphys_vars

 implicit none
 private
 public:: chemistry_get_tend


!Interface between the chemistry schemes and the non-hydrostatic dynamical core.
!
! subroutines in mpas_atmchem_todynamics:
! ---------------------------------------
! chemistry_get_tend: add and mass-weigh tendencies before being added to dynamics tendencies.
! tend_toEdges      : interpolate wind-tendencies from centers to edges of grid-cells.
! ----------------------------------------

 contains
 

!=================================================================================================================
 subroutine chemistry_get_tend( mesh, state, tend, tend_chemistry )
!=================================================================================================================
 use mpas_atm_dimensions
!input variables:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: tend_chemistry
!inout variables:
 type(mpas_pool_type),intent(inout):: tend

!local variables:
 integer        :: i,iCell,k,n
 integer,pointer:: nCells,nCellsSolve

 real(kind=RKIND),dimension(:,:),  pointer :: mass          ! time level 2 rho_zz
 real(kind=RKIND),dimension(:,:,:),pointer :: tend_chem

!=================================================================================================================
 call mpas_pool_get_dimension(mesh, 'nCells',         nCells        )
 call mpas_pool_get_dimension(mesh, 'nCellsSolve',    nCellsSolve   )
 call mpas_pool_get_array(state,'rho_zz',   mass, 2)
 call mpas_pool_get_array(tend ,'chem_tend',tend_chem)

 tend_chem(:,:,:) = 0._RKIND
 call chemistry_get_tend_work(state,nCells,nCellsSolve, &
                           tend_chemistry, &
                           mass,tend_chem  &
                           )

 end subroutine chemistry_get_tend

 !==================================================================================================
 subroutine chemistry_get_tend_work(state,nCells,nCellsSolve, &
                           tend_chemistry, &
                           mass,tend_chem  &
                           )

!==================================================================================================

    use mpas_atm_dimensions
 
    implicit none
    type(mpas_pool_type),intent(in) :: state
    type(mpas_pool_type),intent(in) :: tend_chemistry
    integer, intent(in) :: nCells, nCellsSolve
    real(kind=RKIND),dimension(nVertLevels,nCells+1), intent(in) :: mass
    real(kind=RKIND),dimension(num_chem,nVertLevels,nCells+1), intent(inout) :: tend_chem
! local variables
    real(kind=RKIND),dimension(:,:,:),pointer :: tend_emis
    integer :: i, k, h

    call mpas_pool_get_array(tend_chemistry,'tend_emis'     ,tend_emis     )

    do h = 1, num_chem
    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_chem(h,k,i) = tend_chem(h,k,i) + tend_emis(h,k,i)*mass(k,i)
    enddo
    enddo
    enddo


 end subroutine chemistry_get_tend_work
 

!=================================================================================================================
 end module mpas_atmchem_todynamics
!=================================================================================================================

!************************************************************************
! This computer software was prepared by University of Science and Technology of
! China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Some chemistry utilities.
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
!     MPAS-v5.2:src/core_atmosphere/physics/mpas_atmphys_utilities.F
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************

module mpas_atmchem_utilities

 use mpas_kind_types
! use mpas_abort, only : mpas_dmpar_global_abort
 use mpas_log,only : mpas_log_write
 use mpas_derived_types, only : MPAS_LOG_ERR, MPAS_LOG_CRIT



 implicit none
 private
 public:: chemistry_error_fatal, &
          chemistry_error_message, &
          chemistry_message

 character(len=StrKIND),public:: mpas_err_message

!jwfeng from AEMOL.
!Simple utility subroutines to broadcast messages or abort chemistry parameterizations.



 contains


!=================================================================================================================
 subroutine chemistry_message(str)
!=================================================================================================================

!input arguments:
 character(len=*),intent(in):: str

!-----------------------------------------------------------------------------------------------------------------

 call mpas_log_write(trim(str))

 end subroutine chemistry_message

!=================================================================================================================
 subroutine chemistry_error_message(str)
!=================================================================================================================

!input arguments:
 character(len=*),intent(in):: str

!--------------------------------------------------------------------------------------------------

 call mpas_log_write('   ',                                                                        messageType=MPAS_LOG_ERR)
 call mpas_log_write('--------------------------MPAS core_chemistry ERROR-----------------------', messageType=MPAS_LOG_ERR)
 call mpas_log_write(trim(str),                                                                    messageType=MPAS_LOG_ERR)

 end subroutine chemistry_error_message


!=================================================================================================================
 subroutine chemistry_error_fatal(str)
!=================================================================================================================

!input arguments:
 character(len=*),intent(in):: str

!--------------------------------------------------------------------------------------------------

!call mpas_dmpar_global_abort(' ',                                                                          deferredAbort=.true.)
!call mpas_dmpar_global_abort('------------------------------ FATAL CALLED ------------------------------', deferredAbort=.true.)
!call mpas_dmpar_global_abort(trim(str),                                                                    deferredAbort=.true.)
!call mpas_dmpar_global_abort('ERROR: MPAS core_chemistry abort')

 call mpas_log_write('   ',                                                                        messageType=MPAS_LOG_ERR)
 call mpas_log_write('------------------------------ FATAL CALLED ------------------------------', messageType=MPAS_LOG_ERR)
 call mpas_log_write(trim(str),                                                                    messageType=MPAS_LOG_ERR)
 call mpas_log_write('MPAS core_chemistry abort', messageType=MPAS_LOG_CRIT)
 
 end subroutine chemistry_error_fatal

!=================================================================================================================
 end module mpas_atmchem_utilities
!=================================================================================================================

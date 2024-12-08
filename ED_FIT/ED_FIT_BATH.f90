MODULE ED_FIT_CHI2
  USE SF_CONSTANTS
  USE SF_OPTIMIZE, only:fmin_cg,fmin_cgminimize
  USE SF_LINALG,   only:eye,zeye,inv,inv_her,trace,operator(.x.) !BLAS xgemm operator overloading
  USE SF_IOTOOLS,  only:reg,free_unit,str
  USE SF_ARRAYS,   only:arange
  USE SF_MISC,     only:assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_BATH
  USE ED_BATH_FUNCTIONS
  USE ED_FIT_REPLICA
  USE ED_FIT_GENERAL


  implicit none
  private

  interface ed_chi2_fitgf
     module procedure chi2_fitgf_normal
     module procedure chi2_fitgf_superc
  end interface ed_chi2_fitgf


  public :: ed_chi2_fitgf


contains

  subroutine chi2_fitgf_normal(fg,bath)
    complex(8),dimension(:,:,:,:,:)    :: fg ![Nspin,Nspin,Nimp,Nimp][Lmats]
    real(8),dimension(:),intent(inout) :: bath
    integer                            :: L
    !
#ifdef _MPI    
    if(check_MPI())call ed_set_MpiComm()
#endif
    !
    select case(cg_method)
    case default;stop "ED Error: cg_method > 1"
    case (0);if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"\Chi2 fit: CG-nr, CG-weight: ",cg_weight," on: ",cg_scheme
    case (1);if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"\Chi2 fit: CG-minimize, CG-weight: ",cg_weight," on: ",cg_scheme
    end select
    !
    if(MpiMaster)then
       !
       L = size(fg,5)
       call assert_shape(fg,[Nspin,Nspin,Nimp,Nimp,L],"chi2_fitgf_normal","fg")
       !
       select case(bath_type)
       case('replica');call chi2_fitgf_replica(fg,bath)
       case('general');call chi2_fitgf_general(fg,bath)
       end select
    end if
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,bath)
       if(.not.MpiMaster)write(LOGfile,"(A)")"Bath received from master node"
    endif
#endif
    !
    !set trim_state_list to true after the first fit has been done: this 
    !marks the ends of the cycle of the 1st DMFT loop.
    trim_state_list=.true.
    !DELETE THE LOCAL MPI COMMUNICATOR:
#ifdef _MPI    
    if(check_MPI())call ed_del_MpiComm()
#endif   
    !
  end subroutine chi2_fitgf_normal



  subroutine chi2_fitgf_superc(fg,ff,bath)
    complex(8),dimension(:,:,:,:,:)    :: fg ![Nspin,Nspin,Nimp,Nimp][Lmats]
    complex(8),dimension(:,:,:,:,:)    :: ff ![Nspin,Nspin,Nimp,Nimp][Lmats]
    real(8),dimension(:),intent(inout) :: bath
    integer                            :: L
    !
#ifdef _MPI    
    if(check_MPI())call ed_set_MpiComm()
#endif
    !
    select case(cg_method)
    case default;stop "ED Error: cg_method > 1"
    case (0);if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"\Chi2 fit: CG-nr, CG-weight: ",cg_weight," on: ",cg_scheme
    case (1);if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"\Chi2 fit: CG-minimize, CG-weight: ",cg_weight," on: ",cg_scheme
    end select
    if(MpiMaster)then
       !
       L = size(fg,5)
       call assert_shape(fg,[Nspin,Nspin,Nimp,Nimp,L],"chi2_fitgf_superc","fg")
       call assert_shape(ff,[Nspin,Nspin,Nimp,Nimp,L],"chi2_fitgf_superc","ff")
       !
       select case(bath_type)
       case('replica');call chi2_fitgf_replica(fg,ff,bath)
       case('general');call chi2_fitgf_general(fg,ff,bath)
       end select
    end if
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,bath)
       if(.not.MpiMaster)write(LOGfile,"(A)")"Bath received from master node"
    endif
#endif
    !
    !set trim_state_list to true after the first fit has been done: this 
    !marks the ends of the cycle of the 1st DMFT loop.
    trim_state_list=.true.
    !DELETE THE LOCAL MPI COMMUNICATOR:
#ifdef _MPI    
    if(check_MPI())call ed_del_MpiComm()
#endif   
  end subroutine chi2_fitgf_superc


end MODULE ED_FIT_CHI2

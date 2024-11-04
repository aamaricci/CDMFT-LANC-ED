MODULE ED_BATH_FIT
!Contains routines that fit the Impurity model bath
  USE SF_CONSTANTS
  USE SF_OPTIMIZE, only:fmin_cg,fmin_cgplus,fmin_cgminimize
  USE SF_LINALG,   only:eye,zeye,inv,inv_her,operator(.x.)
  USE SF_IOTOOLS,  only:reg,free_unit,txtfy
  USE SF_ARRAYS,   only:arange
  USE SF_MISC,     only:assert_shape
  !
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL  
  USE ED_AUX_FUNX
  USE ED_BATH
  USE ED_FIT_COMMON
  USE ED_FIT_NORMAL
  USE ED_FIT_HYBRID
  USE ED_FIT_REPLICA
  USE ED_FIT_GENERAL
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif


  implicit none
  private

  interface ed_chi2_fitgf
  !This subroutine realizes the chisquare fit of the Weiss field or hybridization function via
  !an impurity model non-interacting Green's function. The bath levels (levels/internal structure
  !and hybridization strength) are supplied by the user in the :f:var:`bath` array
  !and are the parameters of the fit.
  !The function(s) to fit can have different shapes:
  !
  !  * [:f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`, :f:var:`lfit` ]  
  !  * [:f:var:`nlat` :math:`\cdot` :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nlat` :math:`\cdot` :f:var:`nspin` 
  !    :math:`\cdot` :f:var:`norb`, :f:var:`lfit`  ]  
  !  * [:f:var:`nlat`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`lfit` ] 
  !  * [:f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :f:var:`lfit` ]
  !  * [:f:var:`nlat`, :f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :f:var:`lfit` ]
  !
  !where :f:var:`nlat` is the number of impurity sites in real-space DMFT. Accordingly, the bath array or arrays have rank 2 or 3.
  !Some global variables directly influence the way the fit is performed and can be modified in the input file. See :f:mod:`ed_input_vars`
  !for the description of :f:var:`lfit`, :f:var:`cg_method` , :f:var:`cg_grad`, :f:var:`cg_ftol`, :f:var:`cg_stop` , :f:var:`cg_niter` ,
  !:f:var:`cg_weight` , :f:var:`cg_scheme` , :f:var:`cg_pow` , :f:var:`cg_minimize_ver` , :f:var:`cg_minimize_hh` .      
  !
     module procedure chi2_fitgf_single_normal_n3
     module procedure chi2_fitgf_single_normal_n5
     module procedure chi2_fitgf_single_superc_n3
     module procedure chi2_fitgf_single_superc_n5
     module procedure chi2_fitgf_lattice_normal_n3
     module procedure chi2_fitgf_lattice_normal_n4
     module procedure chi2_fitgf_lattice_normal_n6
     module procedure chi2_fitgf_lattice_superc_n3
     module procedure chi2_fitgf_lattice_superc_n4
     module procedure chi2_fitgf_lattice_superc_n6
  end interface ed_chi2_fitgf

  public :: ed_chi2_fitgf


contains


  !+----------------------------------------------------------------------+
  !PURPOSE  : Chi^2 fit of the G0/Delta 
  !
  ! - CHI2_FITGF_GENERIC_NORMAL interface for the normal case 
  !   * CHI2_FITGF_GENERIC_NORMAL_NOSPIN interface to fixed spin input
  !+----------------------------------------------------------------------+
  subroutine chi2_fitgf_single_normal_n3(g,bath,ispin,iorb,fmpi)
    complex(8),dimension(:,:,:)                       :: g !normal Weiss field or hybridization function to fit
    real(8),dimension(:)                              :: bath !bath parameters array
    integer,optional                                  :: ispin !spin component to be fitted (default = :code:`1` ). Only used if :f:var:`ed_mode` = :code:`normal` and :f:var:`bath_type` = :code:`normal, hybrid`
    integer,optional                                  :: iorb !orbital to be fitted
    logical,optional                                  :: fmpi !flag to automatically broadcast the fit over the MPI communicator (default = :code:`.true.` )
    !
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: fg
    integer                                           :: ispin_,Liw
    logical                                           :: fmpi_
    !
    write(Logfile,"(A)")""
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG chi2_fitgf_generic_normal_mpi: Start Chi**2 fit"
#endif
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    fmpi_=.true.;if(present(fmpi))fmpi_=fmpi
    !
#ifdef _MPI    
    if(check_MPI().AND.fmpi_)call ed_set_MpiComm()
#endif
    !
    select case(cg_method)
    case default
       stop "ED Error: cg_method > 1"
    case (0)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"Chi^2 fit with CG-nr and CG-weight: ",cg_weight," on: ",cg_scheme
    case (1)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"Chi^2 fit with CG-minimize and CG-weight: ",cg_weight," on: ",cg_scheme
    end select
    !
    if(MpiMaster)then
       !
       call assert_shape(g,[Nspin*Norb,Nspin*Norb,Lmats],"chi2_fitgf_generic_normal","g")
       fg = so2nn_reshape(g(1:Nspin*Norb,1:Nspin*Norb,1:Lmats),Nspin,Norb,Lmats)
       !
       select case(bath_type)
       case default
          select case(ed_mode)
          case ("normal")
             if(present(iorb))then
                call chi2_fitgf_normal_normal(fg(ispin_,ispin_,:,:,:),bath,ispin_,iorb)
             else
                call chi2_fitgf_normal_normal(fg(ispin_,ispin_,:,:,:),bath,ispin_)
             endif
          case ("nonsu2")
             call chi2_fitgf_normal_nonsu2(fg(:,:,:,:,:),bath)
          case default
             stop "chi2_fitgf ERROR: ed_mode!=normal/nonsu2 but only NORMAL component is provided"
          end select
          !
       case ("hybrid")
          select case(ed_mode)
          case ("normal")
             call chi2_fitgf_hybrid_normal(fg(ispin_,ispin_,:,:,:),bath,ispin_)
          case ("nonsu2")
             call chi2_fitgf_hybrid_nonsu2(fg(:,:,:,:,:),bath)
          case default
             stop "chi2_fitgf ERROR: ed_mode!=normal/nonsu2 but only NORMAL component is provided" 
          end select
          !
       case ("replica")
          select case(ed_mode)
          case ("normal","nonsu2")
             call chi2_fitgf_replica(fg,bath)
          case default
             stop "chi2_fitgf ERROR: ed_mode!=normal/nonsu2 but only NORMAL component is provided" 
          end select
       case ("general")
          select case(ed_mode)
          case ("normal","nonsu2")
             call chi2_fitgf_general(fg,bath)
          case default
             stop "chi2_fitgf ERROR: ed_mode!=normal/nonsu2 but only NORMAL component is provided" 
          end select
       end select
    endif
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
    if(check_MPI().AND.fmpi_)call ed_del_MpiComm()
#endif   
#ifdef _DEBUG
    write(Logfile,"(A)")""
#endif
  end subroutine chi2_fitgf_single_normal_n3



  subroutine chi2_fitgf_single_normal_n5(g,bath,ispin,iorb,fmpi)
    complex(8),dimension(:,:,:,:,:)                   :: g
    real(8),dimension(:)                              :: bath    
    integer,optional                                  :: ispin,iorb
    logical,optional                                  :: fmpi
    !
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: fg
    integer                                           :: ispin_,Liw
    logical                                           :: fmpi_
    !
    write(Logfile,"(A)")""
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG chi2_fitgf_generic_normal_mpi: Start Chi**2 fit"
#endif
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    fmpi_=.true.;if(present(fmpi))fmpi_=fmpi
    !
#ifdef _MPI    
    if(check_MPI().AND.fmpi_)call ed_set_MpiComm()
#endif
    !
    select case(cg_method)
    case default
       stop "ED Error: cg_method > 1"
    case (0)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"Chi^2 fit with CG-nr and CG-weight: ",cg_weight," on: ",cg_scheme
    case (1)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"Chi^2 fit with CG-minimize and CG-weight: ",cg_weight," on: ",cg_scheme
    end select
    !
    if(MpiMaster)then
       !
       call assert_shape(g,[Nspin,Nspin,Norb,Norb,Lmats],"chi2_fitgf_generic_normal","g")
       fg = g(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats)
       !
       select case(bath_type)
       case default
          select case(ed_mode)
          case ("normal")
             if(present(iorb))then
                call chi2_fitgf_normal_normal(fg(ispin_,ispin_,:,:,:),bath,ispin_,iorb)
             else
                call chi2_fitgf_normal_normal(fg(ispin_,ispin_,:,:,:),bath,ispin_)
             endif
          case ("nonsu2")
             call chi2_fitgf_normal_nonsu2(fg(:,:,:,:,:),bath)
          case default
             stop "chi2_fitgf ERROR: ed_mode!=normal/nonsu2 but only NORMAL component is provided"
          end select
          !
       case ("hybrid")
          select case(ed_mode)
          case ("normal")
             call chi2_fitgf_hybrid_normal(fg(ispin_,ispin_,:,:,:),bath,ispin_)
          case ("nonsu2")
             call chi2_fitgf_hybrid_nonsu2(fg(:,:,:,:,:),bath)
          case default
             stop "chi2_fitgf ERROR: ed_mode!=normal/nonsu2 but only NORMAL component is provided" 
          end select
          !
       case ("replica")
          select case(ed_mode)
          case ("normal","nonsu2")
             call chi2_fitgf_replica(fg,bath)
          case default
             stop "chi2_fitgf ERROR: ed_mode!=normal/nonsu2 but only NORMAL component is provided" 
          end select
       case ("general")
          select case(ed_mode)
          case ("normal","nonsu2")
             call chi2_fitgf_general(fg,bath)
          case default
             stop "chi2_fitgf ERROR: ed_mode!=normal/nonsu2 but only NORMAL component is provided" 
          end select
       end select
    endif
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
    if(check_MPI().AND.fmpi_)call ed_del_MpiComm()
#endif   
#ifdef _DEBUG
    write(Logfile,"(A)")""
#endif
  end subroutine chi2_fitgf_single_normal_n5







  subroutine chi2_fitgf_single_superc_n3(g,f,bath,ispin,iorb,fmpi)
    complex(8),dimension(:,:,:)                         :: g
    complex(8),dimension(:,:,:)                         :: f !anomalous Weiss field or hybridibazion function to fit (only if :f:var:`ed_mode` = :code:`superc` )
    real(8),dimension(:)                                :: bath
    integer,optional                                    :: ispin,iorb
    logical,optional                                    :: fmpi
    !
    complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Lmats) :: fg
    integer                                             :: ispin_
    logical                                             :: fmpi_
    !
    write(Logfile,"(A)")""
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG chi2_fitgf_generic_superc: Start Chi**2 fit"
#endif
    !
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    fmpi_=.true.;if(present(fmpi))fmpi_=fmpi
    !
#ifdef _MPI    
    if(check_MPI().AND.fmpi_)call ed_set_MpiComm()
#endif
    !
    select case(cg_method)
    case default
       stop "ED Error: cg_method > 1"
    case (0)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"master: Chi^2 fit with CG-nr and CG-weight: ",cg_weight," on: ",cg_scheme
    case (1)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"master: Chi^2 fit with CG-minimize and CG-weight: ",cg_weight," on: ",cg_scheme
    end select
    !
    if(MpiMaster)then
       !
       call assert_shape(g,[Nspin*Norb,Nspin*Norb,Lmats],"chi2_fitgf_generic_superc","g")
       fg(1,:,:,:,:,:) = so2nn_reshape(g(1:Nspin*Norb,1:Nspin*Norb,1:Lmats),Nspin,Norb,Lmats)
       call assert_shape(f,[Nspin*Norb,Nspin*Norb,Lmats],"chi2_fitgf_generic_superc","f")
       fg(2,:,:,:,:,:) = so2nn_reshape(f(1:Nspin*Norb,1:Nspin*Norb,1:Lmats),Nspin,Norb,Lmats)
       !       
       select case(bath_type)
       case default
          select case(ed_mode)
          case ("superc")
             if(present(iorb))then
                call chi2_fitgf_normal_superc(fg(:,ispin_,ispin_,:,:,:),bath,ispin_,iorb)
             else
                call chi2_fitgf_normal_superc(fg(:,ispin_,ispin_,:,:,:),bath,ispin_)
             endif
          case default
             write(LOGfile,"(A)") "chi2_fitgf WARNING: ed_mode=normal/nonsu2 but NORMAL & ANOMAL components provided."
             call chi2_fitgf_normal_normal(fg(1,ispin_,ispin_,:,:,:),bath,ispin_)          
          end select
       case ("hybrid")
          select case(ed_mode)
          case ("superc")
             call chi2_fitgf_hybrid_superc(fg(:,ispin_,ispin_,:,:,:),bath,ispin_)
          case default
             write(LOGfile,"(A)") "chi2_fitgf WARNING: ed_mode=normal/nonsu2 but NORMAL & ANOMAL components provided."
             call chi2_fitgf_hybrid_normal(fg(1,ispin_,ispin_,:,:,:),bath,ispin_)       
          end select
       case ("replica")
          select case(ed_mode)
          case ("superc")
             call chi2_fitgf_replica_superc(fg(:,:,:,:,:,:),bath)
          case default
             write(LOGfile,"(A)") "chi2_fitgf WARNING: ed_mode=normal/nonsu2 but NORMAL & ANOMAL components provided."
             call chi2_fitgf_replica(fg(1,:,:,:,:,:),bath)
          end select
       case ("general")
          select case(ed_mode)
          case ("superc")
             call chi2_fitgf_general_superc(fg(:,:,:,:,:,:),bath)
          case default
             write(LOGfile,"(A)") "chi2_fitgf WARNING: ed_mode=normal/nonsu2 but NORMAL & ANOMAL components provided."
             call chi2_fitgf_general(fg(1,:,:,:,:,:),bath)
          end select
       end select
    endif
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
#ifdef _MPI    
    if(check_MPI().AND.fmpi_)call ed_del_MpiComm()
#endif   
#ifdef _DEBUG
    write(Logfile,"(A)")""
#endif
  end subroutine chi2_fitgf_single_superc_n3


  subroutine chi2_fitgf_single_superc_n5(g,f,bath,ispin,iorb,fmpi)
    complex(8),dimension(:,:,:,:,:)                     :: g,f
    real(8),dimension(:)                                :: bath
    integer,optional                                    :: ispin,iorb
    logical,optional                                    :: fmpi
    !
    complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Lmats) :: fg
    integer                                             :: ispin_
    logical                                             :: fmpi_
    !
    write(Logfile,"(A)")""
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG chi2_fitgf_generic_superc: Start Chi**2 fit"
#endif
    !
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    fmpi_=.true.;if(present(fmpi))fmpi_=fmpi
    !
#ifdef _MPI    
    if(check_MPI().AND.fmpi_)call ed_set_MpiComm()
#endif
    !
    select case(cg_method)
    case default
       stop "ED Error: cg_method > 1"
    case (0)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"master: Chi^2 fit with CG-nr and CG-weight: ",cg_weight," on: ",cg_scheme
    case (1)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"master: Chi^2 fit with CG-minimize and CG-weight: ",cg_weight," on: ",cg_scheme
    end select
    !
    if(MpiMaster)then
       !
       call assert_shape(g,[Nspin,Nspin,Norb,Norb,Lmats],"chi2_fitgf_generic_superc","g")
       fg(1,:,:,:,:,:) = g(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats)
       call assert_shape(f,[Nspin,Nspin,Norb,Norb,Lmats],"chi2_fitgf_generic_superc","f")
       fg(2,:,:,:,:,:) = f(1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats)
       !       
       select case(bath_type)
       case default
          select case(ed_mode)
          case ("superc")
             if(present(iorb))then
                call chi2_fitgf_normal_superc(fg(:,ispin_,ispin_,:,:,:),bath,ispin_,iorb)
             else
                call chi2_fitgf_normal_superc(fg(:,ispin_,ispin_,:,:,:),bath,ispin_)
             endif
          case default
             write(LOGfile,"(A)") "chi2_fitgf WARNING: ed_mode=normal/nonsu2 but NORMAL & ANOMAL components provided."
             call chi2_fitgf_normal_normal(fg(1,ispin_,ispin_,:,:,:),bath,ispin_)          
          end select
       case ("hybrid")
          select case(ed_mode)
          case ("superc")
             call chi2_fitgf_hybrid_superc(fg(:,ispin_,ispin_,:,:,:),bath,ispin_)
          case default
             write(LOGfile,"(A)") "chi2_fitgf WARNING: ed_mode=normal/nonsu2 but NORMAL & ANOMAL components provided."
             call chi2_fitgf_hybrid_normal(fg(1,ispin_,ispin_,:,:,:),bath,ispin_)       
          end select
       case ("replica")
          select case(ed_mode)
          case ("superc")
             call chi2_fitgf_replica_superc(fg(:,:,:,:,:,:),bath)
          case default
             write(LOGfile,"(A)") "chi2_fitgf WARNING: ed_mode=normal/nonsu2 but NORMAL & ANOMAL components provided."
             call chi2_fitgf_replica(fg(1,:,:,:,:,:),bath)
          end select
       case ("general")
          select case(ed_mode)
          case ("superc")
             call chi2_fitgf_general_superc(fg(:,:,:,:,:,:),bath)
          case default
             write(LOGfile,"(A)") "chi2_fitgf WARNING: ed_mode=normal/nonsu2 but NORMAL & ANOMAL components provided."
             call chi2_fitgf_general(fg(1,:,:,:,:,:),bath)
          end select
       end select
    endif
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
#ifdef _MPI    
    if(check_MPI().AND.fmpi_)call ed_del_MpiComm()
#endif   
#ifdef _DEBUG
    write(Logfile,"(A)")""
#endif
  end subroutine chi2_fitgf_single_superc_n5





  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################



  subroutine chi2_fitgf_lattice_normal_n3(g,bath,ispin)
    real(8),intent(inout)       :: bath(:,:)
    complex(8),dimension(:,:,:) :: g
    integer,optional            :: ispin
    !
    complex(8)                  :: fg(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    real(8)                     :: bath_tmp(size(bath,1),size(bath,2))
    integer                     :: i,ispin_
    integer                     :: ilat
    integer                     :: iorb,is,io
    integer                     :: jorb,js,jo
    integer                     :: Nsites
    logical                     :: check_dim
    character(len=5)            :: tmp_suffix
    integer                     :: MPI_ID=0
    integer                     :: MPI_SIZE=1
    logical                     :: MPI_MASTER=.true.
    !
    write(Logfile,"(A)")""
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
    !
#ifdef _MPI    
    if(check_MPI())then
       MPI_ID     = get_Rank_MPI()
       MPI_SIZE   = get_Size_MPI()
       MPI_MASTER = get_Master_MPI()
    endif
#endif
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    do ilat=1+MPI_ID,Nsites,MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    fg = zero
    !
    if(ed_verbose>1)write(Logfile,"(A)")"Chi**2 get G with rank:"//str(rank(g))
    call assert_shape(g,[Nsites*Nspin*Norb,Nsites*Nspin*Norb,Lmats],'chi2_fitgf_generic_normal','g')
    fg = lso2nnn_reshape(g(1:Nsites*Nspin*Norb,1:Nsites*Nspin*Norb,1:Lmats),Nsites,Nspin,Norb,Lmats)
    !
    !
    bath_tmp=0d0
    do ilat = 1+MPI_ID,Nsites,MPI_SIZE
       if(ed_verbose>1)write(Logfile,"(A)")"Start Chi**2 fit for site:"//str(ilat)
       bath_tmp(ilat,:)=bath(ilat,:)
       call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       !
       if(present(ispin))then
          call ed_chi2_fitgf(fg(ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin=ispin_,fmpi=.false.)
       else
          call ed_chi2_fitgf(fg(ilat,:,:,:,:,:),bath_tmp(ilat,:),fmpi=.false.)
       end if
    end do
    !
#ifdef _MPI
    if(check_MPI())then
       bath=0d0
       call AllReduce_MPI(MPI_COMM_WORLD,bath_tmp,bath)
    else
       bath = bath_tmp
    endif
#else
    bath = bath_tmp
#endif
    !
    call ed_reset_suffix
  end subroutine chi2_fitgf_lattice_normal_n3

  subroutine chi2_fitgf_lattice_normal_n4(g,bath,ispin)
    real(8),intent(inout)         :: bath(:,:)
    complex(8),dimension(:,:,:,:) :: g
    integer,optional              :: ispin
    !
    complex(8)                    :: fg(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    real(8)                       :: bath_tmp(size(bath,1),size(bath,2))
    integer                       :: i,ispin_
    integer                       :: ilat
    integer                       :: iorb,is,io
    integer                       :: jorb,js,jo
    integer                       :: Nsites
    logical                       :: check_dim
    character(len=5)              :: tmp_suffix
    integer                       :: MPI_ID=0
    integer                       :: MPI_SIZE=1
    logical                       :: MPI_MASTER=.true.
    !
    write(Logfile,"(A)")""
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
    !
#ifdef _MPI    
    if(check_MPI())then
       MPI_ID     = get_Rank_MPI()
       MPI_SIZE   = get_Size_MPI()
       MPI_MASTER = get_Master_MPI()
    endif
#endif
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    do ilat=1+MPI_ID,Nsites,MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    fg = zero
    !
    if(ed_verbose>1)write(Logfile,"(A)")"Chi**2 get G with rank:"//str(rank(g))
    call assert_shape(g,[Nsites,Nspin*Norb,Nspin*Norb,Lmats],'chi2_fitgf_generic_normal','g')
    do ilat=1,Nsites
       fg(ilat,:,:,:,:,:)  = so2nn_reshape(g(ilat,1:Nspin*Norb,1:Nspin*Norb,1:Lmats),Nspin,Norb,Lmats)
    enddo
    !
    !
    bath_tmp=0d0
    do ilat = 1+MPI_ID,Nsites,MPI_SIZE
       if(ed_verbose>1)write(Logfile,"(A)")"Start Chi**2 fit for site:"//str(ilat)
       bath_tmp(ilat,:)=bath(ilat,:)
       call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       !
       if(present(ispin))then
          call ed_chi2_fitgf(fg(ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin=ispin_,fmpi=.false.)
       else
          call ed_chi2_fitgf(fg(ilat,:,:,:,:,:),bath_tmp(ilat,:),fmpi=.false.)
       end if
    end do
    !
#ifdef _MPI
    if(check_MPI())then
       bath=0d0
       call AllReduce_MPI(MPI_COMM_WORLD,bath_tmp,bath)
    else
       bath = bath_tmp
    endif
#else
    bath = bath_tmp
#endif
    !
    call ed_reset_suffix
  end subroutine chi2_fitgf_lattice_normal_n4

  subroutine chi2_fitgf_lattice_normal_n6(g,bath,ispin)
    real(8),intent(inout)             :: bath(:,:)
    complex(8),dimension(:,:,:,:,:,:) :: g
    integer,optional                  :: ispin
    !
    complex(8)                        :: fg(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    real(8)                           :: bath_tmp(size(bath,1),size(bath,2))
    integer                           :: i,ispin_
    integer                           :: ilat
    integer                           :: iorb,is,io
    integer                           :: jorb,js,jo
    integer                           :: Nsites
    logical                           :: check_dim
    character(len=5)                  :: tmp_suffix
    integer                           :: MPI_ID=0
    integer                           :: MPI_SIZE=1
    logical                           :: MPI_MASTER=.true.
    !
    write(Logfile,"(A)")""
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
    !
#ifdef _MPI    
    if(check_MPI())then
       MPI_ID     = get_Rank_MPI()
       MPI_SIZE   = get_Size_MPI()
       MPI_MASTER = get_Master_MPI()
    endif
#endif
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    do ilat=1+MPI_ID,Nsites,MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    fg = zero
    !
    if(ed_verbose>1)write(Logfile,"(A)")"Chi**2 get G with rank:"//str(rank(g))
    call assert_shape(g,[Nsites,Nspin,Nspin,Norb,Norb,Lmats],'chi2_fitgf_generic_normal','g')
    fg = g(1:Nsites,1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats)
    !
    bath_tmp=0d0
    do ilat = 1+MPI_ID,Nsites,MPI_SIZE
       if(ed_verbose>1)write(Logfile,"(A)")"Start Chi**2 fit for site:"//str(ilat)
       bath_tmp(ilat,:)=bath(ilat,:)
       call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       !
       if(present(ispin))then
          call ed_chi2_fitgf(fg(ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin=ispin_,fmpi=.false.)
       else
          call ed_chi2_fitgf(fg(ilat,:,:,:,:,:),bath_tmp(ilat,:),fmpi=.false.)
       end if
    end do
    !
#ifdef _MPI
    if(check_MPI())then
       bath=0d0
       call AllReduce_MPI(MPI_COMM_WORLD,bath_tmp,bath)
    else
       bath = bath_tmp
    endif
#else
    bath = bath_tmp
#endif
    !
    call ed_reset_suffix
  end subroutine chi2_fitgf_lattice_normal_n6





  subroutine chi2_fitgf_lattice_superc_n3(g,f,bath,ispin)
    real(8),intent(inout)       :: bath(:,:)
    complex(8),dimension(:,:,:) :: g,f
    integer,optional            :: ispin
    !
    complex(8)                  :: fg(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    real(8)                     :: bath_tmp(size(bath,1),size(bath,2))
    integer                     :: ilat,i,iorb,ispin_
    integer                     :: Nsites
    logical                     :: check_dim
    character(len=5)            :: tmp_suffix
    integer                     :: MPI_ID=0
    integer                     :: MPI_SIZE=1
    logical                     :: MPI_MASTER=.true.
    !
    write(Logfile,"(A)")""
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
    !
#ifdef _MPI    
    if(check_MPI())then
       MPI_ID     = get_Rank_MPI()
       MPI_SIZE   = get_Size_MPI()
       MPI_MASTER = get_Master_MPI()
    endif
#endif
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    do ilat = 1 + MPI_ID, Nsites, MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    call assert_shape(g,[Nsites*Nspin*Norb,Nsites*Nspin*Norb,Lmats],'chi2_fitgf_generic_superc','g')
    fg(1,:,:,:,:,:,:) = lso2nnn_reshape(g(1:Nsites*Nspin*Norb,1:Nsites*Nspin*Norb,1:Lmats),Nsites,Nspin,Norb,Lmats)
    !
    call assert_shape(f,[Nsites*Nspin*Norb,Nsites*Nspin*Norb,Lmats],'chi2_fitgf_generic_superc','f')
    fg(2,:,:,:,:,:,:) = lso2nnn_reshape(f(1:Nsites*Nspin*Norb,1:Nsites*Nspin*Norb,1:Lmats),Nsites,Nspin,Norb,Lmats)
    !
    !
    bath_tmp=0.d0
    !
    do ilat= 1 + MPI_ID, Nsites, MPI_SIZE
       write(Logfile,"(A)")"ed_fit_bath_sites_superc: Start Chi**2 fit for site:"//str(ilat)
       bath_tmp(ilat,:) = bath(ilat,:)
       call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       !
       if(present(ispin))then
          call ed_chi2_fitgf(fg(1,ilat,:,:,:,:,:),fg(2,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_,fmpi=.false.)
       else
          do ispin_=1,Nspin
             call ed_chi2_fitgf(fg(1,ilat,:,:,:,:,:),fg(2,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_,fmpi=.false.)
          enddo
       endif
    end do
    !
#ifdef _MPI
    if(check_MPI())then
       bath=0d0
       call AllReduce_MPI(MPI_COMM_WORLD,bath_tmp,bath)
    else
       bath = bath_tmp
    endif
#else
    bath = bath_tmp
#endif
    !
    call ed_reset_suffix
  end subroutine chi2_fitgf_lattice_superc_n3






  subroutine chi2_fitgf_lattice_superc_n4(g,f,bath,ispin)
    real(8),intent(inout)         :: bath(:,:)
    complex(8),dimension(:,:,:,:) :: g,f
    integer,optional              :: ispin
    !
    complex(8)                    :: fg(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    real(8)                       :: bath_tmp(size(bath,1),size(bath,2))
    integer                       :: ilat,i,iorb,ispin_
    integer                       :: Nsites
    logical                       :: check_dim
    character(len=5)              :: tmp_suffix
    integer                       :: MPI_ID=0
    integer                       :: MPI_SIZE=1
    logical                       :: MPI_MASTER=.true.
    !
    write(Logfile,"(A)")""
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
    !
#ifdef _MPI    
    if(check_MPI())then
       MPI_ID     = get_Rank_MPI()
       MPI_SIZE   = get_Size_MPI()
       MPI_MASTER = get_Master_MPI()
    endif
#endif
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    do ilat = 1 + MPI_ID, Nsites, MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    call assert_shape(g,[Nsites,Nspin*Norb,Nspin*Norb,Lmats],'chi2_fitgf_generic_superc','g')
    do ilat=1,Nsites
       fg(1,ilat,:,:,:,:,:)  = so2nn_reshape(g(ilat,1:Nspin*Norb,1:Nspin*Norb,1:Lmats),Nspin,Norb,Lmats)
    enddo
    call assert_shape(f,[Nsites,Nspin*Norb,Nspin*Norb,Lmats],'chi2_fitgf_generic_superc','f')
    do ilat=1,Nsites
       fg(2,ilat,:,:,:,:,:)  = so2nn_reshape(f(ilat,1:Nspin*Norb,1:Nspin*Norb,1:Lmats),Nspin,Norb,Lmats)
    enddo
    !
    bath_tmp=0.d0
    !
    do ilat= 1 + MPI_ID, Nsites, MPI_SIZE
       write(Logfile,"(A)")"ed_fit_bath_sites_superc: Start Chi**2 fit for site:"//str(ilat)
       bath_tmp(ilat,:) = bath(ilat,:)
       call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       !
       if(present(ispin))then
          call ed_chi2_fitgf(fg(1,ilat,:,:,:,:,:),fg(2,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_,fmpi=.false.)
       else
          do ispin_=1,Nspin
             call ed_chi2_fitgf(fg(1,ilat,:,:,:,:,:),fg(2,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_,fmpi=.false.)
          enddo
       endif
    end do
    !
#ifdef _MPI
    if(check_MPI())then
       bath=0d0
       call AllReduce_MPI(MPI_COMM_WORLD,bath_tmp,bath)
    else
       bath = bath_tmp
    endif
#else
    bath = bath_tmp
#endif
    !
    call ed_reset_suffix
  end subroutine chi2_fitgf_lattice_superc_n4







  subroutine chi2_fitgf_lattice_superc_n6(g,f,bath,ispin)
    real(8),intent(inout)             :: bath(:,:)
    complex(8),dimension(:,:,:,:,:,:) :: g,f
    integer,optional                  :: ispin
    !
    complex(8)                        :: fg(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    real(8)                           :: bath_tmp(size(bath,1),size(bath,2))
    integer                           :: ilat,i,iorb,ispin_
    integer                           :: Nsites
    logical                           :: check_dim
    character(len=5)                  :: tmp_suffix
    integer                           :: MPI_ID=0
    integer                           :: MPI_SIZE=1
    logical                           :: MPI_MASTER=.true.
    !
    write(Logfile,"(A)")""
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
    !
#ifdef _MPI    
    if(check_MPI())then
       MPI_ID     = get_Rank_MPI()
       MPI_SIZE   = get_Size_MPI()
       MPI_MASTER = get_Master_MPI()
    endif
#endif
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    do ilat = 1 + MPI_ID, Nsites, MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    call assert_shape(g,[Nsites,Nspin,Nspin,Norb,Norb,Lmats],'chi2_fitgf_generic_superc','g')
    fg(1,:,:,:,:,:,:) = g(1:Nsites,1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats)
    !
    call assert_shape(f,[Nsites,Nspin,Nspin,Norb,Norb,Lmats],'chi2_fitgf_generic_superc','f')
    fg(2,:,:,:,:,:,:) = f(1:Nsites,1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats)    
    !
    bath_tmp=0.d0
    !
    do ilat= 1 + MPI_ID, Nsites, MPI_SIZE
       write(Logfile,"(A)")"ed_fit_bath_sites_superc: Start Chi**2 fit for site:"//str(ilat)
       bath_tmp(ilat,:) = bath(ilat,:)
       call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       !
       if(present(ispin))then
          call ed_chi2_fitgf(fg(1,ilat,:,:,:,:,:),fg(2,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_,fmpi=.false.)
       else
          do ispin_=1,Nspin
             call ed_chi2_fitgf(fg(1,ilat,:,:,:,:,:),fg(2,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_,fmpi=.false.)
          enddo
       endif
    end do
    !
#ifdef _MPI
    if(check_MPI())then
       bath=0d0
       call AllReduce_MPI(MPI_COMM_WORLD,bath_tmp,bath)
    else
       bath = bath_tmp
    endif
#else
    bath = bath_tmp
#endif
    !
    call ed_reset_suffix
  end subroutine chi2_fitgf_lattice_superc_n6




end MODULE ED_BATH_FIT

















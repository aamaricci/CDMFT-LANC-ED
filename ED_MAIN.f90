module ED_MAIN
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE, only: state_list,es_delete_espace,delete_eigenspace
  USE ED_AUX_FUNX
  USE ED_SETUP
  USE ED_BATH
  USE ED_HAMILTONIAN
  USE ED_GREENS_FUNCTIONS
  USE ED_OBSERVABLES
  USE ED_DIAG
  USE SF_IOTOOLS, only: str,reg
  USE SF_TIMER,only: start_timer,stop_timer
  implicit none
  private



  interface ed_init_solver
     !
     !Initialize the Exact Diagonalization solver of `CDMFT`. This procedure reserves and allocates all the  
     !memory required by the solver, performs all the consistency check and initializes the bath instance guessing or reading from a file.      
     !It requires as an input a double precision array of rank-1 [ :f:var:`nb` ] for the single-impurity case or
     !or rank-2 [ :f:var:`nb` , :f:var:`nlat` ] for the Real space DMFT case. :f:var:`nlat` is the number of inequivalent impurity sites,
     !while :f:var:`nb` depends on the bath size and geometry and can be obtained from :f:func:`get_bath_dimension` .
     !
     module procedure :: ed_init_solver_single
  end interface ed_init_solver

  interface ed_solve
     !
     !Launch the Exact Diagonalizaton solver for the single-site and multiple-site (R-DMFT) quantum impurity problem.
     !It requires as an input a double precision array of rank-1 [ :f:var:`nb` ] for the single-impurity case or
     !or rank-2 [ :f:var:`nb` , :f:var:`nlat` ] for the Real space DMFT case. :f:var:`nlat` is the number of inequivalent impurity sites,
     !while :f:var:`nb` depends on the bath size and geometry and can be obtained from :f:func:`get_bath_dimension` .
     !
     ! The solution is achieved in this sequence:
     !
     !  #. setup the MPI environment, if any 
     !  #. Set the internal bath instance :f:var:`dmft_bath` copying from the user provided input :f:var:`bath`
     !  #. Get the low energy spectrum: call :f:func:`diagonalize_impurity`
     !  #. Get the impurity Green's functions: call :f:func:`buildgf_impurity` (if :f:var:`sflag` = :code:`.true.` )
     !  #. Get the impurity susceptibilities, if any: call :f:func:`buildchi_impurity` (if :f:var:`sflag` = :code:`.true.` )
     !  #. Get the impurity observables: call :f:func:`observables_impurity`
     !  #. Get the impurity local energy: call :f:func:`local_energy_impurity`
     !  #. Delete MPI environment and deallocate used structures :f:var:`state_list` and :f:var:`dmft_bath`
     !
     module procedure :: ed_solve_single
  end interface ed_solve



  !> FINALIZE SOLVER AND CLEANUP ENVIRONMENT
  interface ed_finalize_solver
     ! 
     ! Finalize the Exact Diagonalization solver, clean up all the allocated memory. 
     !
     module procedure :: ed_finalize_solver_single
  end interface ed_finalize_solver


  public :: ed_init_solver
  public :: ed_solve
  public :: ed_finalize_solver


  logical,save :: isetup=.true. !Allow :f:func:`init_ed_structure` and :f:func:`setup_global` to be called. Set to :f:code:`.false.` by :f:func:`ed_init_solver`, reset by :f:func:`ed_finalize_solver` . Default :code:`.true.`



contains





  ! PURPOSE: allocate and initialize one or multiple baths -+!
  subroutine ed_init_solver_single(bath)
    real(8),dimension(:),intent(inout) :: bath
    logical                            :: check 
    logical,save                       :: isetup=.true.
    integer                            :: i
    !
    !SET THE MPI FRAMEWORK:
#ifdef _MPI
    if(check_MPI())call ed_set_MpiComm()
#endif
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(ed_file_suffix)
    !
    !Init ED Structure & memory
    if(isetup)call init_ed_structure()
    !
    !Init bath:
    !call set_Hloc(Hloc)
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "init_ed_solver_single error: wrong bath dimensions"
    !
    bath = 0d0
    !
    call allocate_dmft_bath(dmft_bath)
    call init_dmft_bath(dmft_bath)
    call get_dmft_bath(dmft_bath,bath)
    !
    if(isetup)call setup_global
    call deallocate_dmft_bath(dmft_bath)
    isetup=.false.
    !
    !DELETE THE MPI FRAMEWORK:
#ifdef _MPI
    if(check_MPI())call ed_del_MpiComm()
#endif
    !
  end subroutine ed_init_solver_single






  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################





  !+-----------------------------------------------------------------------------+!
  !PURPOSE: solve the impurity problems for a single or many independent
  ! lattice site using ED. 
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_single(bath,sflag)
    real(8),dimension(:),intent(in) :: bath
    logical,optional                :: sflag !flag to calculate ( :code:`.true.` ) or not ( :code:`.false.` ) Green's functions and susceptibilities. Default :code:`.true.` . 
    logical                         :: check,iflag
    !
    iflag=.true. ;if(present(sflag))iflag=sflag
    !
#ifdef _MPI    
    if(check_MPI())call ed_set_MpiComm()
#endif
    !
    if(.not.allocated(impHloc))stop "ED_SOLVE ERROR: impHloc not allocated. Please call ed_set_Hloc first."
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "ED_SOLVE error: wrong bath dimensions"
    !
    if(MpiMaster)call save_input_file(str(ed_input_file))
    !
    call allocate_dmft_bath(dmft_bath)
    call set_dmft_bath(bath,dmft_bath)
    call write_dmft_bath(dmft_bath)
    if(MpiMaster)call save_dmft_bath(dmft_bath,used=.true.)
    !
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_impurity()                   !find target states by digonalization of Hamiltonian
    if(iflag)call buildgf_impurity()            !build the one-particle impurity Green's functions & Self-energy
    ! if(iflag)call get_custom_observables()      !obtain custom user-defined observables (if initialized)
    !if(iflag)call buildchi_impurity()         !build the local susceptibilities (todo)
    call observables_impurity()                   !obtain impurity observables as thermal averages.
    call local_energy_impurity()                  !obtain the local energy of the effective impurity problem
    !
    !OPTIONAL (HEAVY) CALCULATIONS 
    if(dm_flag) call density_matrix_impurity()    !build the cluster density matrix (\rho_IMP = Tr_BATH(\rho))
    !
    call deallocate_dmft_bath(dmft_bath)
    call es_delete_espace(state_list)
    !
    !DELETE THE LOCAL MPI COMMUNICATOR:
#ifdef _MPI    
    if(check_MPI())call ed_del_MpiComm()
#endif    
    !
    nullify(spHtimesV_p)
    write(Logfile,"(A)")""
  end subroutine ed_solve_single






  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: deallocate and finalize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_finalize_solver_single()
    logical                            :: check
    integer                            :: i
    !
    !SET THE MPI FRAMEWORK:
#ifdef _MPI
    if(check_MPI())call ed_set_MpiComm()
#endif
    !
    write(LOGfile,"(A)")"FINALIZE SOLVER "
    !
    !just in case deallocate some objects
    call deallocate_dmft_bath(dmft_bath)
    call es_delete_espace(state_list)
    call deallocate_GFmatrix(impGmatrix)
    call deallocate_grids
    nullify(spHtimesV_p)
    Hbath_status=.false.
    !
    !Delete ED Structure & memory
    call delete_ed_structure()
    !
    !Ready to be setup again
    isetup=.true.
    !
    !DELETE THE MPI FRAMEWORK:
#ifdef _MPI
    if(check_MPI())call ed_del_MpiComm()
#endif
    !
  end subroutine ed_finalize_solver_single







end module ED_MAIN









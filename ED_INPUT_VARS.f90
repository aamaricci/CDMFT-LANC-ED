MODULE ED_INPUT_VARS
  USE SF_VERSION
  USE SF_PARSE_INPUT
  USE SF_IOTOOLS, only:str,free_unit
  USE ED_VERSION
  use iso_c_binding
  implicit none


  !input variables
  !=========================================================
  integer(c_int), bind(c, name="Nlat")             :: Nlat                !Number of cluster sites
  integer(c_int), bind(c, name="Norb")             :: Norb                !Number of impurity orbitals
  integer(c_int), bind(c, name="Nspin")            :: Nspin               !Number spin degeneracy (max 2)
  integer(c_int), bind(c, name="Nbath")            :: Nbath               !Number of bath sites (per orbital or not depending on bath_type)
  integer(c_int), bind(c, name="Nambu")            :: Nambu               !Global Nambu factor for SC calculations

  integer(c_int), bind(c, name="Nloop")            :: Nloop               !max dmft loop variables
  real(c_double),dimension(5),bind(c, name="Uloc") :: Uloc                !local interactions
  real(c_double),bind(c, name="Ust")               :: Ust                 !intra-orbitals interactions
  real(c_double),bind(c, name="Jh")                :: Jh                  !J_Hund: Hunds' coupling constant 
  real(c_double),bind(c, name="Jx")                :: Jx                  !J_X: coupling constant for the spin-eXchange interaction term
  real(c_double),bind(c, name="Jp")                :: Jp                  !J_P: coupling constant for the Pair-hopping interaction term 
  real(c_double),bind(c, name="xmu")               :: xmu                 !chemical potential
  real(c_double),bind(c, name="beta")              :: beta                !inverse temperature

  integer(c_int), bind(c, name="Nsuccess")         :: Nsuccess            !Number of repeated success to fall below convergence threshold  
  real(c_double),bind(c, name="dmft_error")        :: dmft_error          !dmft convergence threshold
  real(c_double),bind(c, name="eps")               :: eps                 !broadening
  real(c_double),bind(c, name="wini")              :: wini                !frequency range min
  real(c_double),bind(c, name="wfin")              :: wfin                !frequency range max
  real(c_double),bind(c, name="sb_field")          :: sb_field            !symmetry breaking field
  real(c_double),bind(c, name="nread")             :: nread               !fixed density. if 0.d0 fixed chemical potential calculation.
  logical(c_bool),bind(c, name="ed_twin")          :: ed_twin             !flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry.
  logical                                          :: ed_twin_


  logical                                          :: chiflag             !
  logical                                          :: gf_flag             !flag to evaluate Green's functions (and related quantities)
  logical                                          :: dm_flag             !flag to evaluate the cluster density matrix \rho_IMP = Tr_BATH(\rho)) 
  logical                                          :: HFmode              !flag for HF interaction form U(n-1/2)(n-1/2) VS Unn
  real(8)                                          :: cutoff              !cutoff for spectral summation
  real(8)                                          :: gs_threshold        !Energy threshold for ground state degeneracy loop up
  real(8)                                          :: deltasc             !breaking symmetry field

  character(len=7)                                 :: ed_mode             !flag to set ed symmetry type: normal=normal (default), superc=superconductive, nonsu2=broken SU(2)
  logical                                          :: ed_finite_temp      !flag to select finite temperature method. note that if T then lanc_nstates_total must be > 1 
  logical                                          :: ed_sparse_H         !flag to select  storage of sparse matrix H (mem--, cpu++) if TRUE, or direct on-the-fly H*v product (mem++, cpu--
  logical                                          :: ed_gf_symmetric     !flag to assume G_ij = G_ji
  logical                                          :: ed_print_Sigma      !flag to print impurity Self-energies
  logical                                          :: ed_print_G          !flag to print impurity Green`s functions
  logical                                          :: ed_print_G0         !flag to print impurity non-interacting Green`s functions
  logical                                          :: ed_sectors          !flag to reduce sector scan for the spectrum to specific sectors +/- ed_sectors_shift
  integer                                          :: ed_sectors_shift    !shift to the ed_sectors scan
  integer                                          :: ed_verbose          !
  real(8)                                          :: ed_offset_bath      !half-bandwidth for the bath initialization: flat in -hwband:hwband
  real(8)                                          :: ed_hw_band              !half-bandwidth for the bath initialization (diagonal part in Hbath)

  character(len=12)                                :: lanc_method         !select the lanczos method to be used in the determination of the spectrum. ARPACK (default), LANCZOS (T=0 only) 
  real(8)                                          :: lanc_tolerance      !Tolerance for the Lanczos iterations as used in Arpack and plain lanczos. 
  integer                                          :: lanc_niter          !Max number of Lanczos iterations
  integer                                          :: lanc_ngfiter        !Max number of iteration in resolvant tri-diagonalization
  integer                                          :: lanc_ncv_factor     !Set the size of the block used in Lanczos-Arpack by multiplying the required Neigen (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
  integer                                          :: lanc_ncv_add        !Adds up to the size of the block to prevent it to become too small (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
  integer                                          :: lanc_nstates_sector !Max number of required eigenvalues per sector
  integer                                          :: lanc_nstates_total  !Max number of states hold in the finite T calculation
  integer                                          :: lanc_nstates_step   !Number of states added at each step to determine the optimal spectrum size at finite T
  integer                                          :: lanc_dim_threshold  !Min dimension threshold to use Lanczos determination of the spectrum rather than Lapack based exact diagonalization.
  !
  character(len=5)                                 :: cg_scheme           !fit scheme: delta (default) or weiss for G0and

  integer                                          :: cg_method           !fit routine type:0=CGnr (default), 1=minimize (old f77)
  integer                                          :: cg_grad             !gradient evaluation: 0=analytic, 1=numeric
  integer                                          :: cg_niter            !Max number of iteration in the fit
  real(8)                                          :: cg_ftol             !Tolerance in the cg fit
  integer                                          :: cg_stop             !fit stop condition:0-3, 0=C1.AND.C2, 1=C1, 2=C2 with C1=|F_n-1 -F_n|<tol*(1+F_n), C2=||x_n-1 -x_n||<tol*(1+||x_n||).
  integer                                          :: cg_matrix           !fit weights for matrix elements: for now flat or 'spectral' (normalized according to \sum_iw A(iw))
  integer                                          :: cg_weight           !fit weights on matsubara axis: 0=1, 1=1/n , 2=1/w_n
  integer                                          :: cg_pow              !fit power to generalize the distance as |G0 - G0and|**cg_pow
  character(len=9)                                 :: cg_norm             !fit norm: elemental (default) or frobenius norm to evaluate distances |G0 - G0and|
  logical                                          :: cg_minimize_ver     !flag to pick old (Krauth) or new (Lichtenstein) version of the minimize CG routine
  real(8)                                          :: cg_minimize_hh      !unknown parameter used in the CG minimize procedure.  
  !
  logical                                          :: finiteT             !flag for finite temperature calculation (UNIMPLEMENTED: hard-linked to lanc_nstates_total=1)
  character(len=7)                                 :: bath_type           !bath representation choice (here either "replica" or "general")
  !
  real(8)                                          :: nerr                !fix density threshold. a loop over from 1.d-1 to required nerr is performed
  real(8)                                          :: ndelta              !initial chemical potential step
  real(8)                                          :: ncoeff              !multiplier for the initial ndelta read from a file (ndelta-->ndelta*ncoeff)
  integer                                          :: niter               !

  !Some parameters for function dimension:
  !=========================================================
  !Some parameters for function dimension:
  integer(c_int),bind(c, name="Lmats")             :: Lmats !Number of Matsubara frequencies
  integer(c_int),bind(c, name="Lreal")             :: Lreal !Number of real-axis frequencies
  integer(c_int),bind(c, name="Lfit")              :: Lfit  !Number of frequencies for bath fitting
  integer(c_int),bind(c, name="Ltau")              :: Ltau  !Number of imaginary time points


  !LOG AND Hamiltonian UNITS
  !=========================================================
  character(len=100)                               :: Hfile  !File where to retrieve/store the bath parameters.
  character(len=100)                               :: HLOCfile !File read the input local H
  character(len=100)                               :: SectorFile !File where to retrieve/store the sectors contributing to the spectrum
  integer(c_int),bind(c, name="LOGfile"),save      :: LOGfile  !Logfile unit

  !THIS IS JUST A RELOCATED GLOBAL VARIABLE
  character(len=200)                               :: ed_input_file="" !Name of input file



contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : READ THE INPUT FILE AND SETUP GLOBAL VARIABLES
  !+-------------------------------------------------------------------+
  subroutine ed_read_input(INPUTunit,comm)
#ifdef _MPI
    USE MPI
    USE SF_MPI
#endif
    character(len=*) :: INPUTunit
    integer,optional :: comm
    logical          :: master=.true.,bool
    integer          :: i,rank=0
    integer          :: unit_xmu
#ifdef _MPI
    if(present(comm))then
       master=get_Master_MPI(comm)
       rank  =get_Rank_MPI(comm)
    endif
#endif
    !
    !Store the name of the input file:
    ed_input_file=str(INPUTunit)
    !
    !DEFAULT VALUES OF THE PARAMETERS:
    call parse_input_variable(Nlat,"NLAT",INPUTunit,default=1,comment="Number of cluster sites")
    call parse_input_variable(Norb,"NORB",INPUTunit,default=1,comment="Number of impurity orbitals (max 5).")
    call parse_input_variable(Nspin,"NSPIN",INPUTunit,default=1,comment="Number of spin degeneracy (max 2)")
    call parse_input_variable(Nbath,"NBATH",INPUTunit,default=4,comment="Number of bath clusters (replicas or generalizations)")
    ! call parse_input_variable(Nbath,"NAMBU",INPUTunit,default=1,comment="Global Nambu factor for SC calculations")
    call parse_input_variable(bath_type,"BATH_TYPE",INPUTunit,default='replica',comment="flag to set bath type: 'replica' or 'general'")
    
    call parse_input_variable(uloc,"ULOC",INPUTunit,default=[2d0,0d0,0d0,0d0,0d0],comment="Values of the local interaction per orbital (max 5)")
    call parse_input_variable(ust,"UST",INPUTunit,default=0.d0,comment="Value of the inter-orbital interaction term")
    call parse_input_variable(Jh,"JH",INPUTunit,default=0.d0,comment="Hunds coupling")
    call parse_input_variable(Jx,"JX",INPUTunit,default=0.d0,comment="S-E coupling")
    call parse_input_variable(Jp,"JP",INPUTunit,default=0.d0,comment="P-H coupling")
    call parse_input_variable(beta,"BETA",INPUTunit,default=1000.d0,comment="Inverse temperature, at T=0 is used as a IR cut-off.")
    call parse_input_variable(xmu,"XMU",INPUTunit,default=0.d0,comment="Chemical potential. If HFMODE=T, xmu=0 indicates half-filling condition.")
    call parse_input_variable(nloop,"NLOOP",INPUTunit,default=100,comment="Max number of DMFT iterations.")
    call parse_input_variable(dmft_error,"DMFT_ERROR",INPUTunit,default=0.00001d0,comment="Error threshold for DMFT convergence")
    call parse_input_variable(sb_field,"SB_FIELD",INPUTunit,default=0.1d0,comment="Value of a symmetry breaking field for magnetic solutions.")
    call parse_input_variable(gf_flag,"GF_FLAG",INPUTunit,default=.true.,comment="flag to evaluate GFs and related quantities.")
    call parse_input_variable(dm_flag,"DM_FLAG",INPUTunit,default=.false.,comment="flag to evaluate the cluster density matrix \rho_IMP = Tr_BATH(\rho))")
    !
    call parse_input_variable(ed_mode,"ED_MODE",INPUTunit,default='normal',comment="Flag to set ED type: normal=normal, superc=superconductive, nonsu2=broken SU(2)")
    Nambu = 1; if(ed_mode=="superc")Nambu = 2
    !
    call parse_input_variable(ed_twin_,"ED_TWIN",INPUTunit,default=.false.,comment="flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry.")
    ed_twin = ed_twin_
    !
    call parse_input_variable(ed_sectors,"ED_SECTORS",INPUTunit,default=.false.,comment="flag to reduce sector scan for the spectrum to specific sectors +/- ed_sectors_shift.")
    call parse_input_variable(ed_sectors_shift,"ED_SECTORS_SHIFT",INPUTunit,1,comment="shift to ed_sectors")
    call parse_input_variable(ed_sparse_H,"ED_SPARSE_H",INPUTunit,default=.true.,comment="flag to select  storage of sparse matrix H (mem--, cpu++) if TRUE, or direct on-the-fly H*v product (mem++, cpu--) if FALSE ")
    call parse_input_variable(ed_gf_symmetric,"ED_GF_SYMMETRIC",INPUTunit,default=.false.,comment="flag to assume Gij = Gji")
    call parse_input_variable(ed_print_Sigma,"ED_PRINT_SIGMA",INPUTunit,default=.true.,comment="flag to print impurity Self-energies")
    call parse_input_variable(ed_print_G,"ED_PRINT_G",INPUTunit,default=.true.,comment="flag to print impurity Greens function")
    call parse_input_variable(ed_print_G0,"ED_PRINT_G0",INPUTunit,default=.true.,comment="flag to print non-interacting impurity Greens function")
    call parse_input_variable(ed_verbose,"ED_VERBOSE",INPUTunit,default=3,comment="Verbosity level: 0=almost nothing --> 5:all. Really: all")
    call parse_input_variable(ed_hw_band,"ED_HW_BAND",INPUTunit,default=2d0,comment="half-bandwidth for the bath initialization: flat in -ed_hw_bath:ed_hw_bath")
    call parse_input_variable(ed_offset_bath,"ED_OFFSET_BATH",INPUTunit,default=1d-1,comment="offset for the initialization of diagonal terms in replica/general bath: -offset:offset")

    call parse_input_variable(nsuccess,"NSUCCESS",INPUTunit,default=1,comment="Number of successive iterations below threshold for convergence")
    call parse_input_variable(Lmats,"LMATS",INPUTunit,default=5000,comment="Number of Matsubara frequencies.")
    call parse_input_variable(Lreal,"LREAL",INPUTunit,default=5000,comment="Number of real-axis frequencies.")
    call parse_input_variable(Ltau,"LTAU",INPUTunit,default=1000,comment="Number of imaginary time points.")
    call parse_input_variable(Lfit,"LFIT",INPUTunit,default=1000,comment="Number of Matsubara frequencies used in the \Chi2 fit.")
    call parse_input_variable(nread,"NREAD",INPUTunit,default=0.d0,comment="Objective density for fixed density calculations.")
    call parse_input_variable(nerr,"NERR",INPUTunit,default=1.d-4,comment="Error threshold for fixed density calculations.")
    call parse_input_variable(ndelta,"NDELTA",INPUTunit,default=0.1d0,comment="Initial step for fixed density calculations.")
    call parse_input_variable(ncoeff,"NCOEFF",INPUTunit,default=1d0,comment="multiplier for the initial ndelta read from a file (ndelta-->ndelta*ncoeff).")    
    call parse_input_variable(wini,"WINI",INPUTunit,default=-5.d0,comment="Smallest real-axis frequency")
    call parse_input_variable(wfin,"WFIN",INPUTunit,default=5.d0,comment="Largest real-axis frequency")
    ![todo] call parse_input_variable(chiflag,"CHIFLAG",INPUTunit,default=.false.,comment="Flag to activate spin susceptibility calculation.")
    call parse_input_variable(hfmode,"HFMODE",INPUTunit,default=.true.,comment="Flag to set the Hartree form of the interaction (n-1/2). see xmu.")
    call parse_input_variable(eps,"EPS",INPUTunit,default=0.01d0,comment="Broadening on the real-axis.")
    call parse_input_variable(cutoff,"CUTOFF",INPUTunit,default=1.d-9,comment="Spectrum cut-off, used to determine the number states to be retained.")
    call parse_input_variable(gs_threshold,"GS_THRESHOLD",INPUTunit,default=1.d-9,comment="Energy threshold for ground state degeneracy loop up")
    !    
    call parse_input_variable(lanc_method,"LANC_METHOD",INPUTunit,default="arpack",comment="select the lanczos method to be used in the determination of the spectrum. ARPACK (default), LANCZOS (T=0 only), DVDSON (no MPI)")
    call parse_input_variable(lanc_nstates_sector,"LANC_NSTATES_SECTOR",INPUTunit,default=2,comment="Initial number of states per sector to be determined.")
    call parse_input_variable(lanc_nstates_total,"LANC_NSTATES_TOTAL",INPUTunit,default=1,comment="Initial number of total states to be determined.")
    call parse_input_variable(lanc_nstates_step,"LANC_NSTATES_STEP",INPUTunit,default=2,comment="Number of states added to the spectrum at each step.")
    call parse_input_variable(lanc_ncv_factor,"LANC_NCV_FACTOR",INPUTunit,default=10,comment="Set the size of the block used in Lanczos-Arpack by multiplying the required Neigen (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_ncv_add,"LANC_NCV_ADD",INPUTunit,default=0,comment="Adds up to the size of the block to prevent it to become too small (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_niter,"LANC_NITER",INPUTunit,default=512,comment="Number of Lanczos iteration in spectrum determination.")
    call parse_input_variable(lanc_ngfiter,"LANC_NGFITER",INPUTunit,default=200,comment="Number of Lanczos iteration in GF determination. Number of momenta.")
    call parse_input_variable(lanc_tolerance,"LANC_TOLERANCE",INPUTunit,default=1d-18,comment="Tolerance for the Lanczos iterations as used in Arpack and plain lanczos.")
    call parse_input_variable(lanc_dim_threshold,"LANC_DIM_THRESHOLD",INPUTunit,default=1024,comment="Min dimension threshold to use Lanczos determination of the spectrum rather than Lapack based exact diagonalization.")
    !
    call parse_input_variable(cg_method,"CG_METHOD",INPUTunit,default=1,comment="Conjugate-Gradient method: 0=NumericalRecipes, 1=minimize.")
    call parse_input_variable(cg_grad,"CG_GRAD",INPUTunit,default=1,comment="Gradient evaluation: 0=analytic (default), 1=numeric.")
    call parse_input_variable(cg_ftol,"CG_FTOL",INPUTunit,default=0.00001d0,comment="Conjugate-Gradient tolerance.")
    call parse_input_variable(cg_stop,"CG_STOP",INPUTunit,default=0,comment="Conjugate-Gradient stopping condition: 0-2, 0=C1.AND.C2, 1=C1, 2=C2 with C1=|F_n-1 -F_n|<tol*(1+F_n), C2=||x_n-1 -x_n||<tol*(1+||x_n||).")
    call parse_input_variable(cg_niter,"CG_NITER",INPUTunit,default=500,comment="Max. number of Conjugate-Gradient iterations.")
    call parse_input_variable(cg_weight,"CG_WEIGHT",INPUTunit,default=1,comment="Conjugate-Gradient weight form, on matsubara axis: 1=1.0, 2=1/n , 3=1/w_n.")
    call parse_input_variable(cg_matrix,"CG_MATRIX",INPUTunit,default=1,comment="Conjugate-Gradient matricial weights: 1=spectral, 0=flat. [WORK IN PROGRESS]")
    call parse_input_variable(cg_scheme,"CG_SCHEME",INPUTunit,default='weiss',comment="Conjugate-Gradient fit scheme: delta or weiss.")
    call parse_input_variable(cg_norm,"CG_NORM",INPUTunit,default='elemental',comment="Conjugate-Gradient norm definition: elemental (default) or frobenius.")
    call parse_input_variable(cg_pow,"CG_POW",INPUTunit,default=2,comment="Fit power for the calculation of the generalized distance as |G0 - G0and|**cg_pow")
    call parse_input_variable(cg_minimize_ver,"CG_MINIMIZE_VER",INPUTunit,default=.false.,comment="Flag to pick old/.false. (Krauth) or new/.true. (Lichtenstein) version of the minimize CG routine")
    call parse_input_variable(cg_minimize_hh,"CG_MINIMIZE_HH",INPUTunit,default=1d-4,comment="Unknown parameter used in the CG minimize procedure.")

    call parse_input_variable(SectorFile,"SectorFile",INPUTunit,default="sectors",comment="File where to retrieve/store the sectors contributing to the spectrum.")
    call parse_input_variable(Hfile,"HFILE",INPUTunit,default="hamiltonian",comment="File where to retrieve/store the bath parameters.")
    call parse_input_variable(HLOCfile,"impHfile",INPUTunit,default="inputHLOC.in",comment="File read the input local H.")
    call parse_input_variable(LOGfile,"LOGFILE",INPUTunit,default=6,comment="LOG unit.")


#ifdef _MPI
    if(present(comm))then
       if(.not.master)then
          LOGfile=1000-rank
          open(LOGfile,file="stdOUT.rank"//str(rank)//".ed")
          do i=1,get_Size_MPI(comm)
             if(i==rank)write(*,"(A,I0,A,I0)")"Rank ",rank," writing to unit: ",LOGfile
          enddo
       endif
    endif
#endif
    !
    !
    Ltau=max(int(beta),Ltau)
    if(master)then
       call print_input()
       call save_input(INPUTunit)
       call scifor_version()
       call code_version(version)
    endif
    !
    if(nread .ne. 0d0) then
       inquire(file="xmu.restart",EXIST=bool)
       if(bool)then
          open(free_unit(unit_xmu),file="xmu.restart")
          read(unit_xmu,*)xmu,ndelta
          ndelta=abs(ndelta)*ncoeff
          close(unit_xmu)
          write(*,"(A,F9.7,A)")"Adjusting XMU to ",xmu," as per provided xmu.restart ",LOGfile
       endif
    endif
    !Act on the input variable only after printing.
    !In the new parser variables are hard-linked into the list:
    !any change to the variable is immediately copied into the list... (if you delete .ed it won't be printed out)
    call substring_delete(Hfile,".restart")
    call substring_delete(Hfile,".ed")
  end subroutine ed_read_input




  subroutine substring_delete (s,sub)
    !! S_S_DELETE2 recursively removes a substring from a string.
    !    The remainder is left justified and padded with blanks.
    !    The substitution is recursive, so
    !    that, for example, removing all occurrences of "ab" from
    !    "aaaaabbbbbQ" results in "Q".
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, character ( len = * ) SUB, the substring to be removed.
    !    Output, integer ( kind = 4 ) IREP, the number of occurrences of
    !    the substring.
    integer          :: ihi
    integer          :: irep
    integer          :: loc
    integer          :: nsub
    character(len=*) ::  s
    integer          :: s_length
    character(len=*) :: sub
    s_length = len ( s )
    nsub = len ( sub )
    irep = 0
    ihi = s_length
    do while ( 0 < ihi )
       loc = index ( s(1:ihi), sub )
       if ( loc == 0 ) then
          return
       end if
       irep = irep + 1
       call s_chop ( s, loc, loc+nsub-1 )
       ihi = ihi - nsub
    end do
    return
  end subroutine substring_delete

  subroutine s_chop ( s, ilo, ihi )
    !! S_CHOP "chops out" a portion of a string, and closes up the hole.
    !  Example:
    !    S = 'Fred is not a jerk!'
    !    call s_chop ( S, 9, 12 )
    !    S = 'Fred is a jerk!    '
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, integer ( kind = 4 ) ILO, IHI, the locations of the first and last
    !    characters to be removed.
    integer               ::ihi
    integer               ::ihi2
    integer               ::ilo
    integer               ::ilo2
    character ( len = * ) :: s
    integer               ::s_length
    s_length = len ( s )
    ilo2 = max ( ilo, 1 )
    ihi2 = min ( ihi, s_length )
    if ( ihi2 < ilo2 ) then
       return
    end if
    s(ilo2:s_length+ilo2-ihi2-1) = s(ihi2+1:s_length)
    s(s_length+ilo2-ihi2:s_length) = ' '
    return
  end subroutine s_chop


END MODULE ED_INPUT_VARS

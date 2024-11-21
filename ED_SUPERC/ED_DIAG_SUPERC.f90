!########################################################################
!PURPOSE  : Diagonalize the Effective Impurity Problem
!|{ImpUP1,...,ImpUPN},BathUP>|{ImpDW1,...,ImpDWN},BathDW>
!########################################################################
module ED_DIAG_SUPERC
  USE SF_CONSTANTS
  USE SF_LINALG, only: eigh
  USE SF_TIMER,  only: start_timer,stop_timer,eta
  USE SF_IOTOOLS, only:reg,free_unit,file_length
  USE SF_STAT
  USE SF_SP_LINALG
  !
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_AUX_FUNX
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN_SUPERC
  implicit none
  private


  public :: diagonalize_impurity_superc


contains



  !+-------------------------------------------------------------------+
  !PURPOSE  : Setup the Hilbert space, create the Hamiltonian, get the
  ! GS, build the Green's functions calling all the necessary routines
  !+------------------------------------------------------------------+
  subroutine diagonalize_impurity_superc()
    !
    ! This procedure performs the diagonalization of the Hamiltonian in each symmetry sector for :f:var:`ed_mode` = :f:var:`superc` , i.e. for quantum number :math:`\vec{Q}=S_{z}`.
    !
    ! The diagonalization proceeds in three steps.
    !
    ! #. Setup the diagonalization of selected sectors only if the file :f:var:`sectorfile` exists.
    !
    ! #. Perform a cycle over all the symmetry sectors :math:`{\cal S}(\vec{Q})`
    !
    !   * Decide if the sector should be diagonalized, i.e. if this is a twin sector and :f:var:`ed_twin` is True.
    !
    !   * Setup the Arpack/Lanczos parameters :f:var:`neigen` , :f:var:`nblock` , :f:var:`nitermax`
    !
    !   * Construct the sector Hamiltonian and/or just associate the correct matrix-vector product in :f:var:`build_hv_sector_normal`
    !
    !   * call :f:var:`sp_eigh` (Arpack) or :f:var:`sp_lanc_eigh` (Lanczos) or :f:var:`eigh` (Lapack) according to the dimension of the sector and the value of :f:var:`lanc_method`
    !
    !   * Retrieve the eigen-states and save them to :f:var:`state_list`
    !
    ! #. Analyze the :f:var:`state_list` , find the overall groundstate :math:`|\psi_0\rangle` , trim the list using the :f:var:`cutoff` :math:`\epsilon` to that :math:`e^{-\beta E_{max}} < \epsilon` , create an histogram of the states requested to each sector and use it to increase or decrease those number according to the contribution of each sector to the spectrum.
    !
    call ed_pre_diag
    call ed_diag_c
    call ed_post_diag
  end subroutine diagonalize_impurity_superc






  !+-------------------------------------------------------------------+
  !PURPOSE  : diagonalize the Hamiltonian in each sector and find the 
  ! spectrum DOUBLE PRECISION
  !+------------------------------------------------------------------+
  subroutine ed_diag_d
    integer                :: isector,Dim,istate
    integer                :: isect,izero,sz,nt
    integer                :: i,j,iter,unit,vecDim
    integer                :: Nitermax,Neigen,Nblock
    real(8)                :: oldzero,enemin,Ei
    real(8),allocatable    :: eig_values(:)
    complex(8),allocatable :: eig_basis(:,:),eig_basis_tmp(:,:)
    logical                :: lanc_solve,Tflag,lanc_verbose,bool
    !
    if(state_list%status)call es_delete_espace(state_list)
    state_list=es_init_espace()
    oldzero=1000.d0
    if(MPIMASTER)then
       write(LOGfile,"(A)")"Diagonalize impurity H:"
       call start_timer()
    endif
    !
    lanc_verbose=.false.
    if(ed_verbose>2)lanc_verbose=.true.
    !
    iter=0
    sector: do isector=1,Nsectors
       if(.not.sectors_mask(isector))cycle sector
       if(.not.twin_mask(isector))cycle sector !cycle loop if this sector should not be investigated
       iter=iter+1
       Tflag = twin_mask(isector).AND.ed_twin
       Tflag = Tflag.AND.(getsz(isector)/=0)
       !
       Dim      = getdim(isector)
       !
       select case (lanc_method)
       case default       !use P-ARPACK
          Neigen   = min(dim,neigen_sector(isector))
          Nitermax = min(dim,lanc_niter)
          Nblock   = min(dim,lanc_ncv_factor*max(Neigen,lanc_nstates_sector) + lanc_ncv_add)
       case ("lanczos")
          Neigen   = 1
          Nitermax = min(dim,lanc_niter)
          Nblock   = 1
       end select
       !
       lanc_solve  = .true.
       if(Neigen==dim)lanc_solve=.false.
       if(dim<=max(lanc_dim_threshold,MPISIZE))lanc_solve=.false.
       !
       if(MPIMASTER)then
          if(ed_verbose>2)then
             sz   = getsz(isector)
             if(lanc_solve)then
                write(LOGfile,"(1X,I9,A,I9,A6,I4,A6,I20,A12,3I6)")&
                     iter,"-Solving sector:",isector," sz:",sz," dim=",&
                     getdim(isector),", Lanc Info:",Neigen,Nitermax,Nblock
             else
                write(LOGfile,"(1X,I9,A,I9,A6,I4,A6,I20)")&
                     iter,"-Solving sector:",isector," sz:",sz," dim=",&
                     getdim(isector)
             endif
          elseif(ed_verbose<=2)then
             call eta(iter,count(twin_mask),LOGfile)
          endif
       endif
       !
       !
       if(allocated(eig_values))deallocate(eig_values)
       if(allocated(eig_basis))deallocate(eig_basis)
       !
       if(ed_verbose>=3.AND.MPIMASTER)call start_timer()
       if(lanc_solve)then
          !
          allocate(eig_values(Neigen)) ; eig_values=0d0
          !
          call build_Hv_sector(isector)
          !
          vecDim = vecDim_Hv_sector_superc(isector)
          allocate(eig_basis(vecDim,Neigen))
          eig_basis=zero
          !
          select case (lanc_method)
          case default       !use P-ARPACK
#ifdef _MPI
             if(MpiStatus)then
                call sp_eigh(MpiComm,spHtimesV_p,eig_values,eig_basis,&
                     Nblock,&
                     Nitermax,&
                     tol=lanc_tolerance,&
                     iverbose=(ed_verbose>4))
             else
                call sp_eigh(spHtimesV_p,eig_values,eig_basis,&
                     Nblock,&
                     Nitermax,&
                     tol=lanc_tolerance,&
                     iverbose=(ed_verbose>4))
             endif
#else
             call sp_eigh(spHtimesV_p,eig_values,eig_basis,&
                  Nblock,&
                  Nitermax,&
                  tol=lanc_tolerance,&
                  iverbose=(ed_verbose>4))
#endif             
             !
             !
          case ("lanczos")   !use Simple Lanczos
#ifdef _MPI
             if(MpiStatus)then
                call sp_lanc_eigh(MpiComm,spHtimesV_p,eig_values(1),eig_basis(:,1),Nitermax,&
                     iverbose=(ed_verbose>4),threshold=lanc_tolerance)
             else
                call sp_lanc_eigh(spHtimesV_p,eig_values(1),eig_basis(:,1),Nitermax,&
                     iverbose=(ed_verbose>4),threshold=lanc_tolerance)
             endif
#else
             call sp_lanc_eigh(spHtimesV_p,eig_values(1),eig_basis(:,1),Nitermax,&
                  iverbose=(ed_verbose>4),threshold=lanc_tolerance)
#endif
          end select
          !
          !
          if(MpiMaster.AND.ed_verbose>3)write(LOGfile,*)""
          call delete_Hv_sector_superc()
#ifdef _MPI
          if(MpiStatus)call Bcast_MPI(MpiComm,eig_values)
#endif
          !
          !
       else                     !else LAPACK_SOLVE
          !
          !
          allocate(eig_values(Dim)) ; eig_values=0d0
          allocate(eig_basis_tmp(Dim,Dim)) ; eig_basis_tmp=zero
          call build_Hv_sector_superc(isector,eig_basis_tmp)
          !
          if(MpiMaster)call eigh(eig_basis_tmp,eig_values)
          if(dim==1)eig_basis_tmp(dim,dim)=one
          !
          call delete_Hv_sector_superc()
#ifdef _MPI
          if(MpiStatus)then
             call Bcast_MPI(MpiComm,eig_values)
             vecDim = vecDim_Hv_sector(isector)
             allocate(eig_basis(vecDim,Neigen)) ; eig_basis=zero
             call scatter_basis_MPI(MpiComm,eig_basis_tmp,eig_basis)
          else
             allocate(eig_basis(Dim,Neigen)) ; eig_basis=zero
             eig_basis = eig_basis_tmp(:,1:Neigen)
          endif
#else
          allocate(eig_basis(Dim,Neigen)) ; eig_basis=zero
          eig_basis = eig_basis_tmp(:,1:Neigen)
#endif
          !
       endif
       !
       if(ed_verbose>=3.AND.MPIMASTER)call stop_timer
       if(ed_verbose>=4)write(LOGfile,*)"EigValues: ",eig_values(:Neigen)
       if(ed_verbose>2)write(LOGfile,*)""
       if(ed_verbose>2)write(LOGfile,*)""
       !
       if(finiteT)then
          do i=1,Neigen
             call es_add_state(state_list,eig_values(i),eig_basis(:,i),isector,twin=Tflag,size=lanc_nstates_total)
          enddo
       else
          do i=1,Neigen
             enemin = eig_values(i)
             if (enemin < oldzero-10.d0*gs_threshold)then
                oldzero=enemin
                call es_free_espace(state_list)
                call es_add_state(state_list,enemin,eig_basis(:,i),isector,twin=Tflag)
             elseif(abs(enemin-oldzero) <= gs_threshold)then
                oldzero=min(oldzero,enemin)
                call es_add_state(state_list,enemin,eig_basis(:,i),isector,twin=Tflag)
             endif
          enddo
       endif
       !
       if(MPIMASTER)then
          unit=free_unit()
          open(unit,file="eigenvalues_list"//reg(ed_file_suffix)//".ed",position='append',action='write')
          call print_eigenvalues_list(isector,eig_values(1:Neigen),unit,lanc_solve,mpiAllThreads)
          close(unit)
       endif
       !
       if(allocated(eig_values))deallocate(eig_values)
       if(allocated(eig_basis_tmp))deallocate(eig_basis_tmp)
       if(allocated(eig_basis))deallocate(eig_basis)
       !
    enddo sector
    if(MPIMASTER)call stop_timer()
  end subroutine ed_diag_d










  !###################################################################################################
  !
  !    PRE-PROCESSING ROUTINES
  !
  !###################################################################################################
  subroutine ed_pre_diag
    integer                          :: Nup,Ndw,Mup,Mdw
    integer                          :: Sz,N,Rz,M,twoJz,twoIz
    integer                          :: i,iud,iorb
    integer                          :: isector,jsector
    integer                          :: unit,unit2,status,istate,ishift,isign
    logical                          :: IOfile
    integer                          :: list_len
    integer,dimension(:),allocatable :: list_sector
    !
    sectors_mask=.true.
    !
    if(ed_sectors)then
       inquire(file="sectors_list"//reg(ed_file_suffix)//".restart",exist=IOfile)
       if(IOfile)then
          sectors_mask=.false.
          write(LOGfile,"(A)")"Analysing sectors_list to reduce sectors scan:"
          list_len=file_length("sectors_list"//reg(ed_file_suffix)//".restart")
          !
          open(free_unit(unit),file=reg("sectors_list"//reg(ed_file_suffix)//".restart",status="old")
          open(free_unit(unit2),file="list_of_sectors"//reg(ed_file_suffix)//".ed")
          do istate=1,list_len
             read(unit,*,iostat=status)sz
             isector = getSector(Sz,1)
             sectors_mask(isector)=.true.
             write(unit2,*)isector,sectors_mask(isector),sz
             !
             do ishift=1,ed_sectors_shift
                do isign=-1,1,2
                   Rz = Sz + isign*ishift
                   jsector = getSector(Rz,1)
                   sectors_mask(jsector)=.true.
                   write(unit2,*)jsector,sectors_mask(jsector),Rz
                enddo
             enddo
          enddo
          close(unit)
          close(unit2)
          !
       endif
    endif
    !
  end subroutine ed_pre_diag




  !###################################################################################################
  !
  !    POST-PROCESSING ROUTINES
  !
  !###################################################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : analyse the spectrum and print some information after 
  !lanczos diagonalization. 
  !+------------------------------------------------------------------+
  subroutine ed_post_diag()
    integer             :: nup,ndw,sz,n,isector,dim
    integer             :: istate
    integer             :: i,unit
    integer             :: nups(Ns_Ud),ndws(Ns_Ud),Indices(2*Ns_Ud)
    integer             :: Nsize,NtoBremoved,nstates_below_cutoff
    integer             :: numgs
    real(8)             :: Egs,Ei,Ec,Etmp
    type(histogram)     :: hist
    real(8)             :: hist_a,hist_b,hist_w
    integer             :: hist_n
    integer,allocatable :: list_sector(:),count_sector(:)    
    !POST PROCESSING:
    if(MPIMASTER)then
       open(free_unit(unit),file="state_list"//reg(ed_file_suffix)//".ed")
       call save_state_list(unit)
       close(unit)
    endif
    if(ed_verbose>=2)call print_state_list(LOGfile)
    !
    zeta_function=0d0
    Egs = state_list%emin
    if(finiteT)then
       do i=1,state_list%size
          ei   = es_return_energy(state_list,i)
          zeta_function = zeta_function + exp(-beta*(Ei-Egs))
       enddo
    else
       zeta_function=real(state_list%size,8)
    end if
    !
    !
    numgs=es_return_gs_degeneracy(state_list,gs_threshold)
    if(numgs>Nsectors)stop "ed_diag: too many gs"
    if(MPIMASTER.AND.ed_verbose>=2)then
       do istate=1,numgs
          isector = es_return_sector(state_list,istate)
          Egs     = es_return_energy(state_list,istate)
          sz  = getsz(isector)
          write(LOGfile,"(A,F20.12,I4)")'Egs =',Egs,sz
       enddo
       write(LOGfile,"(A,F20.12)")'Z   =',zeta_function
    endif
    !
    !
    !
    if(.not.finiteT)then
       !generate a sector_list to be reused in case we want to reduce sectors scan
       open(free_unit(unit),file="sectors_list"//reg(ed_file_suffix)//".restart")       
       do istate=1,state_list%size
          isector = es_return_sector(state_list,istate)
          sz  = getsz(isector)
          write(unit,*)sz
       enddo
       close(unit)
    else
       !get histogram distribution of the sector contributing to the evaluated spectrum:
       !go through states list and update the neigen_sector(isector) sector-by-sector
       if(MPIMASTER)then
          unit=free_unit()
          open(unit,file="histogram_states"//reg(ed_file_suffix)//".ed",position='append')
          hist_n = Nsectors
          hist_a = 1d0
          hist_b = dble(Nsectors)
          hist_w = 1d0
          hist = histogram_allocate(hist_n)
          call histogram_set_range_uniform(hist,hist_a,hist_b)
          do i=1,state_list%size
             isector = es_return_sector(state_list,i)
             call histogram_accumulate(hist,dble(isector),hist_w)
          enddo
          call histogram_print(hist,unit)
          write(unit,*)""
          close(unit)
       endif
       !
       !
       !
       allocate(list_sector(state_list%size),count_sector(Nsectors))
       !get the list of actual sectors contributing to the list
       do i=1,state_list%size
          list_sector(i) = es_return_sector(state_list,i)
       enddo
       !count how many times a sector appears in the list
       do i=1,Nsectors
          count_sector(i) = count(list_sector==i)
       enddo
       !adapt the number of required Neig for each sector based on how many
       !appeared in the list.
       do i=1,Nsectors
          if(any(list_sector==i))then !if(count_sector(i)>1)then
             neigen_sector(i)=neigen_sector(i)+1
          else
             neigen_sector(i)=neigen_sector(i)-1
          endif
          !prevent Neig(i) from growing unbounded but 
          !try to put another state in the list from sector i
          if(neigen_sector(i) > count_sector(i))neigen_sector(i)=count_sector(i)+1
          if(neigen_sector(i) <= 0)neigen_sector(i)=1
       enddo
       !check if the number of states is enough to reach the required accuracy:
       !the condition to fullfill is:
       ! exp(-beta(Ec-Egs)) < \epsilon_c
       ! if this condition is violated then required number of states is increased
       ! if number of states is larger than those required to fullfill the cutoff: 
       ! trim the list and number of states.
       Egs  = state_list%emin
       Ec   = state_list%emax
       Nsize= state_list%size
       if(exp(-beta*(Ec-Egs)) > cutoff)then
          lanc_nstates_total=lanc_nstates_total + lanc_nstates_step
          if(MPIMASTER)write(LOGfile,"(A,I4)")"Increasing lanc_nstates_total:",lanc_nstates_total
       else
          ! !Find the energy level beyond which cutoff condition is verified & cut the list to that size
          write(LOGfile,*)
          isector = es_return_sector(state_list,state_list%size)
          Ei      = es_return_energy(state_list,state_list%size)
          do while ( exp(-beta*(Ei-Egs)) <= cutoff )
             if(ed_verbose>=1.AND.MPIMASTER)write(LOGfile,"(A,I4,I5)")"Trimming state:",isector,state_list%size
             call es_pop_state(state_list)
             isector = es_return_sector(state_list,state_list%size)
             Ei      = es_return_energy(state_list,state_list%size)
          enddo
          if(ed_verbose>=1.AND.MPIMASTER)then
             write(LOGfile,*)"Trimmed state list:"          
             call print_state_list(LOGfile)
          endif
          !
          lanc_nstates_total=max(state_list%size,lanc_nstates_step)+lanc_nstates_step
          write(LOGfile,"(A,I4)")"Adjusting lanc_nstates_total to:",lanc_nstates_total
          !
       endif
    endif
  end subroutine ed_post_diag


  subroutine print_state_list(unit)
    integer :: nup,ndw,sz,n,isector
    integer :: istate
    integer :: unit
    real(8) :: Estate
    if(MPIMASTER)then
       write(unit,"(A)")"# i       E_i           exp(-(E-E0)/T)       Sz     Sect     Dim"
       do istate=1,state_list%size
          Estate  = es_return_energy(state_list,istate)
          isector = es_return_sector(state_list,istate)
          sz   = getsz(isector)
          write(unit,"(i6,f18.12,2x,ES19.12,1x,i3,3x,i3,i10)")&
               istate,Estate,exp(-beta*(Estate-state_list%emin)),sz,isector,getdim(isector)
       enddo
    endif
  end subroutine print_state_list


  subroutine save_state_list(unit)
    integer :: indices(2*Ns_Ud),isector
    integer :: istate
    integer :: unit
    if(MPIMASTER)then
       do istate=1,state_list%size
          isector = es_return_sector(state_list,istate)
          sz   = getsz(isector)
          write(unit,"(i8,i12,i8)")istate,isector,sz
       enddo
    endif
  end subroutine save_state_list


  subroutine print_eigenvalues_list(isector,eig_values,unit,lanc,allt)    
    integer              :: isector
    real(8),dimension(:) :: eig_values
    integer              :: unit,i
    logical              :: lanc,allt
    if(MPIMASTER)then
       if(lanc)then
          if(allt)then
             write(unit,"(A9,A4)")" # Sector","Sz"
          else
             write(unit,"(A10,A4)")" #T Sector","Sz"
          endif
       else
          write(unit,"(A10,A4)")" #X Sector","Sz"
       endif
       write(unit,"(I9,I4)")isector,getsz(isector)
       do i=1,size(eig_values)
          write(unit,*)eig_values(i)
       enddo
       write(unit,*)""
    endif
  end subroutine print_eigenvalues_list





end MODULE ED_DIAG_SUPERC










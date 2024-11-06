MODULE ED_SETUP
  !
  !Contains procedures to set up the Exact Diagonalization calculation, executing all internal consistency checks and allocation of the global memory.
  !
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_TIMER
  USE SF_IOTOOLS, only: free_unit,reg,file_length
  USE SF_MISC,    only: assert_shape
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private

  public :: init_ed_structure
  public :: delete_ed_structure
  public :: setup_global
  public :: get_normal_sector_dimension
  public :: get_superc_sector_dimension


contains

  subroutine ed_checks_global
    !
    if(Lfit>Lmats)Lfit=Lmats
    if(Nspin>2)stop "ED ERROR: Nspin > 2 is currently not supported"
    if(Norb>5)stop "ED ERROR: Norb > 5 is currently not supported"
    !
    if(ed_mode=="superc")then
       if(Nspin>1)stop "ED ERROR: SC + Magnetism is not yet supported"
       ! if(bath_type=="replica")stop "ED ERROR: ed_mode=SUPERC + bath_type=replica is not supported"
    endif
    !
    if(Nspin>1.AND.ed_twin.eqv..true.)then
       write(LOGfile,"(A)")"WARNING: using twin_sector with Nspin>1"
       call sleep(1)
    end if
    !
    if(lanc_method=="lanczos")then
       if(lanc_nstates_total>1)stop "ED ERROR: lanc_method==lanczos available only for lanc_nstates_total==1, T=0"
       if(lanc_nstates_sector>1)stop "ED ERROR: lanc_method==lanczos available only for lanc_nstates_sector==1, T=0"
    endif
    if(ed_finite_temp)then
       if(lanc_nstates_total==1)stop "ED ERROR: ed_diag_type==lanc + ed_finite_temp=T *but* lanc_nstates_total==1 => T=0. Increase lanc_nstates_total"
    else
       if(lanc_nstates_total>1)print*, "ED WARNING: ed_diag_type==lanc + ed_finite_temp=F, T=0 *AND* lanc_nstates_total>1. re-Set lanc_nstates_total=1"
       lanc_nstates_total=1
    endif
    !
  end subroutine ed_checks_global


  !+------------------------------------------------------------------+
  !PURPOSE  : Setup Dimensions of the problem
  ! Norb    = # of impurity orbitals
  ! Nbath   = # of bath levels (depending on bath_type)
  ! Ns      = # of levels (per spin)
  ! Nlevels = 2*Ns = Total # of levels (counting spin degeneracy 2) 
  !+------------------------------------------------------------------+
  subroutine ed_setup_dimensions()
    !
    Nimp = Nlat*Norb    !Total number of levels in the impurity cluster
    Ns = Nimp*(Nbath+1) !Total number of levels per spin
    !
    Ns_Orb = Ns
    Ns_Ud  = 1
    !
    select case(ed_mode)
    case default
       Nsectors = ((Ns_Orb+1)*(Ns_Orb+1))**Ns_Ud
    case ("superc")
       Nsectors = Nlevels+1     !sz=-Ns:Ns=2*Ns+1=Nlevels+1
    end select

  end subroutine ed_setup_dimensions



  !+------------------------------------------------------------------+
  !PURPOSE  : Init ED structure and calculation
  !+------------------------------------------------------------------+
  subroutine init_ed_structure()
    ! Initialize the pool of variables and data structures of the ED calculation.
    ! Performs all the checks calling :f:func:`ed_checks_global`, set up the dimensions in :f:func:`ed_setup_dimensions` given the variables :f:var:`ns`, :f:var:`norb`, :f:var:`nspin`, :f:var:`nbath`, :f:var:`bath_type`. Allocate all the dynamic memory which will be stored in the memory till the calculation will be finalized. 
    logical                          :: control
    integer                          :: i,iud,iorb,jorb,ispin,jspin
    integer,dimension(:),allocatable :: DimUps,DimDws
    !
    call ed_checks_global
    !
    call ed_setup_dimensions
    !
    !

    select case(ed_mode)
    case default
       allocate(DimUps(Ns_Ud))
       allocate(DimDws(Ns_Ud))
       do iud=1,Ns_Ud
          DimUps(iud) = get_normal_sector_dimension(Ns_Orb,Ns_Orb/2)
          DimDws(iud) = get_normal_sector_dimension(Ns_Orb,Ns_Orb-Ns_Orb/2)
       enddo
    case ("superc")
       dim_sector_max=get_superc_sector_dimension(0)
    end select

    if(MpiMaster)then
       write(LOGfile,"(A)")"Summary:"
       write(LOGfile,"(A)")"--------------------------------------------"
       write(LOGfile,"(A,I15)")'# of levels/spin      = ',Ns
       write(LOGfile,"(A,I15)")'Total size            = ',2*Ns
       write(LOGfile,"(A,I15)")'# of sites            = ',Nlat
       write(LOGfile,"(A,I15)")'# of impurities       = ',Norb
       write(LOGfile,"(A,I15)")'# of bath/impurity    = ',Nbath
       write(LOGfile,"(A,I15)")'Number of sectors     = ',Nsectors
       select case(ed_mode)
       case default
          write(LOGfile,"(A,"//str(Ns_Ud)//"I8,2X,"//str(Ns_Ud)//"I8,I20)")&
               'Largest Sector(s)     = ',DimUps,DimDws,product(DimUps)*product(DimDws)*DimPh
       case("superc","nonsu2")
          write(LOGfile,"(A,I15)")'Largest Sector(s)    = ',dim_sector_max
       end select
       write(LOGfile,"(A)")"--------------------------------------------"
    endif
    !
    !
    allocate(spH0ups(Ns_Ud))
    allocate(spH0dws(Ns_Ud))
    !
    !Allocate indexing arrays
    allocate(getCsector(Ns_Ud,2,Nsectors))  ;getCsector  =0
    allocate(getCDGsector(Ns_Ud,2,Nsectors));getCDGsector=0
    !
    select case(ed_mode)
    case default
       allocate(getSector(0,0))
    case ("superc")
       allocate(getSector(-Ns:Ns,1))
    end select
    getSector=0

    allocate(getDim(Nsectors));getDim=0
    allocate(getSz(Nsectors));getSz=0
    !
    allocate(getBathStride(Nlat,Norb,Nbath));getBathStride=0
    allocate(twin_mask(Nsectors))
    allocate(sectors_mask(Nsectors))
    allocate(neigen_sector(Nsectors))
    !
    finiteT = ed_finite_temp
    if(lanc_nstates_total==1)then     !is you only want to keep 1 state
       finiteT=.false.                !set to do zero temperature calculations
       write(LOGfile,"(A)")"Required Lanc_nstates_total=1 => set T=0 calculation"
    endif
    !
    !check whether lanc_nstates_sector and lanc_states are even (we do want to keep doublet among states)
    if(finiteT)then
       if(mod(lanc_nstates_sector,2)/=0)then
          lanc_nstates_sector=lanc_nstates_sector+1
          write(LOGfile,"(A,I10)")"Increased Lanc_nstates_sector:",lanc_nstates_sector
       endif
       if(mod(lanc_nstates_total,2)/=0)then
          lanc_nstates_total=lanc_nstates_total+1
          write(LOGfile,"(A,I10)")"Increased Lanc_nstates_total:",lanc_nstates_total
       endif
       write(LOGfile,"(A)")"Lanczos FINITE temperature calculation:"
       !
       write(LOGfile,"(A,I3)")"Nstates x Sector = ", lanc_nstates_sector
       write(LOGfile,"(A,I3)")"Nstates   Total  = ", lanc_nstates_total
       call sleep(1)
    else
       write(LOGfile,"(A)")"Lanczos ZERO temperature calculation:"
       call sleep(1)
    endif
    !
    Jhflag=.FALSE.
    if(Norb>1.AND.(Jx/=0d0.OR.Jp/=0d0))Jhflag=.TRUE.
    !
    !
    offdiag_gf_flag=.true.
    !
    if(nread/=0.d0)then
       i=abs(floor(log10(abs(nerr)))) !modulus of the order of magnitude of nerror
       niter=nloop/3
    endif

    !ALLOCATE impHloc
    if(.not.allocated(impHloc))then
       allocate(impHloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
       impHloc=zero
    else
       call assert_shape(impHloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"init_ed_structure","impHloc")
    endif
    !
    !
    !allocate functions
    allocate(impSmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impSreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(impSAmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impSAreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    impSAmats=zero
    impSAreal=zero
    impSmats=zero
    impSreal=zero
    !
    allocate(impGmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impGreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(impFmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impFreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    impGmats=zero
    impGreal=zero
    impFmats=zero
    impFreal=zero
    !
    !allocate functions
    allocate(impG0mats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impG0real(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(impF0mats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impF0real(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    impG0mats=zero
    impG0real=zero
    impF0mats=zero
    impF0real=zero
    !
    allocate(impGmatrix(Nlat,Nlat,Nnambu*Nspin,Nnambu*Nspin,Norb,Norb))
    !
    !allocate observables
    allocate(ed_dens(Nlat,Norb))
    allocate(ed_docc(Nlat,Norb))
    allocate(ed_mag(Nlat,Norb))
    allocate(ed_phisc(Nlat,Norb))
    allocate(ed_dens_up(Nlat,Norb),ed_dens_dw(Nlat,Norb))
    ed_dens=0d0
    ed_docc=0d0
    ed_mag=0d0
    ed_phisc=0d0
    ed_dens_up=0d0
    ed_dens_dw=0d0
    !
    !
    allocate(single_particle_density_matrix(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    single_particle_density_matrix=zero
    !
    allocate(cluster_density_matrix(4**Nimp,4**Nimp))
    cluster_density_matrix=zero
    !
  end subroutine init_ed_structure



  !+------------------------------------------------------------------+
  !PURPOSE  : Deallocate ED structure and reset environment
  !+------------------------------------------------------------------+
  subroutine delete_ed_structure()
    ! Delete the entire memory pool upon finalization of the ED calculation. 
    logical                          :: control
    integer                          :: i,iud
    integer                          :: dim_sector_max,iorb,jorb,ispin,jspin
    integer,dimension(:),allocatable :: DimUps,DimDws
    !
    Ns       = 0
    Ns_Orb   = 0
    Ns_Ud    = 0
    Nsectors = 0
    !
    if(MpiMaster)write(LOGfile,"(A)")"Cleaning ED structure"
    if(allocated(spH0ups))deallocate(spH0ups)
    if(allocated(spH0dws))deallocate(spH0dws)
    if(allocated(getCsector))deallocate(getCsector)
    if(allocated(getCDGsector))deallocate(getCDGsector)    !
    if(allocated(getSector))deallocate(getSector)
    if(allocated(getDim))deallocate(getDim)
    if(allocated(getSz))deallocate(getSz)
    if(allocated(getN))deallocate(getN)
    if(allocated(getBathStride))deallocate(getBathStride)
    if(allocated(twin_mask))deallocate(twin_mask)
    if(allocated(sectors_mask))deallocate(sectors_mask)
    if(allocated(neigen_sector))deallocate(neigen_sector)
    if(allocated(impHloc))deallocate(impHloc)
    if(allocated(impSmats))deallocate(impSmats)
    if(allocated(impSreal))deallocate(impSreal)
    if(allocated(impSAmats))deallocate(impSAmats)
    if(allocated(impSAreal))deallocate(impSAreal)
    if(allocated(impGmats))deallocate(impGmats)
    if(allocated(impGreal))deallocate(impGreal)
    if(allocated(impFmats))deallocate(impFmats)
    if(allocated(impFreal))deallocate(impFreal)
    if(allocated(impG0mats))deallocate(impG0mats)
    if(allocated(impG0real))deallocate(impG0real)
    if(allocated(impF0mats))deallocate(impF0mats)
    if(allocated(impF0real))deallocate(impF0real)
    if(allocated(impGmatrix))deallocate(impGmatrix)
    if(allocated(ed_dens))deallocate(ed_dens)
    if(allocated(ed_docc))deallocate(ed_docc)
    if(allocated(ed_phisc))deallocate(ed_phisc)
    if(allocated(ed_dens_up))deallocate(ed_dens_up)
    if(allocated(ed_dens_dw))deallocate(ed_dens_dw)
    if(allocated(ed_mag))deallocate(ed_mag)
  end subroutine delete_ed_structure





  !+------------------------------------------------------------------+
  !PURPOSE: SETUP THE GLOBAL POINTERS FOR THE ED CALCULAIONS.
  !+------------------------------------------------------------------+
  subroutine setup_global
    ! Setup the all the dimensions and the local maps according to a given symmetry of the Hamiltonian problem calling the correct procedure for a given :f:var:`ed_mode`.
    !
    ! Setup the local Fock space maps used in the ED calculation for the normal operative mode.
    ! All sectors dimensions, quantum numbers :math:`\{\vec{N_\uparrow},\vec{N_\downarrow}\}`, :math:`S_z`, :math:`N_{tot}`, 
    ! twin sectors and list of requested eigensolutions for each sectors are defined here.
    ! Identify Bath positions stride for a given value of :f:var:`bath_type`.
    ! Determines the sector indices for :math:`\pm` 1-particle with either spin orientations.
    select case(ed_mode)
    case default
       call setup_global_normal()
    case ("superc")
       call setup_global_superc()
    end select
  end subroutine setup_global


  !+------------------------------------------------------------------+
  !PURPOSE: SETUP THE GLOBAL POINTERS FOR THE ED CALCULAIONS.
  !+------------------------------------------------------------------+
  subroutine setup_global_normal
    integer                          :: DimUp,DimDw
    integer                          :: DimUps(Ns_Ud),DimDws(Ns_Ud)
    integer                          :: Indices(2*Ns_Ud),Jndices(2*Ns_Ud)
    integer                          :: Nups(Ns_ud),Ndws(Ns_ud)
    integer                          :: Jups(Ns_ud),Jdws(Ns_ud)
    integer                          :: i,iud,iorb,ilat,stride,ibath
    integer                          :: isector,jsector
    integer                          :: unit,status,istate,ishift,isign
    logical                          :: IOfile
    integer                          :: list_len
    integer,dimension(:),allocatable :: list_sector
    !
    !Store full dimension of the sectors:
    do isector=1,Nsectors
       call get_DimUp(isector,DimUps)
       call get_DimDw(isector,DimDws)
       DimUp = product(DimUps)
       DimDw = product(DimDws)       
       getDim(isector)  = DimUp*DimDw
    enddo
    !
    !
    inquire(file="state_list"//reg(ed_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       write(LOGfile,"(A)")"Restarting from a state_list file:"
       list_len=file_length("state_list"//reg(ed_file_suffix)//".restart")
       allocate(list_sector(list_len))
       !
       open(free_unit(unit),file="state_list"//reg(ed_file_suffix)//".restart",status="old")
       status=0
       do while(status>=0)
          read(unit,*,iostat=status)istate,isector,indices
          list_sector(istate)=isector
          call get_Nup(isector,Nups)
          call get_Ndw(isector,Ndws)
          if(any(Indices /= [Nups,Ndws]))&
               stop "setup_global error: nups!=nups(isector).OR.ndws!=ndws(isector)"
       enddo
       close(unit)
       !
       lanc_nstates_total = list_len
       do isector=1,Nsectors
          neigen_sector(isector) = max(1,count(list_sector==isector))
       enddo
    else
       do isector=1,Nsectors
          neigen_sector(isector) = min(getDim(isector),lanc_nstates_sector)
       enddo
    endif
    !
    twin_mask=.true.
    if(ed_twin)then
       do isector=1,Nsectors
          call get_Nup(isector,Nups)
          call get_Ndw(isector,Ndws)
          if(any(Nups .ne. Ndws))then
             call get_Sector([Ndws,Nups],Ns_Orb,jsector)
             if (twin_mask(jsector))twin_mask(isector)=.false.
          endif
       enddo
       write(LOGfile,"(A,I6,A,I9)")"Looking into ",count(twin_mask)," sectors out of ",Nsectors
    endif
    !
    stride=Nlat*Norb
    do ibath=1,Nbath
       do ilat=1,Nlat
          do iorb=1,Norb
             getBathStride(ilat,iorb,ibath) = stride + imp_state_index(ilat,iorb) + (ibath-1)*Nlat*Norb
          enddo
       enddo
    enddo
    !
    getCsector  = 0
    getCDGsector= 0
    do isector=1,Nsectors
       call get_Nup(isector,Nups)
       call get_Ndw(isector,Ndws)
       !
       !UPs:
       !DEL:
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws 
          Jups(iud)=Jups(iud)-1; if(Jups(iud) < 0)cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCsector(iud,1,isector)=jsector
       enddo
       !ADD
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws
          Jups(iud)=Jups(iud)+1; if(Jups(iud) > Ns_Orb)cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCDGsector(iud,1,isector)=jsector
       enddo
       !
       !DWs:
       !DEL
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws 
          Jdws(iud)=Jdws(iud)-1; if(Jdws(iud) < 0)cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCsector(iud,2,isector)=jsector
       enddo
       !DEL
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws 
          Jdws(iud)=Jdws(iud)+1; if(Jdws(iud) > Ns_Orb)cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCDGsector(iud,2,isector)=jsector
       enddo
    enddo
    return
  end subroutine setup_global_normal





  !SUPERCONDUCTING
  subroutine setup_global_superc
    !Setup the local Fock space maps used in the ED calculation for the **superc** operative mode. All sectors dimensions, quantum numbers, twin sectors and list of requested eigensolutions for each sectors are defined here. Identify Bath positions stride for a given value of :code:`bath_type`. Determines the sector indices for :math:`\pm` -particle with :math:`\sigma=\uparrow,\downarrow`.
    integer                                           :: i,isz,in,dim,isector,jsector
    integer                                           :: sz,iorb,jsz
    integer                                           :: unit,status,istate
    logical                                           :: IOfile
    integer                                           :: anint
    real(8)                                           :: adouble
    integer                                           :: list_len
    integer,dimension(:),allocatable                  :: list_sector
    !
    isector=0
    do isz=-Ns,Ns
       isector=isector+1
       sz=abs(isz)
       getSector(isz,1)=isector
       getSz(isector)  = isz
       dim = get_superc_sector_dimension(isz)
       getDim(isector) = dim
    enddo
    !
    inquire(file="state_list"//reg(ed_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       write(LOGfile,"(A)")"Restarting from a state_list file:"
       list_len=file_length("state_list"//reg(ed_file_suffix)//".restart")
       allocate(list_sector(list_len))
       !
       open(free_unit(unit),file="state_list"//reg(ed_file_suffix)//".restart",status="old")
       status=0
       do while(status>=0)
          read(unit,*,iostat=status) istate,adouble,adouble,sz,isector,anint
          list_sector(istate)=isector
          if(sz/=getsz(isector))stop "setup_pointers_superc error: sz!=getsz(isector)."
       enddo
       close(unit)
       !
       lanc_nstates_total = list_len
       do isector=1,Nsectors
          neigen_sector(isector) = max(1,count(list_sector==isector))
       enddo
    else
       do isector=1,Nsectors
          neigen_sector(isector) = min(getDim(isector),lanc_nstates_sector)
       enddo
    endif
    !
    twin_mask=.true.
    if(ed_twin)then
       write(LOGfile,*)"USE WITH CAUTION: TWIN STATES IN SC CHANNEL!!"
       do isector=1,Nsectors
          sz=getsz(isector)
          if(sz>0)twin_mask(isector)=.false.
       enddo
       write(LOGfile,"(A,I4,A,I4)")"Looking into ",count(twin_mask)," sectors out of ",Nsectors
    endif
    !
    stride=Nlat*Norb
    do ibath=1,Nbath
       do ilat=1,Nlat
          do iorb=1,Norb
             getBathStride(ilat,iorb,ibath) = stride + imp_state_index(ilat,iorb) + (ibath-1)*Nlat*Norb
          enddo
       enddo
    enddo
    !    
    getCsector=0
    !c_up
    do isector=1,Nsectors
       isz=getsz(isector);if(isz==-Ns)cycle
       jsz=isz-1
       jsector=getsector(jsz,1)
       getCsector(1,1,isector)=jsector
    enddo
    !c_dw
    do isector=1,Nsectors
       isz=getsz(isector);if(isz==Ns)cycle
       jsz=isz+1
       jsector=getsector(jsz,1)
       getCsector(1,2,isector)=jsector
    enddo
    !
    getCDGsector=0
    !cdg_up
    do isector=1,Nsectors
       isz=getsz(isector);if(isz==Ns)cycle
       jsz=isz+1
       jsector=getsector(jsz,1)
       getCDGsector(1,1,isector)=jsector
    enddo
    !cdg_dw
    do isector=1,Nsectors
       isz=getsz(isector);if(isz==-Ns)cycle
       jsz=isz-1
       jsector=getsector(jsz,1)
       getCDGsector(1,2,isector)=jsector
    enddo
    return
  end subroutine setup_global_superc









  !##################################################################
  !##################################################################
  !AUXILIARY PROCEDURES - Sectors,Nup,Ndw,DimUp,DimDw,...
  !##################################################################
  !##################################################################
  function get_normal_sector_dimension(n,m) result(dim)
    !
    !Returns the dimension of the symmetry sector per orbital and spin with quantum numbers :math:`\vec{Q}=[\vec{N}_\uparrow,\vec{N}_\downarrow]`. 
    !
    !:f:var:`dim` = :math:`\binom{n}{m}`
    !
    integer,intent(in) :: n,m
    integer            :: dim
    dim = binomial(n,m)
  end function get_normal_sector_dimension

  function get_superc_sector_dimension(mz) result(dim)
    !
    !Returns the dimension of the symmetry sector with quantum numbers :math:`\vec{Q}=S_z=N_\uparrow-N_\downarrow`
    !
    !:f:var:`dim` = :math:`\sum_i 2^{N-mz-2i}\binom{N}{N-mz-2i}\binom{mz+2i}{i}`
    !
    integer :: mz
    integer :: i,dim,Nb
    dim=0
    Nb=Ns-mz
    do i=0,Nb/2 
       dim=dim + 2**(Nb-2*i)*binomial(ns,Nb-2*i)*binomial(ns-Nb+2*i,i)
    enddo
  end function get_superc_sector_dimension





  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the binomial factor n1 over n2
  !+------------------------------------------------------------------+
  elemental function binomial(n1,n2) result(nchoos)
    integer,intent(in) :: n1,n2
    real(8)            :: xh
    integer            :: i
    integer nchoos
    xh = 1.d0
    if(n2<0) then
       nchoos = 0
       return
    endif
    if(n2==0) then
       nchoos = 1
       return
    endif
    do i = 1,n2
       xh = xh*dble(n1+1-i)/dble(i)
    enddo
    nchoos = int(xh + 0.5d0)
  end function binomial



end MODULE ED_SETUP













!   ! SPIN-RESOLVED ALGORITHM: two map objects Hsigma, with a *sparse* 
!   !                          structure capable to store separately
!   !                          the impurity and the bath spin-states
!   subroutine build_sector_spin(isector,Hup,Hdw)
!     integer                            :: isector
!     type(sector_map)                   :: HUP !Map for the UPs
!     type(sector_map)                   :: HDW !Map for the Dws
!     integer,dimension(Ns_Ud)           :: Nups,Ndws
!     integer,dimension(Ns_Ud)           :: DimUps,DimDws
!     integer                            :: DimUp,DimDw,impDIM
!     integer                            :: iup,idw
!     integer                            :: nup_,ndw_
!     integer                            :: imap,iud
!     integer                            :: iIMP,iBATH
!     !
!     impDIM = 2**(Nimp/Ns_ud) !Number of states for the impurity
!     !
!     !Init UP sub-sector:
!     ! > allocate UP-map to binomial(Ns_Orb Nup)
!     call get_Nup(isector,Nups)
!     call get_DimUp(isector,DimUps); DimUp = product(DimUps)
!     call map_allocate(HUP,DimUp,impDIM)
!     !
!     !Init DW sub-sector:
!     ! > allocate DW-map to binomial(Ns_Orb Ndw)
!     call get_Ndw(isector,Ndws)
!     call get_DimDw(isector,DimDws); DimDw = product(DimDws)
!     call map_allocate(HDW,DimDw,impDIM)
!     !
!     !Formally dealing with ed_total_ud == .false.
!     !->iud=1:Ns_Ud (formally because Ns_Ud = 1 in CDMFT code...)
!     ![STILL PROBLEMS WITH Ns_Ud, don't know how to deal with it in density_matrix_impurity()] 
!     do iud=1,Ns_Ud 
! #ifdef _DEBUG
!        if(ed_verbose>3)write(Logfile,"(A)")&
!             "  DEBUG build_sector_spin(): working UP sub-sector"
! #endif
!        imap=0
!        do iup=0,2**Ns_Orb-1
!           nup_ = popcnt(iup) !equivalent to sum(binary_decomposition(iup))
!           if(nup_ /= Nups(iud))cycle !the state does not have the required number of UPs
!           imap = imap+1
!           !HUP(iud)%map(imap) = iup
!           HUP%map(imap) = iup
!           !
!           iIMP  = ibits(iup,0,Nimp)
!           iBATH = ibits(iup,Nimp,Nimp*Nbath) !check: Nimp+Nimp*Nbath=Nimp(1*Nbath)=Ns
!           !call sp_insert_state(HUP(iud)%sp,iIMP,iBATH,imap)
! #ifdef _DEBUG
!           if(ed_verbose>4)then 
!              write(Logfile,"(A)")&
!                   "    DEBUG build_sector_spin(): inserting UP state"
!              write(Logfile,"(A)")&
!                   "      Iup: "//str(iup)//" | IimpUp: "//str(iIMP)//" | IbathUp: "//str(iBATH)
!           endif
! #endif
!           call sp_insert_state(HUP%sp,iIMP,iBATH,imap) 
!           !
!        enddo
!        !
! #ifdef _DEBUG
!        if(ed_verbose>3)write(Logfile,"(A)")&
!             "  DEBUG build_sector_spin(): working DW sub-sector"
! #endif
!        imap=0
!        do idw=0,2**Ns_Orb-1
!           ndw_=popcnt(idw)   !equivalent to sum(binary_decomposition(idw))
!           if(ndw_ /= Ndws(iud))cycle !the state does not have the required number of DWs
!           imap = imap+1
!           !HDW(iud)%map(imap) = idw
!           HDW%map(imap) = idw
!           !
!           iIMP  = ibits(idw,0,Nimp)
!           iBATH = ibits(idw,Nimp,Nimp*Nbath) !check: Nimp+Nimp*Nbath=Nimp(1*Nbath)=Ns
!           !call sp_insert_state(HDW(iud)%sp,iIMP,iBATH,imap)
! #ifdef _DEBUG
!           if(ed_verbose>4)then 
!              write(LOGfile,"(A)")&
!                   "    DEBUG build_sector_spin(): inserting DW state"
!              write(LOGfile,"(A)")&
!                   "      Idw: "//str(idw)//" | IimpDw: "//str(iIMP)//" | IbathDw: "//str(iBATH)
!           endif
! #endif
!           call sp_insert_state(HDW%sp,iIMP,iBATH,imap)
!           !
!        enddo
!     enddo
!     !
!   end subroutine build_sector_spin




! getCsector=0
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup-1; jdw=ndw; if(jup < 0)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCsector(1,1,isector)=jsector
! enddo
! !
! !
! !
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup;jdw=ndw-1;if(jdw < 0)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCsector(1,2,isector)=jsector
! enddo
! !
! !
! !
! getCDGsector=0
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup+1;jdw=ndw;if(jup > Ns)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCDGsector(1,1,isector)=jsector
! enddo
! !
! !
! !
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup;jdw=ndw+1;if(jdw > Ns)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCDGsector(1,2,isector)=jsector
! enddo

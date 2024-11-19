MODULE ED_VARS_GLOBAL
  !
  !Contains all variables, arrays and derived types instances shared throughout the code.
  !Specifically, it contains definitions of the :f:var:`effective_bath`, the :f:var:`gfmatrix` and the :f:var:`sector` data structures. 
  !
  USE SF_CONSTANTS
  USE ED_SPARSE_MATRIX
  USE ED_SPARSE_MAP
  USE ED_INPUT_VARS
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none

  !-------------------- H EXPANSION STRUCTURE ----------------------!
  type H_operator
     ! The matrix storing in the basis [:f:var:`nambu` , :f:var:`nambu` , :f:var:`nspin` , :f:var:`nspin` , :f:var:`nimp` , :f:var:`nimp` ] each element of the Matrix basis decomposing the replica/general bath Hamiltonian :math:`H_p=\sum_{i=1}^{N_{basis}} \lambda_i(p) O_i`, where :math:`N_{basis}` is the dimension of the user defined basis.  
     complex(8),dimension(:,:,:,:,:,:),allocatable           :: O  !Matrix Basis for Bath Hamiltonian decomposition
  end type H_operator


  !-------------------- EFFECTIVE BATH STRUCTURE ----------------------!
  type effective_bath_component
     ! Effective bath component for the replica/general bath. Each istance of this type defines the parameters :math:`\vec{\lambda}` and the amplitudes :math:`\vec{V}`. The first is used to decompose the Hamiltonian of each element of the bath :math:`H_p=\sum_{i=1}^{N_{basis}} \lambda_i(p) O_i`, the latter describes the hopping from/to the impurity.
     real(8),dimension(:),allocatable                        :: v      ![1] for "replica" bath; [Nbath] for "general" bath
     real(8),dimension(:),allocatable                        :: lambda ![:f:var:`nsym`]
  end type effective_bath_component

  type effective_bath
     ! This structure describes the (effective) discretized bath used in the contruct the Hamiltonian of the quantum impurity system. Each element of this structure is allocated and used according the value of :f:var:`ed_mode` = :code:`normal,superc,nonsu2` and :f:var:`bath_type` = :code:`normal,hybrid,replica,general`.  
     integer                                                 :: Nbasis  !The replica/general Matrix basis dimension     
     type(effective_bath_component),dimension(:),allocatable :: item    ![ :f:var:`nbath` ] Replica/General bath components, V included
     logical                                                 :: status=.false.
  end type effective_bath


  !-------------------- CUSTOM OBSERVABLE STRUCTURE ----------------------!
  type observable
     complex(8),dimension(:,:,:),allocatable                 :: sij
     character(len=32)                                       :: o_name
     real(8)                                                 :: o_value
  end type observable

  type custom_observables
     type(observable),dimension(:),allocatable               :: item
     complex(8),dimension(:,:,:),allocatable                 :: Hk
     integer                                                 :: N_asked
     integer                                                 :: N_filled
     logical                                                 :: init=.false.
  end type custom_observables


  !---------------- SECTOR-TO-FOCK SPACE STRUCTURE -------------------!
  type sector_map
     integer,dimension(:),allocatable                        :: map
     type(sparse_map)                                        :: sp
     logical                                                 :: status=.false.
  end type sector_map


  type sector
     integer                                                 :: index       !
     type(sector_map),dimension(:),allocatable               :: H
     integer,dimension(:),allocatable                        :: DimUps
     integer,dimension(:),allocatable                        :: DimDws
     integer                                                 :: DimUp
     integer                                                 :: DimDw
     integer                                                 :: Dim
     integer,dimension(:),allocatable                        :: Nups
     integer,dimension(:),allocatable                        :: Ndws
     integer                                                 :: Nup
     integer                                                 :: Ndw
     integer                                                 :: Sz
     integer                                                 :: Nlanc
     logical                                                 :: status=.false.
  end type sector


  !-------------- GMATRIX FOR FAST EVALUATION OF GF ------------------!

  type GFspectrum
     !The contributions to the GF Kallen-Lehmann sum are stored as :math:`G_{ab,sr}.state.channel.[w,e]`. A couple of weight,poles :math:`[w,e]` is stored for each *channel*, corresponding to c,cdg or any their combination thereof as well as for any state :math:`|n\rangle` of the spectrum such that :math:`G(z) = \sum \frac{w}{z-e}`
     complex(8),dimension(:),allocatable                     :: weight
     complex(8),dimension(:),allocatable                     :: poles
  end type GFspectrum



  type GFchannel
     type(GFspectrum),dimension(:),allocatable               :: channel !N_channel = 2 (c,cdag), 4 (c,cdag,c pm cdag)
  end type GFchannel


  type GFmatrix
     type(GFchannel),dimension(:),allocatable                :: state !state_list%size = # of state in the spectrum 
  end type GFmatrix



  interface allocate_GFmatrix
     module procedure :: allocate_GFmatrix_Nstate
     module procedure :: allocate_GFmatrix_Nchan
     module procedure :: allocate_GFmatrix_Nexc
  end interface allocate_GFmatrix


  interface deallocate_GFmatrix
     module procedure :: deallocate_GFmatrix_single
     module procedure :: deallocate_GFmatrix_all1
     module procedure :: deallocate_GFmatrix_all2
     module procedure :: deallocate_GFmatrix_all3
     module procedure :: deallocate_GFmatrix_all4
     module procedure :: deallocate_GFmatrix_all5
     module procedure :: deallocate_GFmatrix_all6
  end interface deallocate_GFmatrix

  interface write_GFmatrix
     module procedure :: write_GFmatrix_single
     module procedure :: write_GFmatrix_all1
     module procedure :: write_GFmatrix_all2
     module procedure :: write_GFmatrix_all3
     module procedure :: write_GFmatrix_all4
     module procedure :: write_GFmatrix_all5
     module procedure :: write_GFmatrix_all6
  end interface write_GFmatrix

  interface read_GFmatrix
     module procedure :: read_GFmatrix_single
     module procedure :: read_GFmatrix_all1
     module procedure :: read_GFmatrix_all2
     module procedure :: read_GFmatrix_all3
     module procedure :: read_GFmatrix_all4
     module procedure :: read_GFmatrix_all5
     module procedure :: read_GFmatrix_all6
  end interface read_GFmatrix






  !------------------ ABTRACT INTERFACES PROCEDURES ------------------!
  !SPARSE MATRIX-VECTOR PRODUCTS USED IN ED_MATVEC
  !dbleMat*dbleVec
  abstract interface
     subroutine cc_sparse_HxV(Nloc,v,Hv)
       integer                    :: Nloc
       complex(8),dimension(Nloc) :: v
       complex(8),dimension(Nloc) :: Hv
     end subroutine cc_sparse_HxV
  end interface





  !-------------------------- ED  VARIABLES --------------------------!

  !SIZE OF THE PROBLEM
  !=========================================================
  integer,save                                           :: Ns       !Number of levels per spin
  integer,save                                           :: Nsectors !Number of sectors
  integer,save                                           :: Ns_orb
  integer,save                                           :: Ns_ud
  !
  integer                                                :: Nimp     !Total number of levels in the impurity cluster: Nlat*Norb
  integer                                                :: Nambu   !Global Nambu factor for SC calculations
  

  !Some maps between sectors and full Hilbert space (pointers)
  !PRIVATE:
  !=========================================================
  integer,allocatable,dimension(:)                       :: getDim             ! [Nsectors]
  integer,allocatable,dimension(:,:,:)                   :: getCsector         ! [1/Norb,2,NSectors]
  integer,allocatable,dimension(:,:,:)                   :: getCDGsector       ! [1/Norb,2,NSectors]
  integer,allocatable,dimension(:,:)                     :: getBathStride
  logical,allocatable,dimension(:)                       :: twin_mask
  logical,allocatable,dimension(:)                       :: sectors_mask
  integer,allocatable,dimension(:,:)                     :: getSector
  integer,allocatable,dimension(:)                       :: getSz


  !Effective Bath used in the ED code (this is opaque to user)
  !PRIVATE
  !=========================================================
  type(effective_bath)                                   :: dmft_bath !instance of :f:var:`effective_bath` used to store the quantum impurity effective bath in the rest of the code 





  type(H_operator),dimension(:),allocatable              :: Hbath_basis  ![Nsym]
  real(8),dimension(:,:),allocatable                     :: Hbath_lambda ![Nbath,Nsym]
  logical                                                :: Hbath_status=.false.


  !local part of the Hamiltonian
  !INTERNAL USE (accessed thru functions)
  !=========================================================
  complex(8),dimension(:,:,:,:),allocatable              :: impHloc           !local hamiltonian [Nspin][Nspin][Nimp][Nimp]



  !Variables for DIAGONALIZATION
  !PRIVATE
  !=========================================================  
  type(sparse_matrix_csr)                                :: spH0d !diagonal part
  type(sparse_matrix_csr)                                :: spH0nd !non-diagonal part
  type(sparse_matrix_csr),dimension(:),allocatable       :: spH0ups,spH0dws !reduced UP and DW parts
  !
  procedure(cc_sparse_HxV),pointer                       :: spHtimesV_p=>null()


  !Variables for DIAGONALIZATION
  !PRIVATE
  !=========================================================  
  integer,allocatable,dimension(:)                       :: neigen_sector
  logical                                                :: trim_state_list=.false.


  !Partition function
  !PRIVATE
  !=========================================================
  real(8)                                                :: zeta_function
  real(8)                                                :: gs_energy
  real(8)                                                :: max_exc





  !Impurity Green's function and Self-Energies: (Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,:)
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  ! complex(8),allocatable,dimension(:,:,:,:,:,:,:)        :: impGmats
  ! complex(8),allocatable,dimension(:,:,:,:,:,:,:)        :: impGreal
  ! complex(8),allocatable,dimension(:,:,:,:,:,:,:)        :: impG0mats
  ! complex(8),allocatable,dimension(:,:,:,:,:,:,:)        :: impG0real
  ! complex(8),allocatable,dimension(:,:,:,:,:,:,:)        :: impSmats 
  ! complex(8),allocatable,dimension(:,:,:,:,:,:,:)        :: impSreal
  !
  type(GFmatrix),allocatable,dimension(:,:,:,:,:,:)      :: impGmatrix





  !Density and double occupancy
  !Local energies and generalized double occupancies
  !PRIVATE accessible thru routines)
  !=========================================================
  real(8),dimension(:),allocatable                       :: ed_dens
  real(8),dimension(:),allocatable                       :: ed_dens_up,ed_dens_dw
  real(8),dimension(:),allocatable                       :: ed_docc
  real(8),dimension(:),allocatable                       :: ed_mag
  real(8)                                                :: ed_Epot
  real(8)                                                :: ed_Eint
  real(8)                                                :: ed_Ehartree
  real(8)                                                :: ed_Eknot
  real(8)                                                :: ed_Dust
  real(8)                                                :: ed_Dund
  real(8)                                                :: ed_Dse
  real(8)                                                :: ed_Dph
  type(custom_observables)                               :: custom_o


  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable                       :: wm,tau,wr,vm,vr




  !Impurity operators
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:,:)          :: single_particle_density_matrix ![Nambu,Nambu,Nspin,Nspin,Nimp,Nimp]
  complex(8),allocatable,dimension(:,:)                  :: cluster_density_matrix         ![4**(Nimp),4**(Nimp)]






  !File suffixes for printing fine tuning.
  !=========================================================
  character(len=32)                                      :: ed_file_suffix=""
  character(len=10)                                      :: indx_suffix="_indx"
  integer                                                :: indx_padding=4
  logical                                                :: Jhflag  




  !This is the internal Mpi Communicator and variables.
  !=========================================================
#ifdef _MPI
  integer                                                :: MpiComm_Global=MPI_COMM_NULL
  integer                                                :: MpiComm=MPI_COMM_NULL
#endif
  integer                                                :: MpiGroup_Global=MPI_GROUP_NULL
  integer                                                :: MpiGroup=MPI_GROUP_NULL
  logical                                                :: MpiStatus=.false.
  logical                                                :: MpiMaster=.true.
  integer                                                :: MpiRank=0
  integer                                                :: MpiSize=1
  integer,allocatable,dimension(:)                       :: MpiMembers
  integer                                                :: mpiQup=0
  integer                                                :: mpiRup=0
  integer                                                :: mpiQdw=0
  integer                                                :: mpiRdw=0
  integer                                                :: mpiQ=0
  integer                                                :: mpiR=0
  integer                                                :: mpiIstart
  integer                                                :: mpiIend
  integer                                                :: mpiIshift
  logical                                                :: mpiAllThreads=.true.



contains





  !=========================================================
  subroutine ed_set_MpiComm()
#ifdef _MPI
    integer :: comm,ierr
    ! call MPI_Comm_dup(Comm,MpiComm_Global,ierr)
    ! call MPI_Comm_dup(Comm,MpiComm,ierr)
    MpiComm_Global = MPI_COMM_WORLD
    MpiComm        = MPI_COMM_WORLD
    call Mpi_Comm_group(MpiComm_Global,MpiGroup_Global,ierr)
    MpiStatus      = .true.
    MpiSize        = get_Size_MPI(MpiComm_Global)
    MpiRank        = get_Rank_MPI(MpiComm_Global)
    MpiMaster      = get_Master_MPI(MpiComm_Global)
#endif
  end subroutine ed_set_MpiComm

  subroutine ed_del_MpiComm()
#ifdef _MPI    
    MpiComm_Global = MPI_UNDEFINED
    MpiComm        = MPI_UNDEFINED
    MpiGroup_Global= MPI_GROUP_NULL
    MpiStatus      = .false.
    MpiSize        = 1
    MpiRank        = 0
    MpiMaster      = .true.
#endif
  end subroutine ed_del_MpiComm
  !=========================================================







  !=========================================================
  !Allocate the channels in GFmatrix structure
  subroutine allocate_gfmatrix_Nstate(self,Nstate)
    type(GFmatrix) :: self
    integer        :: Nstate
    if(allocated(self%state))deallocate(self%state)
    allocate(self%state(Nstate))
  end subroutine allocate_gfmatrix_Nstate

  subroutine allocate_gfmatrix_Nchan(self,istate,Nchan)
    type(GFmatrix) :: self
    integer        :: istate,Nchan
    if(allocated(self%state(istate)%channel))deallocate(self%state(istate)%channel)
    allocate(self%state(istate)%channel(Nchan))
  end subroutine allocate_gfmatrix_Nchan

  !Allocate the Excitations spectrum at a given channel
  subroutine allocate_gfmatrix_Nexc(self,istate,ichan,Nexc)
    type(GFmatrix) :: self
    integer        :: istate,ichan
    integer        :: Nexc
    if(allocated(self%state(istate)%channel(ichan)%weight))deallocate(self%state(istate)%channel(ichan)%weight)
    if(allocated(self%state(istate)%channel(ichan)%poles))deallocate(self%state(istate)%channel(ichan)%poles)
    allocate(self%state(istate)%channel(ichan)%weight(Nexc))
    allocate(self%state(istate)%channel(ichan)%poles(Nexc))
  end subroutine allocate_gfmatrix_Nexc





  subroutine deallocate_gfmatrix_single(self)
    type(GFmatrix) :: self
    integer        :: istate,ichan
    if(self%status)then
       do istate=1,size(self%state)
          if(allocated(self%state(istate)%channel))then
             do ichan=1,size(self%state(istate)%channel)
                if(allocated(self%state(istate)%channel(ichan)%weight))&
                     deallocate(self%state(istate)%channel(ichan)%weight)
                !
                if(allocated(self%state(istate)%channel(ichan)%poles))&
                     deallocate(self%state(istate)%channel(ichan)%poles)
             enddo
             deallocate(self%state(istate)%channel)
          endif
       enddo
       deallocate(self%state)       
    endif
    self%status=.false.
  end subroutine deallocate_gfmatrix_single

  subroutine deallocate_gfmatrix_all1(self)
    type(GFmatrix),dimension(:) :: self
    integer                     :: i1
    do i1=1,size(self)
       call deallocate_gfmatrix_single(self(i1))
    enddo
  end subroutine deallocate_gfmatrix_all1

  subroutine deallocate_gfmatrix_all2(self)
    type(GFmatrix),dimension(:,:) :: self
    integer                       :: i1,i2
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          call deallocate_gfmatrix_single(self(i1,i2))
       enddo
    enddo
  end subroutine deallocate_gfmatrix_all2

  subroutine deallocate_gfmatrix_all3(self)
    type(GFmatrix),dimension(:,:,:) :: self
    integer                         :: i1,i2,i3
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             call deallocate_gfmatrix_single(self(i1,i2,i3))
          enddo
       enddo
    enddo
  end subroutine deallocate_gfmatrix_all3

  subroutine deallocate_gfmatrix_all4(self)
    type(GFmatrix),dimension(:,:,:,:) :: self
    integer                           :: i1,i2,i3,i4
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             do i4=1,size(self,4)
                call deallocate_gfmatrix_single(self(i1,i2,i3,i4))
             enddo
          enddo
       enddo
    enddo
  end subroutine deallocate_gfmatrix_all4

  subroutine deallocate_gfmatrix_all5(self)
    type(GFmatrix),dimension(:,:,:,:,:) :: self
    integer                           :: i1,i2,i3,i4,i5
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             do i4=1,size(self,4)
                do i5=1,size(self,5)
                   call deallocate_gfmatrix_single(self(i1,i2,i3,i4i5))
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine deallocate_gfmatrix_all5

  subroutine deallocate_gfmatrix_all6(self)
    type(GFmatrix),dimension(:,:,:,:,:,:) :: self
    integer                           :: i1,i2,i3,i4,i5,i6
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             do i4=1,size(self,4)
                do i5=1,size(self,5)
                   do i6=1,size(self,6)
                      call deallocate_gfmatrix_single(self(i1,i2,i3,i4,i5,i6))
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine deallocate_gfmatrix_all6





  !+-------------------------------------------------------------------+
  !PURPOSE  : WRITE GFmatrix to file
  !+-------------------------------------------------------------------+
  subroutine write_gfmatrix_single(self,file)
    class(GFmatrix)    :: self
    character(len=*)   :: file
    integer            :: unit_
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG write_gfmatrix_single: write self"
#endif
    unit_=free_unit()
    open(unit_,file=str(file))
    call write_formatted_gfmatrix(self,unit_)
    close(unit_)
  end subroutine write_gfmatrix_single

  subroutine write_gfmatrix_all1(self,file)
    class(GFmatrix)  :: self(:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1
    unit_=free_unit()
    open(unit_,file=str(file))
    do i1=1,size(self)
       call write_formatted_gfmatrix(self(i1),unit_)
    enddo
    close(unit_)
  end subroutine write_gfmatrix_all1

  subroutine write_gfmatrix_all2(self,file)
    class(GFmatrix)  :: self(:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2
    unit_=free_unit()
    open(unit_,file=str(file))
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          call write_formatted_gfmatrix(self(i1,i2),unit_)
       enddo
    enddo
    close(unit_)
  end subroutine write_gfmatrix_all2

  subroutine write_gfmatrix_all3(self,file)
    class(GFmatrix)  :: self(:,:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2,i3
    unit_=free_unit()
    open(unit_,file=str(file))
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             call write_formatted_gfmatrix(self(i1,i2,i3),unit_)
          enddo
       enddo
    enddo
    close(unit_)
  end subroutine write_gfmatrix_all3

  subroutine write_gfmatrix_all4(self,file)
    class(GFmatrix)  :: self(:,:,:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2,i3,i4
    unit_=free_unit()
    open(unit_,file=str(file))
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             do i4=1,size(self,4)
                call write_formatted_gfmatrix(self(i1,i2,i3,i4),unit_)
             enddo
          enddo
       enddo
    enddo
    close(unit_)
  end subroutine write_gfmatrix_all4

  subroutine write_gfmatrix_all5(self,file)
    class(GFmatrix)  :: self(:,:,:,:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2,i3,i4,i5
    unit_=free_unit()
    open(unit_,file=str(file))
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             do i4=1,size(self,4)
                do i5=1,size(self,5)
                   call write_formatted_gfmatrix(self(i1,i2,i3,i4,i5),unit_)
                enddo
             enddo
          enddo
       enddo
    enddo
    close(unit_)
  end subroutine write_gfmatrix_all5

  subroutine write_gfmatrix_all6(self,file)
    class(GFmatrix)  :: self(:,:,:,:,:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2,i3,i4,i5,i6
    unit_=free_unit()
    open(unit_,file=str(file))
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             do i4=1,size(self,4)
                do i5=1,size(self,5)
                   do i6=1,size(self,6)
                      call write_formatted_gfmatrix(self(i1,i2,i3,i4,i5,i6),unit_)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    close(unit_)
  end subroutine write_gfmatrix_all6




  !+-------------------------------------------------------------------+
  !PURPOSE  : Read cluster GF from file
  !+-------------------------------------------------------------------+
  subroutine read_gfmatrix_single(self,file)
    class(GFmatrix)  :: self
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2,i3,i4
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG reading_gfmatrix_single: reading self"
#endif
    call deallocate_GFmatrix(self)
    unit_=free_unit()
    open(unit_,file=str(file))
    rewind(unit_)
    call read_formatted_gfmatrix(self,unit_)
    close(unit_)
  end subroutine read_gfmatrix_single

  subroutine read_gfmatrix_all1(self,file)
    class(GFmatrix)  :: self(:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1
    call deallocate_GFmatrix(self)
    unit_=free_unit()
    open(unit_,file=str(file))
    rewind(unit_)
    do i1=1,size(self)
       call read_formatted_gfmatrix(self(i1),unit_)
    enddo
    close(unit_)
  end subroutine read_gfmatrix_all1

  subroutine read_gfmatrix_all2(self,file)
    class(GFmatrix)  :: self(:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2
    call deallocate_GFmatrix(self)
    unit_=free_unit()
    open(unit_,file=str(file))
    rewind(unit_)
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          call read_formatted_gfmatrix(self(i1,i2),unit_)
       enddo
    enddo
    close(unit_)
  end subroutine read_gfmatrix_all2

  subroutine read_gfmatrix_all3(self,file)
    class(GFmatrix)  :: self(:,:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2,i3
    call deallocate_GFmatrix(self)
    unit_=free_unit()
    open(unit_,file=str(file))
    rewind(unit_)
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             call read_formatted_gfmatrix(self(i1,i2,i3),unit_)
          enddo
       enddo
    enddo
    close(unit_)
  end subroutine read_gfmatrix_all3

  subroutine read_gfmatrix_all4(self,file)
    class(GFmatrix)  :: self(:,:,:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2,i3,i4
    call deallocate_GFmatrix(self)
    unit_=free_unit()
    open(unit_,file=str(file))
    rewind(unit_)
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             do i4=1,size(self,4)
                call read_formatted_gfmatrix(self(i1,i2,i3,i4),unit_)
             enddo
          enddo
       enddo
    enddo
    close(unit_)
  end subroutine read_gfmatrix_all4

  subroutine read_gfmatrix_all5(self,file)
    class(GFmatrix)  :: self(:,:,:,:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2,i3,i4,i5
    call deallocate_GFmatrix(self)
    unit_=free_unit()
    open(unit_,file=str(file))
    rewind(unit_)
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             do i4=1,size(self,4)
                do i5=1,size(self,5)
                   call read_formatted_gfmatrix(self(i1,i2,i3,i4,i5),unit_)
                enddo
             enddo
          enddo
       enddo
    enddo
    close(unit_)
  end subroutine read_gfmatrix_all5

  subroutine read_gfmatrix_all6(self,file)
    class(GFmatrix)  :: self(:,:,:,:,:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2,i3,i4,i5,i6
    call deallocate_GFmatrix(self)
    unit_=free_unit()
    open(unit_,file=str(file))
    rewind(unit_)
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             do i4=1,size(self,4)
                do i5=1,size(self,5)
                   do i6=1,size(self,6)
                      call read_formatted_gfmatrix(self(i1,i2,i3,i4,i5,i6),unit_)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    close(unit_)
  end subroutine read_gfmatrix_all6


  !+-------------------------------------------------------------------+
  !PURPOSE  : write overload for GFmatrix type (formatted)
  !+-------------------------------------------------------------------+
  subroutine write_formatted_gfmatrix(dtv, unit)
    class(GFmatrix), intent(in)         :: dtv
    integer, intent(in)                 :: unit
    integer                             :: iexc,Ichan,istate
    integer                             :: Nexc,Nchan,Nstates
    write(unit,*) dtv%status
    if(.not.dtv%status)return
    Nstates = size(dtv%state)
    write(unit,*) Nstates
    do istate=1,Nstates
       Nchan = size(dtv%state(istate)%channel)
       write(unit,*)Nchan
       do ichan=1,Nchan
          write(unit,*) size(dtv%state(istate)%channel(ichan)%poles)
          write(unit,*) dtv%state(istate)%channel(ichan)%poles
          write(unit,*) dtv%state(istate)%channel(ichan)%weight
       enddo
    enddo
    write(unit,*)""
  end subroutine write_formatted_gfmatrix

  !+-------------------------------------------------------------------+
  !PURPOSE  : read overload for GFmatrix type (formatted)
  !+-------------------------------------------------------------------+
  subroutine read_formatted_gfmatrix(dtv, unit)
    class(GFmatrix), intent(inout)                :: dtv
    integer, intent(in)                           :: unit
    logical                                       :: alloc
    integer                                       :: ichan,Nchan,Nlanc,istate,Nstates
    !
    read(unit,*) alloc
    if(.not.alloc)return
    read(unit,*)Nstates
    call allocate_GFmatrix(dtv,Nstate=Nstates)
    do istate=1,Nstates
       read(unit,*)Nchan
       call allocate_GFmatrix(dtv,istate=istate,Nchan=Nchan)
       do ichan=1,Nchan
          read(unit,*)Nlanc
          call allocate_GFmatrix(dtv,istate=istate,ichan=ichan,Nexc=Nlanc)
          read(unit,*) dtv%state(istate)%channel(ichan)%poles
          read(unit,*) dtv%state(istate)%channel(ichan)%weight
       enddo
    enddo
    !
  end subroutine read_formatted_gfmatrix

END MODULE ED_VARS_GLOBAL

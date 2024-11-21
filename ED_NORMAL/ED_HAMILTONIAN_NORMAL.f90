MODULE ED_HAMILTONIAN_NORMAL
  !Setup and build the sector Hamiltonian, returns the correct dimension of the vectors in the Arpack/Lanczos procedure in each thread and provides an interface to Tri-Diagonalize the Hamiltonian on a Krylov basis given a starting vector.
  !
  USE ED_HAMILTONIAN_NORMAL_COMMON
  USE ED_HAMILTONIAN_NORMAL_STORED_HxV
  USE ED_HAMILTONIAN_NORMAL_DIRECT_HxV
  !
  implicit none
  private


  !>Build sparse hamiltonian of the sector
  public  :: build_Hv_sector_normal
  public  :: delete_Hv_sector_normal
  public  :: vecDim_Hv_sector_normal

  !> Tridiag sparse Hamiltonian of the sector
  public  :: tridiag_Hv_sector_normal



contains




  !####################################################################
  !                 MAIN ROUTINES: BUILD/DELETE SECTOR
  !####################################################################
  subroutine build_Hv_sector_normal(isector,Hmat)
    !
    ! Builds the matrix-vector product :math:`H\times \vec{v}` in the current sector.
    !
    !   #. Building the sector through :f:func:`build_sector` for :f:var:`isector`
    !   #. Retrieve all dimensions of the sectors, setup the MPI split in parallel mode.
    !   #. If total sector dimension is < :f:var:`lanc_dim_threshold` then Hamiltonian is stored into dense matrix for Lapack diagonalization
    !   #. Else we proceeds according to the followins scheme:
    !
    !.. list-table                     :: Build Hamiltonian, :math:`H\times\vec{v}` products.
    !    :widths: auto
    !    :header-rows: 1
    !    :stub-columns: 1
    !
    !    * - 
    !      - :f:var:`ed_sparse_H` = :code:`T`
    !      - :f:var:`ed_sparse_H` = :code:`F`
    !
    !    * - :f:var:`ed_total_ud` = :code:`T`
    !      - | call :f:func:`ed_buildh_normal_main`
    !        | serial: :f:func:`sphtimesv_p` :code:`=>` :f:func:`spmatvec_normal_main` 
    !        | MPI:    :f:func:`sphtimesv_p` :code:`=>` :f:func:`spmatvec_mpi_normal_main`
    !      - | serial: :f:func:`sphtimesv_p` :code:`=>` :f:func:`directmatvec_normal_main` 
    !        | MPI:    :f:func:`sphtimesv_p` :code:`=>` :f:func:`directmatvec_mpi_normal_main`
    !
    !
    integer                            :: isector !Index of the actual sector to be analyzed
    complex(8),dimension(:,:),optional :: Hmat    !Dense matrix to store the sector Hamiltonian is dim < :f:var:`lanc_dim_threshold`
    integer                            :: irank,ierr
    integer                            :: i,iup,idw
    integer                            :: j,jup,jdw
    !
    Hsector=isector
    Hstatus=.true.
    !
    allocate(Hs(2*Ns_Ud))

    call build_sector(isector,Hsector)
    !
    allocate(DimUps(Ns_Ud))
    allocate(DimDws(Ns_Ud))
    Dim    = Hsector%Dim
    DimUp  = Hsector%DimUp
    DimDw  = Hsector%DimDw
    DimUps = Hsector%DimUps
    DimDws = Hsector%DimDws
    !
    !#################################
    !          MPI SETUP
    !#################################
    mpiAllThreads=.true.
    !>PREAMBLE: check that split of the DW is performed with the minimum #cpu: no idle cpus allowed (with zero elements)
#ifdef _MPI
    if(MpiStatus)then
       if(DimDw < MpiSize)then
          if(MpiMaster.AND.ed_verbose>4)write(*,*)"Reducing N_cpu to DimDw:",DimDw,MpiSize-DimDw
          allocate(MpiMembers(0:DimDw-1))
          forall(irank=0:DimDw-1)MpiMembers(irank)=irank       
          call Mpi_Group_Incl(MpiGroup_Global,DimDw,MpiMembers,MpiGroup,ierr)
          call Mpi_Comm_create(MpiComm_Global,MpiGroup,MpiComm,ierr)
          deallocate(MpiMembers)
          mpiAllThreads=.false.
          call Barrier_MPI(MpiComm_Global)
#ifdef _DEBUG
          if(ed_verbose>4)then
             if(MpiMaster)write(LOGfile,*)&
                  "       mpiRank,   MpiComm, Comm_Global, Comm_World, Comm_Null, Undefined"
             do i=0,MpiSize-1
                call Barrier_MPI(MpiComm_Global)
                if(MpiRank==i)write(*,*)i,MpiComm,MpiComm_Global,Mpi_Comm_World,Mpi_comm_null,Mpi_Undefined
             enddo
             call Barrier_MPI(MpiComm_Global)
          endif
#endif
       endif
       if( MpiComm /= MPI_COMM_NULL )then
          MpiRank = Get_Rank_MPI(MpiComm)
          MpiSize = Get_Size_MPI(MpiComm)
       endif
    endif
#endif
    !
    !Dw split:    
    mpiQdw = DimDw/MpiSize
    mpiRdw = mod(DimDw,MpiSize)
    if(MpiRank < mod(DimDw,MpiSize) ) then
       mpiRdw = 0
       MpiQdw = MpiQdw+1
    endif
    !
    !Total split: split DW \times UP
    mpiQ = DimUp*mpiQdw
    mpiR = DimUp*mpiRdw
    mpiIstart = 1 + MpiRank*mpiQ+mpiR
    mpiIend   = (MpiRank+1)*mpiQ+mpiR
    mpiIshift = MpiRank*mpiQ+mpiR
    !
    !
#ifdef _MPI
#ifdef _DEBUG
    if(MpiStatus.AND.ed_verbose>4.AND.(MpiComm/=Mpi_Comm_Null).AND.MpiSize>=1)then
       if(MpiMaster)write(LOGfile,*)&
            "         mpiRank,   mpi_Q,   mpi_R,      mpi_Qdw,      mpiR_dw,  mpi_Istart,  mpi_Iend,  Iend-Istart,  Comm, Comm_Global"
       do irank=0,MpiSize-1
          call Barrier_MPI(MpiComm)
          if(MpiRank==irank)write(*,*)MpiRank,MpiQ,MpiR,mpiQdw,MpiRdw,MpiIstart,MpiIend,MpiIend-MpiIstart+1,MpiComm,MpiComm_Global
       enddo
       call Barrier_MPI(MpiComm)
    endif
#endif
#endif
    !
    !
    !#################################
    !          HxV SETUP
    !#################################
    if(present(Hmat))then
       spHtimesV_p => null()
       call ed_buildh_normal_main(Hmat)          
       return
    endif
    !
    select case (ed_sparse_H)
    case (.true.)
       spHtimesV_p => spMatVec_normal_main
#ifdef _MPI
       if(MpiStatus)spHtimesV_p => spMatVec_MPI_normal_main
#endif
       call ed_buildh_normal_main()
    case (.false.)
       spHtimesV_p => directMatVec_normal_main
#ifdef _MPI
       if(MpiStatus)spHtimesV_p => directMatVec_MPI_normal_main
#endif
    end select
    !
  end subroutine build_Hv_sector_normal





  subroutine delete_Hv_sector()
    !
    ! Delete the all the memory used to construct the sector Hamiltonian and the corresponding matrix vector products.
    ! The sector is deleted, all the dimensions and MPI splitting variables are reset to zero. All the sparse matrices are deallocated having gone out of scope. The abstract interface pointer :f:var:`spHtimesV_p` for the matrix-vector product is nullified. 
    !
    integer :: iud,ierr,i
    call delete_sector(Hsector)
    deallocate(DimUps)
    deallocate(DimDws)
    Dim    = 0
    DimUp  = 0
    DimDw  = 0
    !
    !
    !There is no difference here between Mpi and serial version, as Local part was removed.
#ifdef _MPI
    if(MpiStatus)then
       call sp_delete_matrix(MpiComm,spH0d)
       if(Jhflag)call sp_delete_matrix(MpiComm,spH0nd)
    else
       call sp_delete_matrix(spH0d)
       if(Jhflag)call sp_delete_matrix(spH0nd)
    endif
#else
    call sp_delete_matrix(spH0d)
    if(Jhflag)call sp_delete_matrix(spH0nd)
#endif
    do iud=1,Ns_Ud
       call sp_delete_matrix(spH0ups(iud))
       call sp_delete_matrix(spH0dws(iud))
    enddo
    !
    spHtimesV_p => null()
    !
#ifdef _MPI
    if(MpiStatus)then
       if(MpiGroup/=Mpi_Group_Null)call Mpi_Group_free(MpiGroup,ierr)
       if(MpiComm/=Mpi_Comm_Null.AND.MpiComm/=Mpi_Comm_World)call Mpi_Comm_Free(MpiComm,ierr)
       MpiComm = MpiComm_Global
       MpiSize = get_Size_MPI(MpiComm_Global)
       MpiRank = get_Rank_MPI(MpiComm_Global)
       mpiQup=0
       mpiRup=0
       mpiQdw=0
       mpiRdw=0
       mpiQ=0
       mpiR=0
       mpiIstart=0
       mpiIend=0
       mpiIshift=0
    endif
#endif
    iter=0
    !
  end subroutine delete_Hv_sector






  function vecDim_Hv_sector_normal(isector) result(vecDim)
    !
    ! Returns the dimensions :f:var:`vecdim` of the vectors used in the Arpack/Lanczos produces given the current sector index :f:var:`isector` . If parallel mode is active the returned dimension corresponds to the correct chunk for each thread. 
    !
    integer :: isector          !current sector index
    integer :: vecDim           !vector or vector chunck dimension  
    integer :: mpiQdw
    integer :: DimUps(Ns_Ud),DimUp
    integer :: DimDws(Ns_Ud),DimDw
    !
    call get_DimUp(isector,DimUps) ; DimUp = product(DimUps)
    call get_DimDw(isector,DimDws) ; DimDw = product(DimDws)
    !
#ifdef _MPI
    if(MpiStatus)then
       !Dw split:
       mpiQdw = DimDw/MpiSize
       if(MpiRank < mod(DimDw,MpiSize) ) MpiQdw = MpiQdw+1
    else
       mpiQdw = DimDw
    endif
#else
    mpiQdw = DimDw
#endif
    !
    vecDim=DimUp*mpiQdw
    !
  end function vecDim_Hv_sector_normal



  subroutine tridiag_Hv_sector_normal(isector,vvinit,alanc,blanc,norm2)
    !
    ! Returns the parameters :math:`\vec{\alpha}` and :math:`\vec{\beta}` , respectively :f:var:`alanc` and :f:var:`blanc` , of the partial tridiagonalization of the sector Hamiltonian on a Krylov basis with starting vector :f:var:`vvinit`.
    !
    ! Input:
    !  * :f:var:`isector`
    !  * :f:var:`vvinit`
    !
    ! Output:
    !  * :f:var:`alanc` corresponding to :math:`\vec{\alpha}`
    !  * :f:var:`blanc` corresponding to :math:`\vec{\beta}`
    !  * :f:var:`norm2` the norm of the input vector  :math:`\langle {\rm vvinit}|{\rm vvinit}\rangle` 
    !
    integer                            :: isector !current sector index
    real(8),dimension(:)               :: vvinit  !input vector for the construction of the tridiagonal or Krylov basis
    real(8),dimension(:),allocatable   :: alanc !:math:`\vec{\alpha}` or diagonal parameters of the tridiagonal basis.
    real(8),dimension(:),allocatable   :: blanc !:math:`\vec{\beta}` or sub-/over-diagonal parameters of the tridiagonal basis.
    real(8)                            :: norm2 !norm of the input vector :f:var:`vvinit`
    !
    real(8),dimension(:),allocatable   :: vvloc
    integer                            :: vecDim
    !
    !
    if(MpiMaster)then
       norm2=dot_product(vvinit,vvinit)
       vvinit=vvinit/sqrt(norm2)
    endif
#ifdef _MPI
    if(MpiStatus)call bcast_MPI(MpiComm,norm2)
#endif
    call build_Hv_sector_normal(isector)
    allocate(alanc(Hsector%Nlanc),blanc(Hsector%Nlanc))
    alanc=0d0 ; blanc=0d0
    if(norm2/=0d0)then
#ifdef _MPI
       if(MpiStatus)then
          vecDim = vecDim_Hv_sector_normal(isector)
          allocate(vvloc(vecDim))
          call scatter_vector_MPI(MpiComm,vvinit,vvloc)
          call sp_lanc_tridiag(MpiComm,spHtimesV_p,vvloc,alanc,blanc)
       else
          call sp_lanc_tridiag(spHtimesV_p,vvinit,alanc,blanc)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_p,vvinit,alanc,blanc)
#endif
    endif
    call delete_Hv_sector_normal()
  end subroutine tridiag_Hv_sector_normal


end MODULE ED_HAMILTONIAN_NORMAL

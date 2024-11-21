MODULE ED_HAMILTONIAN_SUPERC
  USE ED_HAMILTONIAN_SUPERC_COMMON
  USE ED_HAMILTONIAN_SUPERC_STORED_HxV
  USE ED_HAMILTONIAN_SUPERC_DIRECT_HxV
  !
  implicit none
  private

  !>Build sparse hamiltonian of the sector
  public  :: build_Hv_sector_superc
  public  :: delete_Hv_sector_superc
  public  :: vecDim_Hv_sector_superc


  !> Tridiag sparse Hamiltonian of the sector
  public  :: tridiag_Hv_sector_superc




contains




  !####################################################################
  !                 MAIN ROUTINES: BUILD/DELETE SECTOR
  !####################################################################
  subroutine build_Hv_sector_superc(isector,Hmat)
    !
    ! Builds the matrix-vector product :math:`H\times \vec{v}` in the current sector.
    !
    !   #. Building the sector through :f:func:`build_sector` for :f:var:`isector`
    !   #. Retrieve all dimensions of the sectors, setup the MPI split in parallel mode.
    !   #. If total sector dimension is < :f:var:`lanc_dim_threshold` then Hamiltonian is stored into dense matrix for Lapack diagonalization
    !   #. Else we proceeds according to the followins scheme:
    !
    !.. list-table:: Build Hamiltonian, :math:`H\times\vec{v}` products.
    !    :widths: auto
    !    :header-rows: 1
    !
    !    * - :f:var:`ed_sparse_h` = :code:`T`
    !      - :f:var:`ed_sparse_h` = :code:`F`
    !
    !    * - | call :f:var:`ed_buildh_superc_main`
    !        | serial: :f:func:`sphtimesv_p` :code:`=>` :f:func:`spmatvec_superc_main` 
    !        | MPI:    :f:func:`sphtimesv_p` :code:`=>` :f:func:`spmatvec_mpi_superc_main`
    !      - | serial: :f:func:`sphtimesv_p` :code:`=>` :f:func:`directmatvec_superc_main` 
    !        | MPI:    :f:func:`sphtimesv_p` :code:`=>` :f:func:`directmatvec_mpi_superc_main`
    !
    !
    integer                            :: isector !Index of the actual sector to be analyzed
    complex(8),dimension(:,:),optional :: Hmat    !Dense matrix to store the sector Hamiltonian is dim < :f:var:`lanc_dim_threshold`
    integer                            :: SectorDim
    integer                            :: irank
    integer                            :: i,j,Dim
    !
    !
    call build_sector(isector,Hsector)
    !
    Dim    = Hsector%Dim
    !
    !#################################
    !          MPI SETUP
    !#################################
    mpiAllThreads=.true.
    MpiR = 0
    if(ed_sparse_h .or. present(Hmat))then
       MpiQ = DimEl/MpiSize
       if(MpiRank==(MpiSize-1))MpiR=mod(DimEl,MpiSize)
    else
       MpiQ = Dim/MpiSize
       if(MpiRank==(MpiSize-1))MpiR=mod(Dim,MpiSize)       
    endif
    !
    MpiIshift = MpiRank*mpiQ
    MpiIstart = MpiRank*mpiQ + 1
    MpiIend   = (MpiRank+1)*mpiQ + mpiR
    !
#ifdef _MPI
#ifdef _DEBUG
    if(MpiStatus.AND.ed_verbose>4)then
       write(LOGfile,*)&
            "         mpiRank,   mpi_Q,   mpi_R,   mpi_Istart,   mpi_Iend,   mpi_Iend-mpi_Istart"
       do irank=0,MpiSize-1
          call Barrier_MPI(MpiComm)
          if(MpiRank==irank)then
             write(LOGfile,*)MpiRank,MpiQ,MpiR,MpiIstart,MpiIend,MpiIend-MpiIstart+1
          endif
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
       call ed_buildh_superc_main(Hmat)
       return
    endif
    !
    select case (ed_sparse_H)
    case (.true.)
       spHtimesV_p => spMatVec_superc_main
#ifdef _MPI
       if(MpiStatus)spHtimesV_p=> spMatVec_MPI_superc_main
#endif
       call ed_buildh_superc_main()
    case (.false.)
       spHtimesV_p => directMatVec_superc_main
#ifdef _MPI
       if(MpiStatus)spHtimesV_p => directMatVec_MPI_superc_main
#endif
    end select
    !
  end subroutine build_Hv_sector_superc





  subroutine delete_Hv_sector_superc()
    !
    ! Delete the all the memory used to construct the sector Hamiltonian and the corresponding matrix vector products.
    ! The sector is deleted, all the dimensions and MPI splitting variables are reset to zero. All the sparse matrices are deallocated having gone out of scope. The abstract interface pointer :f:var:`spHtimesV_p` for the matrix-vector product is nullified. 
    !    
    call delete_sector(Hsector)
    Dim    = 0
    !
    !There is no difference here between Mpi and serial version, as Local part was removed.
#ifdef _MPI
    if(MpiStatus)then
       call sp_delete_matrix(MpiComm,spH0)
    else
       call sp_delete_matrix(spH0)
    endif
#else
    call sp_delete_matrix(spH0)
#endif
    !
    spHtimesV_p => null()
    !
#ifdef _MPI
    if(MpiStatus)then
       MpiComm = MpiComm_Global
       MpiSize = get_Size_MPI(MpiComm_Global)
       MpiRank = get_Rank_MPI(MpiComm_Global)
       mpiQ=0
       mpiR=0
       mpiIstart=0
       mpiIend=0
       mpiIshift=0
    endif
#endif
    iter=0
    !
  end subroutine delete_Hv_sector_superc






  function vecDim_Hv_sector_superc(isector) result(vecDim)
    !
    ! Returns the dimensions :f:var:`vecdim` of the vectors used in the Arpack/Lanczos produces given the current sector index :f:var:`isector` . If parallel mode is active the returned dimension corresponds to the correct chunk for each thread.
    !
    integer :: isector          !current sector index
    integer :: vecDim           !vector or vector chunck dimension  
    integer :: Dim,DimEl
    !
    Dim   = getdim(isector)
    DimEl = Dim/(nph+1)
    !
#ifdef _MPI
    if(MpiStatus)then
       MpiQ = Dim/MpiSize
       MpiR = 0
       if(MpiRank==(MpiSize-1))MpiR=mod(Dim,MpiSize)
    else
       MpiQ = Dim
       MpiR = 0
    endif
#else
    MpiQ = Dim
    MpiR = 0
#endif
    !
    vecDim=MpiQ + MpiR
    !
  end function vecDim_Hv_sector_superc







  subroutine tridiag_Hv_sector_superc(isector,vvinit,alanc,blanc,norm2)
    !
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
    !
    integer                             :: isector
    complex(8),dimension(:)             :: vvinit
    real(8),dimension(:),allocatable    :: alanc,blanc
    real(8)                             :: norm2
    complex(8),dimension(:),allocatable :: vvloc
    integer                             :: vecDim
    !
    !
    if(MpiMaster)then
       norm2=dot_product(vvinit,vvinit)
       vvinit=vvinit/sqrt(norm2)
    endif
#ifdef _MPI
    if(MpiStatus)call bcast_MPI(MpiComm,norm2)
#endif
    call build_Hv_sector_superc(isector)
    allocate(alanc(Hsector%Nlanc),blanc(Hsector%Nlanc))
    alanc=0d0 ; blanc=0d0
    !
    !
    if(norm2/=0d0)then
#ifdef _MPI
       if(MpiStatus)then
          vecDim = vecDim_Hv_sector_superc(isector)
          allocate(vvloc(vecDim))
          call scatter_vector_MPI(MpiComm,vvinit,vvloc)
          call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvloc,alanc,blanc)
       else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alanc,blanc)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_cc,vvinit,alanc,blanc)
#endif
    endif
    call delete_Hv_sector_superc()
  end subroutine tridiag_Hv_sector_superc




end MODULE ED_HAMILTONIAN_SUPERC

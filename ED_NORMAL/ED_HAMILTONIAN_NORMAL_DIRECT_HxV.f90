MODULE ED_HAMILTONIAN_NORMAL_DIRECT_HxV
  !Constructs and applies on-the-fly each term of the sector Hamiltonian to the input
  !vector :math:`\vec{w} = H\times \vec{v}` in a Arpack/Lanczos framework.
  USE ED_HAMILTONIAN_NORMAL_COMMON
  implicit none
  private


  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_normal_main
#ifdef _MPI
  public  :: directMatVec_MPI_normal_main
#endif


  
contains


  subroutine directMatVec_normal_main(Nloc,vin,Hv)
    !
    ! Serial version of the direct, on-the-fly matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in Arpack/Lanczos algorithm for :f:var:`ed_total_ud` = :code:`True` 
    ! This procedures evaluates the non-zero terms of any part of the global Hamiltonian and applies them to the input vector using serial algorithm.  
    !
    integer                                                       :: Nloc !Global dimension of the problem. :code:`size(v)=Nloc=size(Hv)`
    complex(8),dimension(Nloc)                                    :: vin  !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc)                                    :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    complex(8),dimension(:),allocatable                           :: vt,Hvt
    integer,dimension(Ns)                                         :: ibup,ibdw
    integer,dimension(2*Ns_Ud)                                    :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)                               :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)                                         :: Nup,Ndw
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Norb,Norb,Nbath) :: Hbath_tmp
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    if(Nloc/=getdim(isector))stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    !Get diagonal hybridization, bath energy
    if(allocated(diag_hybr))deallocate(diag_hybr)
    if(allocated(bath_diag))deallocate(bath_diag)
    allocate(diag_hybr(Nspin,Nimp,Nbath));diag_hybr=0d0
    allocate(bath_diag(Nspin,Nimp,Nbath));bath_diag=0d0
    do ibath=1,Nbath
       Hbath_tmp(:,:,:,:,:,:,ibath)=Hbath_build(dmft_bath%item(ibath)%lambda)
       do ispin=1,Nspin
          do iorb=1,Nimp
             diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%v(iorb+(ispin-1)*Nimp)
             bath_diag(ispin,iorb,ibath)=dreal(Hbath_tmp(1,1,ispin,ispin,iorb,iorb,ibath))
          enddo
       enddo
    enddo
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "direct/HxV_local.f90"
    !
    !UP HAMILTONIAN TERMS
    include "direct/HxV_up.f90"
    !    
    !DW HAMILTONIAN TERMS
    include "direct/HxV_dw.f90"
    !
    !NON-LOCAL HAMILTONIAN PART: H_non_loc*vin = vout
    if(Jhflag)then
       include "direct/HxV_non_local.f90"
    endif
    !-----------------------------------------------!
    !
    deallocate(diag_hybr,bath_diag)
    return
  end subroutine directMatVec_normal_main




#ifdef _MPI
  subroutine directMatVec_MPI_normal_main(Nloc,vin,Hv)
    !
    ! MPI parallel version of the direct, on-the-fly matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in P-Arpack/P-Lanczos algorithm for :f:var:`ed_total_ud` = :code:`True` 
    ! This procedures evaluates the non-zero terms of any part of the global Hamiltonian and applies them to a part of the vector own by the thread using parallel algorithm.  
    !
    integer                                                       :: Nloc !Local dimension of the vector chunk. :code:`size(v)=Nloc` with :math:`\sum_p` :f:var:`Nloc` = :f:var:`Dim`
    complex(8),dimension(Nloc)                                    :: vin  !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc)                                    :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    complex(8),dimension(:),allocatable                           :: vt,Hvt
    integer,dimension(Ns)                                         :: ibup,ibdw
    integer,dimension(2*Ns_Ud)                                    :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)                               :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)                                         :: Nup,Ndw
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Norb,Norb,Nbath) :: Hbath_tmp
    integer                                                       :: MpiIerr,N
    integer,allocatable,dimension(:)                              :: Counts
    integer,allocatable,dimension(:)                              :: Offset
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    !
    if(MpiComm==MPI_UNDEFINED.OR.MpiComm==Mpi_Comm_Null)&
         stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    ! if(.not.MpiStatus)stop "directMatVec_MPI_cc ERROR: MpiStatus = F"
    !
    !Get diagonal hybridization, bath energy
    if(allocated(diag_hybr))deallocate(diag_hybr)
    if(allocated(bath_diag))deallocate(bath_diag)
    allocate(diag_hybr(Nspin,Nimp,Nbath));diag_hybr=0d0
    allocate(bath_diag(Nspin,Nimp,Nbath));bath_diag=0d0
    do ibath=1,Nbath
       Hbath_tmp(:,:,:,:,:,:,ibath)=Hbath_build(dmft_bath%item(ibath)%lambda)
       do ispin=1,Nspin
          do iorb=1,Nimp
             diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%v(iorb+(ispin-1)*Nimp)
             bath_diag(ispin,iorb,ibath)=dreal(Hbath_tmp(1,1,ispin,ispin,iorb,iorb,ibath))
          enddo
       enddo
    enddo
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "direct_mpi/HxV_local.f90"
    !
    !UP HAMILTONIAN TERMS: MEMORY CONTIGUOUS
    include "direct_mpi/HxV_up.f90"
    !
    !DW HAMILTONIAN TERMS: MEMORY NON-CONTIGUOUS
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    allocate(vt(mpiQup*DimDw)) ;vt=zero
    allocate(Hvt(mpiQup*DimDw));Hvt=zero
    call vector_transpose_MPI(DimUp,MpiQdw,Vin,DimDw,MpiQup,vt) !Vin^T --> Vt
    include "direct_mpi/HxV_dw.f90"
    deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ;vt=zero        !reallocate Vt
    call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt) !Hvt^T --> Vt
    Hv = Hv + Vt
    deallocate(vt)
    !
    !NON-LOCAL HAMILTONIAN PART: H_non_loc*vin = vout
    if(Jhflag)then
       N = 0
       call AllReduce_MPI(MpiComm,Nloc,N)
       !
       allocate(vt(N)) ; vt = zero
       call allgather_vector_MPI(MpiComm,vin,vt)
       !
       include "direct_mpi/HxV_non_local.f90"
       !
       deallocate(Vt)
    endif
    !-----------------------------------------------!
    !
    deallocate(diag_hybr,bath_diag)
    return
  end subroutine directMatVec_MPI_normal_main
  
#endif


END MODULE ED_HAMILTONIAN_NORMAL_DIRECT_HXV

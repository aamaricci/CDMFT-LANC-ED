! > SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT 
MODULE ED_HAMILTONIAN_SUPERC_DIRECT_HxV
  USE ED_HAMILTONIAN_SUPERC_COMMON
  implicit none
  private


  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_superc_main
#ifdef _MPI
  public  :: directMatVec_MPI_superc_main
#endif



contains


  subroutine directMatVec_superc_main(Nloc,vin,Hv)
    !
    ! Serial version of the direct, on-the-fly matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in Arpack/Lanczos algorithm.
    ! This procedures evaluates the non-zero terms of any part of the global Hamiltonian and applies them to the input vector using serial algorithm.  
    !
    integer                                                         :: Nloc  !Global dimension of the problem. :code:`size(v)=Nloc=size(Hv)`
    complex(8),dimension(Nloc)                                      :: vin   !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc)                                      :: Hv    !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    integer,dimension(Nlevels)                                    :: ib
    integer,dimension(Ns)                                         :: ibup,ibdw
    real(8),dimension(Nimp)                                       :: nup,ndw
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,Nbath) :: Hbath_tmp

    !
    !
    if(.not.Hsector%status)stop "ed_buildh_main ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    Dim   = Hsector%Dim
    !
    if(Nloc/=Dim)stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    !Get diagonal hybridization, bath energy
    allocate(diag_hybr(Nspin,Nimp,Nbath));diag_hybr=0d0
    allocate(bath_diag(Nambu,Nspin,Nimp,Nbath));bath_diag=0d0
    do ibath=1,Nbath
       Hbath_tmp(:,:,:,:,:,:,ibath)=Hbath_build(dmft_bath%item(ibath)%lambda)
       do ispin=1,Nspin
          do iorb=1,Nimp
             diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%v(iorb+(ispin-1)*Nimp)
             do in=1,Nambu
                bath_diag(in,ispin,iorb,ibath)=dreal(Hbath_tmp(in,in,ispin,ispin,iorb,iorb,ibath))
             enddo
          enddo
       enddo
    enddo
    !
    Hv=zero
    states: do i=1,Dim
       m    = Hsector%H(1)%map(i)
       ib   = bdecomp(m,2*Ns)
       !
       do iorb=1,Nimp
          nup(iorb)=dble(ib(iorb))
          ndw(iorb)=dble(ib(iorb+Ns))
       enddo
       !
       !
       !-----------------------------------------------!
       !LOCAL HAMILTONIAN TERMS
       include "direct/HxVimp.f90"
       !
       !LOCAL INTERACTION
       include "direct/HxVint.f90"
       !
       !BATH HAMILTONIAN
       include "direct/HxVbath.f90"
       !
       !IMPURITY- BATH HYBRIDIZATION
       include "direct/HxVimp_bath.f90"
       !-----------------------------------------------!
    enddo states
    !
    return
  end subroutine directMatVec_superc_main


#ifdef _MPI
  subroutine directMatVec_MPI_superc_main(Nloc,v,Hv)
    !
    ! MPI parallel version of the direct, on-the-fly matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in P-Arpack/P-Lanczos algorithm.
    ! This procedures evaluates the non-zero terms of any part of the global Hamiltonian and applies them to a part of the vector own by the thread using parallel algorithm.  
    !
    integer                                                       :: Nloc !Local dimension of the vector chunk. :code:`size(v)=Nloc` with :math:`\sum_p` :f:var:`Nloc` = :f:var:`Dim`
    complex(8),dimension(Nloc)                                    :: v    !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc)                                    :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    integer                                                       :: N
    complex(8),dimension(:),allocatable                           :: vin
    integer,allocatable,dimension(:)                              :: Counts,Offset
    integer                                                       :: isector
    integer,dimension(Nlevels)                                    :: ib
    integer,dimension(Ns)                                         :: ibup,ibdw
    real(8),dimension(Nimp)                                       :: nup,ndw
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,Nbath) :: Hbath_tmp
    integer                                                       :: first_state,last_state
    integer                                                       :: first_state_up,last_state_up
    integer                                                       :: first_state_dw,last_state_dw
    integer                                                       :: mpiIerr
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    Dim   = Hsector%Dim
    !
    !
    if(MpiComm==MPI_UNDEFINED.OR.MpiComm==Mpi_Comm_Null)&
         stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    ! if(.not.MpiStatus)stop "directMatVec_MPI_cc ERROR: MpiStatus = F"
    !
    allocate(diag_hybr(Nspin,Nimp,Nbath));diag_hybr=0d0
    allocate(bath_diag(Nambu,Nspin,Nimp,Nbath));bath_diag=0d0
    do ibath=1,Nbath
       Hbath_tmp(:,:,:,:,:,:,ibath)=Hbath_build(dmft_bath%item(ibath)%lambda)
       do ispin=1,Nspin
          do iorb=1,Nimp
             diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%v(iorb+(ispin-1)*Nimp)
             do in=1,Nambu
                bath_diag(in,ispin,iorb,ibath)=dreal(Hbath_tmp(in,in,ispin,ispin,iorb,iorb,ibath))
             enddo
          enddo
       enddo
    enddo
    !

    N=0
    call AllReduce_MPI(MpiComm,Nloc,N)
    !
    !Reconstruct Vin and get the displacements for AllGatherV call
    allocate(Counts(0:MpiSize-1)) ; Counts(0:)=0
    allocate(Offset(0:MpiSize-1)) ; Offset(0:)=0
    !
    Counts(0:)        = N/MpiSize
    Counts(MpiSize-1) = N/MpiSize+mod(N,MpiSize)
    !
    do i=1,MpiSize-1
       Offset(i) = Counts(i-1) + Offset(i-1)
    enddo
    !
    allocate(vin(N)); vin  = zero
    call MPI_Allgatherv(&
         v(1:Nloc),Nloc,MPI_Double_Complex,&
         vin,Counts,Offset,MPI_Double_Complex,&
         MpiComm,MpiIerr)
    !
    Hv=zero
    !
    states: do i=MpiIstart,MpiIend
       m  = Hsector%H(1)%map(i)
       ib = bdecomp(m,2*Ns)
       !
       do iorb=1,Nimp
          nup(iorb)=dble(ib(iorb))
          ndw(iorb)=dble(ib(iorb+Ns))
       enddo
       !

       !IMPURITY  HAMILTONIAN
       include "direct/HxVimp.f90"
       !
       !LOCAL INTERACTION
       include "direct/HxVint.f90"
       !
       !BATH HAMILTONIAN
       include "direct/HxVbath.f90"
       !
       !IMPURITY- BATH HYBRIDIZATION
       include "direct/HxVimp_bath.f90"

    enddo states
    !
    return
  end subroutine directMatVec_MPI_superc_main
#endif


END MODULE ED_HAMILTONIAN_DIRECT_HXV

! > BUILD SPARSE HAMILTONIAN of the SECTOR
MODULE ED_HAMILTONIAN_SUPERC_STORED_HxV
  USE ED_HAMILTONIAN_SUPERC_COMMON
  implicit none
  private

  !>Sparse Matric constructors
  public :: ed_buildH_superc_main


  !>Sparse Mat-Vec product using stored sparse matrix 
  public  :: spMatVec_superc_main
#ifdef _MPI
  public  :: spMatVec_MPI_superc_main
#endif


contains



  !####################################################################
  !             BUILD SPARSE HAMILTONIAN of the SECTOR
  !####################################################################
  subroutine ed_buildH_superc_main(Hmat)
    !
    !
    ! Builds the sector Hamiltonian :math:`H` and save each term in a suitable sparse matrix instance for :f:var:`ed_total_ed` = :code:`True`. If the dimension :f:var:`dim` of the sector are smaller than :f:var:`lanc_dim_threshold` the global matrix is dumped to the optional variable :f:var:`hmat`.
    !
    ! All the different electronic terms  are collected in the same sparse matrix, possibly using rows splitting and local / non-local blocks according to the :f:var:`MPI_Allgatherv` algorithm: 
    !  * :math:`H_{\rm int} \rightarrow` :f:var:`sph0` : interaction part of the electronic Hamiltonian
    !  * :math:`H_{\rm imp} \rightarrow` :f:var:`sph0` : impurity part of the eletronic Hamiltonian 
    !  * :math:`H_{\rm bath} \rightarrow` :f:var:`sph0`: bath levels part of the eletronic Hamiltonian
    !  * :math:`H_{\rm hyb} \rightarrow` :f:var:`sph0` : impurity - bath coupling part of the eletronic Hamiltonian
    !
    integer                                                       :: isector
    complex(8),dimension(:,:),optional                            :: Hmat  !optional dense matrix
    complex(8),dimension(:,:),allocatable                         :: Htmp_e,Hrdx
    integer,dimension(Nlevels)                                    :: ib
    integer,dimension(Ns)                                         :: ibup,ibdw
    real(8),dimension(Nimp)                                       :: nup,ndw
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,Nbath) :: Hbath_tmp

    !   
#ifdef _MPI
    if(Mpistatus .AND. MpiComm == MPI_COMM_NULL)return
#endif
    !
    if(.not.Hsector%status)stop "ed_buildh_main ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    Dim   = Hsector%Dim
    !
    if(present(Hmat))&
         call assert_shape(Hmat,[getdim(isector), getdim(isector)],"ed_buildh_main","Hmat")
    !
    if(allocated(diag_hybr))deallocate(diag_hybr)
    if(allocated(bath_diag))deallocate(bath_diag)
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
#ifdef _MPI
    if(MpiStatus)then
       call sp_set_mpi_matrix(MpiComm,spH0,mpiIstart,mpiIend,mpiIshift)
    else
       call sp_init_matrix(spH0,Dim)
    endif
#else
    call sp_init_matrix(spH0,Dim)
#endif
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN TERMS
    include "stored/Himp.f90"
    !
    !LOCAL INTERACTION
    include "stored/Hint.f90"
    !
    !BATH HAMILTONIAN
    include "stored/Hbath.f90"
    !
    !IMPURITY- BATH HYBRIDIZATION
    include "stored/Himp_bath.f90"
    !-----------------------------------------------!
    !
    if(present(Hmat))then
       Hmat = zero
       allocate(Htmp_e(DimEl,DimEl)); Htmp_e=0.d0
       !
#ifdef _MPI
       if(MpiStatus)then
          call sp_dump_matrix(MpiComm,spH0,Htmp_e)
       else
          call sp_dump_matrix(spH0,Htmp_e)
       endif
#else
       call sp_dump_matrix(spH0,Htmp_e)
#endif
       !
       Hmat = Htmp_e
       !
       deallocate(Htmp_e)
    endif
    !
    deallocate(diag_hybr,bath_diag)
    !
    return
    !
  end subroutine ed_buildH_superc_main






  !####################################################################
  !        SPARSE MAT-VEC PRODUCT USING STORED SPARSE MATRIX 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Perform the matrix-vector product H*v used in the
  ! - serial
  ! - MPI
  !+------------------------------------------------------------------+
  subroutine spMatVec_superc_main(Nloc,v,Hv)
    !
    ! Serial version of the matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in Arpack/Lanczos algorithm.
    ! This procedures applies one by one each term of the global Hamiltonian to an input vector using the stored sparse matrices.  
    !
    integer                    :: Nloc !Global dimension of the problem. :code:`size(v)=Nloc=size(Hv)`
    complex(8),dimension(Nloc) :: v    !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc) :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    complex(8)                 :: val
    integer                    :: i,iup,idw,j,jup,jdw,jj,i_el,j_el
    !
    !
    Hv=zero
    ! Electron Part
    do i=1,Nloc
       matmul: do jj=1, spH0%row(i)%Size
          val = spH0%row(i)%cvals(jj)
          j   = spH0%row(i)%cols(jj)
          Hv(i) = Hv(i) + val*v(j)
       end do matmul
    end do
    !
  end subroutine spMatVec_superc_main

#ifdef _MPI
  subroutine spMatVec_mpi_superc_main(Nloc,v,Hv)
    !
    !
    ! MPI parallel version of the matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in P-Arpack/P-Lanczos algorithm.
    ! This procedures applies one by one each term of the global Hamiltonian to a part of the vector own by the thread using the stored sparse matrices.
    !
    integer                             :: Nloc !Local dimension of the vector chunk. :code:`size(v)=Nloc` with :math:`\sum_p` :f:var:`Nloc` = :f:var:`Dim`
    complex(8),dimension(Nloc)          :: v    !input vector part (passed by P-Arpack/P-Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc)          :: Hv   !!output vector (required by P-Arpack/P-Lanczos) :math:`\vec{w}`
    complex(8)                          :: val
    integer                             :: i,j,mpiIerr,iph,jp, i_el,j_el
    integer                             :: N,MpiShift
    complex(8),dimension(:),allocatable :: vin
    integer,allocatable,dimension(:)    :: Counts,Offset
    !
    !
    if(MpiComm==MPI_UNDEFINED)stop "spHtimesV_mpi_cc ERRROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "spMatVec_mpi_cc ERROR: MpiStatus = F"
    !
    MpiRank = get_Rank_MPI(MpiComm)
    MpiSize = get_Size_MPI(MpiComm)
    !
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    !
    !Evaluate the local contribution: Hv_loc = Hloc*v
    MpiShift = spH0%Ishift
    Hv=0d0
    do i=1,Nloc
       local: do j=1,spH0%loc(i)%Size
          Hv(i) = Hv(i) + spH0%loc(i)%cvals(j)*v(spH0%loc(i)%cols(j)-MpiShift)
       end do local
    end do
    !
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
    allocate(vin(N)) ; vin = zero
    call MPI_Allgatherv(&
         v(1:Nloc),Nloc,MPI_Double_Complex,&
         vin      ,Counts,Offset,MPI_Double_Complex,&
         MpiComm,MpiIerr)
    ! Electron Part
    do i=1,Nloc                 !==spH0%Nrow
       matmul: do j=1,spH0%row(i)%Size
          Hv(i) = Hv(i) + spH0%row(i)%cvals(j)*vin(spH0%row(i)%cols(j))
       end do matmul
    end do
    !
  end subroutine spMatVec_mpi_superc_main
#endif



end MODULE ED_HAMILTONIAN_SUPERC_STORED_HXV








MODULE ED_HAMILTONIAN_NORMAL_STORED_HxV
  !Constructs each terms of the sector Hamiltonian storing them into different :f:var:`sparse_matrix` instances, implement the corresponding matrix-vector products using stored sparse matrices. 
  !
  USE ED_HAMILTONIAN_NORMAL_COMMON
  implicit none
  private


  !>Sparse Matric constructors
  public :: ed_buildh_normal_main

  !>Sparse Mat-Vec product using stored sparse matrix
  public  :: spMatVec_normal_main
#ifdef _MPI
  public  :: spMatVec_MPI_normal_main
#endif


contains



  subroutine ed_buildh_normal_main(Hmat)
    !
    ! Builds the sector Hamiltonian :math:`H` and save each term in a suitable sparse matrix instance for :f:var:`ed_total_ed` = :code:`True`. If the dimension :f:var:`dim` of the sector are smaller than :f:var:`lanc_dim_threshold` the global matrix is dumped to the optional variable :f:var:`hmat`.
    !
    ! The sparse matrices are:
    !  * :math:`H_d \rightarrow` :f:var:`sph0d` : diagonal part of the electronic Hamiltonian
    !  * :math:`H_\uparrow \rightarrow` :f:var:`sph0ups` : :math:`\uparrow` spin terms of the eletronic Hamiltonian 
    !  * :math:`H_\downarrow \rightarrow` :f:var:`sph0dws` : :math:`\downarrow`  spin terms of the eletronic Hamiltonian
    !  * :math:`H_{nd} \rightarrow` :f:var:`sph0nd` : non-diagonal part of the eletronic Hamiltonian
    !  * :math:`H_{ph} \rightarrow` :f:var:`sph0_ph` : phonon part of the of the global Hamiltonian
    !  * :math:`H_{e-eph} \rightarrow` :f:var:`sph0e_eph` : electron part of the electron-phonon term of the global Hamiltonian
    !  * :math:`H_{ph_eph} \rightarrow` :f:var:`sph0e_ph` : phonon part of the electron-phonon term of the global Hamiltonian
    !
    complex(8),dimension(:,:),optional                            :: Hmat !optional dense matrix
    integer                                                       :: isector   
    complex(8),dimension(:,:),allocatable                         :: Htmp_up,Htmp_dw,Hrdx
    integer,dimension(Ns)                                         :: ibup,ibdw
    integer,dimension(Ns)                                         :: Nup,Ndw
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,Nbath) :: Hbath_tmp
    !
    nup=zero
    ndw=zero
#ifdef _MPI
    if(Mpistatus .AND. MpiComm == MPI_COMM_NULL)return
#endif
    !
    if(.not.Hsector%status)stop "ed_buildh_main ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    if(present(Hmat))&
         call assert_shape(Hmat,[getdim(isector), getdim(isector)],"ed_buildh_main","Hmat")
    !
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
#ifdef _MPI
    if(MpiStatus)then
       call sp_set_mpi_matrix(MpiComm,spH0d,mpiIstart,mpiIend,mpiIshift)
       call sp_init_matrix(MpiComm,spH0d,Dim)
       if(Jhflag)then
          call sp_set_mpi_matrix(MpiComm,spH0nd,mpiIstart,mpiIend,mpiIshift)
          call sp_init_matrix(MpiComm,spH0nd,Dim)
       endif
    else
       call sp_init_matrix(spH0d,Dim)
       if(Jhflag)call sp_init_matrix(spH0nd,Dim)
    endif
#else
    call sp_init_matrix(spH0d,Dim)
    if(Jhflag)call sp_init_matrix(spH0nd,Dim)
#endif
    call sp_init_matrix(spH0dws(1),DimDw)
    call sp_init_matrix(spH0ups(1),DimUp)
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN TERMS
    print*,"entering Hlocal"
    include "stored/H_local.f90"
    !
    !UP TERMS
    print*,"entering Hup"
    include "stored/H_up.f90"
    !
    !DW TERMS
    print*,"entering Hdw"
    include "stored/H_dw.f90"
    !
    !NON-LOCAL HAMILTONIAN TERMS
    if(jhflag)then
       print*,"entering H_non_local"
       include "stored/H_non_local.f90"
    endif
    !
    !-----------------------------------------------!
    !
    if(present(Hmat))then
       Hmat = zero
       allocate(Htmp_up(DimUp,DimUp));Htmp_up=zero
       allocate(Htmp_dw(DimDw,DimDw));Htmp_dw=zero
       !
#ifdef _MPI
       if(MpiStatus)then
          call sp_dump_matrix(MpiComm,spH0d,Hmat)
       else
          call sp_dump_matrix(spH0d,Hmat)
       endif
#else
       call sp_dump_matrix(spH0d,Hmat)
#endif
       !
       if(Jhflag)then
          allocate(Hrdx(Dim,Dim));Hrdx=zero
#ifdef _MPI
          if(MpiStatus)then
             call sp_dump_matrix(MpiComm,spH0nd,Hrdx)
          else
             call sp_dump_matrix(spH0nd,Hrdx)
          endif
#else
          call sp_dump_matrix(spH0nd,Hrdx)
#endif
          Hmat = Hmat + Hrdx
          deallocate(Hrdx)
       endif
       !
       call sp_dump_matrix(spH0ups(1),Htmp_up)
       call sp_dump_matrix(spH0dws(1),Htmp_dw)
       Hmat = Hmat + kronecker_product(Htmp_dw,one*eye(DimUp))
       Hmat = Hmat + kronecker_product(one*eye(DimDw),Htmp_up)
       !
       deallocate(Htmp_up,Htmp_dw)
    endif
    !
    return
    !
  end subroutine ed_buildh_normal_main






  !####################################################################
  !        SPARSE MAT-VEC PRODUCT USING STORED SPARSE MATRIX 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Perform the matrix-vector product H*v used in the
  ! - serial
  ! - MPI
  !+------------------------------------------------------------------+
  subroutine spMatVec_normal_main(Nloc,v,Hv)
    !
    ! Serial version of the matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in Arpack/Lanczos algorithm for :f:var:`ed_total_ud` = :code:`True` 
    ! This procedures applies one by one each term of the global Hamiltonian to an input vector using the stored sparse matrices.  
    !
    integer                    :: Nloc !Global dimension of the problem. :code:`size(v)=Nloc=size(Hv)`
    complex(8),dimension(Nloc) :: v    !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc) :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    complex(8)                 :: val
    integer                    :: i,iup,idw,j,jup,jdw,jj
    !
    !
    Hv=zero
    !
    !Local:
    do i = 1,Nloc
       do j=1,spH0d%row(i)%Size
          Hv(i) = Hv(i) + spH0d%row(i)%vals(j)*v(spH0d%row(i)%cols(j))
       enddo
    enddo
    !
    !DW:
    do iup=1,DimUp
       !
       do idw=1,DimDw
          i = iup + (idw-1)*DimUp
          do jj=1,spH0dws(1)%row(idw)%Size
             jup = iup
             jdw = spH0dws(1)%row(idw)%cols(jj)
             val = spH0dws(1)%row(idw)%vals(jj)
             j     = jup +  (jdw-1)*DimUp
             Hv(i) = Hv(i) + val*V(j)
          enddo
       enddo
       !
    enddo
    !
    !UP:
    do idw=1,DimDw
       !
       do iup=1,DimUp
          i = iup + (idw-1)*DimUp
          do jj=1,spH0ups(1)%row(iup)%Size
             jup = spH0ups(1)%row(iup)%cols(jj)
             jdw = idw
             val = spH0ups(1)%row(iup)%vals(jj)
             j =  jup + (jdw-1)*DimUp
             Hv(i) = Hv(i) + val*V(j)
          enddo
       enddo
       !
    enddo
    !
    !Non-Local:
    if(jhflag)then
       do i = 1,Nloc
          do j=1,spH0nd%row(i)%Size
             val   = spH0nd%row(i)%vals(j)
             jj    = spH0nd%row(i)%cols(j)
             Hv(i) = Hv(i) + val*V(jj)
          enddo
       enddo
    endif
    !
  end subroutine spMatVec_normal_main

#ifdef _MPI
  subroutine spMatVec_mpi_normal_main(Nloc,v,Hv)
    !
    ! MPI parallel version of the matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in P-Arpack/P-Lanczos algorithm for :f:var:`ed_total_ud` = :code:`True`. 
    ! This procedures applies one by one each term of the global Hamiltonian to a part of the vector own by the thread using the stored sparse matrices.
    !
    integer                             :: Nloc !Local dimension of the vector chunk. :code:`size(v)=Nloc` with :math:`\sum_p` :f:var:`Nloc` = :f:var:`Dim`
    complex(8),dimension(Nloc)          :: v    !input vector part (passed by P-Arpack/P-Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc)          :: Hv   !output vector (required by P-Arpack/P-Lanczos) :math:`\vec{w}`
    !
    integer                             :: N
    complex(8),dimension(:),allocatable :: vt,Hvt
    complex(8),dimension(:),allocatable :: vin
    complex(8)                          :: val
    integer                             :: i,iup,idw,j,jup,jdw,jj
    !local MPI
    integer                             :: irank,MpiIerr
    integer,allocatable,dimension(:)    :: Counts
    integer,allocatable,dimension(:)    :: Offset
    !
    ! if(MpiComm==Mpi_Comm_Null)return
    ! if(MpiComm==MPI_UNDEFINED)stop "spMatVec_mpi_cc ERROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "spMatVec_mpi_cc ERROR: MpiStatus = F"
    !
    !Evaluate the local contribution: Hv_loc = Hloc*v
    Hv=zero
    do i=1,Nloc                 !==spH0%Nrow
       do j=1,spH0d%row(i)%Size
          Hv(i) = Hv(i) + spH0d%row(i)%vals(j)*v(i)
       end do
    end do
    !
    !Non-local terms.
    !UP part: contiguous in memory.
    do idw=1,MpiQdw
       do iup=1,DimUp
          i = iup + (idw-1)*DimUp
          hxv_up: do jj=1,spH0ups(1)%row(iup)%Size
             jup = spH0ups(1)%row(iup)%cols(jj)
             jdw = idw
             val = spH0ups(1)%row(iup)%vals(jj)
             j   = jup + (idw-1)*DimUp
             Hv(i) = Hv(i) + val*v(j)
          end do hxv_up
       enddo
    end do
    !
    !DW part: non-contiguous in memory -> MPI transposition
    !Transpose the input vector as a whole:
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    !
    allocate(vt(mpiQup*DimDw)) ;vt=zero
    allocate(Hvt(mpiQup*DimDw));Hvt=zero
    call vector_transpose_MPI(DimUp,MpiQdw,v,DimDw,MpiQup,vt)
    Hvt=0d0    
    do idw=1,MpiQup             !<= Transposed order:  column-wise DW <--> UP  
       do iup=1,DimDw           !<= Transposed order:  column-wise DW <--> UP
          i = iup + (idw-1)*DimDw
          hxv_dw: do jj=1,spH0dws(1)%row(iup)%Size
             jup = spH0dws(1)%row(iup)%cols(jj)
             jdw = idw             
             j   = jup + (jdw-1)*DimDw
             val = spH0dws(1)%row(iup)%vals(jj)
             Hvt(i) = Hvt(i) + val*vt(j)
          end do hxv_dw
       enddo
    end do
    deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ; vt=zero
    call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt)
    Hv = Hv + Vt
    deallocate(vt)
    !
    !
    !Non-Local:
    if(jhflag)then
       N = 0
       call AllReduce_MPI(MpiComm,Nloc,N)
       ! 
       allocate(vt(N)) ; vt = zero
       call allgather_vector_MPI(MpiComm,v,vt)
       !
       do i=1,Nloc
          matmul: do j=1,spH0nd%row(i)%Size
             Hv(i) = Hv(i) + spH0nd%row(i)%vals(j)*Vt(spH0nd%row(i)%cols(j))
          enddo matmul
       enddo
       deallocate(Vt)
    endif
    !
  end subroutine spMatVec_mpi_normal_main
#endif



end MODULE ED_HAMILTONIAN_NORMAL_STORED_HXV








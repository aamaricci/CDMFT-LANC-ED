MODULE ED_BATH_AUX
  !Implements a number of auxiliary procedures used to construct replica/general bath
  !
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,str
  USE SF_LINALG, only: eye,inv
  USE SF_ARRAYS, only:linspace
  USE SF_MISC, only: assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  implicit none

  public :: hreplica_build                   !INTERNAL (for effective_bath)
  public :: hreplica_mask                    !INTERNAL (for effective_bath)
  public :: hreplica_site                    !INTERNAL (for effective_bath)
  !
  public :: hgeneral_build                   !INTERNAL (for effective_bath)
  public :: hgeneral_mask                    !INTERNAL (for effective_bath)
  public :: hgeneral_site                    !INTERNAL (for effective_bath)





  interface get_Whyb_matrix
     ! This subroutine build up the hybridization matrix :math:`W_{\sigma\sigma^{'}}` used in the :f:var:`ed_mode` = nonsu2 with :f:var:`bath_type` = hybrid. The input can have different shape and type:
     !
     !   * :f:var:`u` , :f:var:`v` with dimensions [ |Nspin| , |Norb| ]  
     !   * :f:var:`u` , :f:var:`v` with dimensions [ |Nspin| , |Norb| , :f:var:`nbath`]  
     !   * :f:var:`dmft_bath` 
     !
     module procedure get_Whyb_matrix_1orb
     module procedure get_Whyb_matrix_Aorb
     module procedure get_Whyb_matrix_dmft_bath
  end interface get_Whyb_matrix


  interface is_identity
     !
     ! This subroutine checks if a matrix :math:`\hat{O}`  in the basis of the :code:`replica` or :code:`general` baths is the identity.
     ! 
     ! The input matrix can have different shapes:
     !    *  [ |Nnambu| . |Nspin| . |Norb| , |Nnambu| . |Nspin| . |Norb| ]
     !    *  [ |Nnambu| . |Nspin| , |Nnambu| . |Nspin| , |Norb| , |Norb| ]
     !
     module procedure ::  is_identity_so
     module procedure ::  is_identity_nn
  end interface is_identity

  interface is_diagonal
     !
     ! This subroutine checks if a matrix :math:`\hat{O}`  in the basis of the :code:`replica` or :code:`general` baths is diagonal.
     ! 
     ! The input matrix can have different shapes:
     !    *  [ |Nnambu| . |Nspin| . |Norb| , |Nnambu| . |Nspin| . |Norb| ]
     !    *  [ |Nnambu| . |Nspin| , |Nnambu| . |Nspin| , |Norb| , |Norb| ]
     !
     module procedure ::  is_diagonal_so
     module procedure ::  is_diagonal_nn
  end interface is_diagonal




contains



  subroutine Hreplica_site(site)
    integer :: site
    if(site<1.OR.site>size(Hreplica_lambda_ineq,1))stop "ERROR Hreplica_site: site not in [1,Nlat]"
    if(.not.allocated(Hreplica_lambda_ineq))stop "ERROR Hreplica_site: Hreplica_lambda_ineq not allocated"
    Hreplica_lambda(:,:)  = Hreplica_lambda_ineq(site,:,:)
  end subroutine Hreplica_site



  function Hreplica_build(lambdavec) result(H)
    !
    !This function is used to reconstruct the local bath Hamiltonian from basis expansion given the vector of :math:`\vec{\lambda}` parameters :math:`h^p=\sum_i \lambda^p_i O_i`. The resulting Hamiltonian has dimensions [ |Nspin| , |Nspin| , |Norb| , |Norb| ]
    !
    real(8),dimension(:),optional                             :: lambdavec !the input vector of bath parameters
    real(8),dimension(:),allocatable                          :: lambda
    integer                                                   :: isym
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: H
    !
    if(.not.Hreplica_status)STOP "ERROR Hreplica_build: Hreplica_basis is not setup"
    allocate(lambda(size(Hreplica_basis)));lambda=1d0
    if(present(lambdavec))then
       if(size(lambdavec)/=size(Hreplica_basis)) STOP "ERROR Hreplica_build: Wrong coefficient vector size"
       lambda = lambdavec
    endif
    H=zero
    do isym=1,size(lambda)
       H=H+lambda(isym)*Hreplica_basis(isym)%O
    enddo
  end function Hreplica_build

  function Hreplica_mask(wdiag,uplo) result(Hmask)
    logical,optional                                          :: wdiag,uplo
    logical                                                   :: wdiag_,uplo_
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: H
    logical,dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb)    :: Hmask
    integer                                                   :: iorb,jorb,ispin,jspin,io,jo
    !
    wdiag_=.false.;if(present(wdiag))wdiag_=wdiag
    uplo_ =.false.;if(present(uplo))  uplo_=uplo
    !
    H = Hreplica_build(Hreplica_lambda(Nbath,:)) 
    Hmask=.false.
    where(abs(H)>1d-6)Hmask=.true.
    !
    !
    if(wdiag_)then
       do ispin=1,Nnambu*Nspin
          do iorb=1,Norb
             Hmask(ispin,ispin,iorb,iorb)=.true.
          enddo
       enddo
    endif
    !
    if(uplo_)then
       do ispin=1,Nnambu*Nspin
          do jspin=1,Nnambu*Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = index_stride_so(ispin,iorb)
                   jo = index_stride_so(jspin,jorb)
                   if(io>jo)Hmask(ispin,jspin,iorb,jorb)=.false.
                enddo
             enddo
          enddo
       enddo
    endif
    !
  end function Hreplica_mask





  subroutine Hgeneral_site(site)
    integer :: site
    if(site<1.OR.site>size(Hgeneral_lambda_ineq,1))stop "ERROR Hgeneral_site: site not in [1,Nlat]"
    if(.not.allocated(Hgeneral_lambda_ineq))stop "ERROR Hgeneral_site: Hgeneral_lambda_ineq not allocated"
    Hgeneral_lambda(:,:)  = Hgeneral_lambda_ineq(site,:,:)
  end subroutine Hgeneral_site

  function Hgeneral_build(lambdavec) result(H)
    !
    !This function is used to reconstruct the local bath Hamiltonian from basis expansion given the vector of :math:`\vec{\lambda}` parameters :math:`h^p=\sum_i \lambda^p_i O_i`. The resulting Hamiltonian has dimensions [ |Nspin| , |Nspin| , |Norb| , |Norb| ]
    !
    real(8),dimension(:),optional                             :: lambdavec  !the input vector of bath parameters
    real(8),dimension(:),allocatable                          :: lambda
    integer                                                   :: isym
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: H
    !
    if(.not.Hgeneral_status)STOP "ERROR Hgeneral_build: Hgeneral_basis is not setup"
    allocate(lambda(size(Hgeneral_basis)));lambda=1d0
    if(present(lambdavec))then
       if(size(lambdavec)/=size(Hgeneral_basis)) STOP "ERROR Hgeneral_build: Wrong coefficient vector size"
       lambda = lambdavec
    endif
    H=zero
    do isym=1,size(lambdavec)
       H=H+lambdavec(isym)*Hgeneral_basis(isym)%O
    enddo
  end function Hgeneral_build

  function Hgeneral_mask(wdiag,uplo) result(Hmask)
    logical,optional                                          :: wdiag,uplo
    logical                                                   :: wdiag_,uplo_
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: H
    logical,dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb)    :: Hmask
    integer                                                   :: iorb,jorb,ispin,jspin,io,jo
    !
    wdiag_=.false.;if(present(wdiag))wdiag_=wdiag
    uplo_ =.false.;if(present(uplo))  uplo_=uplo
    !
    H = Hgeneral_build(Hgeneral_lambda(Nbath,:))
    Hmask=.false.
    where(abs(H)>1d-6)Hmask=.true.
    !
    !
    if(wdiag_)then
       do ispin=1,Nnambu*Nspin
          do iorb=1,Norb
             Hmask(ispin,ispin,iorb,iorb)=.true.
          enddo
       enddo
    endif
    !
    if(uplo_)then
       do ispin=1,Nnambu*Nspin
          do jspin=1,Nnambu*Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = index_stride_so(ispin,iorb)
                   jo = index_stride_so(jspin,jorb)
                   if(io>jo)Hmask(ispin,jspin,iorb,jorb)=.false.
                enddo
             enddo
          enddo
       enddo
    endif
    !
  end function Hgeneral_mask











  function get_Whyb_matrix_1orb(v,u) result(w)
    real(8),dimension(Nspin,Nbath)       :: v,u
    real(8),dimension(Nspin,Nspin,Nbath) :: w
    integer                              :: ispin
    do ispin=1,Nspin
       w(ispin,ispin,:) = v(ispin,:)
    enddo
    w(1,Nspin,:) = u(1,:)
    w(Nspin,1,:) = u(2,:)
  end function get_Whyb_matrix_1orb

  function get_Whyb_matrix_Aorb(v,u) result(w)
    real(8),dimension(Nspin,Norb,Nbath)       :: v,u
    real(8),dimension(Nspin,Nspin,Norb,Nbath) :: w
    integer                                   :: ispin
    do ispin=1,Nspin
       w(ispin,ispin,:,:) = v(ispin,:,:)
    enddo
    w(1,Nspin,:,:) = u(1,:,:)
    w(Nspin,1,:,:) = u(2,:,:)
  end function get_Whyb_matrix_Aorb

  function get_Whyb_matrix_dmft_bath(dmft_bath_) result(w)
    type(effective_bath)                      :: dmft_bath_
    real(8),dimension(Nspin,Nspin,Norb,Nbath) :: w
    integer                                   :: ispin
    !
    do ispin=1,Nspin
       w(ispin,ispin,:,:) = dmft_bath_%v(ispin,:,:)
    enddo
    w(1,Nspin,:,:) = dmft_bath_%u(1,:,:)
    w(Nspin,1,:,:) = dmft_bath_%u(2,:,:)
  end function get_Whyb_matrix_dmft_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if a matrix is the identity
  !+-------------------------------------------------------------------+
  function is_identity_nn(mnnn) result(flag)
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: mnnn
    real(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mtmp
    integer                                  :: i,j
    logical                                  :: flag
    !
    flag=.true.
    !
    mtmp=nn2so_reshape(mnnn,Nnambu*Nspin,Norb)
    !
    do i=1,Nnambu*Nspin*Norb-1
       if((mtmp(i,i).ne.mtmp(i+1,i+1)).or.(mtmp(i,i).lt.1.d-6))flag=.false.
    enddo
    !
    do i=1,Nnambu*Nspin*Norb
       do j=1,Nnambu*Nspin*Norb
          if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_identity_nn

  function is_identity_so(mlso) result(flag)
    complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mlso
    real(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mtmp
    integer                                  :: i,j
    logical                                  :: flag
    !
    flag=.true.
    !
    mtmp=mlso
    !
    do i=1,Nnambu*Nspin*Norb-1
       if((mtmp(i,i).ne.mtmp(i+1,i+1)).or.(mtmp(i,i).lt.1.d-6))flag=.false.
    enddo
    !
    do i=1,Nnambu*Nspin*Norb
       do j=1,Nnambu*Nspin*Norb
          if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_identity_so



  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if a matrix is diagonal
  !+-------------------------------------------------------------------+
  function is_diagonal_nn(mnnn) result(flag)
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: mnnn
    complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mtmp
    integer                                     :: i,j
    logical                                     :: flag
    !
    flag=.true.
    !
    mtmp=abs((nn2so_reshape(mnnn,Nnambu*Nspin,Norb)))
    !
    do i=1,Nnambu*Nspin*Norb
       do j=1,Nnambu*Nspin*Norb
          if((i.ne.j).and.(abs(mtmp(i,j)).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_diagonal_nn

  function is_diagonal_so(mlso) result(flag)
    complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mlso
    complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mtmp
    integer                                     :: i,j
    logical                                     :: flag
    !
    flag=.true.
    !
    mtmp=abs((mlso))
    !
    do i=1,Nnambu*Nspin*Norb
       do j=1,Nnambu*Nspin*Norb
          if((i.ne.j).and.(abs(mtmp(i,j)).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_diagonal_so






  function check_herm(A,N,error) result(bool)
    integer,intent(in)                   :: N
    complex(8),dimension(N,N),intent(in) :: A
    logical                              :: bool
    real(8),optional                     :: error
    real(8)                              :: error_
    error_ = 1d-6 ; if(present(error))error_=error
    bool   = all(abs(A - conjg(transpose(A)))<error_)
  end function check_herm


  function check_nambu(A,N,error) result(bool)
    integer,intent(in)                       :: N
    complex(8),dimension(2*N,2*N),intent(in) :: A
    complex(8),dimension(N,N)                :: h11,h22
    logical                                  :: bool
    real(8),optional                         :: error
    real(8)                                  :: error_
    error_ = 1d-6 ; if(present(error))error_=error
    h11    = A(1:N    ,1:N)
    h22    = A(N+1:2*N,N+1:2*N)
    bool   = check_herm(A,2*N,error_) !this checks also for F = A_12, s.t. A_21=herm(A_12)
    bool   = bool.AND.( all(abs(h22 + conjg(h11))<error_) )
  end function check_nambu





END MODULE ED_BATH_AUX


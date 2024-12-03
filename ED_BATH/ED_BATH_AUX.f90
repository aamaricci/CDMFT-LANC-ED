MODULE ED_BATH_AUX
  !Implements a number of auxiliary procedures used to construct replica/general bath
  !
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,str
  USE SF_LINALG, only: eye,inv,trace
  USE SF_MISC, only: assert_shape
  USE SF_ARRAYS, only: linspace
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  implicit none

  private




  interface is_identity
     !
     ! This subroutine checks if a matrix :math:`\hat{O}`  in the basis of the :code:`replica` or :code:`general` baths is the identity.
     ! 
     ! The input matrix can have different shapes:
     !    *  [ |Nlat| . |Nnambu| . |Nspin| . |Norb| , |Nlat| . |Nnambu| . |Nspin| . |Norb| ]
     !    *  [ |Nlat| , |Nlat| , |Nnambu| . |Nspin| , |Nnambu| . |Nspin| , |Norb| , |Norb| ]
     !
     module procedure ::  is_identity_d6
     module procedure ::  is_identity_d4
     module procedure ::  is_identity_d2
  end interface is_identity

  interface is_diagonal
     !
     ! This subroutine checks if a matrix :math:`\hat{O}`  in the basis of the :code:`replica` or :code:`general` baths is diagonal.
     !
     ! The input matrix can have different shapes:
     !    *  [ |Nlat| . |Nnambu| . |Nspin| . |Norb| , |Nlat| . |Nnambu| . |Nspin| . |Norb| ]
     !    *  [ |Nlat| , |Nlat| , |Nnambu| . |Nspin| , |Nnambu| . |Nspin| , |Norb| , |Norb| ]
     !
     module procedure ::  is_diagonal_d6
     module procedure ::  is_diagonal_d4
     module procedure ::  is_diagonal_d2
  end interface is_diagonal


  interface check_herm
     module procedure :: check_herm_d2
     module procedure :: check_herm_d4
     module procedure :: check_herm_d6
  end interface check_herm

  interface check_nambu
     module procedure :: check_nambu_d2
     module procedure :: check_nambu_d6
  end interface check_nambu


  public :: is_identity
  public :: is_diagonal
  public :: check_herm
  public :: check_nambu


  integer :: inam,jnam
  integer :: ilat,jlat
  integer :: iorb,jorb
  integer :: ispin,jspin
  integer :: is,js
  integer :: io,jo
  integer :: i,j


contains




  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if a matrix is diagonal
  !+-------------------------------------------------------------------+

  function is_diagonal_d6(M6) result(flag)
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp) :: M6
    real(8),dimension(Nambu*Nspin*Nimp,Nambu*Nspin*Nimp)    :: M2
    logical                                                 :: flag
    !
    flag=.true.
    !
    if ( any(abs(dimag(M6)) > 1d-6) ) flag=.false.
    !
    do concurrent(inam=1:Nambu,jnam=1:Nambu,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Nimp,jorb=1:Nimp)
       i = iorb + (ispin-1)*Nimp + (inam-1)*Nspin*Nimp
       j = jorb + (jspin-1)*Nimp + (jnam-1)*Nspin*Nimp
       M2(i,j) = abs(dreal(M6(inam,jnam,ispin,jspin,iorb,jorb)))
    enddo
    !
    do i=1,Nambu*Nspin*Nimp
       do j=i+1,Nambu*Nspin*Nimp
          if((M2(i,j)>1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_diagonal_d6

  function is_diagonal_d4(M4) result(flag)
    complex(8),dimension(Nspin,Nspin,Nimp,Nimp) :: M4
    real(8),dimension(Nspin*Nimp,Nspin*Nimp)    :: M2
    logical                                     :: flag
    !
    flag=.true.
    !
    if ( any(abs(dimag(M4)) > 1d-6) ) flag=.false.
    !
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Nimp,jorb=1:Nimp)
       i = iorb + (ispin-1)*Nimp
       j = jorb + (jspin-1)*Nimp
       M2(i,j) = abs(dreal(M4(ispin,jspin,iorb,jorb)))
    enddo
    !
    do i=1,Nspin*Nimp
       do j=i+1,Nspin*Nimp
          if((M2(i,j)>1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_diagonal_d4

  function is_diagonal_d2(M2) result(flag)
    complex(8),dimension(Nspin*Nimp,Nspin*Nimp) :: M2
    logical                                     :: flag
    !
    flag=.true.
    !
    if ( any(abs(dimag(M2)) > 1d-6) ) flag=.false.
    !
    do i=1,Nspin*Nimp
       do j=i+1,Nspin*Nimp
          if( abs(dreal(M2(i,j)))>1.d-6 )flag=.false.
       enddo
    enddo
    !
  end function is_diagonal_d2





  function is_identity_d6(M6) result(flag)
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp) :: M6
    real(8),dimension(Nambu*Nspin*Nimp,Nambu*Nspin*Nimp)    :: M2
    logical                                                 :: flag
    !
    flag=.true.
    !
    if ( any(abs(dimag(M6)) > 1d-6) ) flag=.false.
    !
    do concurrent(inam=1:Nambu,jnam=1:Nambu,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Nimp,jorb=1:Nimp)
       i = iorb + (ispin-1)*Nimp + (inam-1)*Nspin*Nimp
       j = jorb + (jspin-1)*Nimp + (jnam-1)*Nspin*Nimp
       M2(i,j) = abs(dreal(M6(inam,jnam,ispin,jspin,iorb,jorb)))
    enddo
    !
    do i=1,Nambu*Nspin*Nimp-1
       if((M2(i,i)/=M2(i+1,i+1)).OR.(M2(i,i)<1.d-6))flag=.false.
    enddo
    !
    do i=1,Nambu*Nspin*Nimp
       do j=i+1,Nambu*Nspin*Nimp
          if(M2(i,j)>1.d-6)flag=.false.
       enddo
    enddo
    !
  end function is_identity_d6


  function is_identity_d4(M4) result(flag)
    complex(8),dimension(Nspin,Nspin,Nimp,Nimp) :: M4
    real(8),dimension(Nspin*Nimp,Nspin*Nimp)    :: M2
    logical                                     :: flag
    !
    flag=.true.
    !
    if ( any(abs(dimag(M4)) > 1d-6) ) flag=.false.
    !
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Nimp,jorb=1:Nimp)
       i = iorb + (ispin-1)*Nimp
       j = jorb + (jspin-1)*Nimp
       M2(i,j) = abs(dreal(M4(ispin,jspin,iorb,jorb)))
    enddo
    !
    do i=1,Nspin*Nimp-1
       if((M2(i,i)/=M2(i+1,i+1)).OR.(M2(i,i)<1.d-6))flag=.false.
    enddo
    !
    do i=1,Nspin*Nimp
       do j=i+1,Nspin*Nimp
          if(M2(i,j)>1.d-6)flag=.false.
       enddo
    enddo
    !
  end function is_identity_d4



  function is_identity_d2(M2) result(flag)
    complex(8),dimension(Nspin*Nimp,Nspin*Nimp) :: M2
    real(8),dimension(Nspin*Nimp,Nspin*Nimp) :: M2_
    logical                                     :: flag
    !
    flag=.true.
    !
    if ( ANY( abs(dimag(M2)) .gt. 1d-6 ) ) flag=.false.
    !
    M2_ = abs(dreal(M2))
    do i=1,Nspin*Nimp-1
       if(  (M2_(i,i)/=M2_(i+1,i+1)) .OR. &
            (M2_(i,i)<1.d-6) )flag=.false.
    enddo
    !
    do i=1,Nspin*Nimp
       do j=i+1,Nspin*Nimp
          if( M2_(i,j)>1.d-6 )flag=.false.
       enddo
    enddo
    !
  end function is_identity_d2








  !A[N,N]
  function check_herm_d2(A,error) result(bool)
    complex(8),dimension(:,:) :: A
    real(8),optional          :: error
    logical                   :: bool
    real(8)                   :: error_
    error_ = 1d-6 ; if(present(error))error_=error
    bool   = all(abs(A - conjg(transpose(A)))<error_)
  end function check_herm_d2

  !A[Nspin,Nspin,Nimp,Nimp]
  function check_herm_d4(A,error) result(bool)
    complex(8),dimension(Nspin,Nspin,Nimp,Nimp) :: A
    real(8),optional                            :: error
    complex(8),dimension(Nspin*Nimp,Nspin*Nimp) :: A_
    logical                                     :: bool
    real(8)                                     :: error_
    error_ = 1d-6 ; if(present(error))error_=error
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Nimp,jorb=1:Nimp)
       i = iorb + (ispin-1)*Nimp
       j = jorb + (jspin-1)*Nimp
       A_(i,j) = A(ispin,jspin,iorb,jorb)
    enddo
    bool = check_herm_d2(A_,error_)
  end function check_herm_d4

  !A[Nambu,Nambu,Nspin,Nspin,Nimp,Nimp]
  function check_herm_d6(A,error) result(bool)
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp) :: A
    real(8),optional                                        :: error
    complex(8),dimension(Nambu*Nspin*Nimp,Nambu*Nspin*Nimp) :: A_
    logical                                                 :: bool
    real(8)                                                 :: error_
    error_ = 1d-6 ; if(present(error))error_=error
    do concurrent(inam=1:Nambu,jnam=1:Nambu,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Nimp,jorb=1:Nimp)
       i = iorb + (ispin-1)*Nimp + (inam-1)*Nspin*Nimp
       j = jorb + (jspin-1)*Nimp + (jnam-1)*Nspin*Nimp
       A_(i,j) = A(inam,jnam,ispin,jspin,iorb,jorb)
    enddo
    bool = check_herm_d2(A_,error_)
  end function check_herm_d6






  function check_nambu_d2(A,N,error) result(bool)
    integer,intent(in)                       :: N
    complex(8),dimension(2*N,2*N),intent(in) :: A
    complex(8),dimension(N,N)                :: h11,h22
    logical                                  :: bool
    real(8),optional                         :: error
    real(8)                                  :: error_
    error_ = 1d-6 ; if(present(error))error_=error
    h11    = A(1:N    ,1:N)
    h22    = A(N+1:2*N,N+1:2*N)
    bool   = check_herm_d2(A,error_) !this checks also for F = A_12, s.t. A_21=herm(A_12)
    bool   = bool.AND.( all(abs(h22 + conjg(h11))<error_) )
  end function check_nambu_d2

  !A[Nambu,Nambu,Nspin,Nspin,Nimp,Nimp]
  function check_nambu_d6(A,error) result(bool)
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp) :: A
    real(8),optional                                        :: error
    complex(8),dimension(Nambu*Nspin*Nimp,Nambu*Nspin*Nimp) :: A_
    logical                                                 :: bool
    real(8)                                                 :: error_
    error_ = 1d-6 ; if(present(error))error_=error
    do concurrent(inam=1:Nambu,jnam=1:Nambu,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Nimp,jorb=1:Nimp)
       i = iorb + (ispin-1)*Nimp + (inam-1)*Nspin*Nimp
       j = jorb + (jspin-1)*Nimp + (jnam-1)*Nspin*Nimp
       A_(i,j) = A(inam,jnam,ispin,jspin,iorb,jorb)
    enddo
    bool = check_nambu_d2(A_,Nspin*Nimp,error_)
  end function check_nambu_d6


END MODULE ED_BATH_AUX

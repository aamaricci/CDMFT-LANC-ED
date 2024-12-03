MODULE ED_BATH_DIM
  !Returns or check the dimensions to which the user should allocate the bath array.
  !
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,str
  USE SF_LINALG, only: eye,inv,trace
  USE SF_MISC, only: assert_shape
  USE SF_ARRAYS, only: linspace
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  !
  USE ED_BATH_AUX
  USE ED_BATH_REPLICA
  implicit none

  private


  !##################################################################
  !
  !     BATH DIMENSION ROUTINES:
  !
  !##################################################################
  interface get_bath_dimension
     !
     ! Returns the dimension :f:var:`bath_size` to which the user should allocate the user bath array to contains all the parameters according to the provided input variables. The value is obtained counting all the electronic levels of the system compatible with the operational mode :f:var:`ed_mode`, the bath topology specified by :f:var:`bath_type`, the values of :f:var:`norb`, :f:var:`nbath` and :f:var:`nspin`.
     !
     ! If :f:var:`bath_type` is replica/general then a input matrix :f:var:`h_nn` can be used to count the number of parameters, corresponding to its non-zero elements. In alternative the bath size can be estimated by the number of parameters in the linear decomposition of the bath local Hamiltonian :f:var:`nsym` such that :math:`h^p=\sum_{i=1}^{N_{sym}}\lambda^{p}_{i} O_{i}`. 
     !
     !
     module procedure ::  get_bath_dimension_symmetries
     module procedure ::  get_bath_dimension_direct_d2 !only *normal* deprecated
     module procedure ::  get_bath_dimension_direct_d4 !only *normal* deprecated
     module procedure ::  get_bath_dimension_direct_d6 !only *normal* deprecated
  end interface get_bath_dimension

  public :: get_bath_dimension
  public :: check_bath_dimension



  integer :: ibath
  integer :: inam,jnam
  integer :: ilat,jlat
  integer :: iorb,jorb
  integer :: ispin,jspin
  integer :: is,js
  integer :: io,jo
  integer :: i,j
  integer :: in,jn


contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : Inquire the correct bath size to allocate the
  ! the bath array in the calling program.
  !+-------------------------------------------------------------------+
  function get_bath_dimension_symmetries(Nsym) result(bath_size)
    integer :: Nsym !Number of symmetries (for :f:var:`ed_mode` = :code:`replica, general` )
    integer :: bath_size
    integer :: ndx,isym
    !
    select case(bath_type)
    case("replica","general")
       if(.not.Hbath_status)STOP "get_bath_dimension_symmetries: Hbath_basis  not allocated"
       if(Nsym/=size(Hbath_lambda,2))stop "ERROR get_bath_dimension_symmetries:  size(Hbath_basis) != size(Hbath_lambda,2)"
    case default
       stop "ERROR get_bath_dimension_symmetris wiht bath_type!=replica/general"
    end select
    !
    ndx = Nsym
    !
    !for each replica we also store Nsym itself
    ndx = ndx+1
    !
    !number of replicas
    ndx = ndx * Nbath
    !diagonal hybridizations: Vs
    ! select case(bath_type)
    ! case("replica")
    !    ndx = ndx + Nbath
    ! case("general")
    ndx = ndx + Nbath*Nimp*Nspin
    ! end select
    !
    bath_size = ndx
    !
  end function get_bath_dimension_symmetries

  

  function get_bath_dimension_direct_d2(Hloc) result(bath_size)
    complex(8),intent(in)              :: Hloc(:,:) !optional input matrix with dimension [ |Nlso| , |Nlso| ] used to count the number of bath parameters in replica/general values of :f:var:`bath_type`.
    integer                                     :: bath_size
    complex(8),dimension(Nspin,Nspin,Nimp,Nimp) :: H
    integer                                     :: ndx
    !
    if(ed_mode=='superc')stop "get_bath_dimension_direct_d2 ERROR: called with ed_mode=superc"
    call assert_shape(Hloc,[Nspin*Nimp,Nspin*Nimp],"get_bath_dimension_direct_d2","Hloc")
    !
    H = so2nn_reshape(Hloc,Nspin,Nimp)
    !
    ndx=0
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Nimp
             do jorb=1,Nimp
                io = iorb + (ispin-1)*Nimp
                jo = jorb + (jspin-1)*Nimp
                if(io > jo)cycle
                if(dreal(H(ispin,jspin,iorb,jorb)) /= 0d0)ndx=ndx+1
                if(dimag(H(ispin,jspin,iorb,jorb)) /= 0d0)ndx=ndx+1
             enddo
          enddo
       enddo
    enddo
    ndx   = ndx + 1             !we also print Nbasis
    !
    !number of non vanishing elements for each replica
    ndx = ndx * Nbath
    !diagonal hybridizations: Vs
    ! select case(bath_type)
    ! case("replica")
    !    ndx = ndx + Nbath
    ! case("general")
    ndx = ndx + Nbath*Nimp*Nspin
    ! end select
    !
    bath_size = ndx
    !
  end function get_bath_dimension_direct_d2


  function get_bath_dimension_direct_d4(Hloc) result(bath_size)
    complex(8),intent(in)              :: Hloc(:,:,:,:) !optional input matrix with dimension [ |Nspin| , |Nspin| , |Nimp| , |Nimp|] used to count the number of bath parameters in replica/general values of :f:var:`bath_type`.
    integer                                     :: bath_size
    complex(8),dimension(Nspin,Nspin,Nimp,Nimp) :: H
    integer                                     :: ndx
    !
    if(ed_mode=='superc')stop "get_bath_dimension_direct_d2 ERROR: called with ed_mode=superc"
    call assert_shape(Hloc,[Nspin,Nspin,Nimp,Nimp],"get_bath_dimension_direct_d4","Hloc")
    !
    H = Hloc
    !
    ndx=0
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Nimp
             do jorb=1,Nimp
                io = iorb + (ispin-1)*Nimp
                jo = jorb + (jspin-1)*Nimp
                if(io > jo)cycle
                if(dreal(H(ispin,jspin,iorb,jorb)) /= 0d0)ndx=ndx+1
                if(dimag(H(ispin,jspin,iorb,jorb)) /= 0d0)ndx=ndx+1
             enddo
          enddo
       enddo
    enddo
    ndx   = ndx + 1             !we also print Nbasis
    !
    !number of non vanishing elements for each replica
    ndx = ndx * Nbath
    !diagonal hybridizations: Vs
    ! select case(bath_type)
    ! case("replica")
    !    ndx = ndx + Nbath
    ! case("general")
    ndx = ndx + Nbath*Nimp*Nspin
    ! end select
    !
    bath_size = ndx
    !
  end function get_bath_dimension_direct_d4


  function get_bath_dimension_direct_d6(Hloc) result(bath_size)
    complex(8),intent(in)                       :: Hloc(:,:,:,:,:,:) !optional input matrix with dimension [ |Nlat| , |Nlat| , |Nspin| , |Nspin| , |Norb| , |Norb|] used to count the number of bath parameters in replica/general values of :f:var:`bath_type`.
    integer                                     :: bath_size
    complex(8),dimension(Nspin,Nspin,Nimp,Nimp) :: H
    integer                                     :: ndx
    !
    if(ed_mode=='superc')stop "get_bath_dimension_direct_d2 ERROR: called with ed_mode=superc"
    call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Nimp,Nimp],"get_bath_dimension_direct_d6","Hloc")
    !
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       io = iorb+(ilat-1)*Norb
       jo = jorb+(jlat-1)*Norb
       H(ispin,jspin,io,jo) = Hloc(ilat,jlat,ispin,jspin,iorb,jorb)
    end do
    !
    ndx=0
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Nimp
             do jorb=1,Nimp
                io = iorb + (ispin-1)*Nimp
                jo = jorb + (jspin-1)*Nimp
                if(io > jo)cycle
                if(dreal(H(ispin,jspin,iorb,jorb)) /= 0d0)ndx=ndx+1
                if(dimag(H(ispin,jspin,iorb,jorb)) /= 0d0)ndx=ndx+1
             enddo
          enddo
       enddo
    enddo
    ndx   = ndx + 1             !we also print Nbasis
    !
    !number of non vanishing elements for each replica
    ndx = ndx * Nbath
    !diagonal hybridizations: Vs
    ! select case(bath_type)
    ! case("replica")
    !    ndx = ndx + Nbath
    ! case("general")
    ndx = ndx + Nbath*Nimp*Nspin
    ! end select
    !
    bath_size = ndx
    !
  end function get_bath_dimension_direct_d6














  function check_bath_dimension(bath_) result(bool)
    !
    ! Checks the  user bath :f:var:`bath_` on input has the correct dimensions according to the choice of input parameters for the calculations. 
    !
    real(8),dimension(:)           :: bath_
    integer                        :: Ntrue,i
    logical                        :: bool
    !
    if(.not.allocated(Hbath_basis))STOP "check_bath_dimension: Hbasis not allocated"
    !
    select case (bath_type)
    case default;stop "ERROR check_bath_dimension: bath_type!=replica/general"
    case ('replica','general')
       Ntrue   = get_bath_dimension_symmetries(size(Hbath_basis))
    end select
    bool  = ( size(bath_) == Ntrue )
  end function check_bath_dimension




END MODULE ED_BATH_DIM

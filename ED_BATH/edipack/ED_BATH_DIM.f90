MODULE ED_BATH_DIM
  !Returns or check the dimensions to which the user should allocate the bath array.
  !
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,str
  USE SF_LINALG, only: eye,inv
  USE SF_ARRAYS, only:linspace
  USE SF_MISC, only: assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  !
  USE ED_BATH_AUX
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
     module procedure ::  get_bath_dimension_direct
     module procedure ::  get_bath_dimension_symmetries
  end interface get_bath_dimension

  public :: get_bath_dimension
  public :: check_bath_dimension





contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : Inquire the correct bath size to allocate the
  ! the bath array in the calling program.
  !+-------------------------------------------------------------------+
  function get_bath_dimension_direct(H_nn) result(bath_size)
    complex(8),optional,intent(in) :: H_nn(:,:,:,:) !optional input matrix with dimension [ |Nspin| , |Nspin| , |Norb| , |Norb| ] used to count the number of bath parameters in replica/general values of :f:var:`bath_type`. 
    integer                        :: bath_size,ndx,ispin,iorb,jspin,jorb,io,jo,Maxspin
    complex(8),allocatable         :: H(:,:,:,:)
    !
    select case(bath_type)
       !
    case default
       select case(ed_mode)
       case default
          bath_size = Norb*Nbath + Norb*Nbath
          !( e [Nspin][Norb][Nbath] + v [Nspin][Norb][Nbath] )
       case ("superc")
          bath_size = Norb*Nbath + Norb*Nbath + Norb*Nbath
          !( e [Nspin][Norb][Nbath] + d [Nspin][Norb][Nbath] + v [Nspin][Norb][Nbath] )
       case ("nonsu2")
          bath_size = Norb*Nbath + Norb*Nbath + Norb*Nbath
          !( e [Nspin][Norb][Nbath] + v [Nspin][Norb][Nbath] + u [Nspin][Norb][Nbath] )
       end select
       bath_size=Nspin*bath_size
       !
    case('hybrid')
       select case(ed_mode)
       case default
          bath_size = Nbath + Norb*Nbath
          !(e [Nspin][1][Nbath] + v [Nspin][Norb][Nbath] )
       case ("superc")
          bath_size = Nbath + Nbath + Norb*Nbath
          !(e [Nspin][1][Nbath] + d [Nspin][1][Nbath] + v [Nspin][Norb][Nbath] )
       case ("nonsu2")
          bath_size = Nbath + Norb*Nbath + Norb*Nbath
          !(e [Nspin][1][Nbath] + v [Nspin][Norb][Nbath] + u [Nspin][Norb][Nbath] )
       end select
       bath_size=Nspin*bath_size
       !
    case('replica')
       allocate(H(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb))
       if(present(H_nn))then    !User defined H_nn
          H=H_nn
       elseif(Hreplica_status)then !User defined Hreplica_basis
          H=Hreplica_build()!we build with all lambda set to one as we only need to count non-zero elements!Hreplica_lambda(Nbath,:))
       else                        !Error:
          deallocate(H)
          stop "ERROR get_bath_dimension_direct: ed_mode=replica neither H_nn present nor Hreplica_basis defined"
       endif
       !
       !Check Hermiticity:
       ! if( all(abs(nn2so_reshape(H,Nnambu*Nspin,Norb) - conjg(transpose(nn2so_reshape(H,Nnambu*Nspin,Norb))))<1d-6)  )stop "H is not Hermitian"
      
       if( .not. check_herm( nn2so_reshape(H,Nnambu*Nspin,Norb),Nnambu*Nspin*Norb) )stop "H is not Hermitian"
       !
       !Re/Im off-diagonal non-vanishing elements
       ndx=0
       do ispin=1,Nnambu*Nspin
          do jspin=1,Nnambu*Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   if(io > jo)cycle
                   if(dreal(H(ispin,jspin,iorb,jorb)) /= 0d0)ndx=ndx+1
                   if(dimag(H(ispin,jspin,iorb,jorb)) /= 0d0)ndx=ndx+1
                enddo
             enddo
          enddo
       enddo
       !
       ndx = ndx * Nbath !number of non vanishing elements for each replica
       ndx = ndx + Nbath !diagonal hybridizations: Vs (different per spin)
       ndx = ndx + 1     !we also print Nbasis
       bath_size = ndx
    case('general')
       allocate(H(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb))
       if(present(H_nn))then    !User defined H_nn
          H=H_nn
       elseif(Hgeneral_status)then !User defined Hgeneral_basis
          H=Hgeneral_build()!Hgeneral_lambda(Nbath,:))!see above
       else                        !Error:
          deallocate(H)
          stop "ERROR get_bath_dimension_direct: ed_mode=general neither H_nn present nor Hgeneral_basis defined"
       endif
       !
       !Check Hermiticity:
       ! if( all(abs(nn2so_reshape(H,Nnambu*Nspin,Norb) - conjg(transpose(nn2so_reshape(H,Nnambu*Nspin,Norb))))<1d-6)  )stop "H is not Hermitian"
       if( .not. check_herm( nn2so_reshape(H,Nnambu*Nspin,Norb),Nnambu*Nspin*Norb) )stop "H is not Hermitian"
       !
       !Re/Im off-diagonal non-vanishing elements
       ndx=0
       do ispin=1,Nnambu*Nspin
          do jspin=1,Nnambu*Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   if(io > jo)cycle
                   if(dreal(H(ispin,jspin,iorb,jorb)) /= 0d0)ndx=ndx+1
                   if(dimag(H(ispin,jspin,iorb,jorb)) /= 0d0)ndx=ndx+1
                enddo
             enddo
          enddo
       enddo
       !
       ndx = ndx * Nbath            !number of non vanishing elements for each general
       ndx = ndx + Nbath*Norb*Nspin !diagonal hybridizations: Vs (different per spin and orbitals)
       ndx = ndx + 1                !we also print Nbasis
       bath_size = ndx
    end select
  end function get_bath_dimension_direct

  function get_bath_dimension_symmetries(Nsym) result(bath_size)
    integer :: Nsym !Number of symmetries (for :f:var:`ed_mode` = :code:`replica, general` )
    integer :: bath_size,ndx,isym
    !
    select case(bath_type)
    case("replica")
       if(.not.Hreplica_status)STOP "get_bath_dimension_symmetries: H(replica/general)_basis  not allocated"
       if(Nsym/=size(Hreplica_lambda,2))&
            stop "ERROR get_bath_dimension_symmetries:  size(Hreplica_basis) != size(Hreplica_lambda,2)"
    case("general")
       if(.not.Hgeneral_status)STOP "get_bath_dimension_symmetries: H(general/general)_basis  not allocated"
       if(Nsym/=size(Hgeneral_lambda,2))&
            stop "ERROR get_bath_dimension_symmetries:  size(Hgeneral_basis) != size(Hgeneral_lambda,2)"
    case default
       stop "ERROR get_bath_dimension_symmetris wiht bath_type!=replica/general"
    end select
    !
    ndx=Nsym
    !
    !number of replicas
    ndx = ndx * Nbath
    !diagonal hybridizations: Vs
    select case(bath_type)
    case("replica") ! Vk depends only on bath site
       ndx = ndx + Nbath
    case("general")
       ndx = ndx + Nbath*Norb*Nspin
    end select
    !
    !include Nbasis
    ndx=ndx+1
    !
    bath_size = ndx
    !
  end function get_bath_dimension_symmetries



  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if the dimension of the bath array are consistent
  !+-------------------------------------------------------------------+
  function check_bath_dimension(bath_) result(bool)
    !
    ! Checks the  user bath :f:var:`bath_` on input has the correct dimensions according to the choice of input parameters for the calculations. 
    !
    real(8),dimension(:)           :: bath_ !user-accessible bath array
    integer                        :: Ntrue,i
    logical                        :: bool
    select case (bath_type)
    case default
       Ntrue = get_bath_dimension()
    case ('replica')
       Ntrue   = get_bath_dimension_symmetries(size(Hreplica_basis))
    case ('general')
       Ntrue   = get_bath_dimension_symmetries(size(Hgeneral_basis))
    end select
    bool  = ( size(bath_) == Ntrue )
  end function check_bath_dimension

  

END MODULE ED_BATH_DIM














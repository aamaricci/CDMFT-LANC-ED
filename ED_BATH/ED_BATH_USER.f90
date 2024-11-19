MODULE ED_BATH_USER
  !Implements functions the user can use to enforce specific symmetry operations on the bath array.
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
  USE ED_BATH_DIM
  USE ED_BATH_REPLICA
  USE ED_BATH_DMFT
  implicit none

  private

  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  public :: impose_equal_lambda
  public :: impose_bath_offset




contains

  !##################################################################
  !
  !     USER BATH  SYMMETRIES: PREDEFINED AND USER CONTROLLED
  !
  !##################################################################

  subroutine impose_equal_lambda(bath_,ibath,lambdaindex_vec)
    ! Function to impose  :math:`\vec{\lambda}` parameters to be equal to a given average of a subset of values :f:var:`lambdaindex_vec` and for a specific bath element :f:var:`ibath` if :f:var:`bath_type` = :code:`replica` , :code:`general`. 
    !
    !
    real(8),dimension(:) :: bath_      !user bath array
    type(effective_bath) :: dmft_bath_
    real(8)              :: val
    integer,dimension(:) :: lambdaindex_vec!(sub)set of indices :math:`i` of :math:`\lambda_i` to be averaged out
    integer              :: ibath           !index of the bath element
    integer              :: i,N
    !
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    !
    N=size(lambdaindex_vec)
    val=0.d0
    do i=1,N
       val=val+dmft_bath%item(ibath)%lambda(lambdaindex_vec(i))/N
    enddo
    !
    do i=1,N
       dmft_bath%item(ibath)%lambda(lambdaindex_vec(i))=val
    enddo
    !
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine impose_equal_lambda


  subroutine impose_bath_offset(bath_,ibath,offset)
    real(8),dimension(:)    :: bath_
    type(effective_bath) :: dmft_bath_
    real(8)                 :: offset
    integer                 :: isym,N,ibath
    !
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    !
    if(size(Hbath_basis) .NE. dmft_bath_%Nbasis)then
       dmft_bath_%item(ibath)%lambda(dmft_bath_%Nbasis)=offset
    else
       do isym=1,size(Hbath_lambda)
          if(is_identity(Hbath_basis(isym)%O)) dmft_bath%item(ibath)%lambda(isym)=offset
          return
       enddo
    endif
    !
    !
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
    !
  end subroutine impose_bath_offset






END MODULE ED_BATH_USER

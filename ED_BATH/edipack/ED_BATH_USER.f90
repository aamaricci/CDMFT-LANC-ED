MODULE ED_BATH_USER
  !Implements functions the user can use to enforce specific symmetry operations on the bath array.
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
  USE ED_BATH_DIM
  USE ED_BATH_DMFT
  implicit none

  private

  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  !explicit symmetries:
  interface break_symmetry_bath
     !
     ! Function to impose a specific symmetry breaking pattern into the energy levels of the bath. A common case is to find magnetic solution by breaking spin degeneracy of the levels.
     !
     ! * :f:var:`bath` is rank-1 [:] or rank-2 [:, :f:var:`Nlat` ] array type: double precision with assumed size
     !
     !
     module procedure break_symmetry_bath_site
     module procedure break_symmetry_bath_lattice
  end interface break_symmetry_bath

  interface spin_symmetrize_bath
     !
     ! Function to impose a spin symmetry to the parameters of the bath. Enforces a non-magnetic solution
     !
     ! * :f:var:`bath` is rank-1 [:] or rank-2 [:, :f:var:`Nlat` ] array type: double precision with assumed size
     !
     module procedure spin_symmetrize_bath_site
     module procedure spin_symmetrize_bath_lattice
  end interface spin_symmetrize_bath

  interface orb_symmetrize_bath
     !
     ! Function to impose a orbital symmetry to the parameters of the bath. Enforces an orbital non-polarized solution. If two orbital indices :f:var:`orb1` and :f:var:`orb2` are passed symmetry is imposed only among such two orbitals
     !
     !.. warning::
     !
     !   This operation requires the orbital to be degenerate.
     !
     ! * :f:var:`bath` is rank-1 [:] or rank-2 [:, :f:var:`Nlat` ] array type: double precision with assumed size
     !
     module procedure orb_symmetrize_bath_site
     module procedure orb_symmetrize_bath_lattice
     module procedure orb_symmetrize_bath_site_o1o2
     module procedure orb_symmetrize_bath_lattice_o1o2
  end interface orb_symmetrize_bath

  interface orb_equality_bath
     !
     ! Function to impose a orbital equality on the parameters of the bath. 
     !
     ! * :f:var:`bath` is rank-1 [:] or rank-2 [:, :f:var:`Nlat` ] array type: double precision with assumed size
     !
     module procedure orb_equality_bath_site
     module procedure orb_equality_bath_lattice
  end interface orb_equality_bath

  interface ph_symmetrize_bath
     !
     ! Function to impose particle-hole symmetry to the parameters of the bath. 
     !
     ! * :f:var:`bath` is rank-1 [:] or rank-2 [:, :f:var:`Nlat` ] array type: double precision with assumed size
     !
     module procedure ph_symmetrize_bath_site
     module procedure ph_symmetrize_bath_lattice
  end interface ph_symmetrize_bath

  interface ph_trans_bath
     !
     ! Function to perform particle-hole transformation to the parameters of the bath. 
     !
     ! * :f:var:`bath` is rank-1 [:] or rank-2 [:, :f:var:`Nlat` ] array type: double precision with assumed size
     !
     module procedure ph_trans_bath_site
     module procedure ph_trans_bath_lattice
  end interface ph_trans_bath

  interface enforce_normal_bath
     !
     ! Function to impose normal solution to the parameters of the bath, i.e. suppressed superconductivity if any.
     !
     ! * :f:var:`bath` is rank-1 [:] or rank-2 [:, :f:var:`Nlat` ] array type: double precision with assumed size
     !
     module procedure enforce_normal_bath_site
     module procedure enforce_normal_bath_lattice
  end interface enforce_normal_bath

  interface save_array_as_bath
     !
     ! Write the bath parameters to a file following the convention of the internal data structure :f:var:`effective_bath`. 
     !
     ! * :f:var:`bath` is rank-1 [:] or rank-2 [:, :f:var:`Nlat` ] array type: double precision with assumed size
     !
     module procedure save_array_as_bath_site
     module procedure save_array_as_bath_lattice
  end interface save_array_as_bath




  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  public :: break_symmetry_bath
  public :: spin_symmetrize_bath
  public :: orb_symmetrize_bath
  public :: orb_equality_bath
  public :: ph_symmetrize_bath
  public :: ph_trans_bath
  public :: enforce_normal_bath
  public :: impose_equal_lambda
  public :: impose_bath_offset
  public :: save_array_as_bath





contains


  !##################################################################
  !
  !     USER BATH  SYMMETRIES: PREDEFINED AND USER CONTROLLED
  !
  !##################################################################

  !+-------------------------------------------------------------------+
  !PURPOSE  : given a bath array apply a specific transformation or
  ! impose a given symmetry:
  ! - break spin symmetry by applying a symmetry breaking field
  ! - given a bath array set both spin components to have
  !    the same bath, i.e. impose non-magnetic solution
  ! - given a bath array enforces the particle-hole symmetry
  !    by setting the positive energies in modulo identical to the negative
  !    ones.
  ! - given a bath enforce normal (i.e. non superconducting) solution
  ! - given a dmft bath pull/push the components W^{ss'}_\a(l) of the Hybridization
  !    matrix
  ! - given a dmft bath pull/push the nonsu2 components
  !+-------------------------------------------------------------------+
  subroutine impose_equal_lambda(bath_,ibath,lambdaindex_vec)
    !
    ! Function to impose  :math:`\vec{\lambda}` parameters to be equal to a given average of a subset of values :f:var:`lambdaindex_vec` and for a specific bath element :f:var:`ibath` if :f:var:`bath_type` = :code:`replica` , :code:`general`. 
    !
    !
    real(8),dimension(:) :: bath_      !user bath array
    type(effective_bath) :: dmft_bath_
    real(8)              :: val
    integer,dimension(:) :: lambdaindex_vec !(sub)set of indices :math:`i` of :math:`\lambda_i` to be averaged out
    integer              :: ibath           !index of the bath element
    integer              :: i,N
    !
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    !
    N=size(lambdaindex_vec)
    val=0.d0
    do i=1,N
       val=val+dmft_bath_%item(ibath)%lambda(lambdaindex_vec(i))/N
    enddo
    !
    do i=1,N
       dmft_bath_%item(ibath)%lambda(lambdaindex_vec(i))=val
    enddo
    !
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine impose_equal_lambda


  subroutine impose_bath_offset(bath_,ibath,offset)
    real(8),dimension(:) :: bath_
    type(effective_bath) :: dmft_bath_
    real(8)              :: offset
    integer              :: isym,N,ibath
    !
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    !
    select case(bath_type)
    case default
       if(size(Hreplica_basis) .ne. dmft_bath_%Nbasis)then
          dmft_bath_%item(ibath)%lambda(dmft_bath_%Nbasis)=offset
       else
          do isym=1,size(Hreplica_basis)
             if(is_identity(Hreplica_basis(isym)%O)) dmft_bath_%item(ibath)%lambda(isym)=offset
             return
          enddo
       endif
    case("general")
       if(size(Hgeneral_basis) .ne. dmft_bath_%Nbasis)then
          dmft_bath_%item(ibath)%lambda(dmft_bath_%Nbasis)=offset
       else
          do isym=1,size(Hgeneral_basis)
             if(is_identity(Hgeneral_basis(isym)%O)) dmft_bath_%item(ibath)%lambda(isym)=offset
             return
          enddo
       endif
    end select

    !
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
    !
  end subroutine impose_bath_offset


  subroutine break_symmetry_bath_site(bath_,field,sign,save)
    real(8),dimension(:)   :: bath_      !user bath array
    type(effective_bath)   :: dmft_bath_ 
    real(8)                :: field      !tiny symmetry breaking field 
    real(8)                :: sign       !sign of the field
    logical,optional       :: save       !optional flag to save the output bath
    logical                :: save_
    if(bath_type=="replica")stop "break_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "break_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    dmft_bath_%e(1,:,:)    =dmft_bath_%e(1,:,:)      + sign*field
    dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(Nspin,:,:)  - sign*field
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine break_symmetry_bath_site
  !
  subroutine break_symmetry_bath_lattice(bath_,field,sign,save)
    real(8),dimension(:,:) :: bath_
    real(8)                :: field
    real(8)                :: sign
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call break_symmetry_bath_site(bath_(ilat,:),field,sign,save_)
    enddo
    ed_file_suffix=""
  end subroutine break_symmetry_bath_lattice


  !---------------------------------------------------------!


  subroutine spin_symmetrize_bath_site(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    logical,optional       :: save
    logical                :: save_
    integer :: ibath
    if(bath_type=="replica")stop "spin_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "spin_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    if(Nspin==1)then
       write(LOGfile,"(A)")"spin_symmetrize_bath: Nspin=1 nothing to symmetrize"
       return
    endif
    !
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    select case(ed_mode)
    case default
       dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(1,:,:)
       dmft_bath_%v(Nspin,:,:)=dmft_bath_%v(1,:,:)
    case ("superc")
       dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(1,:,:)
       dmft_bath_%v(Nspin,:,:)=dmft_bath_%v(1,:,:)
       dmft_bath_%d(Nspin,:,:)=dmft_bath_%d(1,:,:)
    end select
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine spin_symmetrize_bath_site
  !
  subroutine spin_symmetrize_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call spin_symmetrize_bath_site(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine spin_symmetrize_bath_lattice


  !---------------------------------------------------------!


  subroutine orb_symmetrize_bath_site(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: iorb
    real(8),allocatable    :: lvl(:,:),hyb(:,:)
    if(bath_type=="replica")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    if(Norb==1)then
       write(LOGfile,"(A)")"orb_symmetrize_bath: Norb=1 nothing to symmetrize"
       return
    endif
    !
    call allocate_dmft_bath(dmft_bath_)
    ! if (bath_type=="replica")call init_dmft_bath_mask(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    !
    if(allocated(lvl))deallocate(lvl);allocate(lvl(Nspin,Nbath));lvl=0d0;lvl=sum(dmft_bath_%e,dim=2)/Norb
    if(allocated(hyb))deallocate(hyb);allocate(hyb(Nspin,Nbath));hyb=0d0;hyb=sum(dmft_bath_%v,dim=2)/Norb
    do iorb=1,Norb
       dmft_bath_%e(:,iorb,:)=lvl
       dmft_bath_%v(:,iorb,:)=hyb
    enddo
    !
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine orb_symmetrize_bath_site
  subroutine orb_symmetrize_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    if(bath_type=="replica")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call orb_symmetrize_bath_site(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine orb_symmetrize_bath_lattice


  subroutine orb_symmetrize_bath_site_o1o2(bath_,orb1,orb2,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: iorb,orb1,orb2
    real(8),allocatable    :: lvl(:,:),hyb(:,:)
    if(bath_type=="replica")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    if(Norb==1)then
       write(LOGfile,"(A)")"orb_symmetrize_bath: Norb=1 nothing to symmetrize"
       return
    endif
    !
    call allocate_dmft_bath(dmft_bath_)
    ! if (bath_type=="replica")call init_dmft_bath_mask(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    !
    if(allocated(lvl))deallocate(lvl);allocate(lvl(Nspin,Nbath));lvl=0d0
    if(allocated(hyb))deallocate(hyb);allocate(hyb(Nspin,Nbath));hyb=0d0
    !
    lvl=(dmft_bath_%e(:,orb1,:)+dmft_bath_%e(:,orb2,:))/2d0
    hyb=(dmft_bath_%v(:,orb1,:)+dmft_bath_%v(:,orb2,:))/2d0
    !
    dmft_bath_%e(:,orb1,:)=lvl
    dmft_bath_%v(:,orb1,:)=hyb
    dmft_bath_%e(:,orb2,:)=lvl
    dmft_bath_%v(:,orb2,:)=hyb
    !
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine orb_symmetrize_bath_site_o1o2
  subroutine orb_symmetrize_bath_lattice_o1o2(bath_,orb1,orb2,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat,orb1,orb2
    if(bath_type=="replica")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call orb_symmetrize_bath_site_o1o2(bath_(ilat,:),orb1,orb2,save_)
    enddo
    ed_file_suffix=""
  end subroutine orb_symmetrize_bath_lattice_o1o2

  !---------------------------------------------------------!


  subroutine orb_equality_bath_site(bath_,indx,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    integer,optional       :: indx
    logical,optional       :: save
    integer                :: indx_
    logical                :: save_
    integer                :: iorb
    real(8),allocatable    :: lvl(:,:),hyb(:,:)
    if(bath_type=="replica")stop "orb_equality_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "orb_equality_bath_site ERROR: can not be used with bath_type=general"
    indx_=1     ;if(present(indx))indx_=indx
    save_=.true.;if(present(save))save_=save
    if(Norb==1)then
       write(LOGfile,"(A)")"orb_symmetrize_bath: Norb=1 nothing to symmetrize"
       return
    endif
    !
    call allocate_dmft_bath(dmft_bath_)
    ! if (bath_type=="replica")call init_dmft_bath_mask(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    !
    if(allocated(lvl))deallocate(lvl);allocate(lvl(Nspin,Nbath));lvl=0d0;lvl=dmft_bath_%e(:,indx_,:)
    if(allocated(hyb))deallocate(hyb);allocate(hyb(Nspin,Nbath));hyb=0d0;hyb=dmft_bath_%v(:,indx_,:)
    do iorb=1,Norb
       if(iorb==indx_)cycle
       dmft_bath_%e(:,iorb,:)=lvl
       dmft_bath_%v(:,iorb,:)=hyb
    enddo
    !
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine orb_equality_bath_site
  subroutine orb_equality_bath_lattice(bath_,indx,save)
    real(8),dimension(:,:) :: bath_
    integer,optional       :: indx
    logical,optional       :: save
    integer                :: indx_
    logical                :: save_
    integer                :: iorb
    integer                :: Nsites,ilat
    if(bath_type=="replica")stop "orb_equality_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "orb_equality_bath_site ERROR: can not be used with bath_type=general"
    indx_=1     ;if(present(indx))indx_=indx
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call orb_equality_bath_site(bath_(ilat,:),indx_,save_)
    enddo
    ed_file_suffix=""
  end subroutine orb_equality_bath_lattice



  !---------------------------------------------------------!

  subroutine ph_symmetrize_bath_site(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    integer                :: i
    logical,optional       :: save
    logical                :: save_
    if(bath_type=="replica")stop "ph_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "ph_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    if(Nbath==1)return
    if(mod(Nbath,2)==0)then
       do i=1,Nbath/2
          dmft_bath_%e(:,:,Nbath+1-i)=-dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,Nbath+1-i)= dmft_bath_%v(:,:,i)
          if(ed_mode=="superc")dmft_bath_%d(:,:,Nbath+1-i)=dmft_bath_%d(:,:,i)
       enddo
    else
       do i=1,(Nbath-1)/2
          dmft_bath_%e(:,:,Nbath+1-i)=-dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,Nbath+1-i)= dmft_bath_%v(:,:,i)
          if(ed_mode=="superc")dmft_bath_%d(:,:,Nbath+1-i)=dmft_bath_%d(:,:,i)
       enddo
       dmft_bath_%e(:,:,(Nbath-1)/2+1)=0.d0
    endif
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine ph_symmetrize_bath_site
  subroutine ph_symmetrize_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    if(bath_type=="replica")stop "ph_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "ph_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call ph_symmetrize_bath_site(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine ph_symmetrize_bath_lattice

  !---------------------------------------------------------!

  subroutine ph_trans_bath_site(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    type(effective_bath)   :: tmp_dmft_bath
    integer                :: i
    logical,optional       :: save
    logical                :: save_
    if(bath_type=="replica")stop "ph_trans_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "ph_trans_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    call allocate_dmft_bath(dmft_bath_)
    call allocate_dmft_bath(tmp_dmft_bath)
    call set_dmft_bath(bath_,dmft_bath_)
    if(Nbath==1)return
    do i=1,Nbath
       select case(Norb)
       case default
          ! do nothing
          dmft_bath_%e(:,:,i)= dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,i)= dmft_bath_%v(:,:,i)
       case(1)
          dmft_bath_%e(:,:,i)= -dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,i)=  dmft_bath_%v(:,:,i)
       case(2)
          tmp_dmft_bath%e(:,1,i) = -dmft_bath_%e(:,2,i)
          tmp_dmft_bath%e(:,2,i) = -dmft_bath_%e(:,1,i)
          dmft_bath_%e(:,:,i)    = tmp_dmft_bath%e(:,:,i)
          tmp_dmft_bath%v(:,1,i) = dmft_bath_%v(:,2,i)
          tmp_dmft_bath%v(:,2,i) = dmft_bath_%v(:,1,i)
          dmft_bath_%v(:,:,i)    = tmp_dmft_bath%v(:,:,i)
       end select
    end do
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine ph_trans_bath_site
  subroutine ph_trans_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    if(bath_type=="replica")stop "ph_trans_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "ph_trans_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call ph_trans_bath_site(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine ph_trans_bath_lattice

  !---------------------------------------------------------!

  subroutine enforce_normal_bath_site(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    logical,optional       :: save
    logical                :: save_
    if(bath_type=="replica")stop "enforce_normal_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "enforce_normal_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    if(ed_mode=="superc")dmft_bath_%d(:,:,:)=0.d0
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine enforce_normal_bath_site

  subroutine enforce_normal_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call enforce_normal_bath_site(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine enforce_normal_bath_lattice










  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  check if the specified itype is consistent with the input parameters.
  !+-----------------------------------------------------------------------------+!
  subroutine check_bath_component(type)
    character(len=1) :: type
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          if(type/="e".AND.type/='v')stop "check_bath_component error: type!=e,v"
       case ("superc")
          if(type/="e".AND.type/='v'.AND.type/='d')stop "check_bath_component error: type!=e,v,d"
       case ("nonsu2")
          if(type/="e".AND.type/='v'.AND.type/='u')stop "check_bath_component error: type!=e,v,u"
       end select
    case ("replica","general")
       if(type/="v".AND.type/="l")stop "check_bath_component error: type!=v,l"
    end select
    return
  end subroutine check_bath_component

  !+------------------------------------------------------------------+
  !PURPOSE  : given a array, save it as a bath. Do nothing else.
  !+------------------------------------------------------------------+

  subroutine save_array_as_bath_site(bath_)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    call save_dmft_bath(dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine save_array_as_bath_site

  subroutine save_array_as_bath_lattice(bath_)
    real(8),dimension(:,:) :: bath_
    integer                :: Nsites,ilat
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call save_array_as_bath_site(bath_(ilat,:))
    enddo
    ed_file_suffix=""
  end subroutine save_array_as_bath_lattice




END MODULE ED_BATH_USER














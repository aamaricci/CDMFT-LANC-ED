MODULE ED_BATH_DMFT
  !
  !A class for the the :f:var:`effective_bath` data structure describing the effective bath in the code.
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
  implicit none

  private


  public :: allocate_dmft_bath
  public :: deallocate_dmft_bath
  public :: init_dmft_bath
  public :: write_dmft_bath
  public :: save_dmft_bath
  public :: set_dmft_bath
  public :: get_dmft_bath



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
  !PURPOSE  : Allocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine allocate_dmft_bath(dmft_bath_)
    !
    ! Allocate the  :f:var:`effective_bath` input  :f:var:`dmft_bath_` according to the values of the global input parameters |Nspin| , |Norb|, |Nbath|, |ed_mode| and |bath_type|.
    !
    ! The correct components of the :f:var:`effective_bath` data structure are allocated.   
    !
    ! .. list-table:: allocated components of :f:var:`dmft_bath_`
    !    :widths: auto
    !    :header-rows: 1
    !    :stub-columns: 1
    !
    !    * - 
    !      - :code:`normal`
    !      - :code:`superc`
    !    !
    !    * - :code:`replica`
    !      - :f:var:`nbasis` , :f:var:`item`  , :f:var:`item` % :f:var:`lambda`
    !      - :f:var:`nbasis` , :f:var:`item`  , :f:var:`item` % :f:var:`lambda`
    !
    !    * - :code:`general`
    !      - :f:var:`nbasis` , :f:var:`item`  , :f:var:`item` % :f:var:`lambda` , :f:var:`item` % :f:var:`vg`
    !      - :f:var:`nbasis` , :f:var:`item`  , :f:var:`item` % :f:var:`lambda` , :f:var:`item` % :f:var:`vg`
    !
    !
    type(effective_bath) :: dmft_bath_
    integer              :: isym,Nsym
    !
    if(.not.Hbath_status)stop "ERROR allocate_dmft_bath: Hbath_basis not allocated"
    call deallocate_dmft_bath(dmft_bath_)
    Nsym=size(Hbath_basis)
    !
    allocate(dmft_bath_%item(Nbath))
    dmft_bath_%Nbasis=Nsym
    do ibath=1,Nbath
       allocate(dmft_bath_%item(ibath)%lambda(Nsym))
       allocate(dmft_bath_%item(ibath)%v(Nspin*Nimp)) !allocate all, but V(2:N)=V(1) for Replica, V(:) for General
    enddo
    !
    dmft_bath_%status=.true.
    !
  end subroutine allocate_dmft_bath



  !+-------------------------------------------------------------------+
  !PURPOSE  : Deallocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine deallocate_dmft_bath(dmft_bath_)
    !
    ! Deallocate the  :f:var:`effective_bath` input  :f:var:`dmft_bath_` 
    !
    type(effective_bath) :: dmft_bath_
    integer              :: ibath,isym
    if(.not.dmft_bath_%status)return
    dmft_bath_%Nbasis = 0
    do ibath=1,Nbath
       deallocate(dmft_bath_%item(ibath)%lambda)
       deallocate(dmft_bath_%item(ibath)%v)
    enddo
    deallocate(dmft_bath_%item)
    dmft_bath_%status=.false.
  end subroutine deallocate_dmft_bath




  !+------------------------------------------------------------------+
  !PURPOSE  : Initialize the DMFT loop, building H parameters and/or
  !           reading previous (converged) solution
  !+------------------------------------------------------------------+
  subroutine init_dmft_bath(dmft_bath_,used)
    !
    ! Initialize the :f:var:`effective_bath` input :f:var:`dmft_bath_`.
    ! The replica/general local bath Hamiltonians are built out of the initial values of :math:`\vec{\lambda}`, using the matrix decomposition.
    !
    !.. note::
    !
    !   If the file :f:var:`hfile` with the correct :f:var:`ed_file_suffix` and :f:var:`hsuffix` is found in the running directory, the bath is initilized by reading the from the file.
    !
    !
    type(effective_bath)     :: dmft_bath_
    logical,optional         :: used
    real(8)                  :: re,im
    integer                  :: isym,unit,Nh,Nsym,flen
    integer,dimension(Nbath) :: Nlambdas
    logical                  :: IOfile,used_,diagonal_hsym,all_are_equal
    real(8)                  :: de
    real(8)                  :: rescale(Nbath),offset(Nbath)
    real(8),dimension(Nbath) :: one_lambdaval
    character(len=21)        :: space
    character(len=20)        :: hsuffix
    !
    used_   = .false.   ;if(present(used))used_=used
    hsuffix = ".restart";if(used_)hsuffix=reg(".used")
    if(.not.dmft_bath_%status)stop "ERROR init_dmft_bath error: bath not allocated"
    !
    rescale = 0d0
    if(Nbath>1)rescale = linspace(ed_hw_band/Nbath,ed_hw_band,Nbath)
    offset  = 0d0
    if(Nbath>1) offset = linspace(-ed_offset_bath,ed_offset_bath,Nbath)    
    !
    !
    !BATH V INITIALIZATION
    do ibath=1,Nbath
       dmft_bath_%item(ibath)%v = max(0.1d0,1.d0/sqrt(dble(Nbath)))
    enddo
    !
    !BATH LAMBDAS INITIALIZATION
    Nsym = dmft_bath_%Nbasis
    do isym=1,Nsym
       do ibath=1,Nbath
          dmft_bath_%item(ibath)%lambda(isym) = Hbath_lambda(ibath,isym)
       enddo
       !
       diagonal_hsym = is_diagonal(Hbath_basis(isym)%O)
       one_lambdaval = Hbath_lambda(1,isym)
       all_are_equal = all(Hbath_lambda(:,isym)==one_lambdaval)
       !
       if(diagonal_hsym.AND.all_are_equal.AND.Nbath>1)then
          !CDMFT:
          !> BACK-COMPATIBILITY PATCH: Nbath /degenerate/ lambdas, rescaled internally
          do ibath=1,Nbath
             dmft_bath_%item(ibath)%lambda(isym) = rescale(ibath) * Hbath_lambda(ibath,isym)
          enddo
          !
          write(*,*) "                                                                    "
          write(*,*) "WARNING: some of your lambdasym values have been internally changed "
          write(*,*) "         while calling ed_init_solver. This happens whenever the    "
          write(*,*) "         corresponding Hsym is diagonal and all the generals receive"
          write(*,*) "         the same initial lambda value, due to the deprecated legacy"
          write(*,*) "         practice of defining a unique lambda vector forall generals"
          write(*,*) "         and let the solver decide how to handle these degeneracies."
          write(*,*) "         >>> If you really intend to have a degenerate diagonal term"
          write(*,*) "             in the bath you can define a suitable restart file.    "
          write(*,*) "         >>> If instead this is what you expected please consider to"
          write(*,*) "             move the desired rescaling in your driver, since this  "
          write(*,*) "             funcionality might be removed in a future update.      "
          write(*,*) "                                                                    "
          !
       endif
    enddo
    !
    !Read from file if exist:
    inquire(file=trim(Hfile)//trim(ed_file_suffix)//trim(hsuffix),exist=IOfile)
    if(IOfile)then
       write(LOGfile,"(A)")'Reading bath from file'//trim(Hfile)//trim(ed_file_suffix)//trim(hsuffix)
       unit = free_unit()
       flen = file_length(trim(Hfile)//trim(ed_file_suffix)//trim(hsuffix))
       !
       open(unit,file=trim(Hfile)//trim(ed_file_suffix)//trim(hsuffix))
       read(unit,*)
       read(unit,*)dmft_bath_%Nbasis       
       do ibath=1,Nbath
          select case(bath_type)
          case('replica')
             read(unit,*)dmft_bath_%item(ibath)%v(1),(dmft_bath_%item(ibath)%lambda(io),io=1,dmft_bath_%Nbasis)
             dmft_bath_%item(ibath)%v(2:Nspin*Nimp) = dmft_bath_%item(ibath)%v(1)
          case('general')
             read(unit,*)dmft_bath_%item(ibath)%v(:),(dmft_bath_%item(ibath)%lambda(io),io=1,dmft_bath_%Nbasis)
          end select
       enddo
       close(unit)
       !
    endif
  end subroutine init_dmft_bath








  !+-------------------------------------------------------------------+
  !PURPOSE  : write out the bath to a given unit [temporary format]
  !+-------------------------------------------------------------------+
  subroutine write_dmft_bath(dmft_bath_,unit)
    !
    ! Write the :f:var:`effective_bath` on input to std.output or to a file associated to the [optional] unit.
    !
    type(effective_bath) :: dmft_bath_
    integer,optional     :: unit
    integer              :: unit_,isym
    complex(8)           :: Ho(Nambu*Nspin*Nimp,Nambu*Nspin*Nimp)
    character(len=64)    :: string_fmt,string_fmt_first
    !
    unit_ = LOGfile;if(present(unit))unit_=unit
    if(.not.dmft_bath_%status)stop "write_dmft_bath error: bath not allocated"
    !
    string_fmt  ="("//str(Nspin*Nimp)//"(A1,F5.2,A1,F5.2,A1,2x))"
    !
    select case(bath_type)
    case('replica')
       write(unit_,"(A1,90(A21,1X))")"#",&
            ("V_i"//reg(str(io)),io=1,1),&
            ("Lambda_i"//reg(str(io)),io=1,dmft_bath_%Nbasis)
    case('general')
       write(unit_,"(A1,90(A21,1X))")"#",&
            ("V_i"//reg(str(io)),io=1,Nspin*Nimp),&
            ("Lambda_i"//reg(str(io)),io=1,dmft_bath_%Nbasis)
    end select
    !
    write(unit_,"(I3)")dmft_bath_%Nbasis
    !
    do i=1,Nbath
       select case(bath_type)
       case('replica')
          write(unit_,"(90(ES21.12,1X))")&
               (dmft_bath_%item(i)%v(io),io=1,1),&
               (dmft_bath_%item(i)%lambda(io),io=1,dmft_bath_%Nbasis)
       case('general')
          write(unit_,"(90(ES21.12,1X))")&
               (dmft_bath_%item(i)%v(io),io=1,Nspin*Nimp),&
               (dmft_bath_%item(i)%lambda(io),io=1,dmft_bath_%Nbasis)
       end select
    enddo
    !
    if(unit_/=LOGfile .AND. Nambu*Nspin*Nimp<10)then
       write(unit_,*)""
       do isym=1,size(Hbath_basis)
          !
          do concurrent(in=1:Nambu,jn=1:Nambu,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Nimp,jorb=1:Nimp)
             i = iorb + (ispin-1)*Nimp + (in-1)*Nspin*Nimp
             j = jorb + (jspin-1)*Nimp + (jn-1)*Nspin*Nimp
             Ho(i,j) = Hbath_basis(isym)%O(in,jn,ispin,jspin,iorb,jorb)
          enddo
          !
          do i=1,Nambu*Nspin*Nimp
             write(unit_,string_fmt)&
                  ('(',dreal(Ho(i,j)),',',dimag(Ho(i,j)),')',j =1,Nambu*Nspin*Nimp)
          enddo
          write(unit_,*)""
       enddo
    endif
  end subroutine write_dmft_bath



  subroutine save_dmft_bath(dmft_bath_,file,used)
    !
    ! Save the :f:var:`effective_bath` to a file with .used or .restart extension according to input.
    !
    type(effective_bath)      :: dmft_bath_
    character(len=*),optional :: file
    character(len=256)        :: file_
    logical,optional          :: used
    logical                   :: used_
    character(len=16)         :: extension
    integer                   :: unit_
    if(.not.dmft_bath_%status)stop "save_dmft_bath error: bath is not allocated"
    used_=.false.;if(present(used))used_=used
    extension=".restart";if(used_)extension=".used"
    file_=str(str(Hfile)//str(ed_file_suffix)//str(extension))
    if(present(file))file_=str(file)
    unit_=free_unit()
    open(unit_,file=str(file_))
    call write_dmft_bath(dmft_bath_,unit_)
    close(unit_)
  end subroutine save_dmft_bath





  subroutine set_dmft_bath(bath_,dmft_bath_)
    !
    ! Set the :f:var:`effective_bath` components from the input user bath :f:var:`bath_` , i.e. it dumps the user bath to the internal data structure. 
    !
    real(8),dimension(:) :: bath_
    type(effective_bath) :: dmft_bath_
    integer              :: stride,ibath,Nmask,isym
    logical              :: check
    !
    if(.not.dmft_bath_%status)stop "get_dmft_bath error: bath not allocated"
    !
    check=check_bath_dimension(bath_)
    if(.not.check)stop "get_dmft_bath error: wrong bath dimensions"
    !
    do ibath=1,Nbath
       dmft_bath_%item(ibath)%v=0d0
       dmft_bath_%item(ibath)%lambda=0d0
    enddo
    !
    stride = 1
    !Get Nbasis
    dmft_bath_%Nbasis = NINT(bath_(stride))
    do ibath=1,Nbath
       !get Vs
       select case(bath_type)
       case('replica')
          dmft_bath_%item(ibath)%v(:) = bath_(stride+1)
          stride = stride + 1
       case('general')
          dmft_bath_%item(ibath)%v(:) = bath_(stride+1:stride+Nspin*Nimp)
          stride = stride + Nspin*Nimp
       end select
       !get Lambdas
       dmft_bath_%item(ibath)%lambda=bath_(stride+1:stride+dmft_bath_%Nbasis)
       stride=stride+dmft_bath_%Nbasis
    enddo
    !
  end subroutine set_dmft_bath





  !+-------------------------------------------------------------------+
  !PURPOSE  : set the bath components from a given user provided
  !           bath-array
  !+-------------------------------------------------------------------+
  subroutine get_dmft_bath(dmft_bath_,bath_)
    !
    ! Set the user input bath :f:var:`bath_` from the components of the :f:var:`effective_bath` :f:var:`dmft_bath_` , i.e. it dumps the internal data structure to the user bath. 
    !
    type(effective_bath)   :: dmft_bath_
    real(8),dimension(:)   :: bath_
    integer                :: stride,Nmask,isym
    logical                :: check
    !
    if(.not.dmft_bath_%status)stop "set_dmft_bath error: bath not allocated"
    !
    check = check_bath_dimension(bath_)
    if(.not.check)stop "set_dmft_bath error: wrong bath dimensions"
    !
    bath_=0d0
    !
    stride = 1
    bath_(stride)=dmft_bath_%Nbasis
    !
    do ibath=1,Nbath
       select case(bath_type)
       case('replica')
          bath_(stride+1)=dmft_bath_%item(ibath)%v(1)
          stride = stride + 1
       case('general')
          bath_(stride+1:stride+Nspin*Nimp)=dmft_bath_%item(ibath)%v(:)
          stride = stride + Nspin*Nimp
       end select
       bath_(stride+1:stride+dmft_bath_%Nbasis)=dmft_bath_%item(ibath)%lambda
       stride=stride+dmft_bath_%Nbasis
    enddo
  end subroutine get_dmft_bath




END MODULE ED_BATH_DMFT

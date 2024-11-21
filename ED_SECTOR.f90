MODULE ED_SECTOR
  !
  !Contains procedures to construct the symmetry sectors corresponding to a given set of quantum numbers :math:`\vec{Q}`, in particular it allocated and build the  :f:var:`sector_map` connecting the states of a given sector with the corresponding Fock ones. 
  !
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_TIMER
  USE SF_IOTOOLS, only: free_unit,reg,file_length
  USE SF_MISC,    only: assert_shape
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private


  interface map_allocate
     !
     ! Allocate the map(s) connecting the states in the symmetry sector to thos in the Fock space. 
     !
     module procedure :: map_allocate_scalar
     module procedure :: map_allocate_vector
  end interface map_allocate

  interface map_deallocate
     module procedure :: map_deallocate_scalar
     module procedure :: map_deallocate_vector
  end interface map_deallocate

  interface flip_state
     module procedure :: flip_state_normal
     module procedure :: flip_state_other
  end interface flip_state

  interface get_Sector
     !
     ! Returns the index of :f:var:`isector` of the symmetry sector given the quantum numbers :f:var:`qn`
     !
     module procedure :: get_Sector_normal
     module procedure :: get_Sector_superc
  end interface get_Sector

  interface get_QuantumNumbers
     !
     !Returns the quantum numbers :f:var:`qn` given the  index of :f:var:`isector` of the symmetry sector
     !
     module procedure :: get_QuantumNumbers_normal
     module procedure :: get_QuantumNumbers_other
  end interface get_QuantumNumbers



  public :: build_sector
  public :: delete_sector
  !
  public :: apply_Cops
  public :: apply_op_C
  public :: apply_op_CDG
  !
  public :: get_Sector
  public :: get_Indices
  public :: get_Nup
  public :: get_Ndw
  public :: get_Sz
  public :: get_DimUp
  public :: get_DimDw
  public :: get_Dim
  !
  public :: indices2state
  public :: state2indices
  public :: iup_index
  public :: idw_index
  !
  !    
  public :: twin_sector_order
  public :: get_twin_sector
  public :: flip_state
  !





contains






  !######################################################################
  !######################################################################
  !BUILD SECTORS: given isector build / delete the corresponding map(s)
  !######################################################################
  !######################################################################
  subroutine build_sector(isector,self,Itrace)
    !
    ! This procedure build the sector :f:type:`sector` :math:`{\cal S}(\vec{Q})` given the index :f:var:`isector` which univocally determines the quantum numbers :math:`\vec{Q}`. All the components of the data structure :f:type:`sector` are determined here. Moreover the map connecting the states of the sector to thos in the Fock space is constructed.
    ! To avoid integer overflow the loop over Fock space states is decomposed in spin up and down integer so that :math:`I = I_\uparrow + 2^{N}I_\downarrow`.  
    !
    integer                  :: isector
    type(sector)             :: self
    logical,optional         :: itrace
    logical                  :: itrace_
    integer,dimension(Ns_Ud) :: Nups,Ndws
    integer,dimension(Ns_Ud) :: DimUps,DimDws
    integer                  :: DimUp,DimDw,impDIM
    integer                  :: iup,idw
    integer                  :: nup_,ndw_
    integer                  :: imap,iud
    integer                  :: iIMP,iBATH
    !
    itrace_=.false. ; if(present(itrace))itrace_=itrace
    if(self%status)call delete_sector(self)
    !
    self%index = isector
    !
    impDim = 2**(Nimp/Ns_ud)
    !
    select case(ed_mode)
    case default
       allocate(self%H(2*Ns_Ud))
       allocate(self%DimUps(Ns_Ud))
       allocate(self%DimDws(Ns_Ud))
       allocate(self%Nups(Ns_Ud))
       allocate(self%Ndws(Ns_Ud))
       !
       call get_Nup(isector,Nups);self%Nup=sum(self%Nups)
       call get_Ndw(isector,Ndws);self%Ndw=sum(self%Ndws)
       call get_DimUp(isector,DimUps);self%DimUp=product(self%DimUps)
       call get_DimDw(isector,DimDws);self%DimDw=product(self%DimDws)
       !
       self%Dim = self%DimUp*self%DimDw
       !
       if(itrace_)then
          call map_allocate(H,[self%DimUps,self%DimDws],impDim)
       else
          call map_allocate(H,[self%DimUps,self%DimDws])
       endif
       !
       do iud=1,Ns_Ud
          !UP    
          imap=0
          do iup=0,2**Ns_Orb-1
             nup_ = popcnt(iup)  
             if(nup_ /= Nups(iud))cycle
             imap  = imap+1
             self%H(iud)%map(imap) = iup
             if(.not.itrace_)cycle
             iIMP  = ibits(iup,0,Nimp)
             iBATH = ibits(iup,Nimp,Nimp*Nbath)
             call sp_insert_state(self%H(iud)%sp,iIMP,iBATH,imap) 
          enddo
          !DW
          imap=0
          do idw=0,2**Ns_Orb-1
             ndw_= popcnt(idw)
             if(ndw_ /= Ndws(iud))cycle
             imap = imap+1
             self%H(iud+Ns_Ud)%map(imap) = idw
             if(.not.itrace_)cycle
             iIMP  = ibits(idw,0,Nimp)
             iBATH = ibits(idw,Nimp,Nimp*Nbath)
             call sp_insert_state(self%H(iud+Ns_Ud)%sp,iIMP,iBATH,imap) 
          enddo
       enddo
       !       
    case ("superc")
       allocate(self%H(1))
       self%Sz  = getSz(isector)
       self%Dim = getDim(isector)
       !
       if(itrace_)then
          call map_allocate(self%H,[self%Dim],impDim)
       else
          call map_allocate(self%H,[self%Dim])
       endif
       !
       imap=0
       do idw=0,2**Ns-1
          ndw_= popcnt(idw)
          do iup=0,2**Ns-1
             nup_ = popcnt(iup)
             sz_  = nup_ - ndw_
             if(sz_ /= self%Sz)cycle
             imap=imap+1
             self%H(1)%map(imap) = iup + idw*2**Ns
             if(.not.itrace_)cycle
             iIMP  = ibits(iup + idw*2**Ns,0,Nimp)
             iBATH = ibits(iup + idw*2**Ns,Nimp,Nimp*Nbath)
             call sp_insert_state(self%H(1)%sp,iIMP,iBATH,imap) 
          enddo
       enddo

    end select
    !
  end subroutine build_sector


  subroutine delete_sector(self)
    type(sector) :: self
    call map_deallocate(self%H)
    if(allocated(self%H))deallocate(self%H)
    if(allocated(self%DimUps))deallocate(self%DimUps)
    if(allocated(self%DimDws))deallocate(self%DimDws)
    if(allocated(self%Nups))deallocate(self%Nups)
    if(allocated(self%Ndws))deallocate(self%Ndws)
    self%index=0
    self%DimUp=0
    self%DimDw=0
    self%Dim=0
    self%Nup=0
    self%Ndw=0
    self%Sz=-1000
    self%Nlanc=0
    self%status=.false.
  end subroutine delete_sector










  !=========================================================
  subroutine map_allocate_scalar(H,N,Nsp)
    type(sector_map) :: H
    integer          :: N
    integer,optional :: Nsp
    if(H%status) call map_deallocate_scalar(H)
    allocate(H%map(N))
    if(present(Nsp))call sp_init_map(H%sp,Nsp)
    H%status=.true.
  end subroutine map_allocate_scalar
  !
  subroutine map_allocate_vector(H,N,Nsp)
    type(sector_map),dimension(:) :: H
    integer,dimension(size(H))    :: N
    integer,optional              :: Nsp
    integer                       :: i
    do i=1,size(H)
       if(present(Nsp))then
          call map_allocate_scalar(H(i),N(i),Nsp)
       else
          call map_allocate_scalar(H(i),N(i))
       endif
    enddo
  end subroutine map_allocate_vector


  !=========================================================
  subroutine map_deallocate_scalar(H)
    type(sector_map) :: H
    if(.not.H%status)then
       write(*,*) "WARNING map_deallocate_scalar: H is not allocated"
       return
    endif
    if(allocated(H%map))deallocate(H%map)
    call sp_delete_map(H%sp)
    H%status=.false.
  end subroutine map_deallocate_scalar
  !
  subroutine map_deallocate_vector(H)
    type(sector_map),dimension(:) :: H
    integer                       :: i
    do i=1,size(H)
       call map_deallocate_scalar(H(i))
    enddo
  end subroutine map_deallocate_vector

  !=========================================================




  function apply_op_C(V,ipos,ispin,jsector,sectorI) result(OV)
    complex(8),dimension(:),intent(in)  :: V
    integer, intent(in)                 :: ipos,ispin,jsector
    type(sector),intent(in)             :: sectorI
    type(sector)                        :: sectorJ
    complex(8),dimension(:),allocatable :: OV
    integer                             :: i,j
    real(8)                             :: sgn
    integer                             :: r
    integer                             :: fi
    integer,dimension(2*Ns_Ud)          :: Indices
    integer,dimension(2*Ns_Ud)          :: Jndices
    integer,dimension(2,Ns_Orb)         :: Nud !Nbits(Ns_Orb)
    integer,dimension(2)                :: Iud
    integer,dimension(2*Ns)             :: ib
    !
    if(MpiMaster)then
       !
       if(size(V)/=sectorI%Dim)stop "apply_op_C ERROR: size(V) != sectorI.Dim"
       !
       call build_sector(jsector,sectorJ)
       !
       if(allocated(OV))deallocate(OV(sectorJ%Dim))
       allocate(OV(sectorJ%Dim)) ; OV=zero
       !
       !       
       if(ed_verbose>2)then
          select case(ed_mode)
          case default
             write(LOGfile,"(2(A,I6,2I4))")&
                  'From:',isector,sectorI%Nups,sectorI%Ndws,&
                  ' -> apply C:',jsector,sectorJ%Nups,sectorJ%Ndws
          case ("superc")
             write(LOGfile,"(2(A,I6,I3))")&
                  'From:',isector,sectorI%Sz,&
                  'apply C:',jsector,sectorJ%Sz
          end select
       endif
       !
       do i=1,sectorI%Dim
          !
          select case(ed_mode)
          case default
             call state2indices(i,[sectorI%DimUps,sectorI%DimDws],Indices)
             iud(1)   = sectorI%H(1)%map(Indices(1))
             iud(2)   = sectorI%H(2)%map(Indices(2))
             Nud(1,:) = Bdecomp(iud(1),Ns_Orb)
             Nud(2,:) = Bdecomp(iud(2),Ns_Orb)
             if(Nud(ispin,ipos)/=1)cycle
             call c(ipos,iud(ispin),r,sgn)
             Jndices        = Indices
             Jndices(ispin) = binary_search(sectorJ%H(ispin)%map,r)
             call indices2state(Jndices,[sectorJ%DimUps,sectorJ%DimDws],j)
          case("superc","nonsu2")
             fi = sectorI%H(1)%map(i)
             ib = bdecomp(fi,2*Ns)
             if(ib(ipos + (ispin-1)*Ns)/=1)cycle
             call c(ipos + (ispin-1)*Ns,fi,r,sgn)
             j  = binary_search(sectorJ%H(1)%map,r)
          end select
          !
          OV(j) = sgn*V(i)
       enddo
       call delete_sector(sectorJ)
    else
       if(allocated(OV))deallocate(OV)
       allocate(OV(1)) ; OV=zero
    end if
    !
  end function apply_op_C


  function apply_op_CDG(V,ipos,ispin,jsector,sectorI) result(OV)
    complex(8),dimension(:),intent(in)  :: V
    integer, intent(in)                 :: ipos,ispin,jsector
    type(sector),intent(in)             :: sectorI
    type(sector)                        :: sectorJ
    complex(8),dimension(:),allocatable :: OV
    integer                             :: i,j
    real(8)                             :: sgn
    integer                             :: r
    integer                             :: fi
    integer,dimension(2*Ns_Ud)          :: Indices
    integer,dimension(2*Ns_Ud)          :: Jndices
    integer,dimension(2,Ns_Orb)         :: Nud !Nbits(Ns_Orb)
    integer,dimension(2)                :: Iud
    integer,dimension(2*Ns)             :: ib
    !
    if(MpiMaster)then
       !
       if(size(V)/=sectorI%Dim)stop "apply_op_CDG ERROR: size(V) != sectorI.Dim"
       !
       call build_sector(jsector,sectorJ)
       !
       if(allocated(OV))deallocate(OV(sectorJ%Dim))
       allocate(OV(sectorJ%Dim)) ; OV=zero
       !
       if(ed_verbose>2)then
          select case(ed_mode)
          case default
             write(LOGfile,"(2(A,I6,2I4))")&
                  'From:',isector,sectorI%Nups,sectorI%Ndws,&
                  ' -> apply C^+:',jsector,sectorJ%Nups,sectorJ%Ndws
          case ("superc")
             write(LOGfile,"(2(A,I6,I3))")&
                  'From:',isector,sectorI%Sz,&
                  'apply C^+:',jsector,sectorJ%Sz
          end select
       endif
       !
       do i=1,sectorI%Dim
          !
          select case(ed_mode)
          case default
             call state2indices(i,[sectorI%DimUps,sectorI%DimDws],Indices)
             iud(1)   = sectorI%H(1)%map(Indices(1))
             iud(2)   = sectorI%H(2)%map(Indices(2))
             Nud(1,:) = Bdecomp(iud(1),Ns_Orb)
             Nud(2,:) = Bdecomp(iud(2),Ns_Orb)
             if(Nud(ispin,ipos)/=0)cycle
             call cdg(ipos,iud(ispin),r,sgn)
             Jndices        = Indices
             Jndices(ispin) = binary_search(sectorJ%H(ispin)%map,r)
             call indices2state(Jndices,[sectorJ%DimUps,sectorJ%DimDws],j)
          case("superc","nonsu2")
             fi = sectorI%H(1)%map(i)
             ib = bdecomp(fi,2*Ns)
             if(ib(ipos + (ispin-1)*Ns)/=0)cycle
             call cdg(ipos + (ispin-1)*Ns,fi,r,sgn)
             j  = binary_search(sectorJ%H(1)%map,r)
          end select
          !
          OV(j) = OV(j) + sgn*V(i)
       enddo
       !
       call delete_sector(sectorJ)
       !
    else
       if(allocated(OV))deallocate(OV)
       allocate(OV(1)) ; OV=zero
    end if
  end function apply_op_CDG



  function apply_COps(V,As,Os,Pos,Spin,jsector,sectorI) result(OV)
    complex(8),dimension(:),intent(in)     :: V
    complex(8),dimension(:),intent(in)     :: As
    integer,dimension(size(As)),intent(in) :: Os
    integer,dimension(size(As)),intent(in) :: Pos,Spin
    integer, intent(in)                    :: jsector
    type(sector),intent(in)                :: sectorI
    type(sector)                           :: sectorJ
    complex(8),dimension(:),allocatable    :: OV
    integer                                :: ipos,ispin,ios
    integer                                :: i,j,Nos,is
    real(8)                                :: sgn
    integer                                :: r
    integer                                :: fi
    integer,dimension(2*Ns_Ud)             :: Indices
    integer,dimension(2*Ns_Ud)             :: Jndices
    integer,dimension(2,Ns_Orb)            :: Nud !Nbits(Ns_Orb)
    integer,dimension(2)                   :: Iud
    integer,dimension(2*Ns)                :: ib
    character(3),dimension(-1:1)           :: Cstr = ["C  "," x ","C^+"]
    !
    if(MpiMaster)then
       !
       if(size(V)/=sectorI%Dim)stop "apply_COps ERROR: size(V) != sectorI.Dim"
       !
       call build_sector(jsector,sectorJ)
       !
       if(allocated(OV))deallocate(OV)
       allocate(OV(sectorJ%Dim)) ; OV=zero
       !
       if(ed_verbose>2)then
          select case(ed_mode)
          case default
             write(LOGfile,"(2(A,I6,2I4))")&
                  'From:',isector,sectorI%Nups,sectorI%Ndws,&
                  ' -> apply C:',jsector,sectorJ%Nups,sectorJ%Ndws
          case ("superc")
             write(LOGfile,"(2(A,I6,I3))")&
                  'From:',isector,sectorI%Sz,&
                  'apply C:',jsector,sectorJ%Sz
          end select
       endif

       do is=1,size(As)
          ipos  = Pos(is)
          ispin = Spin(is)
          ios   = Os(is)
          if(ed_verbose>2)write(LOGfile,"(A)")
          'apply '//str(Cstr(ios))//" l:"//str(ipos)//" s:"//str(ispin)
          !
          do i=1,sectorI%Dim
             select case(ed_mode)
             case default
                call state2indices(i,[sectorI%DimUps,sectorI%DimDws],Indices)
                iud(1)   = sectorI%H(1)%map(Indices(1))
                iud(2)   = sectorI%H(2)%map(Indices(2))
                Nud(1,:) = Bdecomp(iud(1),Ns_Orb)
                Nud(2,:) = Bdecomp(iud(2),Ns_Orb)             
                select case(ios)
                case default;stop "apply_COps ERROR: ios sign not \in [-1,1]"
                case(-1)
                   if(Nud(ispin,ipos)/=1)cycle
                   call c(ipos,iud(ispin),r,sgn)
                case(1)
                   if(Nud(ispin,ipos)/=0)cycle
                   call cdg(ipos,iud(ispin),r,sgn)
                end select
                Jndices        = Indices
                Jndices(ispin) = binary_search(sectorJ%H(ispin)%map,r)
                call indices2state(Jndices,[sectorJ%DimUps,sectorJ%DimDws],j)
             case("superc","nonsu2")
                fi = sectorI%H(1)%map(i)
                ib = bdecomp(fi,2*Ns)
                select case(ios)
                case default;stop "apply_COps ERROR: ios sign not \in [-1,1]"
                case(-1)
                   if(ib(ipos + (ispin-1)*Ns)/=1)cycle
                   call c(ipos + (ispin-1)*Ns,fi,r,sgn)
                case(1)
                   if(ib(ipos + (ispin-1)*Ns)/=0)cycle
                   call cdg(ipos + (ispin-1)*Ns,fi,r,sgn)
                end select
                j  = binary_search(sectorJ%H(1)%map,r)
             end select
             !
             OV(j) = OV(j) + sgn*V(i)*As(is)
             !
          enddo
       enddo
       call delete_sector(sectorJ)
    else
       if(allocated(OV))deallocate(OV)
       allocate(OV(1)) ; OV=zero
    end if
  end function apply_COps






  ! subroutine build_op_Ns(i,Nup,Ndw,sectorI) 
  !   integer, intent(in)             :: i
  !   type(sector),intent(in)         :: sectorI
  !   integer,dimension(Ns)           :: Nup,Ndw  ![Ns]
  !   integer                         :: iph,i_el,ii,iorb
  !   integer,dimension(2*Ns_Ud)      :: Indices
  !   integer,dimension(Ns_Ud,Ns_Orb) :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
  !   integer,dimension(2*Ns)         :: Ib
  !   integer,dimension(2)            :: Iud
  !   !
  !   select case(ed_mode)
  !   case default
  !      iph = (i-1)/(sectorI%DimEl) + 1
  !      i_el = mod(i-1,sectorI%DimEl) + 1
  !      !
  !      call state2indices(i_el,[sectorI%DimUps,sectorI%DimDws],Indices)
  !      do ii=1,Ns_Ud
  !         iud(1) = sectorI%H(ii)%map(Indices(ii))
  !         iud(2) = sectorI%H(ii+Ns_Ud)%map(Indices(ii+Ns_ud))
  !         Nups(ii,:) = Bdecomp(iud(1),Ns_Orb) ![Norb,1+Nbath]
  !         Ndws(ii,:) = Bdecomp(iud(2),Ns_Orb)
  !      enddo
  !      Nup = Breorder(Nups)
  !      Ndw = Breorder(Ndws)
  !      !
  !   case("superc","nonsu2")
  !      ii = sectorI%H(1)%map(i)
  !      Ib = bdecomp(ii,2*Ns)
  !      Nup = Ib(1:Ns)
  !      Ndw = Ib(Ns+1:)
  !      ! do ii=1,Ns
  !      !    Nup(ii)= ib(ii)
  !      !    Ndw(ii)= ib(ii+Ns)
  !      ! enddo
  !   end select
  !   !
  ! end subroutine build_op_Ns



















  subroutine get_Sector_normal(QN,N,isector)
    integer,dimension(:) :: QN
    integer              :: N
    integer              :: isector
    integer              :: i,Nind,factor
    Nind = size(QN)
    Factor = N+1
    isector = 1
    do i=Nind,1,-1
       isector = isector + QN(i)*(Factor)**(Nind-i)
    enddo
  end subroutine get_Sector_normal
  !
  subroutine get_Sector_superc(QN,isector)
    integer :: QN
    integer :: isector
    isector=getSector(QN,1)
  end subroutine get_Sector_superc




  subroutine get_QuantumNumbers_normal(isector,N,QN)
    integer                          :: isector,N
    integer,dimension(:)             :: QN
    integer                          :: i,count,Dim
    integer,dimension(size(QN)) :: QN_
    !
    Dim = size(QN)
    if(mod(Dim,2)/=0)stop "get_Indices_main error: Dim%2 != 0"
    count=isector-1
    do i=1,Dim
       QN_(i) = mod(count,N+1)
       count      = count/(N+1)
    enddo
    QN = QN_(Dim:1:-1)
  end subroutine get_QuantumNumbers_normal

  subroutine get_QuantumNumbers_other(isector,QN)
    integer                     :: isector
    integer                     :: QN
    select case(ed_mode)
    case ("superc")
       QN=getSz(isector)
    case default
       stop "get_QuantumNumbers_other ERROR: invoked with ed_mode=normal"
    end select
  end subroutine get_QuantumNumbers_other







  subroutine get_Nup(isector,Nup)
    integer                   :: isector,Nup(Ns_Ud)
    integer                   :: i,count
    integer,dimension(2*Ns_Ud)  :: indices_
    count=isector-1
    do i=1,2*Ns_Ud
       indices_(i) = mod(count,Ns_Orb+1)
       count      = count/(Ns_Orb+1)
    enddo
    Nup = indices_(2*Ns_Ud:Ns_Ud+1:-1)
  end subroutine get_Nup

  subroutine get_Ndw(isector,Ndw)
    integer                   :: isector,Ndw(Ns_Ud)
    integer                   :: i,count
    integer,dimension(2*Ns_Ud) :: indices_
    count=isector-1
    do i=1,2*Ns_Ud
       indices_(i) = mod(count,Ns_Orb+1)
       count      = count/(Ns_Orb+1)
    enddo
    Ndw = indices_(Ns_Ud:1:-1)
  end subroutine get_Ndw

  subroutine get_Sz(isector,Sz)
    integer :: isector,Sz
    Sz = getSz(isector)
  end subroutine get_Sz







  subroutine  get_DimUp(isector,DimUps)
    integer                :: isector,DimUps(Ns_Ud)
    integer                :: Nups(Ns_Ud),iud
    call get_Nup(isector,Nups)
    do iud=1,Ns_Ud
       DimUps(iud) = binomial(Ns_Orb,Nups(iud))
    enddo
  end subroutine get_DimUp


  subroutine get_DimDw(isector,DimDws)
    integer                :: isector,DimDws(Ns_Ud)
    integer                :: Ndws(Ns_Ud),iud
    call get_Ndw(isector,Ndws)
    do iud=1,Ns_Ud
       DimDws(iud) = binomial(Ns_Orb,Ndws(iud))
    enddo
  end subroutine get_DimDw

  subroutine  get_Dim(isector,Dim)
    !
    !Returns the dimension of the symmetry sector with index :f:var:`isector`
    !
    integer                :: isector,Dim
    Dim=getDim(isector)
  end subroutine get_Dim






  subroutine indices2state(ivec,Nvec,istate)
    integer,dimension(:)          :: ivec
    integer,dimension(size(ivec)) :: Nvec
    integer                       :: istate,i
    istate=ivec(1)
    do i=2,size(ivec)
       istate = istate + (ivec(i)-1)*product(Nvec(1:i-1))
    enddo
  end subroutine indices2state

  subroutine state2indices(istate,Nvec,ivec)
    integer                       :: istate
    integer,dimension(:)          :: Nvec
    integer,dimension(size(Nvec)) :: Ivec
    integer                       :: i,count,N
    count = istate-1
    N     = size(Nvec)
    do i=1,N
       Ivec(i) = mod(count,Nvec(i))+1
       count   = count/Nvec(i)
    enddo
  end subroutine state2indices


  function iup_index(i,DimUp) result(iup)
    integer :: i
    integer :: DimUp
    integer :: iup
    iup = mod(i,DimUp);if(iup==0)iup=DimUp
  end function iup_index


  function idw_index(i,DimUp) result(idw)
    integer :: i
    integer :: DimUp
    integer :: idw
    idw = (i-1)/DimUp+1
  end function idw_index

















  !##################################################################
  !##################################################################
  !TWIN SECTORS ROUTINES:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : Build the re-ordering map to go from sector A(nup,ndw)
  ! to its twin sector B(ndw,nup), with nup!=ndw.
  !
  !- build the map from the A-sector to \HHH
  !- get the list of states in \HHH corresponding to sector B twin of A
  !- return the ordering of B-states in \HHH with respect to those of A
  !+------------------------------------------------------------------+
  subroutine twin_sector_order(isector,order)
    integer                             :: isector
    integer,dimension(:)                :: order
    type(sector)                        :: sectorH
    type(sector_map),dimension(2*Ns_Ud) :: H
    integer,dimension(2*Ns_Ud)          :: Indices,Istates
    integer,dimension(Ns_Ud)            :: DimUps,DimDws
    integer                             :: Dim
    integer                             :: i,iud
    !
    if(size(Order)/=Dim)stop "twin_sector_order error: wrong dimensions of *order* array"
    !
    call build_sector(isector,sectorH)
    select case(ed_mode)
    case ("normal")
       do i=1,sectorH%Dim
          call state2indices(i,[sectorH%DimUps,sectorH%DimDws],Indices)
          forall(iud=1:2*Ns_Ud)Istates(iud) = sectorH%H(iud)%map(Indices(iud))
          Order(i) = flip_state( Istates )
       enddo
    case default
       do i=1,sectorH%dim
          Order(i) = flip_state(sectorH%H(1)%map(i))
       enddo
    end select
    !
    call delete_sector(sectorH)
    call sort_array(Order)
    !
  end subroutine twin_sector_order



  !+------------------------------------------------------------------+
  !PURPOSE  : Flip an Hilbert space state m=|{up}>|{dw}> into:
  !
  ! normal: j=|{dw}>|{up}>  , nup --> ndw
  !+------------------------------------------------------------------+
  function flip_state_normal(istate) result(j)
    integer,dimension(2*Ns_Ud) :: istate
    integer                    :: j
    integer,dimension(Ns_Ud)   :: jups,jdws
    integer,dimension(2*Ns_Ud) :: dims
    !
    jups = istate(Ns_Ud+1:2*Ns_Ud)
    jdws = istate(1:Ns_Ud)
    dims = 2**Ns_Orb
    call indices2state([jups,jdws],Dims,j)
    !
  end function flip_state_normal
  function flip_state_other(istate) result(j)
    integer          :: istate
    integer          :: j
    integer          :: ivec(2*Ns),foo(2*Ns)
    !
    Ivec = bdecomp(istate,2*Ns)
    select case(ed_mode)
    case("superc")  !Invert the overall spin sign: |{up}> <---> |{dw}>
       foo(1:Ns)     = Ivec(Ns+1:2*Ns)
       foo(Ns+1:2*Ns)= Ivec(1:Ns)
    case default    !Exchange UP-config |{up}> with DW-config |{dw}>
       stop "flip_state_other error: called with ed_mode==normal"
    end select
    !
    j = bjoin(foo,2*Ns)
    !
  end function flip_state_other


  !+------------------------------------------------------------------+
  !PURPOSE  : get the twin of a given sector (the one with opposite 
  ! quantum numbers): 
  ! nup,ndw ==> ndw,nup (spin-exchange)
  !+------------------------------------------------------------------+
  function get_twin_sector(isector) result(jsector)
    integer,intent(in)       :: isector
    integer                  :: jsector
    integer,dimension(Ns_Ud) :: Iups,Idws
    integer                  :: Sz
    select case(ed_mode)
    case default
       call get_Nup(isector,iups)
       call get_Ndw(isector,idws)
       call get_Sector([idws,iups],Ns_Orb,jsector)
    case ("superc")
       call get_Sz(isector,Sz)
       Sz = -Sz
       call get_Sector(Sz,jsector)
    end select
  end function get_twin_sector










  !##################################################################
  !##################################################################
  !AUXILIARY COMPUTATIONAL ROUTINES ARE HERE BELOW:
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE : sort array of integer using random algorithm
  !+------------------------------------------------------------------+
  subroutine sort_array(array)
    integer,dimension(:),intent(inout)      :: array
    integer,dimension(size(array))          :: order
    integer                                 :: i
    forall(i=1:size(array))order(i)=i
    call qsort_sort( array, order, 1, size(array) )
    array=order
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      integer, dimension(:)                 :: array
      integer, dimension(:)                 :: order
      integer                               :: left
      integer                               :: right
      integer                               :: i
      integer                               :: last
      if ( left .ge. right ) return
      call qsort_swap( order, left, qsort_rand(left,right) )
      last = left
      do i = left+1, right
         if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
            last = last + 1
            call qsort_swap( order, last, i )
         endif
      enddo
      call qsort_swap( order, left, last )
      call qsort_sort( array, order, left, last-1 )
      call qsort_sort( array, order, last+1, right )
    end subroutine qsort_sort
    !---------------------------------------------!
    subroutine qsort_swap( order, first, second )
      integer, dimension(:)                 :: order
      integer                               :: first, second
      integer                               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine qsort_swap
    !---------------------------------------------!
    function qsort_rand( lower, upper )
      implicit none
      integer                               :: lower, upper
      real(8)                               :: r
      integer                               :: qsort_rand
      call random_number(r)
      qsort_rand =  lower + nint(r * (upper-lower))
    end function qsort_rand
    function compare(f,g)
      integer                               :: f,g
      integer                               :: compare
      compare=1
      if(f<g)compare=-1
    end function compare
  end subroutine sort_array









  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the binomial factor n1 over n2
  !+------------------------------------------------------------------+
  elemental function binomial(n1,n2) result(nchoos)
    integer,intent(in) :: n1,n2
    real(8)            :: xh
    integer            :: i
    integer nchoos
    xh = 1.d0
    if(n2<0) then
       nchoos = 0
       return
    endif
    if(n2==0) then
       nchoos = 1
       return
    endif
    do i = 1,n2
       xh = xh*dble(n1+1-i)/dble(i)
    enddo
    nchoos = int(xh + 0.5d0)
  end function binomial





end MODULE ED_SECTOR













!   ! SPIN-RESOLVED ALGORITHM: two map objects Hsigma, with a *sparse* 
!   !                          structure capable to store separately
!   !                          the impurity and the bath spin-states
!   subroutine build_sector_spin(isector,Hup,Hdw)
!     integer                            :: isector
!     type(sector_map)                   :: HUP !Map for the UPs
!     type(sector_map)                   :: HDW !Map for the Dws
!     integer,dimension(Ns_Ud)           :: Nups,Ndws
!     integer,dimension(Ns_Ud)           :: DimUps,DimDws
!     integer                            :: DimUp,DimDw,impDIM
!     integer                            :: iup,idw
!     integer                            :: nup_,ndw_
!     integer                            :: imap,iud
!     integer                            :: iIMP,iBATH
!     !
!     impDIM = 2**(Nimp/Ns_ud) !Number of states for the impurity
!     !
!     !Init UP sub-sector:
!     ! > allocate UP-map to binomial(Ns_Orb Nup)
!     call get_Nup(isector,Nups)
!     call get_DimUp(isector,DimUps); DimUp = product(DimUps)
!     call map_allocate(HUP,DimUp,impDIM)
!     !
!     !Init DW sub-sector:
!     ! > allocate DW-map to binomial(Ns_Orb Ndw)
!     call get_Ndw(isector,Ndws)
!     call get_DimDw(isector,DimDws); DimDw = product(DimDws)
!     call map_allocate(HDW,DimDw,impDIM)
!     !
!     !Formally dealing with ed_total_ud == .false.
!     !->iud=1:Ns_Ud (formally because Ns_Ud = 1 in CDMFT code...)
!     ![STILL PROBLEMS WITH Ns_Ud, don't know how to deal with it in density_matrix_impurity()] 
!     do iud=1,Ns_Ud 
! #ifdef _DEBUG
!        if(ed_verbose>3)write(Logfile,"(A)")&
!             "  DEBUG build_sector_spin(): working UP sub-sector"
! #endif
!        imap=0
!        do iup=0,2**Ns_Orb-1
!           nup_ = popcnt(iup) !equivalent to sum(binary_decomposition(iup))
!           if(nup_ /= Nups(iud))cycle !the state does not have the required number of UPs
!           imap = imap+1
!           !HUP(iud)%map(imap) = iup
!           HUP%map(imap) = iup
!           !
!           iIMP  = ibits(iup,0,Nimp)
!           iBATH = ibits(iup,Nimp,Nimp*Nbath) !check: Nimp+Nimp*Nbath=Nimp(1*Nbath)=Ns
!           !call sp_insert_state(HUP(iud)%sp,iIMP,iBATH,imap)
! #ifdef _DEBUG
!           if(ed_verbose>4)then 
!              write(Logfile,"(A)")&
!                   "    DEBUG build_sector_spin(): inserting UP state"
!              write(Logfile,"(A)")&
!                   "      Iup: "//str(iup)//" | IimpUp: "//str(iIMP)//" | IbathUp: "//str(iBATH)
!           endif
! #endif
!           call sp_insert_state(HUP%sp,iIMP,iBATH,imap) 
!           !
!        enddo
!        !
! #ifdef _DEBUG
!        if(ed_verbose>3)write(Logfile,"(A)")&
!             "  DEBUG build_sector_spin(): working DW sub-sector"
! #endif
!        imap=0
!        do idw=0,2**Ns_Orb-1
!           ndw_=popcnt(idw)   !equivalent to sum(binary_decomposition(idw))
!           if(ndw_ /= Ndws(iud))cycle !the state does not have the required number of DWs
!           imap = imap+1
!           !HDW(iud)%map(imap) = idw
!           HDW%map(imap) = idw
!           !
!           iIMP  = ibits(idw,0,Nimp)
!           iBATH = ibits(idw,Nimp,Nimp*Nbath) !check: Nimp+Nimp*Nbath=Nimp(1*Nbath)=Ns
!           !call sp_insert_state(HDW(iud)%sp,iIMP,iBATH,imap)
! #ifdef _DEBUG
!           if(ed_verbose>4)then 
!              write(LOGfile,"(A)")&
!                   "    DEBUG build_sector_spin(): inserting DW state"
!              write(LOGfile,"(A)")&
!                   "      Idw: "//str(idw)//" | IimpDw: "//str(iIMP)//" | IbathDw: "//str(iBATH)
!           endif
! #endif
!           call sp_insert_state(HDW%sp,iIMP,iBATH,imap)
!           !
!        enddo
!     enddo
!     !
!   end subroutine build_sector_spin




! getCsector=0
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup-1; jdw=ndw; if(jup < 0)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCsector(1,1,isector)=jsector
! enddo
! !
! !
! !
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup;jdw=ndw-1;if(jdw < 0)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCsector(1,2,isector)=jsector
! enddo
! !
! !
! !
! getCDGsector=0
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup+1;jdw=ndw;if(jup > Ns)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCDGsector(1,1,isector)=jsector
! enddo
! !
! !
! !
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup;jdw=ndw+1;if(jdw > Ns)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCDGsector(1,2,isector)=jsector
! enddo

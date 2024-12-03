MODULE ED_BATH_REPLICA
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
  implicit none


  private



  interface set_Hbath
     !
     ! This function sets the matrix basis :math:`\{ \hat{O}_i \}_{i=1,\dots,N_{sym}}` used to decompose the single bath hamiltonian :math:`h^p`. It also sets the initial values of the variational parameters :math:`\vec{\lambda}` in the :f:var:`replica` ath type.
     !
     ! Input: :f:var:`Hvec` 
     !  * rank-7, matrix basis: [ |Nambu| , |Nambu| , |Nspin| , |Nspin| , |Nimp| , |Nimp| , |Nsym| ]
     !  * rank-5, matrix basis: [ |Nspin| , |Nspin| , |Nimp|  , |Nimp|  , |Nsym| ] only *normal*
     !  * rank-3, matrix basis: [ |Ntot|  , |Ntot|  , |Nsym| ]
     !  * rank-6, Hloc [ |Nlat| , |Nlat|  , |Nspin| , |Nspin| , |Norb| , |Norb| ] only *normal* deprecated
     !  * rank-4, Hloc [ |Nspin| , |Nspin| , |Nimp| , |Nimp| ] only *normal* deprecated
     !  * rank-2, Hloc [ |Nlso| , |Nlso| ] only *normal* deprecated
     !
     ! Input: :f:var:`lambdavec`
     !  * rank-2, dimensions: [ |Nbath| , |Nsym| ]
     !  * rank-1, dimensions: [ |Nsym| ] Legacy support: degenerate parameters. deprecated.
     !
     module procedure init_Hbath_symmetries_d7
     module procedure init_Hbath_symmetries_d5
     module procedure init_Hbath_symmetries_d3
     module procedure init_Hbath_symmetries_LEGACY ! (deprecation-cycle)
     module procedure init_Hbath_direct_lso6         !only *normal* deprecated
     module procedure init_Hbath_direct_lso4         !only *normal* deprecated
     module procedure init_Hbath_direct_lso2         !only *normal* deprecated
  end interface set_Hbath



  public :: Hbath_build
  public :: set_Hbath


  logical :: bool
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


  !allocate GLOBAL basis for H (used for Hbath) and vectors coefficient
  subroutine allocate_Hbath(Nsym)
    integer          :: Nsym
    integer          :: isym
    !
    if(allocated(Hbath_basis))deallocate(Hbath_basis)
    if(allocated(Hbath_lambda))deallocate(Hbath_lambda)
    !
    allocate(Hbath_basis(Nsym))
    allocate(Hbath_lambda(Nbath,Nsym))
    do isym=1,Nsym
       allocate(Hbath_basis(isym)%O(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp))
       Hbath_basis(isym)%O=zero
       Hbath_lambda(:,isym)=0d0
    enddo
    Hbath_status=.true.
  end subroutine allocate_Hbath


  !deallocate GLOBAL basis for H (used for impHloc and bath) and vectors coefficient
  subroutine deallocate_Hbath()
    integer              :: isym
    !
    do isym=1,size(Hbath_basis)
       if(allocated(Hbath_basis(isym)%O))deallocate(Hbath_basis(isym)%O)
    enddo
    if(allocated(Hbath_basis))deallocate(Hbath_basis)
    if(allocated(Hbath_lambda))deallocate(Hbath_lambda)
  end subroutine deallocate_Hbath






  !##################################################################
  !##################################################################
  !##################################################################


  function Hbath_build(lambdavec) result(H)
    !
    !This function is used to reconstruct the local bath Hamiltonian from basis expansion given the vector of :math:`\vec{\lambda}` parameters :math:`h^p=\sum_i \lambda^p_i O_i`. The resulting Hamiltonian has dimensions [ |Nambu| , |Nambu| , |Nspin| , |Nspin| , |Nimp| , |Nimp| ]
    !
    real(8),dimension(:),optional                           :: lambdavec  !the input vector of bath parameters
    real(8),dimension(:),allocatable                        :: lambda
    integer                                                 :: isym
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp) :: H
    !
    if(.not.Hbath_status)STOP "ERROR Hbath_build: Hbath_basis is not setup"
    allocate(lambda(size(Hbath_basis)));lambda=1d0
    !
    if(present(lambdavec))then
       if(size(lambdavec)/=size(Hbath_basis)) STOP "ERROR Hbath_build: Wrong coefficient vector size"
       lambda = lambdavec
    endif
    !
    H=zero
    do isym=1,size(lambda)
       H=H+lambda(isym)*Hbath_basis(isym)%O
    enddo
  end function Hbath_build




  !##################################################################
  !##################################################################
  !##################################################################



  subroutine init_Hbath_symmetries_d7(Hvec,lambdavec)
    complex(8),dimension(:,:,:,:,:,:,:)   :: Hvec       ![Nambu,Nambu,Nspin,Nspin,Nimp,Nimp][Nsym]
    real(8),dimension(:,:)                :: lambdavec ![Nbath,Nsym]
    integer                               :: isym,Nsym,ibath
    complex(8),dimension(:,:),allocatable :: H
    !
    call assert_shape(Hvec,[Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(lambdavec,2)],"init_Hbath_symmetries_d7","Hvec")
    !
    Nsym=size(lambdavec,2)
    if(size(lambdavec(:,1))/=Nbath) stop "init_Hbath_symmetries_d7 ERROR: size(lambdavec,1) != Nbath"
    !
    call allocate_Hbath(Nsym)
    !
    select case(ed_mode)
    case default
       if(Nambu>1)stop "init_Hbath_symmetries_d7 ERROR: Nambu>1 with ed_mode=normal"
       do isym=1,Nsym
          bool = check_herm(Hvec(1,1,:,:,:,:,isym))
          if(.not.bool)stop "init_Hbath_symmetries_d7 ERROR: not Hermitian of replica basis O_"//str(isym)
          Hbath_lambda(:,isym) = lambdavec(:,isym)
          Hbath_basis(isym)%O  = Hvec(:,:,:,:,:,:,isym)
       enddo
    case("superc")
       if(Nambu==1)stop "init_Hbath_symmetries_d7 ERROR: Nambu=1 with ed_mode=superc"
       do isym=1,Nsym
          bool = check_nambu(Hvec(:,:,:,:,:,:,isym))
          if(.not.bool) stop "init_Hbath_symmetries_d7 ERROR: not Nambu of replica basis O_"//str(isym)
          !
          Hbath_lambda(:,isym) = lambdavec(:,isym)
          Hbath_basis(isym)%O  = Hvec(:,:,:,:,:,:,isym)
       enddo
    end select
    !  
    if(ed_verbose>2)then
       do ibath=1,Nbath
          write(*,*) "Hbath #"//str(ibath)//":"
          call print_hloc(Hbath_build(Hbath_lambda(ibath,:)))
       enddo
    endif
    !
  end subroutine init_Hbath_symmetries_d7




  subroutine init_Hbath_symmetries_d5(Hvec,lambdavec)
    complex(8),dimension(:,:,:,:,:)       :: Hvec       ![Nspin,Nspin,Nimp,Nimp][Nsym]
    real(8),dimension(:,:)                :: lambdavec ![Nbath,Nsym]
    integer                               :: isym,Nsym,ibath
    complex(8),dimension(:,:),allocatable :: H
    !
    call assert_shape(Hvec,[Nspin,Nspin,Nimp,Nimp,size(lambdavec,2)],"init_Hbath_symmetries_d5","Hvec")
    !
    Nsym=size(lambdavec,2)
    if(size(lambdavec(:,1))/=Nbath) stop "init_Hbath_symmetries_d5 ERROR: size(lambdavec,1) != Nbath"
    !
    call allocate_Hbath(Nsym)
    !
    select case(ed_mode)
    case default
       do isym=1,Nsym
          bool = check_herm(Hvec(:,:,:,:,isym))
          if(.not.bool) stop "init_Hbath_symmetries_d5 ERROR: not Hermitian of replica basis O_"//str(isym)
          Hbath_lambda(:,isym) = lambdavec(:,isym)
          Hbath_basis(isym)%O(1,1,:,:,:,:)  = Hvec(:,:,:,:,isym)
       enddo
    case("superc")
       stop "init_Hbath_symmetries_d5 ERROR: called with ed_mode=superc but shape(Hvec)=[Nspin,Nspin,Nimp,Nimp]"
    end select
    !  
    if(ed_verbose>2)then
       do ibath=1,Nbath
          write(*,*) "Hbath #"//str(ibath)//":"
          call print_hloc(Hbath_build(Hbath_lambda(ibath,:)))
       enddo
    endif
    !
  end subroutine init_Hbath_symmetries_d5




  subroutine init_Hbath_symmetries_d3(Hvec,lambdavec)
    complex(8),dimension(:,:,:) :: Hvec       ![Ntot,Ntot,Nsym] /Ntot = Nambu*Nspin*Nimp
    real(8),dimension(:,:)      :: lambdavec ![Nbath,Nsym]
    integer                     :: isym,Nsym,ibath
    !
    call assert_shape(Hvec,[Nambu*Nspin*Nimp,Nambu*Nspin*Nimp,size(lambdavec,2)],"init_Hbath_symmetries_d3","Hvec")
    !
    Nsym=size(lambdavec,2)
    if(size(lambdavec(:,1))/=Nbath) stop "init_Hbath_symmetries_d3 ERROR: size(lambdavec,1) != Nbath"
    !
    call allocate_Hbath(Nsym)
    !
    select case(ed_mode)
    case default
       do isym=1,Nsym
          bool = check_herm(Hvec(:,:,isym))
          if(.not.bool) stop "init_Hbath_symmetries_d3 ERROR: not Hermitian of replica basis O_"//str(isym)
          Hbath_lambda(:,isym) = lambdavec(:,isym)          
          Hbath_basis(isym)%O(1,1,:,:,:,:)  = so2nn_reshape(Hvec(:,:,isym),Nspin,Nimp)
       enddo
    case("superc")
       do isym=1,Nsym
          bool = check_nambu(Hvec(:,:,isym),Nspin*Nimp)
          if(.not.bool) stop "init_Hbath_symmetries_d3 ERROR: not Nambu of replica basis O_"//str(isym)
          Hbath_lambda(:,isym) = lambdavec(:,isym)
          do concurrent(in=1:Nambu,jn=1:Nambu,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Nimp,jorb=1:Nimp)
             i = iorb + (ispin-1)*Nimp + (in-1)*Nspin*Nimp
             j = jorb + (jspin-1)*Nimp + (jn-1)*Nspin*Nimp
             Hbath_basis(isym)%O(in,jn,ispin,jspin,iorb,jorb) = Hvec(i,j,isym)
          enddo
       enddo
    end select
    !
    if(ed_verbose>2)then
       do ibath=1,Nbath
          write(*,*) "Hbath #"//str(ibath)//":"
          call print_hloc(Hbath_build(Hbath_lambda(ibath,:)))
       enddo
    endif
    !
  end subroutine init_Hbath_symmetries_d3




  subroutine init_Hbath_symmetries_LEGACY(Hvec,lambdavec)
    complex(8),dimension(:,:,:,:,:,:,:)   :: Hvec      ![Nambu,Nambu,Nspin,Nspin,Nimp,Nimp][Nsym]
    real(8),dimension(:)                  :: lambdavec ![Nsym]
    integer                               :: isym,Nsym,ibath
    complex(8),dimension(:,:),allocatable :: H
    !
    Nsym=size(lambdavec)
    !
    call assert_shape(Hvec,[Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,Nsym],"init_Hbath_symmetries superc","Hvec")
    !
    !
    call allocate_Hbath(Nsym)
    !
    select case(ed_mode)
    case default
       do isym=1,Nsym
          bool = check_herm(Hvec(:,:,:,:,:,:,isym))
          if(.not.bool) stop "init_Hbath_symmetries_d6 ERROR: not Hermitian of replica basis O_"//str(isym)
          Hbath_lambda(:,isym) = lambdavec(isym)
          Hbath_basis(isym)%O  = Hvec(:,:,:,:,:,:,isym)
       enddo
       !
    case("superc")
       do isym=1,Nsym
          bool = check_nambu(Hvec(:,:,:,:,:,:,isym))
          if(.not.bool) stop "init_Hbath_symmetries_d6 ERROR: not Nambu of replica basis O_"//str(isym)
          Hbath_lambda(:,isym) = lambdavec(isym)
          Hbath_basis(isym)%O  = Hvec(:,:,:,:,:,:,isym)
       enddo
    end select
    !    
    ! PRINT DEPRECATION MESSAGE TO LOG
    write(*,*) "                                                                               "
    write(*,*) "WARNING: Passing a single lambdasym vector to ed_set_Hbath is /deprecated/. "
    write(*,*) "         You should instead define a different lambda for each bath component, "
    write(*,*) "         namely passing a [Nbath]x[Nsym] array instead of a [Nsym] vector.     "
    write(*,*) "         Your single lambda vector has been internally copied into the required"
    write(*,*) "         higher-rank array, so giving each replica the same set of lambdas.    "
    write(*,*) "         >>> This back-compatibility patch might be removed in a future update."
    write(*,*) "                                                                               "
    !
    if(ed_verbose>2)then
       do ibath=1,Nbath
          write(*,*) "Hbath #"//str(ibath)//":"
          call print_hloc(Hbath_build(Hbath_lambda(ibath,:)))
       enddo
    endif
    !
  end subroutine init_Hbath_symmetries_LEGACY





  !DEPRECATED
  !+------------------------------------------------------------------+
  !PURPOSE  : Set Hbath from user defined Hloc
  !1: [Nlat,Nlat,Nspin,Nspin,Norb,Norb]
  !2: [Nspin,Nspin,Nimp,Nimp]
  !3: [Nlso,Nlso]
  !+------------------------------------------------------------------+
  subroutine init_Hbath_direct_lso6(Hloc)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hloc
    integer                                               :: counter,Nsym
    !
    if(ed_mode=='superc')stop "init_Hbath_direct_nnn ERROR: called with ed_mode=superc"
    !
    !SPIN DIAGONAL
    counter=0
    do ispin=1,Nspin
       do jspin=1,Nspin
          do ilat=1,Nlat
             do jlat=1,Nlat
                do iorb=1,Norb
                   do jorb=1,Norb
                      io=index_stride_lso(ilat,ispin,iorb)
                      jo=index_stride_lso(jlat,jspin,jorb)
                      if((Hloc(ilat,jlat,ispin,jspin,iorb,jorb) /= zero).AND.(io<jo))then
                         if(DREAL(Hloc(ilat,jlat,ispin,jspin,iorb,jorb))/=0d0)counter=counter+1
                         if(DIMAG(Hloc(ilat,jlat,ispin,jspin,iorb,jorb))/=0d0)counter=counter+1
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    call allocate_Hbath(counter)
    !
    counter=0
    do ispin=1,Nspin
       do jspin=1,Nspin
          do ilat=1,Nlat
             do jlat=1,Nlat
                do iorb=1,Norb
                   do jorb=1,Norb
                      io=index_stride_lso(ilat,ispin,iorb)
                      jo=index_stride_lso(jlat,jspin,jorb)
                      if((Hloc(ilat,jlat,ispin,jspin,iorb,jorb) /= zero).AND.(io<jo))then
                         if(dreal(Hloc(ilat,jlat,ispin,jspin,iorb,jorb))/=0d0)then
                            counter=counter+1
                            i = iorb+(ilat-1)*Norb
                            j = jorb+(jlat-1)*Norb
                            Hbath_basis(counter)%O(1,1,ispin,jspin,i,j)=one
                            Hbath_basis(counter)%O(1,1,ispin,jspin,j,i)=one
                            Hbath_lambda(:,counter)=dreal(Hloc(1,1,ispin,ispin,i,j))
                         endif
                         !
                         if(dimag(Hloc(ilat,jlat,ispin,jspin,iorb,jorb))/=0d0)then
                            counter=counter+1
                            Hbath_basis(counter)%O(1,1,ispin,jspin,i,j)=xi
                            Hbath_basis(counter)%O(1,1,ispin,jspin,j,i)=xi
                            Hbath_lambda(:,counter)=dimag(Hloc(1,1,ispin,ispin,i,j))
                         endif
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine init_Hbath_direct_lso6


  subroutine init_Hbath_direct_lso4(Hloc)
    complex(8),dimension(Nspin,Nspin,Nimp,Nimp) :: Hloc
    integer                                     :: counter,Nsym
    !
    if(ed_mode=='superc')stop "init_Hbath_direct_nnn ERROR: called with ed_mode=superc"
    !
    !SPIN DIAGONAL
    counter=0
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Nimp
             do jorb=1,Nimp
                io=iorb + (ispin-1)*Nimp
                jo=jorb + (jspin-1)*Nimp
                if((Hloc(ispin,jspin,iorb,jorb) /= zero).AND.(io<jo))then
                   if(dreal(Hloc(ispin,jspin,iorb,jorb))/=0d0)counter=counter+1
                   if(dimag(Hloc(ispin,jspin,iorb,jorb))/=0d0)counter=counter+1
                endif
             enddo
          enddo
       enddo
    enddo
    !
    call allocate_Hbath(counter)
    !
    counter=0
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Nimp
             do jorb=1,Nimp
                io=iorb + (ispin-1)*Nimp
                jo=jorb + (jspin-1)*Nimp
                if((Hloc(ispin,jspin,iorb,jorb) /= zero).AND.(io<jo))then
                   if(dreal(Hloc(ispin,jspin,iorb,jorb))/=0d0)then
                      counter=counter+1
                      Hbath_basis(counter)%O(1,1,ispin,jspin,iorb,jorb)=one
                      Hbath_basis(counter)%O(1,1,ispin,jspin,jorb,iorb)=one
                      Hbath_lambda(:,counter)=dreal(Hloc(ispin,ispin,iorb,jorb))
                   endif
                   !
                   if(dimag(Hloc(ispin,jspin,iorb,jorb))/=0d0)then
                      counter=counter+1
                      Hbath_basis(counter)%O(1,1,ispin,jspin,iorb,jorb)=xi
                      Hbath_basis(counter)%O(1,1,ispin,jspin,jorb,iorb)=xi
                      Hbath_lambda(:,counter)=dimag(Hloc(ispin,ispin,iorb,jorb))
                   endif
                endif
             enddo
          enddo
       enddo
    enddo
  end subroutine init_Hbath_direct_lso4

  subroutine init_Hbath_direct_lso2(Hloc)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hloc
    integer                                               :: counter,Nsym
    !
    if(ed_mode=='superc')stop "init_Hbath_direct_lso ERROR: called with ed_mode=superc"
    !
    !SPIN DIAGONAL
    do ispin=1,Nspin
       do jspin=1,Nspin
          do ilat=1,Nlat
             do jlat=1,Nlat
                do iorb=1,Norb
                   do jorb=1,Norb
                      io=index_stride_lso(ilat,ispin,iorb)
                      jo=index_stride_lso(jlat,jspin,jorb)
                      if((Hloc(io,jo)/=zero).AND.(io<jo))then
                         if(dreal(Hloc(io,jo))/=0d0)counter=counter+1
                         if(dimag(Hloc(io,jo))/=0d0)counter=counter+1
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    call allocate_Hbath(counter)
    !
    counter=0
    !
    do ispin=1,Nspin
       do jspin=1,Nspin
          do ilat=1,Nlat
             do jlat=1,Nlat
                do iorb=1,Norb
                   do jorb=1,Norb
                      io=index_stride_lso(ilat,ispin,iorb)
                      jo=index_stride_lso(jlat,jspin,jorb)
                      if((Hloc(io,jo)/=zero).AND.(io<jo))then
                         if(dreal(Hloc(io,jo))/=0d0)then
                            counter=counter+1
                            i = iorb+(ilat-1)*Norb
                            j = jorb+(jlat-1)*Norb
                            Hbath_basis(counter)%O(1,1,ispin,jspin,i,j)=one
                            Hbath_basis(counter)%O(1,1,ispin,jspin,j,i)=one
                            Hbath_lambda(:,counter)=dreal(Hloc(io,jo))
                         endif
                         !
                         if(dimag(Hloc(io,jo))/=0d0)then
                            counter=counter+1
                            Hbath_basis(counter)%O(1,1,ispin,jspin,i,j)=xi
                            Hbath_basis(counter)%O(1,1,ispin,jspin,j,i)=xi
                            Hbath_lambda(:,counter)=dimag(Hloc(io,jo))
                         endif
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine init_Hbath_direct_lso2




  function index_stride_lso(ilat,ispin,iorb) result(io)
    integer :: ilat,ispin,iorb
    integer :: io
    io = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
  end function index_stride_lso


END MODULE ED_BATH_REPLICA









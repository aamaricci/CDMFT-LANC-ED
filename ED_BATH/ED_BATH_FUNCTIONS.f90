MODULE ED_BATH_FUNCTIONS
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,str
  USE SF_LINALG, only: eye,inv,zeye,inv_her
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  implicit none

  private

  !\DELTA, THE HYBRIDIZATION FUNCTION
  interface delta_bath_function
     !
     ! Evaluates the normal hybridization function :math:`\Delta(x)`.
     !
     ! Output:
     !   * :f:var:`delta` : complex rank-5 array with dimension [ |Nlat| , |Nlat| , |Nns| , |Nns| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure delta_bath_array
  end interface delta_bath_function


  !NON-INTERACTING GREEN'S FUNCTION 
  interface g0and_bath_function
     !
     ! Evaluates the normal non-interacting Green's function :math:`G_0(x)`.
     !
     ! Output:
     !   * :f:var:`g0and` : complex rank-5 array with dimension [ |Nlat| , |Nlat| , |Nns| , |Nns| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure g0and_bath_array
  end interface g0and_bath_function

  !INVERSE NON-INTERACTING GREEN'S FUNCTION 
  interface invg0_bath_function
     !
     ! Evaluates the inverse of the normal non-interacting Green's function :math:`G^{-1}_0(x)`.
     !
     ! Output:
     !   * :f:var:`g0and` : complex rank-5 array with dimension [ |Nlat| , |Nlat| , |Nns| , |Nns| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure invg0_bath_array
  end interface invg0_bath_function




  public :: delta_bath_function
  public :: g0and_bath_function
  public :: invg0_bath_function





  integer :: ibath
  integer :: inam,jnam
  integer :: ilat,jlat
  integer :: iorb,jorb
  integer :: ispin,jspin
  integer :: is,js
  integer :: io,jo
  integer :: i,j
  integer :: L,Ntot

contains



  function delta_bath_array(x,dmft_bath_,axis) result(Delta)
    complex(8),dimension(:),intent(in)                              :: x          !complex  array for the frequency
    type(effective_bath)                                            :: dmft_bath_ !the current :f:var:`effective_bath` instance
    character(len=*),optional                                       :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(x)) :: Delta
    ! complex(8),dimension(Nambu*Nspin*Nimp,Nambu*Nspin*Nimp)         :: JJ
    complex(8),dimension(Nambu*Nspin*Nimp,Nambu*Nspin*Nimp)         :: Vk
    complex(8),dimension(Nambu*Nspin*Nimp,Nambu*Nspin*Nimp)         :: Hk,invHk
    character(len=4)                                                :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    !
    Delta=zero
    !
    L    = size(x)
    Ntot = Nambu*Nspin*Nimp
    !
    select case(bath_type)
    case default
       if(Nambu/=1)stop "delta_bath_array ERROR: Nambu != 1 with ed_mode=normal"
       do ibath=1,Nbath
          Vk = one*diag(dmft_bath%item(ibath)%v(:))
          Hk = nnn2nso_reshape(Hbath_build(dmft_Bath%item(ibath)%lambda)) !rank6 -> rank2
          do i=1,L
             invHk = zeye(Ntot)*x(i) - Hk
             call inv(invHk)
             invHk = matmul(matmul(Vk,invHk),Vk)
             Delta(:,:,:,:,:,:,i)=Delta(:,:,:,:,:,:,i) + nso2nnn_reshape(invHk)
          enddo
       enddo
    case ("superc")             !treat all the 2x2 Nambu components at once
       if(Nambu/=2)stop "delta_bath_array ERROR: Nambu != 2 with ed_mode=superc"
       ! select case(axis_)
       ! case default ; JJ=zeye(Ntot)
       ! case ('real'); JJ=kron(pauli_sigma_z,zeye(Nspin*Nimp))
       ! end select
       do ibath=1,Nbath
          Vk = kron(pauli_sigma_z,one*diag(dmft_bath_%item(ibath)%v(:)))
          Hk = nnn2nso_reshape(Hbath_build(dmft_Bath%item(ibath)%lambda)) !rank6 -> rank2          
          do i=1,L
             invHk   = diag(zeta_superc(x,0d0,axis_)) - Hk !x(i)*JJ - Hk
             call inv(invHk)
             invHk   = matmul(matmul(Vk,invHk),Vk)
             Delta(:,:,:,:,:,:,i) = Delta(:,:,:,:,:,:,i) + nso2nnn_reshape(invHk)
          enddo
       enddo
    end select
    !>ENFORCE SYMMETRY HERE?
    !D(2,1,:,:,:,:) = conjg(D(1,2,:,:,:,:)) (or dagger?)
    !D(2,2,:,:,:,:) = -conjg(D(1,1,:,:,:,:/L:1:-1))
  end function delta_bath_array




  !G0 = (z+mu)1 - Hloc - Delta
  function invg0_bath_array(x,dmft_bath_,axis) result(invG0and)
    complex(8),dimension(:),intent(in)                              :: x !complex  array for the frequency
    type(effective_bath)                                            :: dmft_bath_ !the current :f:var:`effective_bath` instance
    character(len=*),optional                                       :: axis !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(x)) :: invG0and,Delta
    character(len=4)                                                :: axis_    
    complex(8),dimension(Nnambu*Nspin*Nimp,Nnambu*Nspin*Nimp)       :: iG0
    complex(8),dimension(Nambu*Nspin*Nimp,size(x))                  :: zeta
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    L = size(x)
    Ntot = Nambu*Nspin*Nimp
    !
    !> Get Delta:
    Delta = delta_bath_array(x,dmft_bath_,axis_)
    !
    select case(bath_type)
    case default
       if(Nambu/=1)stop "invG0_bath_array ERROR: Nambu != 1 with ed_mode=normal"
       do i=1,L
          invG0and(1,1,:,:,:,:,i) = &
               so2nn_reshape((x(i)+xmu)*zeye(Nspin*Nimp),Nspin,Nimp) - impHloc - Delta(1,1,:,:,:,:,i)
       enddo
    case ("superc")             !treat all the 2x2 Nambu components at once
       if(Nambu/=2)stop "delta_bath_array ERROR: Nambu != 2 with ed_mode=superc"
       zeta = zeta_superc(x,mu,axis_)
       do i=1,L
          invG0and(:,:,:,:,:,:,i) = nso2nnn_reshape(diag(zeta(:,i))) -  hloc_superc(impHloc) - Delta(:,:,:,:,:,:,i)
       enddo
    end select
    !>ENFORCE SYMMETRY HERE?
    !D(2,1,:,:,:,:) = conjg(D(1,2,:,:,:,:)) (or dagger?)
    !D(2,2,:,:,:,:) = -conjg(D(1,1,:,:,:,:/L:1:-1))
  end function invg0_bath_array





  function g0and_bath_array(x,dmft_bath_,axis) result(G0and)
    complex(8),dimension(:),intent(in)                              :: x !complex  array for the frequency
    type(effective_bath)                                            :: dmft_bath_ !the current :f:var:`effective_bath` instance
    character(len=*),optional                                       :: axis !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
    character(len=4)                                                :: axis_    
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(x)) :: G0and
    complex(8),dimension(Nnambu*Nspin*Nimp,Nnambu*Nspin*Nimp)       :: iG0
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    L     = size(x)
    Ntot  = Nambu*Nspin*Nimp
    !
    G0and = invg0_bath(x,dmft_bath_,axis_)
    !
    do i=1,L
       iG0=nnn2nso_reshape(G0and(:,:,:,:,:,:,i))
       call inv(iG0)
       G0and(:,:,:,:,:,:,i)=nso2nnn_reshape(iG0)
    enddo
    !>ENFORCE SYMMETRY HERE?
    !D(2,1,:,:,:,:) = conjg(D(1,2,:,:,:,:)) (or dagger?)
    !D(2,2,:,:,:,:) = -conjg(D(1,1,:,:,:,:/L:1:-1))
  end function g0and_bath_array













  function zeta_superc(x,mu,axis) result(zeta)
    complex(8),dimension(:)                        :: x
    real(8)                                        :: mu
    character(len=*)                               :: axis
    complex(8),dimension(Nambu*Nspin*Nimp,size(x)) :: zeta
    integer                                        :: N,L
    N = Nspin*Nimp
    L = size(x)
    do concurrent(iorb=1:N)
       select case(axis_)
       case default
          zeta(iorb ,1:L)   = x(1:L) + mu
          zeta(N+iorb,1:L)  = x(1:L) - mu
       case ('real')
          zeta(iorb ,1:L)   = x(1:L) + mu
          zeta(N+iorb:,1:L) = -conjg(x(L:1:-1) + mu)
       end select
    enddo
  end function zeta_superc


  function hloc_superc(h) result(zh)
    complex(8),dimension(Nspin,Nspin,Nimp,Nimp)             :: h
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp) :: zh
    if(Nambu/=2)stop "hloc_superc ERROR: called with Nambu!=2"
    zh              = zero 
    zh(1,1,:,:,:,:) = h
    zh(2,2,:,:,:,:) = -conjg(h)
  end function hloc_superc


END MODULE ED_BATH_FUNCTIONS

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

  interface Fdelta_bath_function
     !
     ! Evaluates the anomalouse hybridization function :math:`\Theta(x)`.
     !
     ! Output:
     !   * :f:var:`fdelta` : complex rank-5 array with dimension [ |Nlat| , |Nlat| , |Nns| , |Nns| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure Fdelta_bath_array
  end interface Fdelta_bath_function




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

  interface f0and_bath_function
     !
     ! Evaluates the anomalous non-interacting Green's function :math:`F_0(x)`.
     !
     ! Output:
     !   * :f:var:`f0and` : complex rank-5 array with dimension [ |Nlat| , |Nlat| , |Nns| , |Nns| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure f0and_bath_array
  end interface f0and_bath_function



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

  interface invf0_bath_function
     !
     ! Evaluates the inverse of the anomalous non-interacting Green's function :math:`F^{-1}_0(x)`.
     !
     ! Output:
     !   * :f:var:`f0and` : complex rank-5 array with dimension [ |Nlat| , |Nlat| , |Nns| , |Nns| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure invf0_bath_array
  end interface invf0_bath_function


  public :: delta_bath_function,fdelta_bath_function
  public :: g0and_bath_function,f0and_bath_function
  public :: invg0_bath_function,invf0_bath_function



contains



  function delta_bath_array(x,dmft_bath_,axis) result(Delta)
    complex(8),dimension(:),intent(in)                            :: x          !complex  array for the frequency
    type(effective_bath)                                          :: dmft_bath_ !the current :f:var:`effective_bath` instance
    character(len=*),optional                                     :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    integer                                                       :: i,ih,L
    integer                                                       :: ilat,jlat,iorb,jorb,ispin,jspin,ibath
    integer                                                       :: io,jo
    !
    complex(8),dimension(Nlat*Nns*Norb,Nlat*Nns*Norb,size(x))     :: zeta
    complex(8),dimension(Nlat*Nns*Norb,Nlat*Nns*Norb)             :: Vk
    complex(8),dimension(Nlat*Nns*Norb,Nlat*Nns*Norb)             :: invH_k
    complex(8),dimension(Nlat,Nlat,Nns,Nns,Norb,Norb)             :: invH_knn
    character(len=4)                                              :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    !
    Delta=zero
    !
    L    = size(x)
    !
    select case(bath_type)
    case default
       invH_k=zero
       do i=1,L
          do ibath=1,Nbath
             Vk       = dzdiag(dmft_bath%item(ibath)%v(:))
             invH_knn = Hbath_build(dmft_Bath%item(ibath)%lambda)
             invH_k   = nnn2lso_reshape(invH_knn,Nlat,Nspin,Norb)
             invH_k   = zeye(Nlat*Nspin*Norb)*x(i) - invH_k
             call inv(invH_k)
             invH_k   = matmul(matmul(Vk,invH_k),Vk)
             invH_knn = lso2nnn_reshape(invH_k,Nlat,Nspin,Norb)
             Delta(:,:,:,:,:,:,i)=Delta(:,:,:,:,:,:,i) + invH_knn
          enddo
       enddo
    case ("superc")
       JJ=kron(pauli_sigma_z,zeye(Nlat*Norb))
       do i=1,L
          select case(axis_)
          case default ; zeta(:,:,i) = x(i)*zeye(Nlat*Nns*Norb)
          case ('real'); zeta(:,:,i) = x(i)*JJ
          end select
          do ibath=1,Nbath
             Vk       = kron(pauli_sigma_z,dzdiag(dmft_bath_%item(ibath)%v(:)))
             invH_knn = Hbath_build(dmft_bath_%item(ibath)%lambda)
             invH_k   = nnn2lso_reshape(invH_knn,Nlat,Nns,Norb)
             invH_k   = zeta(:,:,i) - invH_k
             call inv(invH_k)
             invH_k   = matmul(matmul(Vk,invH_k),Vk)
             invH_knn = lso2nnn_reshape(invH_k,Nlat,Nns,Norb)
             Delta(:,:,:,:,:,:,i)=Delta(:,:,:,:,:,:,i) + invH_knn(:,:,1,1,:,:) !get component normal 11
          enddo
       enddo
    end select
    !
  end function delta_bath_array



  function fdelta_bath_array(x,dmft_bath_,axis) result(Fdelta)
    complex(8),dimension(:),intent(in)                                          :: x          !complex  array for the frequency
    type(effective_bath)                                                        :: dmft_bath_ !the current :f:var:`effective_bath` instance
    character(len=*),optional                                                   :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x))               :: Fdelta
    integer                                                                     :: i,ih,L
    integer                                                                     :: ilat,jlat,iorb,jorb,ispin,jspin,ibath
    integer                                                                     :: io,jo
    !
    complex(8),dimension(Nlat*Nns*Norb,Nlat*Nns*Norb,size(x)) :: zeta
    complex(8),dimension(Nlat*Nns*Norb,Nlat*Nns*Norb)         :: Vk
    complex(8),dimension(Nlat*Nns*Norb,Nlat*Nns*Norb)         :: invH_k
    complex(8),dimension(Nlat,Nlat,Nns,Nns,Norb,Norb)         :: invH_knn
    character(len=4)                                                            :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    !
    Delta=zero
    !
    L    = size(x)
    !
    select case(ed_mode)
    case default
       stop "Fdelta_bath_array ERROR: called with ed_mode=normal."
    case ("superc")
       JJ=kron(pauli_sigma_z,zeye(Nlat*Norb))
       do i=1,L
          select case(axis_)
          case default ; zeta(:,:,i) = x(i)*zeye(Nlat*Nns*Norb)
          case ('real'); zeta(:,:,i) = x(i)*JJ
          end select
          do ibath=1,Nbath
             Vk       = kron(pauli_sigma_z,dzdiag(dmft_bath_%item(ibath)%v(:)))
             invH_knn = Hbath_build(dmft_bath_%item(ibath)%lambda)
             invH_k   = nnn2lso_reshape(invH_knn,Nlat,Nns,Norb)
             invH_k   = zeta(:,:,i) - invH_k
             call inv(invH_k)
             invH_k   = matmul(matmul(Vk,invH_k),Vk)
             invH_knn = lso2nnn_reshape(invH_k,Nlat,Nns,Norb)
             Fdelta(:,:,:,:,:,:,i)=Fdelta(:,:,:,:,:,:,i) + invH_knn(:,:,1,2,:,:) !get the 12 anomalous component
          enddo
       enddo
    end select
    !
  end function fdelta_bath_array







  function g0and_bath(g0and_bath_array(x,dmft_bath_,axis) result(G0and)
    complex(8),dimension(:),intent(in)                            :: x !complex  array for the frequency
    type(effective_bath)                                          :: dmft_bath_ !the current :f:var:`effective_bath` instance
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta
    integer                                                       :: ilat,jlat,iorb,jorb,ispin,jspin,io,jo,Nso,i,L
    real(8),dimension(size(x))                                    :: det
    complex(8),dimension(:,:),allocatable                         :: fg,zeta
    character(len=*),optional                                     :: axis !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
    character(len=4)                                              :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    !
    L=size(x)
    !
    G0and = zero
    !
    select case(ed_mode)
    case default
       allocate(fg(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
       G0and = invg0_bath(x,dmft_bath_,axis_)
       do i=1,L
          fg=nnn2lso_reshape(G0and(:,:,:,:,:,:,i),Nlat,Nspin,Norb)
          call inv(fg)
          G0and(:,:,:,:,:,:,i)=lso2nnn_reshape(fg,Nlat,Nspin,Norb)
       enddo
       deallocate(fg)
       !
    case ("superc")
       allocate(fgorb(Nlat*Nns*Norb,Nlat*Nns*Norb),zeta(Nlat*Nns*Norb,Nlat*Nns*Norb))
       Delta =  delta_bath_array(x,dmft_bath_,axis_)
       Fdelta= fdelta_bath_array(x,dmft_bath_,axis_)

       do ispin=1,Nspin   !==1
          do i=1,L
             Nso = Nlat*Norb
             zeta = zero
             fgorb= zero
             select case(axis_)
             case default
                do concurrent(ilat=1:Nlat,iorb=1:Norb)
                   io = iorb + (ilat-1)*Norb
                   zeta(io,io)           = x(i) + xmu
                   zeta(io+Nso,iorb+Nso) = x(i) - xmu
                enddo
                do concurrent(ilat=1:Nlat,jlat=1:Nlat,iorb=1:Norb,jorb=1:Norb)
                   io = iorb + (ilat-1)*Norb
                   jo = jorb + (jlat-1)*Norb
                   fgorb(io,jo)         = zeta(io,jo)         - impHloc(ilat,jlat,ispin,ispin,iorb,jorb)        - Delta(ilat,jlat,ispin,ispin,iorb,jorb,i)
                   fgorb(io,jo+Norb)    = zeta(io,jo+Nso)                                                       - Fdelta(ilat,jlat,ispin,ispin,iorb,jorb,i)
                   fgorb(io+Nso,jo)     = zeta(io+Nso,jo)                                                       - conjg(Fdelta(ilat,jlat,ispin,ispin,iorb,jorb,i))
                   fgorb(io+Nso,jo+Nso) = zeta(io+Nso,jo+Nso) + conjg(impHloc(ilat,jlat,ispin,ispin,iorb,jorb)) + conjg( Delta(ilat,jlat,ispin,ispin,iorb,jorb,i) )
                enddo
             case("real")
                do iorb=1,Norb
                   zeta(iorb,iorb)           =        x(i)     + xmu
                   zeta(iorb+Norb,iorb+Norb) = -conjg(x(L-i+1) + xmu) !as above this is == -x(i)-mu == -w-i*eta - mu
                enddo
                do iorb=1,Norb
                   do jorb=1,Norb
                      fgorb(iorb,jorb)           = zeta(iorb,jorb)           - impHloc(ilat,jlat,ispin,ispin,iorb,jorb)         - Delta(ilat,jlat,ispin,ispin,iorb,jorb,i)
                      fgorb(iorb,jorb+Norb)      = zeta(iorb,jorb+Norb)                                                         - Fdelta(ilat,jlat,ispin,ispin,iorb,jorb,i)
                      fgorb(iorb+Norb,jorb)      = zeta(iorb+Norb,jorb)                                                         - conjg(Fdelta(ilat,jlat,ispin,ispin,iorb,jorb,i))
                      fgorb(iorb+Norb,jorb+Norb) = zeta(iorb+Norb,jorb+Norb) + conjg(impHloc(ilat,jlat,ispin,ispin,iorb,jorb))  + conjg( Delta(ilat,jlat,ispin,ispin,iorb,jorb,L-i+1) )
                   enddo
                enddo
             end select
             !
             call inv(fgorb)
             G0and(ispin,ispin,:,:,i) = fgorb(1:Norb,1:Norb)
             !
             deallocate(fgorb,zeta)
          enddo
       enddo
    end select





    !
  end function g0and_bath







  function invg0_bath_array(x,dmft_bath_,axis) result(G0and)
    complex(8),dimension(:),intent(in)                                          :: x !complex  array for the frequency
    type(effective_bath)                                                        :: dmft_bath_ !the current :f:var:`effective_bath` instance
    complex(8),dimension(Nlat,Nlat,Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,size(x)) :: G0and,Delta
    integer                                                                     :: i,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,Nso,L
    complex(8),dimension(Nlat*Nnambu*Nspin*Norb,Nlat*Nnambu*Nspin*Norb)         :: zeta
    character(len=*),optional                                                   :: axis !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
    character(len=4)                                                            :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    L = size(x)
    !
    Delta = delta_bath_array(x,dmft_bath_,axis_)
    G0and = zero
    do i=1,L
       zeta = (x(i)+xmu)*zeye(Nlat*Nnambu*Nspin*Norb)
       G0and(:,:,:,:,:,:,i) = lso2nnn_reshape(zeta,Nlat,Nnambu*Nspin,Norb)-impHloc-Delta(:,:,:,:,:,:,i)
    enddo
    !
  end function invg0_bath_array




  function invf0_bath_array(x,dmft_bath_,axis) result(F0and)
    complex(8),dimension(:),intent(in)                                          :: x !complex  array for the frequency
    type(effective_bath)                                                        :: dmft_bath_ !the current :f:var:`effective_bath` instance
    complex(8),dimension(Nlat,Nlat,Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,size(x)) :: F0and,Fdelta
    integer                                                                     :: i,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,Nso,L
    character(len=*),optional                                                   :: axis !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
    character(len=4)                                                            :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    L = size(x)
    !
    Fdelta = fdelta_bath_array(x,dmft_bath_,axis_)
    F0and  = zero
    do i=1,L
       F0and(:,:,:,:,:,:,i) = -Fdelta(:,:,:,:,:,:,i)
    enddo
    !
  end function invf0_bath_array





  ! Auxiliary ℝ -> ℂ vector-to-diagonal-matrix constructor
  function dzdiag(x) result(A)
    real(8),dimension(:)                   :: x
    complex(8),dimension(:,:),allocatable  :: A
    integer                                :: N,i
    N=size(x,1)
    allocate(A(N,N))
    A=0.d0
    do i=1,N
       A(i,i)=x(i)
    enddo
  end function dzdiag

END MODULE ED_BATH_FUNCTIONS

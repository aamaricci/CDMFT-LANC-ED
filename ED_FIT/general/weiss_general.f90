!######################################################################
! WEISS : NORMAL
!######################################################################
function g0and_general_normal(a) result(G0and)
  real(8),dimension(:)                               :: a
  complex(8),dimension(Nspin,Nspin,Nimp,Nimp,Ldelta) :: G0and,Delta
  complex(8),dimension(Nimp*Nspin,Nimp*Nspin)        :: fgorb
  integer                                            :: i
  complex(8)                                         :: z
  !
  Delta = delta_general_normal(a)
  !
  do i=1,Ldelta
     Z     = xi*Xdelta(i)+xmu
     FGorb = Z*zeye(Nimp*Nspin) - nn2so_reshape(impHloc + Delta(:,:,:,:,i),Nspin,Nimp)
     call inv(FGorb)
     G0and(:,:,:,:,i) = so2nn_reshape(FGorb,Nspin,Nimp)
  enddo
  !
end function g0and_general_normal




!######################################################################
! WEISS : SUPERC
!######################################################################
function g0and_general_superc(a) result(G0and)
  real(8),dimension(:)                                           :: a
  complex(8),dimension(2,Nspin,Nspin,Nimp,Nimp,Ldelta)           :: G0and
  complex(8),dimension(Nambu*Nspin*Nimp,Nambu*Nspin*Nimp)        :: FG,Haux,Htmp,Z,Zloc,Vaux
  complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp)        :: g0
  complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,Ldelta) :: delta
  integer                                                        :: i,stride
  real(8),dimension(Nbath,Nspin*Nimp)                            :: Vk
  real(8),dimension(Nbath,Nlambdas)                              :: Lk
  !!Get Hs
  stride = 0
  do ibath=1,Nbath
     !Get Vs
     Vk(ibath,:)= a(stride+1:stride+Nspin*Nimp)
     stride     = stride + Nspin*Norb
     !Get Lambdas
     Lk(ibath,:)= a(stride+1:stride+Nlambdas)
     stride     = stride + Nlambdas
  enddo
  !
  !> get all Delta components:
  Delta=zero
  do ibath=1,Nbath
     Vaux = kron(pauli_sigma_z,one*diag(Vk(ibath,:)))
     Htmp = nnn2nso_reshape(Hbath_build(Lk(ibath,:)))
     do i=1,Ldelta
        Z  = xi*Xdelta(i)*zeye(Nambu*Nspin*Nimp)
        Haux = Z - Htmp
        call inv(Haux)
        Delta(:,:,:,:,:,:,i)=Delta(:,:,:,:,:,:,i) + nso2nnn_reshape( matmul(matmul(Vaux,Haux),Vaux) )
     enddo
  enddo
  !
  Zloc = xmu*kron(pauli_tau_z,zeye(Nspin*Nimp)) - &
       kron(pauli_tau_z,nn2so_reshape(impHloc,Nspin,Nimp))
  do i=1,Ldelta
     Z    = xi*Xdelta(i)*zeye(Nambu*Nspin*Nimp) + Zloc
     FG = Z - nnn2nso_reshape(Delta(:,:,:,:,:,:,i))
     call inv(FG)
     g0 = nso2nnn_reshape(FG)
     G0and(1,:,:,:,:,i) = g0(1,1,:,:,:,:)
     G0and(2,:,:,:,:,i) = g0(1,2,:,:,:,:)
  enddo
  !
end function g0and_general_superc






!######################################################################
! WEISS GRAD : NORMAL
!######################################################################
function grad_g0and_general_normal(a) result(dG0and)
  real(8),dimension(:)                                       :: a
  complex(8),dimension(Nspin,Nspin,Nimp,Nimp,Ldelta,size(a)) :: dG0and,dDelta
  complex(8),dimension(Nspin,Nspin,Nimp,Nimp,Ldelta)         :: G0and
  complex(8),dimension(Nspin*Nimp,Nspin*Nimp)                :: dDelta_lso,dG0and_lso,G0and_lso
  integer                                                    :: ispin,iorb,jorb
  integer                                                    :: ik,l
  !
  G0and  = g0and_general_normal(a)
  dDelta = grad_delta_general_normal(a)
  !
  dG0and = zero
  !
  do l=1,Ldelta
     G0and_lso=nn2so_reshape(g0and(:,:,:,:,l),Nspin,Nimp)
     do ik=1,size(a)
        dDelta_lso=nn2so_reshape(dDelta(:,:,:,:,l,ik),Nspin,Nimp)
        dG0and_lso = (G0and_lso .x. dDelta_lso) .x. G0and_lso
        dG0and(:,:,:,:,l,ik)=so2nn_reshape(dG0and_lso,Nspin,Nimp)
     enddo
  enddo
  !
end function grad_g0and_general_normal






!######################################################################
! WEISS GRAD : SUPERC ??
!######################################################################

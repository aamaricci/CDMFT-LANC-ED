!######################################################################
! DELTA : NORMAL
!######################################################################
function delta_replica_normal(a) result(Delta)
  real(8),dimension(:)                                    :: a
  complex(8),dimension(Nspin,Nspin,Nimp,Nimp,Ldelta)      :: Delta
  integer                                                 :: i,ibath,stride
  complex(8),dimension(Nspin*Nimp,Nspin*Nimp)             :: Haux,Htmp
  real(8),dimension(Nbath)                                :: Vk
  real(8),dimension(Nbath,Nlambdas)                       :: Lk
  complex(8)                                              :: iw
  !
  !Get Hs
  stride = 1
  do ibath=1,Nbath
     !Get Vs
     Vk(ibath)  = a(stride)
     stride     = stride + 1
     !Get Lambdas
     Lk(ibath,:)= a(stride+1:stride+Nlambdas)
     stride     = stride + Nlambdas
  enddo
  !
  !
  Delta=zero
  do ibath=1,Nbath
     Htmp   = nnn2nso_reshape(Hbath_build(Lk(ibath,:)))
     do i=1,Ldelta
        iw   = xi*Xdelta(i)
        Haux = zeye(Nspin*Nimp)*iw - Htmp
        call inv(Haux)
        Delta(:,:,:,:,i)=Delta(:,:,:,:,i)+ Vk(ibath)*so2nn_reshape(Haux,Nspin,Nimp)*Vk(ibath)
     enddo
  enddo
end function delta_replica_normal






!######################################################################
! DELTA  : SUPERC
!######################################################################
function delta_replica_superc(a) result(Delta)
  real(8),dimension(:)                                    :: a
  complex(8),dimension(2,Nspin,Nspin,Nimp,Nimp,Ldelta)    :: Delta
  integer                                                 :: ibath
  integer                                                 :: i,io,jo,ndx,stride
  complex(8),dimension(Nspin*Nimp,Nspin*Nimp)             :: Haux,Htmp
  complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp) :: invHk
  real(8),dimension(Nbath)                                :: Vk
  real(8),dimension(Nbath,Nlambdas)                       :: Lk
  complex(8)                                              :: iw
  !
  stride = 1
  do ibath=1,Nbath
     !Get Vs
     Vk(ibath)   = a(stride)
     stride      = stride + 1
     !Get Lambdas
     Lk(ibath,:) = a(stride+1:stride+Nlambdas)
     stride      = stride + Nlambdas
  enddo
  !
  !
  Delta=zero
  do ibath=1,Nbath
     Htmp = nnn2nso_reshape(Hbath_build(Lk(ibath,:)))
     do i=1,Ldelta
        iw   = xi*Xdelta(i)
        Haux = zeye(Nambu*Nspin*Nimp)*iw - Htmp
        call inv(Haux)
        invHk= nso2nnn_reshape(Haux)
        Delta(1,:,:,:,:,i)=Delta(1,:,:,:,:,i)+ Vk(ibath)*invHk(1,1,:,:,:,:)*Vk(ibath)
        Delta(2,:,:,:,:,i)=Delta(2,:,:,:,:,i)+ Vk(ibath)*invHk(1,2,:,:,:,:)*Vk(ibath)
     enddo
  enddo
end function delta_replica_superc






!######################################################################
! DELTA GRAD: NORMAL
!######################################################################
function grad_delta_replica_normal(a) result(dDelta)
  real(8),dimension(:)                                       :: a
  complex(8),dimension(Nspin,Nspin,Nimp,Nimp,Ldelta,size(a)) :: dDelta
  integer                                                    :: ibath
  integer                                                    :: k,l,stride
  complex(8),dimension(Nspin*Nimp,Nspin*Nimp)                :: H_reconstructed, Htmp,O_k
  complex(8),dimension(Nspin*Nimp,Nspin*Nimp,Ldelta)         :: Haux
  real(8),dimension(Nbath)                                   :: Vk !TODO) to extend to Vup != Vdw: 1->NSPIN
  real(8),dimension(Nbath,Nlambdas)                          :: Lk
  !
  dDelta=zero
  !
  !Get Hs
  stride = 1
  do ibath=1,Nbath
     Vk(ibath)  = a(stride)
     stride     = stride + 1
     Lk(ibath,:)= a(stride+1:stride+Nlambdas)
     stride     = stride + Nlambdas
  enddo
  !
  stride = 1
  do ibath=1,Nbath
     !Derivate_Vp
     H_reconstructed = nnn2nso_reshape(Hbath_build(Lk(ibath,:))) !Nambu=1
     do l=1,Ldelta
        Haux(:,:,l) = zeye(Nspin*Nimp)*(xi*Xdelta(l)) - H_reconstructed
        call inv(Haux(:,:,l))
        ! do ispin=1,Nspin
        dDelta(:,:,:,:,l,stride) = 2d0*Vk(ibath)*so2nn_reshape(Haux(:,:,l),Nspin,Nimp)
        ! enddo
     enddo
     stride=stride+1
     !
     !Derivate_lambda_p
     do k=1,Nlambdas
        O_k=nnn2nso_reshape(Hbath_basis(k)%O) !-> [Nso,Nso]
        do l=1,Ldelta
           Htmp = (Haux(:,:,l) .x. O_k) .x. Haux(:,:,l)
           dDelta(:,:,:,:,l,stride)=Vk(ibath)*so2nn_reshape(Htmp,Nspin,Nimp)*Vk(ibath)
        enddo
        stride = stride + 1
     enddo
  enddo
  !
end function grad_delta_replica_normal





!######################################################################
! DELTA GRAD: SUPERC ??
!######################################################################



!##################################################################
!##################################################################
!##################################################################




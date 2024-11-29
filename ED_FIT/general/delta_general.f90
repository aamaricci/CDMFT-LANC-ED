function delta_general_normal(a) result(Delta)
  real(8),dimension(:)                                    :: a
  complex(8),dimension(Nspin,Nspin,Nimp,Nimp,Ldelta)      :: Delta
  integer                                                 :: i,ibath
  integer                                                 :: i,stride
  complex(8),dimension(Nspin*Nimp,Nspin*Nimp)             :: Haux,Vaux
  real(8),dimension(Nbath,Nspin*Nimp)                     :: Vk
  real(8),dimension(Nbath,Nlambdas)                       :: Lk
  complex(8)                                              :: iw
  !
  !Get Hs
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
  !
  Delta=zero
  do i=1,Ldelta
     iw = xi*Xdelta(i)
     do ibath=1,Nbath
        Vaux = one*diag(Vk(ibath,:))
        Haux = zeye(Nspin*Nimp)*iw - nnn2nso_reshape(Hbath_build(Lk(ibath,:))) !Nambu=1
        call inv(Haux)
        Haux = matmul(Vaux,matmul(Haux,Vaux))
        Delta(:,:,:,:,i)=Delta(:,:,:,:,i)+ so2nn_reshape(Haux)
     enddo
  enddo
  !
end function delta_general_normal




function grad_delta_general_normal(a) result(dDelta)
  real(8),dimension(:)                                       :: a
  complex(8),dimension(Nspin,Nspin,Nimp,Nimp,Ldelta,size(a)) :: dDelta
  integer                                                    :: ibath
  integer                                                    :: k,l,stride
  complex(8),dimension(Nspin*Nimp,Nspin*Nimp)                :: H_reconstructed, Htmp,O_k,Vaux
  complex(8),dimension(Nspin*Nimp,Nspin*Nimp,Ldelta)         :: Haux
  real(8),dimension(Nbath,Nspin*Nimp)                        :: Vk
  real(8),dimension(Nbath,Nlambdas)                          :: Lk
  !
  dDelta=zero
  !
  !Get Hs
  stride = 0
  do ibath=1,Nbath
     Vk(ibath,:)= a(stride+1:stride+Nspin*Nimp)
     stride     = stride + Nspin*Nimp
     Lk(ibath,:)= a(stride+1:stride+Nlambdas)
     stride     = stride + Nlambdas
  enddo
  !
  stride = 0
  do ibath=1,Nbath
     !Derivate_Vp
     Vaux = one*diag(Vk(ibath,:))      
     H_reconstructed = nnn2nso_reshape(Hbath_build(Lk(ibath,:))) !Nambu=1
     do l=1,Ldelta
        Haux(:,:,l) = zeye(Nspin*Nimp)*(xi*Xdelta(l)) - H_reconstructed
        call inv(Haux(:,:,l))
        do k = 1,Nlat*Nspin*Norb
           Htmp      = zero
           Htmp(:,k) = Htmp(:,k) + VH_prod(V_k,Haux(:,:,l),k)
           Htmp(k,:) = Htmp(k,:) + HV_prod(V_k,Haux(:,:,l),k)
           dDelta(:,:,:,:,l,stride+k) = so2nn_reshape(Htmp,Nspin,Nimp)
        enddo
     enddo
     !
     !Derivate_lambda_p
     do k=1,Nlambdas
        stride = stride + 1
        O_k=nnn2nso_reshape(Hbath_basis(k)%O) !-> [Nso,Nso]
        do l=1,Ldelta
           Htmp = (Haux(:,:,l) .x. O_k) .x. Haux(:,:,l)
           Htmp = (V_k .x. Htmp) .x. V_k
           dDelta(:,:,:,:,l,stride)=so2nn_reshape(Htmp,Nspin,Nimp)
        enddo
     enddo
  enddo
  !
contains
  !Auxiliary functions for general derivative in V_k
  function HV_prod(vv,HH,ig) result(AA)
    complex(8),dimension(size(vv,1),size(vv,1))             :: HH,vv
    complex(8),dimension(size(vv,1),size(vv,1)),allocatable :: Atmp
    complex(8),dimension(size(vv,1)),allocatable            :: AA
    integer                                                 :: ig
    Atmp=matmul(HH,vv)
    AA  =Atmp(ig,:)
  end function HV_prod
  !
  function VH_prod(vv,HH,ig) result(AA)
    complex(8),dimension(size(vv,1),size(vv,1))             :: HH,vv
    complex(8),dimension(size(vv,1),size(vv,1)),allocatable :: Atmp
    complex(8),dimension(size(vv,1)),allocatable            :: AA
    integer                               :: N,i,ig
    Atmp=matmul(vv,HH)
    AA  =Atmp(:,ig)
  end function VH_prod

end function grad_delta_general_normal




!##################################################################
!##################################################################
!##################################################################





function delta_general_superc(a) result(Delta)
  real(8),dimension(:)                                    :: a
  complex(8),dimension(2,Nspin,Nspin,Nimp,Nimp,Ldelta)    :: Delta
  integer                                                 :: ibath
  integer                                                 :: i,io,jo,ndx,stride
  complex(8),dimension(Nambu*Nspin*Nimp,Nambu*Nspin*Nimp) :: Haux,Vaux,Z
  complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp) :: invHk
  real(8),dimension(Nbath,Nspin*Nimp)                     :: Vk
  real(8),dimension(Nbath,Nlambdas)                       :: Lk
  complex(8)                                              :: iw
  !
  stride = 1
  do ibath=1,Nbath
     !Get Vs
     Vk(ibath,:)= a(stride+1:stride+Nspin*Nimp)
     stride     = stride + Nspin*Norb
     !Get Lambdas
     Lk(ibath,:)= a(stride+1:stride+Nlambdas)
     stride     = stride + Nlambdas
  enddo
  !
  !
  Delta=zero
  do i=1,Ldelta
     Z  = xi*Xdelta(i)*zeye(Nambu*Nspin*Nimp)
     do ibath=1,Nbath
        Vaux = kron(pauli_sigma_z,one*diag(Vk(ibath,:)))
        Haux = Z - nnn2nso_reshape(Hbath_build(Lk(ibath,:)))
        call inv(Haux)
        invHk  = nso2nnn_reshape( matmul(matmul(Vaux,Haux),Vaux) )
        Delta(1,:,:,:,:,i)=Delta(1,:,:,:,:,i)+ invHk(1,1,:,:,:,:)
        Delta(2,:,:,:,:,i)=Delta(2,:,:,:,:,i)+ invHk(1,2,:,:,:,:)
     enddo
  enddo
  !
end function delta_general_superc

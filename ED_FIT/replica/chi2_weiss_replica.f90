!######################################################################
! WEISS CHI2 FIT : NORMAL
!######################################################################
function chi2_weiss_replica_normal(a) result(chi2)
  real(8),dimension(:) :: a
  real(8)              :: chi2
  select case(cg_norm)
  case default;stop "chi2_fitgf_replica_normal error: cg_norm != [elemental,frobenius]"
  case ("elemental");chi2 = chi2_weiss_replica_normal_elemental(a)
  case ("frobenius");chi2 = chi2_weiss_replica_normal_frobenius(a)
  end select
end function chi2_weiss_replica_normal

!> ELEMENTAL NORM: weighted sum over i\omega for each matrix element, then weighted sum over elements
function chi2_weiss_replica_normal_elemental(a) result(chi2)
  real(8),dimension(:)                               :: a
  real(8)                                            :: chi2
  real(8),dimension(Ldelta)                          :: chi2_freq
  real(8),dimension(Nspin,Nspin,Nimp,Nimp)           :: chi2_mtrx,Wmat
  complex(8),dimension(Nspin,Nspin,Nimp,Nimp,Ldelta) :: g0and
  integer                                            :: ilat,jlat,iorb,jorb,ispin,jspin
  !
  g0and = g0and_replica_normal(a)
  !
  do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Nimp,jorb=1:Nimp)
     chi2_freq = abs(g0and(ispin,jspin,iorb,jorb,1:Ldelta) - FGmatrix(ispin,jspin,iorb,jorb,1:Ldelta))
     chi2_mtrx(ispin,jspin,iorb,jorb) = sum(chi2_freq**cg_pow/Wdelta)
     select case(cg_matrix)
     case(0);Wmat(ispin,jspin,iorb,jorb) = 1d0
     case(1);Wmat(ispin,jspin,iorb,jorb) = abs(sum(FGmatrix(ispin,jspin,iorb,jorb,1:Lmats)))/beta
     end select
  enddo
  !
  chi2 = sum(chi2_mtrx / Wmat, Hmask) !Weighted sum over matrix elements
  chi2 = chi2 / Ldelta / count(Hmask) !Normalization over {iw} and Hmask
  !
end function chi2_weiss_replica_normal_elemental

!> FROBENIUS NORM: global \chi^2 for all components, only i\omega are weighted
function chi2_weiss_replica_normal_frobenius(a) result(chi2)
  real(8),dimension(:)                               :: a
  real(8)                                            :: chi2
  real(8),dimension(Ldelta)                          :: chi2_freq
  complex(8),dimension(Nspin,Nspin,Nimp,Nimp,Ldelta) :: g0and
  complex(8),dimension(Nspin*Nimp,Nspin*Nimp)        :: Delta_lso
  integer                                            :: l
  !
  g0and = g0and_replica_normal(a)
  !
  do l=1,Ldelta
     Delta_lso    =  nn2so_reshape(g0and(:,:,:,:,l) - FGmatrix(:,:,:,:,l),Nspin,Nimp)
     chi2_freq(l) =  sqrt(trace(matmul(Delta_lso,conjg(transpose(Delta_lso)))))
  enddo
  !
  chi2 = sum(chi2_freq**cg_pow/Wdelta) !Weighted sum over matsubara frqs
  chi2 = chi2/Ldelta/(Nspin*Nimp) !Normalization over {iw} and Nlso
  !
end function chi2_weiss_replica_normal_frobenius









!######################################################################
! WEISS CHI2 FIT : SUPERC
!######################################################################
function chi2_weiss_replica_superc(a) result(chi2)
  !Evaluate the \chi^2 distance of \Weiss_Anderson function.
  real(8),dimension(:) :: a
  real(8)              :: chi2
  !
  select case(cg_norm)
  case default;stop "chi2_fitgf_replica_normal error: cg_norm != [elemental,frobenius]"
  case ("elemental");chi2 = chi2_weiss_replica_superc_elemental(a)
  case ("frobenius");chi2 = chi2_weiss_replica_superc_frobenius(a)
  end select
end function chi2_weiss_replica_superc

function chi2_weiss_replica_superc_elemental(a) result(chi2)
  real(8),dimension(:)                                 :: a
  real(8)                                              :: chi2
  real(8),dimension(2,Ldelta)                          :: chi2_freq
  real(8),dimension(Nspin,Nspin,Nimp,Nimp)             :: chi2_mtrx,Wmat
  complex(8),dimension(2,Nspin,Nspin,Nimp,Nimp,Ldelta) :: g0and
  integer                                              :: ilat,jlat,iorb,jorb,ispin,jspin
  !
  g0and = g0and_replica_superc(a)
  !
  do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Nimp,jorb=1:Nimp)
     chi2_freq(1,:) = abs(g0and(1,ispin,jspin,iorb,jorb,1:Ldelta) - FGmatrix(ispin,jspin,iorb,jorb,1:Ldelta))
     chi2_freq(2,:) = abs(g0and(2,ispin,jspin,iorb,jorb,1:Ldelta) - FFmatrix(ispin,jspin,iorb,jorb,1:Ldelta))
     chi2_mtrx(ispin,jspin,iorb,jorb) = sum(chi2_freq(1,:)**cg_pow/Wdelta) + sum(chi2_freq(2,:)**cg_pow/Wdelta)
     select case(cg_matrix)
     case(0);Wmat(ispin,jspin,iorb,jorb) = 1d0
     case(1);Wmat(ispin,jspin,iorb,jorb) = abs(sum(FGmatrix(ispin,jspin,iorb,jorb,1:Lmats)))/beta
     end select
  enddo
  !
  chi2 = sum(chi2_mtrx / Wmat, Hmask)
  chi2 = chi2 / Ldelta / count(Hmask)
  !
end function chi2_weiss_replica_superc_elemental

function chi2_weiss_replica_superc_frobenius(a) result(chi2)
  real(8),dimension(:)                                 :: a
  real(8)                                              :: chi2
  real(8),dimension(Ldelta)                          :: chi2_freq
  complex(8),dimension(2,Nspin,Nspin,Nimp,Nimp,Ldelta) :: g0and
  complex(8),dimension(Nspin*Nimp,Nspin*Nimp)        :: D,F
  integer                                              :: l
  !
  g0and = g0and_replica_superc(a)
  !
  do l=1,Ldelta
     D    =  nn2so_reshape(g0and(1,:,:,:,:,l) - FGmatrix(:,:,:,:,l),Nspin,Nimp)
     F    =  nn2so_reshape(g0and(2,:,:,:,:,l) - FFmatrix(:,:,:,:,l),Nspin,Nimp)          
     chi2_freq(l) =  &
          sqrt(trace(matmul(D,conjg(transpose(D)))))  + &
          sqrt(trace(matmul(F,conjg(transpose(F)))))
  enddo
  !
  chi2 = sum(chi2_freq**cg_pow/Wdelta) !Weighted sum over matsubara frqs
  chi2 = chi2/Ldelta/(Nspin*Nimp)      !Normalization over {iw} and Nlso
  !
end function chi2_weiss_replica_superc_frobenius















!######################################################################
! WEISS GRAD chi2: NORMAL
!######################################################################
function grad_chi2_weiss_replica_normal(a) result(dchi2)
  real(8),dimension(:)       :: a
  real(8),dimension(size(a)) :: dchi2
  !
  select case(cg_norm)
  case default;stop "chi2_fitgf_replica_normal error: cg_norm != [elemental,frobenius]"
  case ("elemental");dchi2 = grad_chi2_weiss_replica_normal_elemental(a)
  case ("frobenius");dchi2 = grad_chi2_weiss_replica_normal_frobenius(a)
  end select
  !
end function grad_chi2_weiss_replica_normal

!> ELEMENTAL NORM: weighted sum over i\omega for each matrix element, then weighted sum over elements
function grad_chi2_weiss_replica_normal_elemental(a) result(dchi2)
  real(8),dimension(:)                                       :: a
  real(8),dimension(size(a))                                 :: dchi2
  real(8),dimension(Nspin,Nspin,Nimp,Nimp,size(a))           :: df
  real(8),dimension(Nspin,Nspin,Nimp,Nimp)                   :: Wmat
  complex(8),dimension(Nspin,Nspin,Nimp,Nimp,Ldelta)         :: g0and
  complex(8),dimension(Nspin,Nspin,Nimp,Nimp,Ldelta,size(a)) :: dg0and
  complex(8),dimension(Ldelta)                               :: Ftmp
  real(8),dimension(Ldelta)                                  :: Ctmp
  integer                                                    :: ia
  !
  !
  g0and  = g0and_replica_normal(a)
  dg0and = grad_g0and_replica_normal(a)
  !
  do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Nimp,jorb=1:Nimp)
     Ftmp = g0and(ispin,jspin,iorb,jorb,1:Ldelta) - FGmatrix(ispin,jspin,iorb,jorb,1:Ldelta)
     Ctmp = abs(Ftmp)**(cg_pow-2)
     do ia = 1,size(a)
        df(ispin,jspin,iorb,jorb,ia) = & !Weighted sum over matsubara frqs
             sum( dreal(Ftmp) * dreal(dg0and(ispin,jspin,iorb,jorb,:,ia)) * Ctmp/Wdelta ) + &
             sum( dimag(Ftmp) * dimag(dg0and(ispin,jspin,iorb,jorb,:,ia)) * Ctmp/Wdelta )
     enddo
     select case(cg_matrix)
     case(0) 
        Wmat(ispin,jspin,iorb,jorb) = 1d0
     case(1) 
        Wmat(ispin,jspin,iorb,jorb) = abs(sum(FGmatrix(ispin,jspin,iorb,jorb,1:Lmats)))/beta
     end select
  enddo
  !
  do ia=1,size(a)
     dchi2(ia) = + cg_pow * sum( df(:,:,:,:,ia) / Wmat, Hmask) !Weighted sum over matrix elements
     dchi2(ia) = dchi2(ia) / Ldelta / count(Hmask)             !Normalization over {iw} and Hmask
  enddo
  !
end function grad_chi2_weiss_replica_normal_elemental

function grad_chi2_weiss_replica_normal_frobenius(a) result(dchi2)
  real(8),dimension(:)                                       :: a
  real(8),dimension(size(a))                                 :: dchi2
  real(8),dimension(Ldelta,size(a))                          :: df
  complex(8),dimension(Nspin,Nspin,Nimp,Nimp,Ldelta)         :: g0and
  complex(8),dimension(Nspin,Nspin,Nimp,Nimp,Ldelta,size(a)) :: dg0and
  complex(8),dimension(Ldelta)                               :: Ftmp
  real(8),dimension(Ldelta,size(a))                          :: dChi_freq
  integer                                                    :: i,j,idelta,ilat,jlat,iorb,jorb,ispin,jspin
  !
  !
  print*, "                                                                     "
  print*, "WARNING: the analytic gradient of the Weiss field Frobenius distance "
  print*, "         has been found giving /QUALITATIVELY WRONG/ fitted spectra! "
  print*, "         > the issue has emerged in some Nlat=Nspin=Nimp=1 test runs "
  print*, "         > while the bug is investigated please switch cg_scheme to  "
  print*, "           'delta' or cg_norm to 'elemental' (or give up analytic cg)"
  print*, "                                                                     "
  !
  !
  g0and  = g0and_replica_normal(a)
  dg0and = grad_g0and_replica_normal(a)
  Ftmp=zero
  df=zero
  !
  do idelta=1,Ldelta
     do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Nimp,jorb=1:Nimp)
        !
        Ftmp(idelta) = Ftmp(idelta) + abs(g0and(ispin,jspin,iorb,jorb,idelta)-FGmatrix(ispin,jspin,iorb,jorb,idelta))**2
        !
        do j=1,size(a)
           df(idelta,j) = df(idelta,j) + &
                dreal(g0and(ispin,jspin,iorb,jorb,idelta) - FGmatrix(ispin,jspin,iorb,jorb,idelta))* &
                dreal(dg0and(ispin,jspin,iorb,jorb,idelta,j)) + &
                dimag(g0and(ispin,jspin,iorb,jorb,idelta) - FGmatrix(ispin,jspin,iorb,jorb,idelta)) * &
                dimag(dg0and(ispin,jspin,iorb,jorb,idelta,j))
        enddo
     enddo
     Ftmp(idelta) = cg_pow * (sqrt(Ftmp(idelta))**(cg_pow-2)) / Wdelta(idelta)
     dchi_freq(idelta,:) = Ftmp(idelta) * df(idelta,:)
  enddo
  !
  dchi2 = sum(dchi_freq,1)/Ldelta/(Nspin*Nimp)
  !
end function grad_chi2_weiss_replica_normal_frobenius



!######################################################################
! WEISS GRAD chi2: SUPERC ??
!######################################################################

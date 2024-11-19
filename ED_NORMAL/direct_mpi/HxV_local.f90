  do i=1,Nloc
     iup = iup_index(i+mpiIshift,DimUp)
     idw = idw_index(i+mpiIshift,DimUp)
     !
     mup = Hs(1)%map(iup)
     mdw = Hs(2)%map(idw)
     !
     nup = bdecomp(mup,Ns)
     ndw = bdecomp(mdw,Ns)
     !
     !> HxV_imp: Diagonal Elements, i.e. local part
     htmp = zero
     do iorb=1,Nimp
        htmp = htmp + impHloc(1,1,iorb,iorb)*Nup(iorb)
        htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb)*Ndw(iorb)
        htmp = htmp - xmu*(Nup(iorb)+Ndw(iorb))
     enddo
     if(any(spin_field(:,3)/=0d0))then
        !F_z.S^z:= F_z.(n_up-n_dw)
        do iorb=1,Nimp
           htmp = htmp + spin_field(iorb,3)*(nup(iorb)-ndw(iorb))
        enddo
     endif
     !
     !
     !
     !> H_Int: Kanamori interaction part. non-local S-E and P-H terms commented below.
     !
     ! density-density interaction: same orbital, opposite spins:
     !  = \sum_\a U_\a*(n_{\a,up}*n_{\a,dw})
     do iorb=1,Nimp
        htmp = htmp + Uloc(iorb)*Nup(iorb)*Ndw(iorb)
     enddo
     if(Nimp>1)then
        !density-density interaction: different orbitals, opposite spins:
        ! =   U'   *     sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
        ! =  (Uloc-2*Jh)*sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
        do iorb=1,Nimp
           do jorb=iorb+1,Nimp
              htmp = htmp + Ust*(Nup(iorb)*Ndw(jorb) + Nup(jorb)*Ndw(iorb))
           enddo
        enddo
        !density-density interaction: different orbitals, parallel spins
        ! = \sum_{i<j}    U''     *[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
        ! = \sum_{i<j} (Uloc-3*Jh)*[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
        do iorb=1,Nimp
           do jorb=iorb+1,Nimp
              htmp = htmp + (Ust-Jh)*(Nup(iorb)*Nup(jorb) + Ndw(iorb)*Ndw(jorb))
           enddo
        enddo
     endif
     !if using the Hartree-shifted chemical potential: mu=0 for half-filling
     !sum up the contributions of hartree terms:
     if(hfmode)then
        do iorb=1,Nimp
           htmp = htmp - 0.5d0*Uloc(iorb)*(Nup(iorb)+Ndw(iorb)) !+ 0.25d0*uloc(iorb)
        enddo
        if(Nimp>1)then
           do iorb=1,Nimp
              do jorb=iorb+1,Nimp
                 htmp=htmp-0.5d0*Ust*(Nup(iorb)+Ndw(iorb)+Nup(jorb)+Ndw(jorb))!+0.25d0*Ust
                 htmp=htmp-0.5d0*(Ust-Jh)*(Nup(iorb)+Ndw(iorb)+Nup(jorb)+Ndw(jorb))!+0.25d0*(Ust-Jh)
              enddo
           enddo
        endif
     endif
     !
     !
     !> H_Bath: local bath energy contribution.
     !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Nimp\sum_{l=1,Nbath}\e^a_l n^a_l
     do iorb=1,size(bath_diag,2)
        do ibath=1,Nbath
           ialfa = getBathStride(iorb,ibath)
           htmp =htmp + bath_diag(1    ,iorb,ibath)*Nup(ialfa) !UP
           htmp =htmp + bath_diag(Nspin,iorb,ibath)*Ndw(ialfa) !DW
        enddo
     enddo
     !
     hv(i) = hv(i) + htmp*vin(i)
     !
  enddo
  



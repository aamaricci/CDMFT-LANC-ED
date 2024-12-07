  do i=MpiIstart,MpiIend
     iup = iup_index(i,DimUp)
     idw = idw_index(i,DimUp)
     !
     mup = Hsector%H(1)%map(iup)
     mdw = Hsector%H(2)%map(idw)
     !
     Nup = bdecomp(mup,Ns)
     Ndw = bdecomp(mdw,Ns)
     !
     !> H_Imp: Diagonal Elements, i.e. local part
     htmp = zero
     do iorb=1,Nimp
        htmp = htmp + impHloc(1,1,iorb,iorb)*Nup(iorb)
        htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb)*Ndw(iorb)
        htmp = htmp - xmu*( Nup(iorb)+Ndw(iorb) )
     enddo
     ! if(any(spin_field(:,3)/=0d0))then
     !    !F_z.S^z:= F_z.(n_up-n_dw)
     !    do iorb=1,Nimp
     !       htmp = htmp + spin_field(iorb,3)*(nup(iorb)-ndw(iorb))
     !    enddo
     ! endif
     !
     !
     !> H_Int: Kanamori interaction part. non-local S-E and P-H terms commented below.
     !
     !density-density interaction: same orbital, opposite spins:
     ! = \sum_\a U_\a*(n_{\a,up}*n_{\a,dw})
     do ilat=1,Nlat
        do iorb=1,Norb
           io = iorb + (ilat-1)*Norb
           htmp = htmp + Uloc(iorb)*nup(io)*ndw(io)
        enddo
     enddo
     if(Norb>1)then
        !density-density interaction: different orbitals, opposite spins:
        ! =   U'   *     sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
        ! =  (Uloc-2*Jh)*sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
        do ilat=1,Nlat
           do iorb=1,Norb
              do jorb=iorb+1,Norb
                 io = iorb + (ilat-1)*Norb
                 jo = jorb + (ilat-1)*Norb
                 htmp = htmp + Ust*(nup(io)*ndw(jo) + nup(jo)*ndw(io))
              enddo
           enddo
        enddo
        !density-density interaction: different orbitals, parallel spins
        ! = \sum_{i<j}    U''     *[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
        ! = \sum_{i<j} (Uloc-3*Jh)*[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
        do ilat=1,Nlat
           do iorb=1,Norb
              do jorb=iorb+1,Norb
                 io = iorb + (ilat-1)*Norb
                 jo = jorb + (ilat-1)*Norb
                 htmp = htmp + (Ust-Jh)*(nup(io)*nup(jo) + ndw(io)*ndw(jo))
              enddo
           enddo
        enddo
     endif
     !if using the Hartree-shifted chemical potential: mu=0 for half-filling
     !sum up the contributions of hartree terms:
     if(hfmode)then
        do ilat=1,Nlat
           do iorb=1,Norb
              io   = iorb + (ilat-1)*Norb
              htmp = htmp - 0.5d0*Uloc(iorb)*(nup(io)+ndw(io)) + 0.25d0*Uloc(iorb)
           enddo
        enddo
        if(Norb>1)then
           do ilat=1,Nlat
              do iorb=1,Norb
                 do jorb=iorb+1,Norb
                    io = iorb + (ilat-1)*Norb
                    jo = jorb + (ilat-1)*Norb
                    htmp=htmp-0.5d0*Ust*(nup(io)+ndw(io)+nup(jo)+ndw(jo)) +0.25d0*Ust
                    htmp=htmp-0.5d0*(Ust-Jh)*(nup(io)+ndw(io)+nup(jo)+ndw(jo))+0.25d0*(Ust-Jh)
                 enddo
              enddo
           enddo
        endif
     endif
     !
     !
     !> H_Bath: local bath energy contribution.
     !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
     do iorb=1,Nimp
        do ibath=1,Nbath
           ialfa = getBathStride(iorb,ibath)
           htmp = htmp + bath_diag(1    ,iorb,ibath)*Nup(ialfa) !UP
           htmp = htmp + bath_diag(Nspin,iorb,ibath)*Ndw(ialfa) !DW
        enddo
     enddo
     !
     select case(MpiStatus)
     case (.true.)
        call sp_insert_element(MpiComm,spH0d,htmp,i,i)
     case (.false.)
        call sp_insert_element(spH0d,htmp,i,i)
     end select
     !
  enddo


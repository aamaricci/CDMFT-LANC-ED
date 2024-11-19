  !We build the transposed H_non_local here (symmetric)
  !to comply with the MPI decomposition of the matrix.
  !A better MPI handling might be necessary here...
  do i=1,Nloc
     iup = iup_index(i+mpiIshift,DimUp)
     idw = idw_index(i+mpiIshift,DimUp)
     !
     mup = Hs(1)%map(iup)
     mdw = Hs(2)%map(idw)
     !
     Nuo = bdecomp(mup,Ns)
     Ndw = bdecomp(mdw,Ns)
     !
     ! SPIN-EXCHANGE (S-E) and PAIR-HOPPING TERMS
     !    S-E: J c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up  (i.ne.j) 
     !    S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
     if(Jhflag.AND.Jx/=0d0)then
        do iorb=1,Nimp
           do jorb=1,Nimp
              Jcondition=(&
                   (iorb/=jorb).AND.&
                   (nup(jorb)==1).AND.&
                   (ndw(iorb)==1).AND.&
                   (ndw(jorb)==0).AND.&
                   (nup(iorb)==0))
              if(Jcondition)then
                 call c(iorb,mdw,k1,sg1)  !DW
                 call cdg(jorb,k1,k2,sg2) !DW
                 jdw=binary_search(Hs(2)%map,k2)
                 call c(jorb,mup,k3,sg3)  !UP
                 call cdg(iorb,k3,k4,sg4) !UP
                 jup=binary_search(Hs(1)%map,k4)
                 htmp = Jx*sg1*sg2*sg3*sg4
                 j = jup + (jdw-1)*dimup
                 !
                 Hv(i) = Hv(i) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
     endif
     !
     ! PAIR-HOPPING (P-H) TERMS
     !    P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
     !    P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
     if(Jhflag.AND.Jp/=0d0)then
        do iorb=1,Nimp
           do jorb=1,Nimp
              Jcondition=(&
                   (nup(jorb)==1).AND.&
                   (ndw(jorb)==1).AND.&
                   (ndw(iorb)==0).AND.&
                   (nup(iorb)==0))
              if(Jcondition)then
                 call c(jorb,mdw,k1,sg1)       !c_jorb_dw
                 call cdg(iorb,k1,k2,sg2)      !c^+_iorb_dw
                 jdw = binary_search(Hs(2)%map,k2)
                 call c(jorb,mup,k3,sg3)       !c_jorb_up
                 call cdg(iorb,k3,k4,sg4)      !c^+_iorb_up
                 jup = binary_search(Hs(1)%map,k4)
                 htmp = Jp*sg1*sg2*sg3*sg4
                 j = jup + (jdw-1)*dimup
                 !
                 Hv(i) = Hv(i) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
     endif
     !
  enddo

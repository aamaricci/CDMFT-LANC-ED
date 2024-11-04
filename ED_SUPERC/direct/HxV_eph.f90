  ! We build the electronic part of the electron-phonon interaction:
  ! Diagonal terms: Sum_iorb g_iorb,iorb n_iorb*(bdg + b)
  htmp=zero
  do iorb=1,Norb
     htmp = htmp + g_ph(iorb,iorb)*(nup(iorb)+ndw(iorb)  - 1.d0) !electron part
  enddo
  !
  if(iph<DimPh)then !bdg
     j = i_el + (iph  )*DimEl
     Hv(j-MpiIshift) = Hv(j-MpiIshift) + htmp*vin(i)*sqrt(dble(iph))
  endif
  if(iph>1)then !b
     j = i_el + (iph-2)*DimEl
     Hv(j-MpiIshift) = Hv(j-MpiIshift) + htmp*vin(i)*sqrt(dble(iph-1))
  endif
  !
  ! Off-Diagonal terms: Sum_iorb,jorb g_iorb,jorb cdg_iorb*c_jorb*(bdg+b)
  ! UP spin
  ! remark: iorb=jorb can't have simultaneously n=0 and n=1 (Jcondition)
  do iorb=1,Norb
     do jorb=1,Norb
        Jcondition = (g_ph(iorb,jorb)/=zero) .AND. &
             (nup(jorb)==1) .AND. (nup(iorb)==0)
        if(Jcondition)then
           call c(jorb,m,k1,sg1)
           call cdg(iorb,k1,k2,sg2)
           j_el  = binary_search(Hsector%H(1)%map,k2)
           htmp  = g_ph(iorb,jorb)*sg1*sg2
           !
           if(iph<DimPh)then !bdg
              j     = j_el + (iph)*DimEl
              hv(j-MpiIshift) = hv(j-MpiIshift) + htmp*vin(i)*sqrt(dble(iph))
           endif
           if(iph>1)then !b
              j     = j_el + (iph-2)*DimEl
              hv(j-MpiIshift) = hv(j-MpiIshift) + htmp*vin(i)*sqrt(dble(iph-1))
           endif
        endif
     enddo
  enddo
  ! DW spin
  ! remark: iorb=jorb can't have simultaneously n=0 and n=1 (Jcondition)
  do iorb=1,Norb
     do jorb=1,Norb
        Jcondition = &
             (g_ph(iorb,jorb)/=zero) .AND. &
             (Ndw(jorb)==1) .AND. (Ndw(iorb)==0)
        if(Jcondition)then
           call c(jorb+Ns,m,k1,sg1)
           call cdg(iorb+Ns,k1,k2,sg2)
           j_el  = binary_search(Hsector%H(1)%map,k2)
           htmp = g_ph(iorb,jorb)*sg1*sg2
           !
           if(iph<DimPh)then !bdg
              j     = j_el + (iph)*DimEl
              hv(j-MpiIshift) = hv(j-MpiIshift) + htmp*vin(i)*sqrt(dble(iph))
           endif
           if(iph>1)then !b
              j     = j_el + (iph-2)*DimEl
              hv(j-MpiIshift) = hv(j-MpiIshift) + htmp*vin(i)*sqrt(dble(iph-1))
           endif
        endif
        !              
     enddo
  enddo

 !We build the electronic part of the electron-phonon interaction:
  do i=MpiIstart,MpiIend
     !
     m  = Hsector%H(1)%map(i)
     ib = bdecomp(m,2*Ns)
     !
     do iorb=1,Norb
        nup(iorb)=dble(ib(iorb))
        ndw(iorb)=dble(ib(iorb+Ns))
     enddo
     ! Diagonal terms: Sum_iorb g_iorb,iorb n_iorb
     htmp = zero
     do iorb=1,Norb
        htmp = htmp + g_ph(iorb,iorb)*(nup(iorb)+ndw(iorb) -1.d0)
     enddo
     !
     select case(MpiStatus)
     case (.true.)
        call sp_insert_element(MpiComm,spH0e_eph,htmp,i,i)
     case (.false.)
        call sp_insert_element(spH0e_eph,htmp,i,i)
     end select
     !
     ! Off-Diagonal terms: Sum_iorb,jorb g_iorb,jorb cdg_iorb*c_jorb
     ! UP spin
     ! remark: iorb=jorb can't have simultaneously n=0 and n=1 (Jcondition)
     do iorb=1,Norb
        do jorb=1,Norb
           Jcondition = &
                (g_ph(iorb,jorb)/=zero) .AND. &
                (ib(jorb)==1) .AND. (ib(iorb)==0)
           if(Jcondition)then
              call c(jorb,m,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              j = binary_search(Hsector%H(1)%map,k2)
              htmp = conjg(g_ph(iorb,jorb))*sg1*sg2
              !
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0e_eph,htmp,i,j)
              case (.false.)
                 call sp_insert_element(spH0e_eph,htmp,i,j)
              end select
              
           endif
        enddo
     enddo
     ! DW spin
     ! remark: iorb=jorb can't have simultaneously n=0 and n=1 (Jcondition)
     do iorb=1,Norb
        do jorb=1,Norb
           Jcondition = &
                (g_ph(iorb,jorb)/=zero) .AND. &
                (ib(jorb+Ns)==1) .AND. (ib(iorb+Ns)==0)
           if (Jcondition) then
              call c(jorb+Ns,m,k1,sg1)
              call cdg(iorb+Ns,k1,k2,sg2)
              j  = binary_search(Hsector%H(1)%map,k2)
              htmp = conjg(g_ph(iorb,jorb))*sg1*sg2
              !
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0e_eph,htmp,i,j)
              case (.false.)
                 call sp_insert_element(spH0e_eph,htmp,i,j)
              end select
              !
           endif
              
        enddo
     enddo
     !
  enddo
  ! Here we build the phoninc part of the electron-phonon interaction: (b^+ + b)
  htmp = zero
  do iph=1,DimPh
     if(iph < DimPh) then ! bdg|iph>
        htmp = sqrt(dble(iph))
        call sp_insert_element(spH0ph_eph,htmp,iph,iph+1)
     end if
     if(iph>1) then ! b|iph>
        htmp = sqrt(dble(iph - 1))
        call sp_insert_element(spH0ph_eph,htmp,iph,iph-1)
     end if
  end do

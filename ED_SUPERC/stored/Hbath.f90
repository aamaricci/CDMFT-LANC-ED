  do i=MpiIstart,MpiIend
     m  = Hsector%H(1)%map(i)
     ib = bdecomp(m,2*Ns)
     !
     do iorb=1,Norb
        nup(iorb)=dble(ib(iorb))
        ndw(iorb)=dble(ib(iorb+Ns))
     enddo

     !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
     htmp=zero
     do kp=1,Nbath
        do iorb=1,Nimp
           ialfa = getBathStride(iorb,kp)
           htmp = htmp + bath_diag(1,          1,iorb,kp)*ib(ialfa)     !UP
           htmp = htmp + bath_diag(Nnambu,     1,iorb,kp)*ib(ialfa))    !UP CHANGE due to Nambu repr
           htmp = htmp + bath_diag(1,      Nspin,iorb,kp)*ib(ialfa+Ns)  !DW
           htmp = htmp + bath_diag(Nnambu, Nspin,iorb,kp)*ib(ialfa+Ns)) !DW CHANGE due to Nambu repr
        enddo
     enddo
     !
     select case(MpiStatus)
     case (.true.)
        call sp_insert_element(MpiComm,spH0,htmp,i,i)
     case (.false.)
        call sp_insert_element(spH0,htmp,i,i)
     end select
     !
     !off-diagonal elements
     do kp=1,Nbath
        do iorb=1,Nimp
           do jorb=1,Nimp
              !Nambu=1.1
              !UP
              ialfa = getBathStride(iorb,kp)
              ibeta = getBathStride(jorb,kp)
              Jcondition = &
                   (hbath_tmp(1,1,1,1,iorb,jorb,kp)/=zero)   .AND. &
                   (ib(ibeta)==1) .AND. (ib(ialfa)==0)
              if (Jcondition)then
                 call c(ibeta,m,k1,sg1)
                 call cdg(ialfa,k1,k2,sg2)
                 j    = binary_search(Hsector%H(1)%map,k2)
                 htmp = conjg( hbath_tmp(1,1,1,1,iorb,jorb,kp) )*sg1*sg2
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
              endif
              ! !DW
              ! ialfa = getBathStride(iorb,kp) + Ns
              ! ibeta = getBathStride(jorb,kp) + Ns
              ! Jcondition = &
              !      (hbath_tmp(1,1,Nspin,Nspin,iorb,jorb,kp)/=zero) .AND. &
              !      (ib(ibeta)==0).AND. (ib(ialfa)==1)
              ! if (Jcondition)then
              !    call cdg(ibeta,m,k1,sg1)
              !    call c(ialfa,k1,k2,sg2)
              !    j = binary_search(Hsector%H(1)%map,k2)
              !    htmp = conjg(hbath_tmp(1,1,Nspin,Nspin,iorb,jorb,kp) )*sg1*sg2
              !    select case(MpiStatus)
              !    case (.true.)
              !       call sp_insert_element(MpiComm,spH0,htmp,i,j)
              !    case (.false.)
              !       call sp_insert_element(spH0,htmp,i,j)
              !    end select
              ! endif
              !
              !
              !Nambu=2.2
              !UP
              ! ialfa = getBathStride(iorb,kp)
              ! ibeta = getBathStride(jorb,kp)
              ! Jcondition = &
              !      (hbath_tmp(Nambu,Nambu,1,1,iorb,jorb,kp)/=zero)   .AND. &
              !      (ib(ibeta)==1) .AND. (ib(ialfa)==0)
              ! if (Jcondition)then
              !    call c(ibeta,m,k1,sg1)
              !    call cdg(ialfa,k1,k2,sg2)
              !    j    = binary_search(Hsector%H(1)%map,k2)
              !    htmp = conjg( hbath_tmp(Nambu,Nambu,1,1,iorb,jorb,kp) )*sg1*sg2
              !    select case(MpiStatus)
              !    case (.true.)
              !       call sp_insert_element(MpiComm,spH0,htmp,i,j)
              !    case (.false.)
              !       call sp_insert_element(spH0,htmp,i,j)
              !    end select
              ! endif
              !DW
              ialfa = getBathStride(iorb,kp) + Ns
              ibeta = getBathStride(jorb,kp) + Ns
              Jcondition = &
                   (hbath_tmp(Nnambu,Nnambu,Nspin,Nspin,iorb,jorb,kp)/=zero) .AND. &
                   (ib(ibeta)==0).AND. (ib(ialfa)==1)
              if (Jcondition)then
                 call cdg(ibeta,m,k1,sg1)
                 call c(ialfa,k1,k2,sg2)
                 j = binary_search(Hsector%H(1)%map,k2)
                 htmp = -(hbath_tmp(Nambu,Nambu,Nspin,Nspin,iorb,jorb,kp) )*sg1*sg2
                 !htmp = conjg(hbath_tmp(Nambu,Nambu,Nspin,Nspin,iorb,jorb,kp) )*sg1*sg2
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
              endif

           enddo
        enddo
     enddo
     !
     !
     !anomalous pair-creation/destruction
     do kp=1,Nbath
        !
        do iorb=1,Nimp
           do jorb=1,Nimp
              !
              !UP-DW \Delta_l cdg_{\up,a} cdg_{\dw,b}
              ialfa = getBathStride(iorb,kp)
              ibeta = getBathStride(jorb,kp) + Ns
              Jcondition = &
                   (hbath_tmp(1,Nambu,1,1,iorb,jorb,kp)/=zero) .AND. &
                   (ib(ibeta)==0) .AND. (ib(ialfa)==0)
              if(Jcondition)then
                 call cdg(ibeta,m,k1,sg1)
                 call cdg(ialfa,k1,k2,sg2)
                 j=binary_search(Hsector%H(1)%map,k2)
                 htmp=conjg(hbath_tmp(1,Nambu,iorb,jorb,kp))*sg1*sg2
                 !
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
              endif
              !DW-UP \Delta_l c_{\dw,a} c_{\up,b}
              ialfa = getBathStride(iorb,kp) + Ns
              ibeta = getBathStride(jorb,kp)
              Jcondition = &
                   (hbath_tmp(Nambu,1,1,1,iorb,jorb,kp)/=zero) .AND. &
                   (ib(ibeta)==1) .AND. (ib(ialfa)==1)
              if(Jcondition)then
                 call c(ibeta,m,k1,sg1)
                 call c(ialfa,k1,k2,sg2)
                 j=binary_search(Hsector%H(1)%map,k2)
                 htmp=conjg(hbath_tmp(Nnambu,1,iorb,jorb,kp))*sg1*sg2
                 !
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
              endif
              !
           enddo
        enddo
        !
     enddo
     !
  enddo

  do i=MpiIstart,MpiIend
     m  = Hsector%H(1)%map(i)
     ib = bdecomp(m,2*Ns)
     !
     do iorb=1,Norb
        nup(iorb)=dble(ib(iorb))
        ndw(iorb)=dble(ib(iorb+Ns))
     enddo

     !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
     select case(bath_type)
     case default 
        htmp=zero
        do iorb=1,size(dmft_bath%e,2)
           do kp=1,Nbath
              ialfa=getBathStride(iorb,kp)
              htmp =htmp + dmft_bath%e(1,iorb,kp)*ib(ialfa)        !UP
              htmp =htmp + dmft_bath%e(Nspin,iorb,kp)*ib(ialfa+Ns) !DW
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
     case ("replica","general")
        htmp=zero
        do kp=1,Nbath
           do iorb=1,Norb
              ialfa = getBathStride(iorb,kp)
              htmp = htmp + bath_diag(1          ,iorb,kp)*ib(ialfa)    !UP
              htmp = htmp -conjg(bath_diag(Nnambu*Nspin,iorb,kp)*ib(ialfa+Ns)) !DW CHANGE due to Nambu repr
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
           do iorb=1,Norb
              do jorb=1,Norb
                 !UP
                 ialfa = getBathStride(iorb,kp)
                 ibeta = getBathStride(jorb,kp)
                 Jcondition = &
                      (hbath_tmp(1,1,iorb,jorb,kp)/=zero)   .AND. &
                      (ib(ibeta)==1) .AND. (ib(ialfa)==0)
                 if (Jcondition)then
                    call c(ibeta,m,k1,sg1)
                    call cdg(ialfa,k1,k2,sg2)
                    j = binary_search(Hsector%H(1)%map,k2)
                    htmp = conjg(hbath_tmp(1,1,iorb,jorb,kp))*sg1*sg2
                    !
                    select case(MpiStatus)
                    case (.true.)
                       call sp_insert_element(MpiComm,spH0,htmp,i,j)
                    case (.false.)
                       call sp_insert_element(spH0,htmp,i,j)
                    end select
                    !
                 endif
                 !DW
                 ialfa = getBathStride(iorb,kp) + Ns
                 ibeta = getBathStride(jorb,kp) + Ns
                 Jcondition = &
                      (hbath_tmp(Nnambu*Nspin,Nnambu*Nspin,iorb,jorb,kp)/=zero) .AND. &
                      (ib(ibeta)==0).AND. (ib(ialfa)==1)
                 if (Jcondition)then
                    call cdg(ibeta,m,k1,sg1)
                    call c(ialfa,k1,k2,sg2)
                    j = binary_search(Hsector%H(1)%map,k2)
                    htmp = -(hbath_tmp(Nnambu*Nspin,Nnambu*Nspin,iorb,jorb,kp))*sg1*sg2
                    select case(MpiStatus)
                    case (.true.)
                       call sp_insert_element(MpiComm,spH0,htmp,i,j)
                    case (.false.)
                       call sp_insert_element(spH0,htmp,i,j)
                    end select
                    !
                 endif
              enddo
           enddo
        enddo
        !
     end select


     !anomalous pair-creation/destruction
     select case(bath_type)
     case default 
        do iorb=1,size(dmft_bath%e,2)
           do kp=1,Nbath
              ms=getBathStride(iorb,kp)
              !\Delta_l c_{\up,ms} c_{\dw,ms}
              if( (dmft_bath%d(1,iorb,kp)/=0d0) .AND. (ib(ms)==1) .AND. (ib(ms+Ns)==1) )then
                 call c(ms,m,k1,sg1)
                 call c(ms+Ns,k1,k2,sg2)
                 j=binary_search(Hsector%H(1)%map,k2)
                 htmp=one*dmft_bath%d(1,iorb,kp)*sg1*sg2
                 !
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
                 !
              endif
              !\Delta_l cdg_{\up,ms} cdg_{\dw,ms}
              if( (dmft_bath%d(1,iorb,kp)/=0d0) .AND. (ib(ms)==0) .AND. (ib(ms+Ns)==0) )then
                 call cdg(ms+Ns,m,k1,sg1)
                 call cdg(ms,k1,k2,sg2)
                 j=binary_search(Hsector%H(1)%map,k2)
                 htmp=one*dmft_bath%d(1,iorb,kp)*sg1*sg2 !
                 !
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
                 !
              endif
           enddo
        enddo

     case ("replica","general")

        do kp=1,Nbath
           !
           do iorb=1,Norb
              do jorb=1,Norb
                 !
                 !UP-DW \Delta_l cdg_{\up,a} cdg_{\dw,b}
                 ialfa = getBathStride(iorb,kp)
                 ibeta = getBathStride(jorb,kp) + Ns
                 Jcondition = &
                      (hbath_tmp(1,Nnambu,iorb,jorb,kp)/=zero) .AND. &
                      (ib(ibeta)==0) .AND. (ib(ialfa)==0)
                 if(Jcondition)then
                    call cdg(ibeta,m,k1,sg1)
                    call cdg(ialfa,k1,k2,sg2)
                    j=binary_search(Hsector%H(1)%map,k2)
                    htmp=conjg(hbath_tmp(1,Nnambu,iorb,jorb,kp))*sg1*sg2
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
                      (hbath_tmp(Nnambu,1,iorb,jorb,kp)/=zero) .AND. &
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
     end select
     !
  enddo

  do jup=1,DimUp
     do jdw=1,DimDw
        mdw  = Hsector%H(2)%map(jdw)
        Ndw  = bdecomp(mdw,Ns)
        j    = jup + (jdw-1)*dimUp
        !
        !
        !> H_imp: Off-diagonal elements, i.e. non-local part. 
        !remark: iorb=jorb cant have simultaneously n=0 and n=1 (Jcondition)
        do iorb=1,Nimp
           do jorb=1,Nimp
              Jcondition = (impHloc(Nspin,Nspin,iorb,jorb)/=0d0)&
                   .AND.(ibdw(jorb)==1).AND.(ibdw(iorb)==0)
              if (Jcondition) then
                 call c(jorb,mdw,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 idw = binary_search(Hsector%H(2)%map,k2)
                 iup = jup
                 i   = iup + (idw-1)*DimUp
                 htmp = impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
                 !
                 Hv(i) = Hv(i) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
        !
        !
        do ibath=1,Nbath
           do iorb=1,Nimp
              do jorb=1,Nimp
                 !
                 ialfa = getBathStride(iorb,ibath)
                 ibeta = getBathStride(jorb,ibath)
                 Jcondition = &
                      (hbath_tmp(1,1,Nspin,Nspin,iorb,jorb,ibath)/=zero) &
                      .AND. (ibdw(ibeta)==1) .AND. (ibdw(ialfa)==0)
                 !
                 if (Jcondition)then
                    call c(ibeta,mdw,k1,sg1)
                    call cdg(ialfa,k1,k2,sg2)
                    idw = binary_search(Hsector%H(2)%map,k2)
                    iup = jup
                    i   = iup + (idw-1)*DimUp
                    htmp = hbath_tmp(1,1,Nspin,Nspin,iorb,jorb,ibath)*sg1*sg2
                    !
                    Hv(i) = Hv(i) + htmp*vin(j)
                    !
                 endif
              enddo
           enddo
        enddo
        !
        !
        !>H_hyb: hopping terms for a given spin (imp <--> bath)
        do iorb=1,Nimp
           do ibath=1,Nbath
              ialfa=getBathStride(iorb,ibath)
              if( (diag_hybr(Nspin,iorb,ibath)/=0d0) &
                   .AND. (ibdw(iorb)==1) .AND. (ibdw(ialfa)==0) )then
                 call c(iorb,mdw,k1,sg1)
                 call cdg(ialfa,k1,k2,sg2)
                 idw = binary_search(Hsector%H(2)%map,k2)
                 iup = jup
                 i   = iup + (idw-1)*DimUp
                 htmp=diag_hybr(Nspin,iorb,ibath)*sg1*sg2
                 !
                 Hv(i) = Hv(i) + htmp*vin(j)
                 !
              endif
              if( (diag_hybr(Nspin,iorb,ibath)/=0d0) &
                   .AND. (ibdw(iorb)==0) .AND. (ibdw(ialfa)==1) )then
                 call c(ialfa,mdw,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 idw = binary_search(Hsector%H(2)%map,k2)
                 iup = jup
                 i   = iup + (idw-1)*DimUp
                 htmp=diag_hybr(Nspin,iorb,ibath)*sg1*sg2
                 !
                 Hv(i) = Hv(i) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
        !
        !

        !F_0. T^0_ab :=   F_0 . (+ C^+_{a,dw}C_{b,dw})
        !F_z. T^z_ab :=   F_z . (- C^+_{a,dw}C_{b,dw})
        ! if(any(exc_field/=0d0))then
        !    do iorb=1,Nimp
        !       do jorb=iorb+1,Nimp
        !          Jcondition = (Ndw(jorb)==1) .AND. (Ndw(iorb)==0)
        !          if (Jcondition) then
        !             call c(jorb,mdw,k1,sg1)
        !             call cdg(iorb,k1,k2,sg2)
        !             idw = binary_search(Hsector%H(2)%map,k2)
        !             iup = jup
        !             i   = iup + (idw-1)*DimUp 
        !             !
        !             htmp = exc_field(1)*sg1*sg2
        !             Hv(i) = Hv(i) + htmp*vin(j)
        !             !
        !             htmp = -exc_field(4)*sg1*sg2
        !             Hv(i) = Hv(i) + htmp*vin(j)
        !          endif
        !       enddo
        !    enddo
        ! endif

     end do
  enddo

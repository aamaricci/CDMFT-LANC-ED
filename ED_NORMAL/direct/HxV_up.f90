  do jdw=1,DimDw
     do jup=1,DimUp
        mup = Hsector%H(1)%map(jup)
        Nup = bdecomp(mup,Ns)
        j   = jup + (jdw-1)*dimUp
        !
        !
        !> H_imp: Off-diagonal elements, i.e. non-local part. 
        !remark: iorb=jorb cant have simultaneously n=0 and n=1 (Jcondition)
        do iorb=1,Nimp
           do jorb=1,Nimp
              Jcondition = (impHloc(1,1,iorb,jorb)/=0d0)&
                   .AND.(ibup(jorb)==1).AND.(ibup(iorb)==0)
              if (Jcondition) then
                 call c(jorb,mup,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 iup = binary_search(Hsector%H(1)%map,k2)
                 idw = jdw 
                 i   = iup + (idw-1)*DimUp
                 htmp = impHloc(1,1,iorb,jorb)*sg1*sg2
                 !
                 Hv(i) = Hv(i) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
        !
        !> H_Bath: inter-orbital bath hopping contribution.
        do ibath=1,Nbath
           do iorb=1,Nimp
              do jorb=1,Nimp
                 !
                 ialfa = getBathStride(iorb,ibath)
                 ibeta = getBathStride(jorb,ibath)
                 Jcondition = &
                      (hbath_tmp(1,1,1,1,iorb,jorb,ibath)/=0d0) .AND. &
                      (ibup(ibeta)==1) .AND. (ibup(ialfa)==0)
                 !
                 if (Jcondition)then
                    call c(ibeta,mup,k1,sg1)
                    call cdg(ialfa,k1,k2,sg2)
                    iup = binary_search(Hsector%H(1)%map,k2)
                    idw = jdw 
                    i   = iup + (idw-1)*DimUp                    
                    htmp = hbath_tmp(1,1,1,1,iorb,jorb,ibath)*sg1*sg2
                    !
                    hv(i) = hv(i) + htmp*vin(j)
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
              if( (diag_hybr(1,iorb,ibath)/=0d0) .AND. &
                   (ibup(iorb)==1) .AND. (ibup(ialfa)==0) )then              
                 call c(iorb,mup,k1,sg1)
                 call cdg(ialfa,k1,k2,sg2)
                 iup = binary_search(Hsector%H(1)%map,k2)
                 idw = jdw 
                 i   = iup + (idw-1)*DimUp
                 htmp = diag_hybr(1,iorb,ibath)*sg1*sg2
                 !
                 hv(i) = hv(i) + htmp*vin(j)
                 !
              endif
              if( (diag_hybr(1,iorb,ibath)/=0d0) .AND. &
                   (ibup(iorb)==0) .AND. (ibup(ialfa)==1) )then
                 call c(ialfa,mup,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 iup = binary_search(Hsector%H(1)%map,k2)
                 idw = jdw 
                 i   = iup + (idw-1)*DimUp
                 htmp = diag_hybr(1,iorb,ibath)*sg1*sg2
                 !
                 hv(i) = hv(i) + htmp*vin(j)
                 !
              endif
           enddo
        enddo


        !F_0. T^0_ab :=   F_0 . C^+_{a,up}C_{b,up}
        !F_z. T^z_ab :=   F_z . C^+_{a,up}C_{b,up}
        ! if(any(exc_field/=0d0))then
        !    do iorb=1,Nimp
        !       do jorb=iorb+1,Nimp
        !          Jcondition = (Nup(jorb)==1) .AND. (Nup(iorb)==0)
        !          if (Jcondition) then
        !             call c(jorb,mup,k1,sg1)
        !             call cdg(iorb,k1,k2,sg2)
        !             iup = binary_search(Hsector%H(1)%map,k2)
        !             idw = jdw 
        !             i   = iup + (idw-1)*DimUp
        !             !
        !             htmp = exc_field(1)*sg1*sg2
        !             hv(i) = hv(i) + htmp*vin(j)
        !             !
        !             htmp = exc_field(4)*sg1*sg2
        !             hv(i) = hv(i) + htmp*vin(j)
        !          endif
        !       enddo
        !    enddo
        ! endif

     enddo
  enddo
  





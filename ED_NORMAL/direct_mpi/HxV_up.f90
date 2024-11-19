  do jdw=1,MpiQdw
     do jup=1,DimUp
        mup = Hs(1)%map(jup)
        nup = bdecomp(mup,Ns)
        !
        j    = jup + (jdw-1)*dimUp
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
                 iup  = binary_search(Hs(1)%map,k2)
                 idw  = jdw
                 i    = iup + (idw-1)*dimUp
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
                      (hbath_Reconstructed(1,1,iorb,jorb,ibath)/=0d0) &
                      .AND. (ibup(ibeta)==1) .AND. (ibup(ialfa)==0)
                 !
                 if (Jcondition)then
                    call c(ibeta,mup,k1,sg1)
                    call cdg(ialfa,k1,k2,sg2)
                    iup  = binary_search(Hs(1)%map,k2)
                    idw  = jdw
                    i    = iup + (idw-1)*dimUp
                    htmp = hbath_Reconstructed(1,1,iorb,jorb,ibath)*sg1*sg2
                    !
                    hv(i) = hv(i) + htmp*vin(j)
                    !
                 endif
              enddo
           enddo
        enddo
        !
        !>H_hyb: hopping terms for a given spin (imp <--> bath)
        do iorb=1,Nimp
           do ibath=1,Nbath
              ialfa=getBathStride(iorb,ibath)
              if( (diag_hybr(1,iorb,ibath)/=0d0) .AND. &
                   (ibup(iorb)==1) .AND. (ibup(ialfa)==0) )then              
                 call c(iorb,mup,k1,sg1)
                 call cdg(ialfa,k1,k2,sg2)
                 iup = binary_search(Hs(1)%map,k2)
                 idw  = jdw
                 i    = iup + (idw-1)*dimUp
                 htmp = diag_hybr(1,iorb,ibath)*sg1*sg2
                 !
                 hv(i) = hv(i) + htmp*vin(j)
                 !
              endif
              if( (diag_hybr(1,iorb,ibath)/=0d0) .AND. &
                   (ibup(iorb)==0) .AND. (ibup(ialfa)==1) )then
                 call c(ialfa,mup,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 iup = binary_search(Hs(1)%map,k2)
                 idw  = jdw
                 i    = iup + (idw-1)*dimUp
                 htmp = diag_hybr(1,iorb,ibath)*sg1*sg2
                 !
                 hv(i) = hv(i) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
     enddo
  enddo






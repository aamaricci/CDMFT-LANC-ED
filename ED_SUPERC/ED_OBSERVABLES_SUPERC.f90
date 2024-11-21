MODULE ED_OBSERVABLES_SUPERC
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_LINALG
  USE SF_TIMER
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_BATH
  USE ED_HAMILTONIAN_SUPERC
  implicit none
  private
  !

  public :: observables_impurity
  public :: local_energy_impurity
  ! public :: density_matrix_impurity


  logical,save                        :: iolegend=.true.
  real(8),dimension(:),allocatable    :: dens,dens_up,dens_dw
  real(8),dimension(:),allocatable    :: docc
  real(8),dimension(:),allocatable    :: magZ
  real(8),dimension(:,:),allocatable  :: phiscAB
  real(8),dimension(:),allocatable    :: phisc
  real(8),dimension(:,:),allocatable  :: sz2,n2
  real(8)                             :: s2tot
  real(8)                             :: Egs
  real(8)                             :: Ei
  !
  integer                             :: iorb,jorb,istate
  integer                             :: ispin,jspin
  integer                             :: isite,jsite
  integer                             :: ibath
  integer                             :: r,m,k,k1,k2,k3,k4
  integer                             :: iup,idw
  integer                             :: jup,jdw
  integer                             :: mup,mdw
  integer                             :: iph,i_el,j_el,isz
  real(8)                             :: sgn,sgn1,sgn2,sg1,sg2,sg3,sg4
  real(8)                             :: gs_weight
  !
  real(8)                             :: peso
  real(8)                             :: norm
  !
  integer                             :: i,j,ii
  integer                             :: isector,jsector
  integer                             :: idim,idimUP,idimDW
  !
  complex(8),dimension(:),allocatable :: vvinit
  complex(8),dimension(:),allocatable :: state_cvec
  logical                             :: Jcondition
  !
  type(sector)                        :: sectorI,sectorJ


contains 


  !+-------------------------------------------------------------------+
  !PURPOSE  : Lanc method
  !+-------------------------------------------------------------------+
  subroutine observables_superc()
    !Calculate the values of the local observables
    integer                 :: val
    integer,dimension(2*Ns) :: ib
    integer,dimension(2,Ns) :: Nud
    integer,dimension(Ns)   :: IbUp,IbDw
    real(8),dimension(Nimp) :: nup,ndw,Sz,nt

    !
    !LOCAL OBSERVABLES:
    allocate(dens(Nimp))
    allocate(dens_up(Nimp))
    allocate(dens_dw(Nimp))
    allocate(docc(Nimp))
    allocate(magz(Nimp))
    allocate(sz2(Nimp,Nimp))
    allocate(n2(Nimp,Nimp))
    allocate(phisc(Nimp))
    allocate(phiscAB(Nimp,Nimp))
    !
    Egs     = state_list%emin
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0
    phisc   = 0.d0
    phiscAB = 0.d0
    magz    = 0.d0
    sz2     = 0.d0
    n2      = 0.d0
    s2tot   = 0.d0
    !
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       !
#ifdef _MPI
       if(MpiStatus)then
          call es_return_cvector(MpiComm,state_list,istate,state_cvec)
       else
          call es_return_cvector(state_list,istate,state_cvec)
       endif
#else
       call es_return_cvector(state_list,istate,state_cvec)
#endif
       !
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          do i = 1,sectorI%Dim
             gs_weight=peso*abs(state_cvec(i))**2
             !
             m    = sectorI%H(1)%map(i)
             ib   = bdecomp(m,2*Ns)
             do iorb=1,Nimp
                nup(iorb)= dble(ib(iorb))
                ndw(iorb)= dble(ib(iorb+Ns))
             enddo
             sz = (nup-ndw)/2d0
             nt =  nup+ndw
             !
             !
             !Evaluate averages of observables:
             do iorb=1,Nimp
                dens(iorb)     = dens(iorb)      +  nt(iorb)*gs_weight
                dens_up(iorb)  = dens_up(iorb)   +  nup(iorb)*gs_weight
                dens_dw(iorb)  = dens_dw(iorb)   +  ndw(iorb)*gs_weight
                docc(iorb)     = docc(iorb)      +  nup(iorb)*ndw(iorb)*gs_weight
                magz(iorb)     = magz(iorb)      +  (nup(iorb)-ndw(iorb))*gs_weight
                sz2(iorb,iorb) = sz2(iorb,iorb)  +  (sz(iorb)*sz(iorb))*gs_weight
                n2(iorb,iorb)  = n2(iorb,iorb)   +  (nt(iorb)*nt(iorb))*gs_weight
                do jorb=iorb+1,Nimp
                   sz2(iorb,jorb) = sz2(iorb,jorb)  +  (sz(iorb)*sz(jorb))*gs_weight
                   sz2(jorb,iorb) = sz2(jorb,iorb)  +  (sz(jorb)*sz(iorb))*gs_weight
                   n2(iorb,jorb)  = n2(iorb,jorb)   +  (nt(iorb)*nt(jorb))*gs_weight
                   n2(jorb,iorb)  = n2(jorb,iorb)   +  (nt(jorb)*nt(iorb))*gs_weight
                enddo
             enddo
             s2tot = s2tot  + (sum(sz))**2*gs_weight
             !            
          enddo


          !GET <(b_up + adg_dw)(bdg_up + a_dw)> = 
          !<b_up*bdg_up> + <adg_dw*a_dw> + <b_up*a_dw> + <adg_dw*bdg_up> = 
          !<n_a,dw> + < 1 - n_b,up> + 2*<PHI>_ab
          ! do ispin=1,Nspin 
          do iorb=1,Nimp !A
             do jorb=1,Nimp !B
                isz = getsz(isector)
                if(isz<Ns)then
                   jsector = getsector(isz+1,1)
                   vvinit  = apply_Cops(v_state,[one,one],[ 1, -1],[iorb,jorb],[1,2],jsector,sectorI)
                   phiscAB(iorb,jorb) = phiscAB(iorb,jorb) + dot_product(vvinit,vvinit)*peso
                   if(allocated(vvinit))deallocate(vvinit)
                endif
             enddo
          enddo
          ! enddo
          call delete_sector(isector,HI)
       endif
       !
       !
       if(allocated(state_cvec))deallocate(state_cvec)
       !
    enddo
    !
    !
    do iorb=1,Nimp
       do jorb=1,Nimp
          phiscAB(iorb,jorb) = 0.5d0*(phiscAB(iorb,jorb) - dens_dw(iorb) - (1.d0-dens_up(jorb)))
       enddo
       phisc(iorb)=phiscAB(iorb,iorb)
    enddo
    !
    !
    if(MPIMASTER)then
       if(iolegend)call write_legend
       call write_observables()
    endif
    write(LOGfile,"(A,50f18.12,f18.12,A)")"dens  "//reg(ed_file_suffix)//"=",((dens(iorb+(ilat-1)*Norb),iorb=1,Norb),ilat=1,Nlat),sum(dens)
    write(LOGfile,"(A,50f18.12,A)")       "docc  "//reg(ed_file_suffix)//"=",((docc(iorb+(ilat-1)*Norb),iorb=1,Norb),ilat=1,Nlat)
    if(Nspin==2)write(LOGfile,"(A,10f18.12,A)")    "magZ"//reg(ed_file_suffix)//"=",((magz(iorb+(ilat-1)*Norb),iorb=1,Norb),ilat=1,Nlat)
    write(LOGfile,"(A,100f18.12,A)")      "phiAB "//reg(ed_file_suffix)//"=",((phiscAB(iorb,jorb),iorb=1,Nimp),jorb=1,Nimp)
    write(LOGfile,"(A,100f18.12,A)")      " | phiAA*Uloc ",(abs(uloc(iorb))*phisc(iorb),iorb=1,Nimp)
    !
    do iorb=1,Nimp
       ed_dens_up(iorb)=dens_up(iorb)
       ed_dens_dw(iorb)=dens_dw(iorb)
       ed_dens(iorb)   =dens(iorb)
       ed_docc(iorb)   =docc(iorb)
       ed_mag(3,iorb)  =magZ(iorb)
       ed_phisc(iorb)  =phisc(iorb)
    enddo
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_dens_up)
       call Bcast_MPI(MpiComm,ed_dens_dw)
       call Bcast_MPI(MpiComm,ed_dens)
       call Bcast_MPI(MpiComm,ed_docc)
       call Bcast_MPI(MpiComm,ed_phisc)
       call Bcast_MPI(MpiComm,ed_mag)
       if(allocated(imp_density_matrix))call Bcast_MPI(MpiComm,imp_density_matrix)
    endif
#endif
    !
    deallocate(dens,docc,phiscAB,phisc,dens_up,dens_dw,magz,sz2,n2)
  end subroutine observables_superc






  !+-------------------------------------------------------------------+
  !PURPOSE  : Get internal energy from the Impurity problem.
  !+-------------------------------------------------------------------+
  subroutine local_energy_superc()
    !Calculate the values of the local observables
    integer,dimension(2*Ns)       :: ib
    integer,dimension(2,Ns)       :: Nud
    integer,dimension(Ns)         :: IbUp,IbDw
    real(8),dimension(Nimp)       :: nup,ndw
    real(8),dimension(Nspin,Nimp) :: eloc
    !
    Egs     = state_list%emin
    ed_Ehartree= 0.d0
    ed_Eknot   = 0.d0
    ed_Epot    = 0.d0
    ed_Dust    = 0.d0
    ed_Dund    = 0.d0
    ed_Dse     = 0.d0
    ed_Dph     = 0.d0
    !
    !
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          call es_return_cvector(MpiComm,state_list,istate,state_cvec)
       else
          call es_return_cvector(state_list,istate,state_cvec)
       endif
#else
       call es_return_cvector(state_list,istate,state_cvec)
#endif
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       !Master:
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          do i=1,sectorI%Dim
             m    = sectorI%H(1)%map(i)
             ib   = bdecomp(m,2*Ns)
             do iorb=1,Norb
                nup(iorb)=dble(ib(iorb))
                ndw(iorb)=dble(ib(iorb+Ns))
             enddo
             gs_weight=peso*abs(state_cvec(i))**2
             !
             !start evaluating the Tr(H_loc) to estimate potential energy
             !> H_Imp: Diagonal Elements, i.e. local part
             do iorb=1,Nimp
                ed_Eknot = ed_Eknot + impHloc(1,1,iorb,iorb)*Nup(iorb)*gs_weight
                ed_Eknot = ed_Eknot + impHloc(Nspin,Nspin,iorb,iorb)*Ndw(iorb)*gs_weight
             enddo
             ! !> H_imp: Off-diagonal elements, i.e. non-local part. 
             do iorb=1,Nimp
                do jorb=1,Nimp
                   !SPIN UP
                   Jcondition = &
                        (impHloc(1,1,iorb,jorb)/=zero) .AND. &
                        (ib(jorb)==1)                  .AND. &
                        (ib(iorb)==0)
                   if (Jcondition) then
                      call c(jorb,m,k1,sg1)
                      call cdg(iorb,k1,k2,sg2)
                      j=binary_search(sectorI%H(1)%map,k2)
                      ed_Eknot = ed_Eknot + impHloc(1,1,iorb,jorb)*sg1*sg2*state_cvec(i)*conjg(state_cvec(j))*peso
                   endif
                   !SPIN DW
                   Jcondition = &
                        (impHloc(Nspin,Nspin,iorb,jorb)/=zero) .AND. &
                        (ib(jorb+Ns)==1)                       .AND. &
                        (ib(iorb+Ns)==0)
                   if (Jcondition) then
                      call c(jorb+Ns,m,k1,sg1)
                      call cdg(iorb+Ns,k1,k2,sg2)
                      j =binary_search(sectorI%H(1)%map,k2)
                      ed_Eknot = ed_Eknot + impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2*state_cvec(i)*conjg(state_cvec(j))*peso
                   endif
                enddo
             enddo
             !
             !
             !DENSITY-DENSITY INTERACTION: SAME ORBITAL, OPPOSITE SPINS
             !Euloc=\sum=i U_i*(n_u*n_d)_i
             !ed_Epot = ed_Epot + dot_product(uloc,nup*ndw)*gs_weight
             do iorb=1,Nimp
                ed_Epot = ed_Epot + Uloc(iorb)*nup(iorb)*ndw(iorb)*gs_weight
             enddo
             !
             !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, OPPOSITE SPINS
             !Eust=\sum_ij Ust*(n_up_i*n_dn_j + n_up_j*n_dn_i)
             !    "="\sum_ij (Uloc - 2*Jh)*(n_up_i*n_dn_j + n_up_j*n_dn_i)
             if(Norb>1)then
                do iorb=1,Nimp
                   do jorb=iorb+1,Nimp
                      ed_Epot = ed_Epot + Ust*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*gs_weight
                      ed_Dust = ed_Dust + (nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*gs_weight
                   enddo
                enddo
             endif
             !
             !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, PARALLEL SPINS
             !Eund = \sum_ij Und*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             !    "="\sum_ij (Ust-Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             !    "="\sum_ij (Uloc-3*Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             if(Norb>1)then
                do iorb=1,Nimp
                   do jorb=iorb+1,Nimp
                      ed_Epot = ed_Epot + (Ust-Jh)*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*gs_weight
                      ed_Dund = ed_Dund + (nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*gs_weight
                   enddo
                enddo
             endif
             !
             !SPIN-EXCHANGE (S-E) TERMS
             !S-E: Jh *( c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up )  (i.ne.j) 
             if(Norb>1.AND.(Jx/=0d0.OR.Jp/=0d0))then
                do iorb=1,Nimp
                   do jorb=1,Nimp
                      Jcondition=((iorb/=jorb).AND.&
                           (ib(jorb)==1)      .AND.&
                           (ib(iorb+Ns)==1)   .AND.&
                           (ib(jorb+Ns)==0)   .AND.&
                           (ib(iorb)==0))
                      if(Jcondition)then
                         call c(jorb,m,k1,sg1)
                         call c(iorb+Ns,k1,k2,sg2)
                         call cdg(jorb+Ns,k2,k3,sg3)
                         call cdg(iorb,k3,k4,sg4)
                         j_el=binary_search(sectorI%H(1)%map,k4)
                         j   = j_el + (iph-1)*sectorI%DimEl
                         ed_Epot = ed_Epot + Jx*sg1*sg2*sg3*sg4*state_cvec(i)*conjg(state_cvec(j))*peso
                         ed_Dse  = ed_Dse  + sg1*sg2*sg3*sg4*state_cvec(i)*conjg(state_cvec(j))*peso
                      endif
                   enddo
                enddo
             endif
             !
             !
             !PAIR-HOPPING (P-H) TERMS
             !P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
             !P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
             if(Norb>1.AND.(Jx/=0d0.OR.Jp/=0d0))then
                do iorb=1,Nimp
                   do jorb=1,Nimp
                      Jcondition=((iorb/=jorb).AND.&
                           (ib(jorb)==1)      .AND.&
                           (ib(jorb+Ns)==1)   .AND.&
                           (ib(iorb+Ns)==0)   .AND.&
                           (ib(iorb)==0))
                      if(Jcondition)then
                         call c(jorb,m,k1,sg1)
                         call c(jorb+Ns,k1,k2,sg2)
                         call cdg(iorb+Ns,k2,k3,sg3)
                         call cdg(iorb,k3,k4,sg4)
                         j_el=binary_search(sectorI%H(1)%map,k4)
                         j   = j_el + (iph-1)*sectorI%DimEl
                         ed_Epot = ed_Epot + Jp*sg1*sg2*sg3*sg4*state_cvec(i)*conjg(state_cvec(j))*peso
                         ed_Dph  = ed_Dph  + sg1*sg2*sg3*sg4*state_cvec(i)*conjg(state_cvec(j))*peso
                      endif
                   enddo
                enddo
             endif
             !
             !HARTREE-TERMS CONTRIBUTION:
             if(hfmode)then               
                do iorb=1,Nimp
                   ed_Ehartree=ed_Ehartree - 0.5d0*uloc(iorb)*(nup(iorb)+ndw(iorb))*gs_weight + 0.25d0*uloc(iorb)*gs_weight
                enddo
                if(Norb>1)then
                   do iorb=1,Nimp
                      do jorb=iorb+1,Nimp
                         ed_Ehartree=ed_Ehartree - 0.5d0*Ust*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*gs_weight + 0.5d0*Ust*gs_weight
                         ed_Ehartree=ed_Ehartree - 0.5d0*(Ust-Jh)*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*gs_weight + 0.5d0*(Ust-Jh)*gs_weight
                      enddo
                   enddo
                endif
             endif
          enddo
          call delete_sector(isector,H)         
       endif
       !
       if(allocated(state_cvec))deallocate(state_cvec)
       !
    enddo
    !
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_Epot)
       call Bcast_MPI(MpiComm,ed_Eknot)
       call Bcast_MPI(MpiComm,ed_Ehartree)
       call Bcast_MPI(MpiComm,ed_Dust)
       call Bcast_MPI(MpiComm,ed_Dund)
    endif
#endif
    !
    ed_Epot = ed_Epot + ed_Ehartree
    !
    if(ed_verbose>=3)then
       write(LOGfile,"(A,10f18.12)")"<Hint>  =",ed_Epot
       write(LOGfile,"(A,10f18.12)")"<V>     =",ed_Epot-ed_Ehartree
       write(LOGfile,"(A,10f18.12)")"<E0>    =",ed_Eknot
       write(LOGfile,"(A,10f18.12)")"<Ehf>   =",ed_Ehartree    
       write(LOGfile,"(A,10f18.12)")"Dust    =",ed_Dust
       write(LOGfile,"(A,10f18.12)")"Dund    =",ed_Dund
       write(LOGfile,"(A)")" "
    endif
    !
    if(MPIMASTER)then
       call write_energy_info()
       call write_energy()
    endif
    !
    !
  end subroutine local_energy_superc









  !   !+-------------------------------------------------------------------+
  !   !PURPOSE  : Build the impurity density matrices
  !   !+-------------------------------------------------------------------+
  !   subroutine density_matrix_impurity()
  !     integer                             :: istate
  !     integer                             :: iUP,iDW,jUP,jDW
  !     integer                             :: IimpUp,IimpDw,JimpUp,JimpDw
  !     integer                             :: IbathUp,IbathDw,bUP,bDW
  !     integer                             :: lenBATHup,lenBATHdw
  !     integer,allocatable                 :: BATHup(:),BATHdw(:)
  !     integer                             :: Nup, Ndw
  !     integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws         ![1,Ns]-[Norb,1+Nbath]
  !     integer,dimension(Ns_Ud)            :: iDimUps,iDimDws   ![1]-[Norb]
  !     integer,dimension(Ns)               :: IbUp,IbDw         ![Ns]
  !     type(sector_map),dimension(2*Ns_Ud) :: HI                ![2]-[2*Norb]
  !     integer,dimension(2*Ns_Ud)          :: Indices,Jndices   ![2]-[2*Norb]
  !     integer                             :: Nud(2,Ns),iud(2),jud(2),is,js,io,jo
  !     !
  !     ! Here we build two different density matrices related to the impurity problem
  !     ! > The CLUSTER-Reduced density matrix: \rho_IMP = Tr_BATH(\rho) 
  !     ! > The SINGLE-PARTICLE-Reduced density matrix: \rho_sp = <C^+_a C_b> 
  !     !
  !     write(LOGfile,"(A)")"Get cluster density matrix:"
  !     !CLUSTER DENSITY MATRIX [\rho_IMP = Tr_BATH(\rho)]
  !     cluster_density_matrix=zero
  !     !
  !     if(MpiMaster) call start_timer()  
  !     do istate=1,state_list%size
  !        isector = es_return_sector(state_list,istate)
  !        Ei      = es_return_energy(state_list,istate)
  ! #ifdef _MPI
  !        if(MpiStatus)then
  !           call es_return_cvector(MpiComm,state_list,istate,state_cvec)
  !        else
  !           call es_return_cvector(state_list,istate,state_cvec)
  !        endif
  ! #else
  !        call es_return_cvector(state_list,istate,state_cvec)
  ! #endif
  !        !Finite temperature weighting factor
  !        peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
  !        peso = peso/zeta_function
  !        !
  !        !Sector sizes (per spin)
  !        call get_DimUp(isector,iDimUps)
  !        call get_DimDw(isector,iDimDws)
  !        iDimUp = product(iDimUps)
  !        iDimDw = product(iDimDws)
  !        !
  !        if(MpiMaster)then
  !           call build_sector(isector,HI,itrace=.true.)
  !           !
  !           do IimpUp=0,2**Nimp-1
  !              do JimpUp=0,2**Nimp-1
  !                 !
  !                 !Finding the unique bath states connecting IimpUp and JimpUp -> BATHup(:)
  !                 call sp_return_intersection(HI(1)%sp,IimpUp,JimpUp,BATHup,lenBATHup)
  !                 if(lenBATHup==0)cycle  !there are no bath states intersecting IimpUp,JimpUp
  !                 !
  !                 do IimpDw=0,2**Nimp-1
  !                    do JimpDw=0,2**Nimp-1
  !                       !
  !                       !Finding the unique bath states connecting IimpDw and JimpDw -> BATHdw(:)
  !                       call sp_return_intersection(HI(2)%sp,IimpDw,JimpDw,BATHdw,lenBATHdw)
  !                       if(lenBATHdw==0)cycle  !there are no bath states intersecting IimpDw,JimpDw
  !                       !------------------------------------------------------------------------------------
  !                       !Arrived here we have finally determined the {IbathUp,IbathDw} states to trace on.
  !                       !We shall use them to build back the Fock space states (IimpSIGMA,IbathSIGMA) and
  !                       !through those search the map to retrieve the formal labels i=(iUP,iDW), j=(jUP,jDW)
  !                       !------------------------------------------------------------------------------------
  !                       !=== >>> TRACE over bath states <<< =================================================
  !                       do bUP=1,lenBATHup
  !                          IbathUp = BATHup(bUP)
  !                          do bDW=1,lenBATHdw
  !                             IbathDw = BATHdw(bDW)
  !                             !-----------------------------------------------------
  !                             !Allowed spin Fock space Istates: 
  !                             !Iup = IimpUp +  2^Nimp * IbathUp
  !                             !Idw = IimpDw +  2^Nimp * IbathDw
  !                             !
  !                             !Corresponding sector indices per spin:
  !                             iUP= binary_search(HI(1)%map,IimpUp + 2**Nimp*IbathUp)
  !                             iDW= binary_search(HI(2)%map,IimpDw + 2**Nimp*IbathDw)
  !                             !
  !                             !Global sector index:
  !                             i  = iUP + (iDW-1)*iDimUp
  !                             !----------------------------------------------------- 
  !                             !Allowed spin Fock space Jstates: 
  !                             !Jup = JimpUp +  2^Nimp * IbathUp
  !                             !Jdw = JimpDw +  2^Nimp * IbathDw
  !                             !
  !                             !Corresponding sector jndices per spin:
  !                             jUP= binary_search(HI(1)%map,JimpUp + 2**Nimp*IbathUp)
  !                             jDW= binary_search(HI(2)%map,JimpDw + 2**Nimp*IbathDw)
  !                             !
  !                             !Global sector jndex:
  !                             j  = jUP + (jDW-1)*iDimUp
  !                             !-----------------------------------------------------
  !                             !Final composition for the impurity:
  !                             io = (IimpUp + 2**Nimp*IimpDw) + 1 
  !                             jo = (JimpUp + 2**Nimp*JimpDw) + 1
  !                             !-----------------------------------------------------------------
  !                             !(i,j)_th contribution to the (io,jo)_th element of \rho_IMP
  !                             cluster_density_matrix(io,jo) = cluster_density_matrix(io,jo) + &
  !                                  state_cvec(i)*conjg(state_cvec(j))*peso
  !                             !-----------------------------------------------------------------
  !                          enddo
  !                       enddo !=============================================================================
  !                       !
  !                    enddo
  !                 enddo
  !                 !
  !              enddo
  !           enddo
  !           !
  !           call delete_sector(isector,HI)       
  !        endif
  !        !
  !        if(allocated(state_cvec))deallocate(state_cvec)
  !        !
  !     enddo
  !     if(MpiMaster) call stop_timer()
  !     !
  !     !
  !     !
  !     !
  !     !
  !     !
  !     !
  !     !
  !     !
  !     !
  !     write(LOGfile,"(A)")"Get single-particle density matrix (no print)"
  !     !SINGLE PARTICLE DENSITY MATRIX [<C^+_a C_b>]
  !     do istate=1,state_list%size
  !        isector = es_return_sector(state_list,istate)
  !        Ei      = es_return_energy(state_list,istate)
  ! #ifdef _MPI
  !        if(MpiStatus)then
  !           call es_return_cvector(MpiComm,state_list,istate,state_cvec)
  !        else
  !           call es_return_cvector(state_list,istate,state_cvec)
  !        endif
  ! #else
  !        call es_return_cvector(state_list,istate,state_cvec)
  ! #endif
  !        !
  !        peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
  !        peso = peso/zeta_function
  !        !
  !        idim  = getdim(isector)
  !        call get_DimUp(isector,iDimUps)
  !        call get_DimDw(isector,iDimDws)
  !        iDimUp = product(iDimUps)
  !        iDimDw = product(iDimDws)
  !        !
  !        if(MpiMaster)then
  !           call build_sector(isector,HI)
  !           do i=1,iDim
  !              call state2indices(i,[iDimUps,iDimDws],Indices)
  !              do ii=1,Ns_Ud
  !                 mup = HI(ii)%map(Indices(ii))
  !                 mdw = HI(ii+Ns_Ud)%map(Indices(ii+Ns_ud))
  !                 Nups(ii,:) = Bdecomp(mup,Ns_Orb) ![Ns_Orb = Ns = Nlat*Norb*(Nbath+1) in CDMFT code]
  !                 Ndws(ii,:) = Bdecomp(mdw,Ns_Orb)
  !              enddo
  !              Nud(1,:) = Breorder(Nups) !Actually, they are already reordered in CDMFT code...
  !              Nud(2,:) = Breorder(Ndws) !...and breorder() corresponds to an identity: look in ED_SETUP!
  !              !
  !              !Diagonal densities
  !              do ilat=1,Nlat
  !                 do ispin=1,Nspin
  !                    do iorb=1,Norb
  !                       is = imp_state_index(ilat,iorb)
  !                       single_particle_density_matrix(ilat,ilat,ispin,ispin,iorb,iorb) = &
  !                            single_particle_density_matrix(ilat,ilat,ispin,ispin,iorb,iorb) + &
  !                            peso*nud(ispin,is)*conjg(state_cvec(i))*state_cvec(i)
  !                    enddo
  !                 enddo
  !              enddo
  !              !
  !              !Off-diagonal
  !              do ispin=1,Nspin
  !                 do ilat=1,Nlat
  !                    do jlat=1,Nlat
  !                       do iorb=1,Norb
  !                          do jorb=1,Norb
  !                             is = imp_state_index(ilat,iorb)
  !                             js = imp_state_index(jlat,jorb)
  !                             if((Nud(ispin,js)==1).and.(Nud(ispin,is)==0))then
  !                                iud(1) = HI(1)%map(Indices(1))
  !                                iud(2) = HI(2)%map(Indices(2))
  !                                call c(js,iud(ispin),r,sgn1)
  !                                call cdg(is,r,k,sgn2)
  !                                Jndices = Indices
  !                                Jndices(1+(ispin-1)*Ns_Ud) = &
  !                                     binary_search(HI(1+(ispin-1)*Ns_Ud)%map,k)
  !                                call indices2state(Jndices,[iDimUps,iDimDws],j)
  !                                !
  !                                single_particle_density_matrix(ilat,jlat,ispin,ispin,iorb,jorb) = &
  !                                     single_particle_density_matrix(ilat,jlat,ispin,ispin,iorb,jorb) + &
  !                                     peso*sgn1*state_cvec(i)*sgn2*conjg(state_cvec(j))
  !                             endif
  !                          enddo
  !                       enddo
  !                    enddo
  !                 enddo
  !              enddo
  !              !endif
  !              !
  !              !
  !           enddo
  !           call delete_sector(isector,HI)         
  !        endif
  !        !
  !        if(allocated(state_cvec))deallocate(state_cvec)
  !        !
  !     enddo
  !     ! 
  !     !
  !     !
  !     !
  !     !
  !     !
  !   end subroutine density_matrix_impurity










  !####################################################################
  !                        PRINTING ROUTINES
  !####################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : write legend, i.e. info about columns 
  !+-------------------------------------------------------------------+
  subroutine write_legend()
    integer :: unit,iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="observables_info.ed")
    write(unit,"(A1,90(A10,6X))")"#",&
         (reg(txtfy(iorb))//"dens_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(Norb+iorb))//"docc_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(2*Norb+iorb))//"nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(3*Norb+iorb))//"ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(4*Norb+iorb))//"mag_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(5*Norb+1))//"s2",&
         reg(txtfy(5*Norb+2))//"egs",&
         ((reg(txtfy(5*Norb+2+(iorb-1)*Norb+jorb))//"sz2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
         ((reg(txtfy((5+Norb)*Norb+2+(iorb-1)*Norb+jorb))//"n2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="parameters_info.ed")
    write(unit,"(A1,90(A14,1X))")"#","1xmu","2beta",&
         (reg(txtfy(2+iorb))//"U_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(2+Norb+1))//"U'",reg(txtfy(2+Norb+2))//"Jh"
    close(unit)
    !
    iolegend=.false.
  end subroutine write_legend

  subroutine write_custom_legend()
    integer :: unit,i
    unit = free_unit()
    open(unit,file="custom_observables_info.ed")
    write(unit,"(A1,90(A10,6X))")"#",(reg(txtfy(i))//reg(custom_o%item(i)%o_name),i=1,custom_o%N_filled)
    close(unit)
  end subroutine write_custom_legend


  subroutine write_energy_info()
    integer :: unit
    unit = free_unit()
    open(unit,file="energy_info.ed")
    write(unit,"(A1,90(A14,1X))")"#",&
         reg(txtfy(1))//"<Hi>",&
         reg(txtfy(2))//"<V>=<Hi-Ehf>",&
         reg(txtfy(3))//"<Eloc>",&
         reg(txtfy(4))//"<Ehf>",&
         reg(txtfy(5))//"<Dst>",&
         reg(txtfy(6))//"<Dnd>"
    close(unit)
  end subroutine write_energy_info


  !+-------------------------------------------------------------------+
  !PURPOSE  : write observables to file
  !+-------------------------------------------------------------------+
  subroutine write_observables()
    integer :: unit
    integer :: iorb,jorb,ispin,ilat
    unit = free_unit()
    do ilat=1,Nlat
       open(unit,file="observables_all"//reg(ed_file_suffix)//"_site"//str(ilat,3)//".ed",position='append')
       write(unit,"(90(F15.9,1X))")&
            (dens(iorb),iorb=1,Nimp),&
            (docc(iorb),iorb=1,Nimp),&
            (dens_up(iorb),iorb=1,Nimp),&
            (dens_dw(iorb),iorb=1,Nimp),&
            (magz(iorb),iorb=1,Nimp),&
            s2tot,egs,&
            ((sz2(iorb,jorb),jorb=1,Nimp),iorb=1,Nimp),&
            ((n2(iorb,jorb),jorb=1,Nimp),iorb=1,Nimp)
       close(unit)
    enddo
    !
    unit = free_unit()
    do ilat=1,Nlat
       open(unit,file="parameters_last"//reg(ed_file_suffix)//".ed")
       write(unit,"(90F15.9)")xmu,beta,(uloc(iorb),iorb=1,Norb),Ust,Jh,Jx,Jp
       close(unit)
    enddo
    !
    unit = free_unit()
    do ilat=1,Nlat
       open(unit,file="observables_last"//reg(ed_file_suffix)//"_site"//str(ilat,3)//".ed")
       write(unit,"(90(F15.9,1X))")&
            (dens(iorb),iorb=1,Nimp),&
            (docc(iorb),iorb=1,Nimp),&
            (dens_up(iorb),iorb=1,Nimp),&
            (dens_dw(iorb),iorb=1,Nimp),&
            (magz(iorb),iorb=1,Nimp),&
            s2tot(ilat),egs,&
            ((sz2(iorb,jorb),jorb=1,Nimp),iorb=1,Nimp),&
            ((n2(iorb,jorb),jorb=1,Nimp),iorb=1,Nimp)
       close(unit)
    enddo

    unit = free_unit()
    open(unit,file="Sz_ab_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(A)")"#a, b, Sz(a,b)"
    do iorb=1,Nimp
       do jorb=1,Nimp
          write(unit,"(4I15,F15.9)")iorb,jorb,sz2(iorb,jorb)
       enddo
    enddo
    close(unit)

    open(unit,file="N2_ab_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(A)")"#a, b, N2(a,b)"
    do iorb=1,Nimp
       do jorb=1,Nimp
          write(unit,"(4I15,F15.9)")iorb,jorb,n2(iorb,jorb)
       enddo
    enddo
    close(unit)
  end subroutine write_observables




  subroutine write_energy()
    integer :: unit
    unit = free_unit()
    open(unit,file="energy_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")ed_Epot,ed_Epot-ed_Ehartree,ed_Eknot,ed_Ehartree,ed_Dust,ed_Dund!,ed_Dse,ed_Dph
    close(unit)
  end subroutine write_energy



END MODULE ED_OBSERVABLES_SUPERC







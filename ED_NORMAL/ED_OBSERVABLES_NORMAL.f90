MODULE ED_OBSERVABLES_NORMAL
  !This module calculates a series of observables, and stores them in aptly named plain-text files. :f:var:`ed_mode` = :code:`normal`
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_INTEGRATE, only: quad
  USE SF_LINALG
  USE SF_TIMER
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_BATH
  USE ED_HAMILTONIAN_NORMAL
  implicit none
  private



  ! interface add_custom_observable
  !    module procedure :: add_custom_observable_local
  !    module procedure :: add_custom_observable_kdep
  ! end interface add_custom_observable
  !   public :: init_custom_observables
  ! public :: add_custom_observable
  ! public :: get_custom_observables
  ! public :: clear_custom_observables


  public :: observables_normal
  public :: local_energy_normal
  public :: density_matrix_normal


  logical,save                            :: iolegend=.true.
  real(8),dimension(:),allocatable        :: dens    ! orbital-resolved charge density
  real(8),dimension(:),allocatable        :: dens_up ! orbital-resolved spin-:math:`\uparrow` electron density
  real(8),dimension(:),allocatable        :: dens_dw ! orbital-resolved spin-:math:`\downarrow` electron density
  real(8),dimension(:),allocatable        :: docc    ! orbital-resolved double occupation
  real(8),dimension(:),allocatable        :: magz    ! orbital-resolved magnetization ( :code:`z` component )
  real(8),dimension(:,:),allocatable      :: n2      ! :math:`\langle n_{i} n_{j} \rangle` for i,j orbitals
  real(8),dimension(:,:),allocatable      :: sz2     ! :math:`\langle S^{z}_{i} S^{z}_{j} \rangle` for i,j orbitals
  real(8)                                 :: s2tot   ! :math:`\langle S_{z}^{2} \rangle`
  real(8)                                 :: Egs     ! Ground-state energy
  real(8)                                 :: Ei
  real(8)                                 :: integrationR
  complex(8),allocatable,dimension(:,:,:) :: sij
  complex(8),allocatable,dimension(:,:,:) :: Hk
  !
  integer                                 :: iorb,jorb,iorb1,jorb1
  integer                                 :: ispin,jspin
  integer                                 :: ibath,jbath
  integer                                 :: r,m,k,k1,k2,k3,k4
  integer                                 :: iup,idw
  integer                                 :: jup,jdw
  integer                                 :: mup,mdw
  real(8)                                 :: sgn,sgn1,sgn2,sg1,sg2,sg3,sg4
  real(8)                                 :: gs_weight
  !
  real(8)                                 :: peso
  real(8)                                 :: norm
  !
  integer                                 :: i,j,ii
  integer                                 :: isector,jsector
  integer                                 :: idim,idimUP,idimDW
  !
  real(8),dimension(:),allocatable        :: vvinit
  complex(8),dimension(:),allocatable     :: state_cvec
  logical                                 :: Jcondition
  !
  type(sector)                            :: sectorI,sectorJ



contains 



  !+-------------------------------------------------------------------+
  !PURPOSE  : Lanc method
  !+-------------------------------------------------------------------+
  subroutine observables_normal()
    integer                             :: istate,Nud(2,Ns),iud(2),jud(2),is,js
    integer,dimension(2*Ns_Ud)          :: Indices
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns_Ud)            :: iDimUps,iDimDws
    integer,dimension(Ns)               :: IbUp,IbDw  ![Ns]
    real(8),dimension(Nimp)             :: nup,ndw,Sz,nt

    !
    !LOCAL OBSERVABLES:
    allocate(dens(Nimp),dens_up(Nimp),dens_dw(Nimp))
    allocate(docc(Nimp))
    allocate(magz(Nimp),sz2(Nimp,Nimp),n2(Nimp,Nimp))
    !
    Egs     = state_list%emin
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0   
    magz    = 0.d0
    sz2     = 0.d0
    n2      = 0.d0
    s2tot   = 0.d0
    !
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       state_cvec  =  es_return_vector(state_list,istate)
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          do i = 1,sectorI%Dim
             !
             gs_weight=peso*abs(state_cvec(i))**2
             !
             call state2indices(i,[sectorI%DimUps,sectorI%DimDws],Indices)
             mup  = sectorI%H(1)%map(Indices(1))
             mdw  = sectorI%H(2)%map(Indices(2))
             IbUp = Bdecomp(mup,Ns_Orb) ![Ns_Orb = Ns = Nlat*Norb*(Nbath+1) in CDMFT code]
             IbDw = Bdecomp(mdw,Ns_Orb)
             nup  = IbUp(1:Nimp)
             ndw  = IbDw(1:Nimp)
             sz   = (nup-ndw)/2d0
             nt   =  nup+ndw
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
             !
          enddo
          call delete_sector(sectorI)
       endif
       !
       if(allocated(state_cvec))deallocate(state_cvec)
       !
    enddo
    !
    !
    if(MPIMASTER)then
       if(iolegend)call write_legend
       call write_observables()
       !
       write(LOGfile,"(A,10f18.12,f18.12)")"dens "//reg(ed_file_suffix)//"=",(dens(iorb),iorb=1,Nimp),sum(dens)
       write(LOGfile,"(A,10f18.12)")       "docc "//reg(ed_file_suffix)//"=",(docc(iorb),iorb=1,Nimp)
       if(Nspin==2)&
            write(LOGfile,"(A,10f18.12)")  "mag  "//reg(ed_file_suffix)//"=",(magz(iorb),iorb=1,Nimp)
       write(LOGfile,"(A)")" "
       !
       ed_dens_up=dens_up
       ed_dens_dw=dens_dw
       ed_dens   =dens
       ed_docc   =docc
    endif
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_dens_up)
       call Bcast_MPI(MpiComm,ed_dens_dw)
       call Bcast_MPI(MpiComm,ed_dens)
       call Bcast_MPI(MpiComm,ed_docc)
    endif
#endif
    !
    deallocate(dens,docc,dens_up,dens_dw,magz,sz2,n2)
  end subroutine observables_normal






  !+-------------------------------------------------------------------+
  !PURPOSE  : Get internal energy from the Impurity problem.
  !+-------------------------------------------------------------------+
  subroutine local_energy_normal()
    !Calculate the value of the local energy components
    integer                             :: istate,iud(2),jud(2),is,js
    integer,dimension(2*Ns_Ud)          :: Indices
    integer,dimension(Ns_Ud)            :: iDimUps,iDimDws
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    real(8),dimension(Ns)               :: Nup,Ndw
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
             call state2indices(i,[sectorI%DimUps,sectorI%DimDws],Indices)
             mup = sectorI%H(1)%map(Indices(1))
             mdw = sectorI%H(2)%map(Indices(2))
             Nups(1,:) = Bdecomp(mup,Ns_Orb) ![Ns_Orb = Ns = Nlat*Norb*(Nbath+1) in CDMFT code]
             Ndws(1,:) = Bdecomp(mdw,Ns_Orb)
             Nup = Breorder(Nups) !Actually, they are already reordered in CDMFT code...
             Ndw = Breorder(Ndws) !...and breorder() corresponds to an identity: look in ED_SETUP!
             !
             gs_weight=peso*abs(state_cvec(i))**2
             !
             !> H_Imp: Diagonal Elements, i.e. local part
             do iorb=1,Nimp
                ed_Eknot = ed_Eknot + impHloc(1,1,iorb,iorb)*Nup(iorb)*gs_weight
                ed_Eknot = ed_Eknot + impHloc(Nspin,Nspin,iorb,iorb)*Ndw(iorb)*gs_weight
             enddo
             ! !> H_imp: Off-diagonal elements, i.e. non-local part. 
             iup = Indices(1)
             idw = Indices(2)
             mup = sectorI%H(1)%map(iup)
             mdw = sectorI%H(2)%map(idw)
             do iorb=1,Nimp
                do jorb=1,Nimp
                   !UP
                   Jcondition = &
                        (impHloc(1,1,iorb,jorb)/=0d0) .AND. &
                        (Nup(jorb)==1) .AND. (Nup(iorb)==0)
                   if (Jcondition) then
                      call c(jorb,mup,k1,sg1)
                      call cdg(iorb,k1,k2,sg2)
                      jup = binary_search(sectorI%H(1)%map,k2)
                      jdw = idw
                      j   = jup + (jdw-1)*iDimUp
                      ed_Eknot = ed_Eknot + &
                           impHloc(1,1,iorb,jorb)*sg1*sg2*state_cvec(i)*conjg(state_cvec(j))*peso
                   endif
                   !
                   !DW
                   Jcondition = &
                        (impHloc(Nspin,Nspin,iorb,jorb)/=0d0) .AND. &
                        (ndw(jorb)==1) .AND. (ndw(iorb)==0)
                   if (Jcondition) then
                      call c(jorb,mdw,k1,sg1)
                      call cdg(iorb,k1,k2,sg2)
                      jdw = binary_search(sectorI%H(2)%map,k2)
                      jup = iup
                      j   = jup + (jdw-1)*iDimUp
                      ed_Eknot = ed_Eknot + &
                           impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2*state_cvec(i)*conjg(state_cvec(j))*peso
                   endif
                enddo
             enddo
             !
             !SPIN-EXCHANGE Jx
             if(Nimp>1.AND.Jx/=0d0)then
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
                         jdw=binary_search(sectorI%H(2)%map,k2)
                         call c(jorb,mup,k3,sg3)  !UP
                         call cdg(iorb,k3,k4,sg4) !UP
                         jup=binary_search(sectorI%H(1)%map,k4)
                         j = jup + (jdw-1)*sectorI%DimUp
                         !
                         ed_Epot = ed_Epot + Jx*sg1*sg2*sg3*sg4*state_cvec(i)*conjg(state_cvec(j))*peso
                         ed_Dse  = ed_Dse + sg1*sg2*sg3*sg4*state_cvec(i)*conjg(state_cvec(j))*peso
                         !
                      endif
                   enddo
                enddo
             endif
             !
             ! PAIR-HOPPING Jp
             if(Norb>1.AND.Jp/=0d0)then
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
                         jdw = binary_search(sectorI%H(2)%map,k2)
                         call c(jorb,mup,k3,sg3)       !c_jorb_up
                         call cdg(iorb,k3,k4,sg4)      !c^+_iorb_up
                         jup = binary_search(sectorI%H(1)%map,k4)
                         j = jup + (jdw-1)*sectorI%DimUp
                         !
                         ed_Epot = ed_Epot + Jp*sg1*sg2*sg3*sg4*state_cvec(i)*conjg(state_cvec(j))*peso
                         ed_Dph = ed_Dph + sg1*sg2*sg3*sg4*state_cvec(i)*conjg(state_cvec(j))*peso
                         !
                      endif
                   enddo
                enddo
             endif
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
             !HARTREE-TERMS CONTRIBUTION:
             if(hfmode)then
                !ed_Ehartree=ed_Ehartree - 0.5d0*dot_product(uloc,nup+ndw)*gs_weight + 0.25d0*sum(uloc)*gs_weight
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
          call delete_sector(sectorI)         
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
  end subroutine local_energy_normal







  !+-------------------------------------------------------------------+
  !PURPOSE  : Build the impurity density matrices
  !+-------------------------------------------------------------------+
  subroutine density_matrix_normal()
    integer                             :: istate
    integer                             :: iUP,iDW,jUP,jDW
    integer                             :: IimpUp,IimpDw,JimpUp,JimpDw
    integer                             :: IbathUp,IbathDw,bUP,bDW
    integer                             :: lenBATHup,lenBATHdw
    integer,allocatable                 :: BATHup(:),BATHdw(:)
    integer                             :: Nup, Ndw
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws         ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns_Ud)            :: iDimUps,iDimDws   ![1]-[Norb]
    integer,dimension(Ns)               :: IbUp,IbDw         ![Ns]
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices   ![2]-[2*Norb]
    integer                             :: Nud(2,Ns),iud(2),jud(2),is,js,io,jo
    !
    ! Here we build two different density matrices related to the impurity problem
    ! > The CLUSTER-Reduced density matrix: \rho_IMP = Tr_BATH(\rho) 
    ! > The SINGLE-PARTICLE-Reduced density matrix: \rho_sp = <C^+_a C_b> 
    !
    write(LOGfile,"(A)")"Get cluster density matrix:"
    !CLUSTER DENSITY MATRIX [\rho_IMP = Tr_BATH(\rho)]
    cluster_density_matrix=zero
    !
    if(MpiMaster) call start_timer()  
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
       !Finite temperature weighting factor
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI,itrace=.true.)
          !
          do IimpUp=0,2**Nimp-1
             do JimpUp=0,2**Nimp-1
                !
                !Finding the unique bath states connecting IimpUp and JimpUp -> BATHup(:)
                call sp_return_intersection(sectorI%H(1)%sp,IimpUp,JimpUp,BATHup,lenBATHup)
                if(lenBATHup==0)cycle  !there are no bath states intersecting IimpUp,JimpUp
                !
                do IimpDw=0,2**Nimp-1
                   do JimpDw=0,2**Nimp-1
                      !
                      !Finding the unique bath states connecting IimpDw and JimpDw -> BATHdw(:)
                      call sp_return_intersection(sectorI%H(2)%sp,IimpDw,JimpDw,BATHdw,lenBATHdw)
                      if(lenBATHdw==0)cycle  !there are no bath states intersecting IimpDw,JimpDw
                      !------------------------------------------------------------------------------------
                      !Arrived here we have finally determined the {IbathUp,IbathDw} states to trace on.
                      !We shall use them to build back the Fock space states (IimpSIGMA,IbathSIGMA) and
                      !through those search the map to retrieve the formal labels i=(iUP,iDW), j=(jUP,jDW)
                      !------------------------------------------------------------------------------------
                      !=== >>> TRACE over bath states <<< =================================================
                      do bUP=1,lenBATHup
                         IbathUp = BATHup(bUP)
                         do bDW=1,lenBATHdw
                            IbathDw = BATHdw(bDW)
                            !-----------------------------------------------------
                            !Allowed spin Fock space Istates: 
                            !Iup = IimpUp +  2^Nimp * IbathUp
                            !Idw = IimpDw +  2^Nimp * IbathDw
                            !
                            !Corresponding sector indices per spin:
                            iUP= binary_search(sectorI%H(1)%map,IimpUp + 2**Nimp*IbathUp)
                            iDW= binary_search(sectorI%H(2)%map,IimpDw + 2**Nimp*IbathDw)
                            !
                            !Global sector index:
                            i  = iUP + (iDW-1)*iDimUp
                            !----------------------------------------------------- 
                            !Allowed spin Fock space Jstates: 
                            !Jup = JimpUp +  2^Nimp * IbathUp
                            !Jdw = JimpDw +  2^Nimp * IbathDw
                            !
                            !Corresponding sector jndices per spin:
                            jUP= binary_search(sectorI%H(1)%map,JimpUp + 2**Nimp*IbathUp)
                            jDW= binary_search(sectorI%H(2)%map,JimpDw + 2**Nimp*IbathDw)
                            !
                            !Global sector jndex:
                            j  = jUP + (jDW-1)*iDimUp
                            !-----------------------------------------------------
                            !Final composition for the impurity:
                            io = (IimpUp + 2**Nimp*IimpDw) + 1 
                            jo = (JimpUp + 2**Nimp*JimpDw) + 1
                            !-----------------------------------------------------------------
                            !(i,j)_th contribution to the (io,jo)_th element of \rho_IMP
                            cluster_density_matrix(io,jo) = cluster_density_matrix(io,jo) + &
                                 state_cvec(i)*conjg(state_cvec(j))*peso
                            !-----------------------------------------------------------------
                         enddo
                      enddo !=============================================================================
                      !
                   enddo
                enddo
                !
             enddo
          enddo
          !
          call delete_sector(sectorI)       
       endif
       !
       if(allocated(state_cvec))deallocate(state_cvec)
       !
    enddo
    if(MpiMaster) call stop_timer()
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    write(LOGfile,"(A)")"Get single-particle density matrix (no print)"
    !SINGLE PARTICLE DENSITY MATRIX [<C^+_a C_b>]
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
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          do i=1,sectorI%Dim
             call state2indices(i,[iDimUps,iDimDws],Indices)
             mup = sectorI%H(1)%map(Indices(1))
             mdw = sectorI%H(2)%map(Indices(2))
             Nups(1,:) = Bdecomp(mup,Ns_Orb) ![Ns_Orb = Ns = Nlat*Norb*(Nbath+1) in CDMFT code]
             Ndws(1,:) = Bdecomp(mdw,Ns_Orb)
             Nud(1,:) = Breorder(Nups) !Actually, they are already reordered in CDMFT code...
             Nud(2,:) = Breorder(Ndws) !...and breorder() corresponds to an identity: look in ED_SETUP!
             !
             !Diagonal densities
             do ispin=1,Nspin
                do iorb=1,Nimp
                   single_particle_density_matrix(1,1,ispin,ispin,iorb,iorb) = &
                        single_particle_density_matrix(1,1,ispin,ispin,iorb,iorb) + &
                        peso*nud(ispin,iorb)*conjg(state_cvec(i))*state_cvec(i)
                enddo
             enddo
             !
             !Off-diagonal
             do ispin=1,Nspin
                do iorb=1,Nimp
                   do jorb=1,Nimp
                      if((Nud(ispin,jorb)==1).and.(Nud(ispin,iorb)==0))then
                         iud(1) = sectorI%H(1)%map(Indices(1))
                         iud(2) = sectorI%H(2)%map(Indices(2))
                         call c(jorb,iud(ispin),r,sgn1)
                         call cdg(iorb,r,k,sgn2)
                         Jndices = Indices
                         Jndices(1+(ispin-1)*Ns_Ud) = &
                              binary_search(sectorI%H(1+(ispin-1)*Ns_Ud)%map,k)
                         call indices2state(Jndices,[iDimUps,iDimDws],j)
                         !
                         single_particle_density_matrix(1,1,ispin,ispin,iorb,jorb) = &
                              single_particle_density_matrix(1,1,ispin,ispin,iorb,jorb) + &
                              peso*sgn1*state_cvec(i)*sgn2*conjg(state_cvec(j))
                      endif
                   enddo
                enddo
             enddo
             !
             !
          enddo
          call delete_sector(sectorI)         
       endif
       !
       if(allocated(state_cvec))deallocate(state_cvec)
       !
    enddo
    ! 
    !
    !
    !
    !
    !
  end subroutine density_matrix_normal











  ! !####################################################################
  ! !                    COMPUTATIONAL ROUTINES
  ! !####################################################################
  ! subroutine init_custom_observables(N,Hk)
  !   integer                      :: N
  !   complex(8),dimension(:,:,:)  :: Hk
  !   !
  !   if(MpiMaster)then
  !      custom_o%N_filled=0
  !      custom_o%N_asked=N
  !      allocate(custom_o%Hk(size(Hk,1),size(Hk,2),size(Hk,3)))
  !      custom_o%Hk=Hk
  !      allocate(custom_o%item(N))
  !      custom_o%init=.true.
  !   endif
  !   !
  ! end subroutine init_custom_observables

  ! subroutine add_custom_observable_local(o_name,sij)
  !   integer                               :: i
  !   complex(8),dimension(:,:)             :: sij
  !   character(len=*)                      :: o_name
  !   !
  !   if(MpiMaster .and. custom_o%init)then
  !      if(custom_o%N_filled .gt. custom_o%N_asked)then
  !         STOP "add_custom_observable: too many observables given"
  !         call clear_custom_observables
  !      endif
  !      !
  !      custom_o%N_filled=custom_o%N_filled+1
  !      custom_o%item(custom_o%N_filled)%o_name=o_name
  !      custom_o%item(custom_o%N_filled)%o_value=0.d0
  !      !
  !      allocate(custom_o%item(custom_o%N_filled)%sij(size(custom_o%Hk,1),size(custom_o%Hk,2),size(custom_o%Hk,3)))
  !      do i=1,size(custom_o%item(custom_o%N_filled)%sij,3)
  !         custom_o%item(custom_o%N_filled)%sij(:,:,i)=sij
  !      enddo
  !   else
  !      STOP "add_custom_observable: custom observables not initialized"
  !   endif
  ! end subroutine add_custom_observable_local


  ! subroutine add_custom_observable_kdep(o_name,sijk)
  !   integer                               :: i
  !   complex(8),dimension(:,:,:)           :: sijk
  !   character(len=*)                      :: o_name
  !   !
  !   if(MpiMaster .and. custom_o%init)then
  !      if(custom_o%N_filled .gt. custom_o%N_asked)then
  !         STOP "add_custom_observable: too many observables given"
  !         call clear_custom_observables
  !      endif
  !      !
  !      custom_o%N_filled=custom_o%N_filled+1
  !      custom_o%item(custom_o%N_filled)%o_name=o_name
  !      custom_o%item(custom_o%N_filled)%o_value=0.d0
  !      !
  !      allocate(custom_o%item(custom_o%N_filled)%sij(size(custom_o%Hk,1),size(custom_o%Hk,2),size(custom_o%Hk,3)))
  !      custom_o%item(custom_o%N_filled)%sij=sijk
  !   else
  !      STOP "add_custom_observable: custom observables not initialized"
  !   endif
  ! end subroutine add_custom_observable_kdep


  ! subroutine get_custom_observables()
  !   integer            :: i
  !   !
  !   if(MpiMaster .and. custom_o%init)then
  !      if(custom_o%N_filled .eq. 0)then
  !         write(Logfile,*)"WARNING! Custom observables initialized but none given."
  !         RETURN
  !      endif
  !      !
  !      write(LOGfile,*)"Calculating custom observables"
  !      !
  !      allocate(sij(size(custom_o%Hk,1),size(custom_o%Hk,2),size(custom_o%Hk,3)))
  !      allocate(Hk(size(custom_o%Hk,1),size(custom_o%Hk,2),size(custom_o%Hk,3)))
  !      sij=zero
  !      Hk=zero
  !      !
  !      Hk=custom_o%Hk
  !      do i=1,custom_o%N_filled
  !         sij=custom_o%item(i)%sij
  !         if(finiteT) then
  !            custom_o%item(i)%o_value=calculate_observable_integral_finite_t()
  !         else
  !            custom_o%item(i)%o_value=calculate_observable_integral_zero_t()
  !         endif
  !         write(LOGfile,"(A,10f18.12,A)")reg(custom_o%item(i)%o_name)//" = ",custom_o%item(i)%o_value
  !      enddo
  !      call write_custom_legend()
  !      call write_custom_observables()
  !      deallocate(sij,Hk)
  !   endif
  !   !
  ! end subroutine get_custom_observables


  ! subroutine clear_custom_observables()
  !   integer                       :: i
  !   if(MpiMaster .and. custom_o%init)then 
  !      do i=1,custom_o%N_filled
  !         deallocate(custom_o%item(i)%sij)
  !         custom_o%item(i)%o_name=""
  !         custom_o%item(i)%o_value=0.d0
  !      enddo
  !      deallocate(custom_o%Hk)
  !      custom_o%N_asked=0
  !      custom_o%N_filled=0
  !      custom_o%init=.false.
  !   endif
  ! end subroutine clear_custom_observables



  ! !+---------------------------------------------------------------------------------+
  ! !PURPOSE  : Evaluate and print out custom observable
  ! !+---------------------------------------------------------------------------------+

  ! !T=0
  ! function calculate_observable_integral_zero_t() result(out_2)
  !   integer                                                   :: inf
  !   real(8)                                                   :: out_2,spin_multiplicity
  !   !
  !   out_2=0.d0
  !   spin_multiplicity=3.d0-Nspin
  !   !
  !   call quad(sum_observable_kmesh,a=0.0d0,inf=1,verbose=(ED_VERBOSE>=3),result=out_2,strict=.false.)
  !   !
  !   out_2=spin_multiplicity*out_2/pi 
  !   return
  ! end function calculate_observable_integral_zero_t

  ! !T finite
  ! function calculate_observable_integral_finite_t() result(out_2)
  !   integer                         :: inf,Nmax,ii
  !   real(8)                         :: out_2,spin_multiplicity,omegamax,integralpart
  !   !
  !   !1) Find the real omegamax
  !   nmax=int(2*(abs(max_exc)+2.d0*hwband)*beta/pi)
  !   if (mod(nmax,2)==0)then
  !      nmax=nmax/2    
  !   else
  !      nmax=(nmax+1)/2
  !   endif
  !   integrationR=2*(nmax+1)*pi/beta
  !   !2) Evaluate discrete sum
  !   !
  !   out_2=0.d0
  !   do ii=0,Nmax
  !      out_2=out_2+dreal(sum_observable_kmesh_complex(xi*(2*ii+1)*pi/beta))
  !   enddo
  !   !
  !   out_2=2.d0*(1/beta)*out_2
  !   !
  !   !3) Evaluate integral part
  !   integralpart=0.d0
  !   call quad(integral_contour,a=-pi,b=pi,verbose=(ED_VERBOSE>=3),key=6,result=integralpart,strict=.false.)
  !   !
  !   !4) Sum all
  !   out_2=out_2+integralpart
  !   !5) Spin trick
  !   spin_multiplicity=3.d0-Nspin
  !   out_2=spin_multiplicity*out_2
  !   return
  ! end function calculate_observable_integral_finite_t


  ! !
  ! !
  ! !+-------------------------------------------------------------------+
  ! !PURPOSE  : complex integration
  ! !+-------------------------------------------------------------------+
  ! function integral_contour(zeta) result(f)
  !   real(8)                 :: zeta,f
  !   complex(8)              :: w,fermi
  !   !
  !   w=integrationR*exp(xi*zeta)
  !   if(dreal((w-XMU)*beta)>= 100)then
  !      fermi=0.d0
  !   else
  !      fermi=(1/(exp(beta*(w-XMU))+1))
  !   endif
  !   !
  !   f=dreal((1.d0/pi)*w*fermi*sum_observable_kmesh_complex(w))
  ! end function integral_contour
  ! !
  ! !
  ! !+-------------------------------------------------------------------+
  ! !PURPOSE  : sum on k-vectors
  ! !+-------------------------------------------------------------------+
  ! !
  ! function sum_observable_kmesh(omega) result(out_1)
  !   integer                                                           :: ii,jj,kk
  !   real(8)                                                           :: omega
  !   real(8)                                                           :: out_1
  !   complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)             :: g,invg0
  !   complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)             :: invg_lso,invg0_lso,sigma,Gk_lso
  !   !
  !   out_1=0.d0
  !   !
  !   g=zero
  !   invg0=zero
  !   invg_lso=zero
  !   invg0_lso=zero
  !   Gk_lso=zero
  !   Sigma=zero
  !   !
  !   !Obtain Sigma(iw)
  !   call ed_gf_cluster(dcmplx(0.d0,omega),g)
  !   invg_lso=nnn2lso_reshape(g,Nlat,Nspin,Norb)
  !   call inv(invg_lso)
  !   invg0=invg0_bath(dcmplx(0.d0,omega))
  !   invg0_lso=nnn2lso_reshape(invg0,Nlat,Nspin,Norb)
  !   Sigma=invg0_lso-invg_lso
  !   !
  !   !Do the k-sum
  !   do ii=1,size(Hk,3)
  !      Gk_lso=(dcmplx(0d0,omega)+xmu)*eye(Nlat*Nspin*Norb) - Hk(:,:,ii) - Sigma    
  !      call inv(Gk_lso)
  !      out_1=out_1+DREAL(trace(matmul(sij(:,:,ii),Gk_lso)) - trace(sij(:,:,ii))/(dcmplx(-1.1d0,omega)))
  !   enddo
  !   !
  !   out_1=out_1/(size(Hk,3))
  !   !
  !   return
  !   !
  ! end function sum_observable_kmesh
  ! !
  ! function sum_observable_kmesh_complex(omega) result(out_1)
  !   integer                                                           :: ii,jj,kk
  !   complex(8)                                                        :: omega
  !   complex(8)                                                        :: out_1
  !   complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)             :: g,invg0
  !   complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)             :: invg_lso,invg0_lso,sigma,Gk_lso
  !   !
  !   out_1=0.d0
  !   !
  !   g=zero
  !   invg0=zero
  !   invg_lso=zero
  !   invg0_lso=zero
  !   Gk_lso=zero
  !   Sigma=zero
  !   !
  !   !    
  !   call ed_gf_cluster(omega,g)
  !   invg_lso=nnn2lso_reshape(g,Nlat,Nspin,Norb)
  !   call inv(invg_lso)
  !   invg0=invg0_bath(omega)
  !   invg0_lso=nnn2lso_reshape(invg0,Nlat,Nspin,Norb)
  !   Sigma=invg0_lso-invg_lso
  !   !
  !   do ii=1,size(Hk,3)
  !      Gk_lso=(xi*omega+xmu)*eye(Nlat*Nspin*Norb) - Hk(:,:,ii)- Sigma    
  !      call inv(Gk_lso)
  !      out_1=out_1+DREAL(trace(matmul(sij(:,:,ii),Gk_lso)))
  !   enddo
  !   out_1=out_1/(size(Hk,3))
  !   !
  !   return
  !   !
  ! end function sum_observable_kmesh_complex





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
         (reg(txtfy(iorb))//"dens_"//reg(txtfy(iorb)),iorb=1,Nimp),&
         (reg(txtfy(Nimp+iorb))//"docc_"//reg(txtfy(iorb)),iorb=1,Nimp),&
         (reg(txtfy(2*Nimp+iorb))//"nup_"//reg(txtfy(iorb)),iorb=1,Nimp),&
         (reg(txtfy(3*Nimp+iorb))//"ndw_"//reg(txtfy(iorb)),iorb=1,Nimp),&
         (reg(txtfy(4*Nimp+iorb))//"mag_"//reg(txtfy(iorb)),iorb=1,Nimp),&
         reg(txtfy(5*Nimp+1))//"s2",&
         reg(txtfy(5*Nimp+2))//"egs",&
         ((reg(txtfy(5*Nimp+2+(iorb-1)*Nimp+jorb))//"sz2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Nimp),iorb=1,Nimp),&
         ((reg(txtfy((5+Nimp)*Nimp+2+(iorb-1)*Nimp+jorb))//"n2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Nimp),iorb=1,Nimp)
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

  ! subroutine write_custom_legend()
  !   integer :: unit,i
  !   unit = free_unit()
  !   open(unit,file="custom_observables_info.ed")
  !   write(unit,"(A1,90(A10,6X))")"#",(reg(txtfy(i))//reg(custom_o%item(i)%o_name),i=1,custom_o%N_filled)
  !   close(unit)
  ! end subroutine write_custom_legend


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
    integer :: iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="observables_all"//reg(ed_file_suffix)//".ed",position='append')
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
    !
    unit = free_unit()
    open(unit,file="parameters_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")xmu,beta,(uloc(iorb),iorb=1,Norb),Ust,Jh,Jx,Jp
    close(unit)
    !
    unit = free_unit()
    open(unit,file="observables_last"//reg(ed_file_suffix)//".ed")
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

    unit = free_unit()
    open(unit,file="Sz_ij_ab_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(A)")"#a, b, Sz(I,J,a,b)"
    do iorb=1,Nimp
       do jorb=1,Nimp
          write(unit,"(2I15,F15.9)")iorb,jorb,sz2(iorb,jorb)
       enddo
    enddo
    close(unit)

    open(unit,file="N2_ij_ab_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(A)")"#a, b, N2(I,J,a,b)"
    do iorb=1,Nimp
       do jorb=1,Nimp
          write(unit,"(2I15,F15.9)")iorb,jorb,n2(iorb,jorb)
       enddo
    enddo
    close(unit)
  end subroutine write_observables



  ! subroutine write_custom_observables()
  !   integer :: i
  !   integer :: unit
  !   unit = free_unit()
  !   open(unit,file="custom_observables_all.ed",position='append')
  !   write(unit,"(90(F15.9,1X))")&
  !        (custom_o%item(i)%o_value,i=1,custom_o%N_filled)
  !   close(unit)
  !   !
  !   unit = free_unit()
  !   open(unit,file="custom_observables_last.ed")
  !   write(unit,"(90(F15.9,1X))")&
  !        (custom_o%item(i)%o_value,i=1,custom_o%N_filled)
  !   close(unit)
  !   !
  ! end subroutine write_custom_observables




  subroutine write_energy()
    integer :: unit
    unit = free_unit()
    open(unit,file="energy_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")ed_Epot,ed_Epot-ed_Ehartree,ed_Eknot,ed_Ehartree,ed_Dust,ed_Dund,ed_Dse,ed_Dph
    close(unit)
  end subroutine write_energy



END MODULE ED_OBSERVABLES_NORMAL







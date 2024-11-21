MODULE ED_GF_SUPERC
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN_SUPERC
  !
  implicit none
  private



  public :: build_gf_superc
  public :: get_Gimp_superc
  public :: get_Sigma_superc


  integer                   :: istate
  integer                   :: isector,jsector
  integer                   :: idim,idimUP,idimDW
  integer                   :: jdim,jdimUP,jdimDW
  complex(8),allocatable    :: vvinit(:)
  real(8),allocatable       :: alfa_(:),beta_(:)
  integer                   :: ialfa,ibeta
  integer                   :: jalfa,jbeta
  integer                   :: r
  integer                   :: i,iup,idw
  integer                   :: j,jup,jdw  
  integer                   :: m,mup,mdw
  real(8)                   :: sgn,norm2,norm0
  integer                   :: Nitermax,Nlanc,vecDim


  !Lanczos shared variables
  !=========================================================
  complex(8),dimension(:),allocatable   :: v_state
  real(8)                               :: e_state



  !AUX GF
  !=========================================================
  complex(8),allocatable,dimension(:,:) :: auxGmats,auxGreal




  integer :: in,jn
  integer :: inam,jnam
  integer :: ilat,jlat
  integer :: iorb,jorb
  integer :: ispin,jspin
  integer :: is,js
  integer :: io,jo
  integer :: i,j

contains



  !+------------------------------------------------------------------+
  !                        SUPERC
  !+------------------------------------------------------------------+
  subroutine build_gf_superc()
    !
    !
    !Evaluates the impurity electrons Green's functions :math:`G(z)` and :math:`F(z)` using dynamical Lanczos method. The result is stored in rank-5 arrays :f:var:`impgmats`, :f:var:`impgreal` , :f:var:`impfmats` , :f:var:`impfreal` of dimensions [ |Nspin| , |Nspin| , |Norb| , |Norb| , :f:var:`Lmats` / :f:var:`Lreal` ] and rank-1 array :f:var:`impdmats`, :f:var:`impdreal`.    
    !
    !The off-diagonal components of :math:`G_{ab}` with :math:`a \neq b` as well as the anomalous Green's functions :math:`F_{ab}(z)\, \forall a,b` are obtained using algebraic manipulation to ensure working with hermitian conjugate operators in the dynamical Lanczos procedure.  
    !
    !The weights and the poles obtained in this procedure are saved in a hierarchical data structure (for every state, every channel (creation or annihilation of excitations, normal or anomalous) and every degree of freedom) :f:var:`impgmatrix` of type :f:var:`gfmatrix`. 
    !
    ! .. _j.cpc.2021.108261: https://doi.org/10.1016/j.cpc.2021.108261
    !
    integer    :: iorb,jorb,ispin,i,isign,in
    !
    !
    if(allocated(impGmatrix))deallocate(impGmatrix)
    allocate(impGmatrix(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp))
    !
    ! 
    max_exc =-huge(1d0)          !find the max excitation
    ispin=1
    !    
    if(MpiMaster)call start_timer
    !
    !get      G^{aa}_{uu,dd} (in impG(1,1,u/d,u/d,a,a)
    do ispin=1,Nspin
       do iorb=1,Nimp
          call GFmatrix_allocate(impGmatrix(1,1,ispin,ispin,iorb,iorb),Nstate=state_list%size)
          call lanc_build_gf_superc_Gdiag(ispin,iorb)
       enddo
    enddo
    !get \bar.G^{aa}_{dd,uu} (in impG(2,2,u/d,u/d,a,a)
    do ispin=1,Nspin
       do iorb=1,Nimp
          call GFmatrix_allocate(impGmatrix(2,2,ispin,ispin,iorb,iorb),Nstate=state_list%size)
          call lanc_build_gf_superc_barGdiag(ispin,iorb)
       enddo
    enddo
    !
    !
    !Get G^{ab}_{ss} --> saved in impG(1,1,s,s,a,b)
    do ispin=1,Nspin
       do iorb=1,Nimp
          do jorb=1,Nimp
             if(iorb==jorb)cycle
             call allocate_GFmatrix(impGmatrix(1,1,ispin,ispin,iorb,jorb),Nstate=state_list%size)
             call lanc_build_gf_superc_Gmix(1,ispin,iorb,jorb)
          enddo
       enddo
    enddo
    !we do not need \bar{G}^{ab} in impG(2,2,s,s,a,b)
    !
    !
    !get F^{ab}_{ud}  --> saved in impG(1,2,s,s,a,b)
    do ispin=1,Nspin
       do iorb=1,Nimp
          do jorb=1,Nimp
             call allocate_GFmatrix(impGmatrix(1,2,ispin,ispin,iorb,jorb),Nstate=state_list%size)
             call lanc_build_gf_superc_Fmix(1,ispin,iorb,jorb)
          enddo
       enddo
    enddo
    !we do not need \bar{F}^{ab} in impG(2,1,s,s,a,b)
    !
    if(MPIMASTER)call stop_timer()
    !
  end subroutine build_gf_superc






  !################################################################
  !################################################################
  !################################################################
  !################################################################






  !in=1 => get     G^{aa}_{u/d} = <CCdg>_u/d + <CdgC>_u/d
  subroutine lanc_build_gf_superc_Gdiag(ispin,iorb)
    integer      :: ispin,iorb
    integer      :: is,in
    type(sector) :: sectorI
    !
    if(MPIMASTER)call start_timer(unit=LOGfile)
    !
    in = 1
    is = ispin
    !
    if(ed_verbose>1)write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(iorb)//"_s"//str(ispin)//"_in"//str(in)
    !
    !
    do istate=1,state_list%size
       call allocate_GFmatrix(impGmatrix(in,in,is,is,iorb,iorb),istate,Nchan=2) !2*[is,is]
       !
       isector  =  es_return_sector(state_list,istate)
       e_state  =  es_return_energy(state_list,istate)
       v_state  =  es_return_vector(state_list,istate)
       !
       if(MpiMaster)call build_sector(isector,sectorI)
       !
       ! c^+_{in,iorb}|v> 
       jsector = getCDGsector(1,is,isector)
       if(jsector/=0)then
          vvinit =  apply_op_CDG(v_state,iorb,is,jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,1,1,istate,impGmatrix(in,in,is,is,iorb,iorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,in,is,is,iorb,iorb),istate,1,Nexc=0)
       endif
       !
       ! c_{up,iorb}|v> 
       jsector = getCsector(1,is,isector)
       if(jsector/=0)then
          vvinit =  apply_op_C(v_state,iorb,is,jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,-1,2,istate,impGmatrix(in,in,is,is,iorb,iorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,in,is,is,iorb,iorb),istate,2,Nexc=0)
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    !
    if(MPIMASTER)call stop_timer()
    !
    return
  end subroutine lanc_build_gf_superc_Gdiag



  !in=2 => get \barG^{aa}_{d/u} = <CdgC>_d/u + <CCdg>_d/u
  subroutine lanc_build_gf_superc_barGdiag(ispin,iorb)
    integer      :: ispin,iorb
    integer      :: in,is
    type(sector) :: sectorI
    !
    if(MPIMASTER)call start_timer(unit=LOGfile)
    !
    in = 2
    is = 3-ispin
    !
    if(ed_verbose>1)write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(iorb)//"_s"//str(ispin)//"_in"//str(in)
    !
    !
    do istate=1,state_list%size
       call allocate_GFmatrix(impGmatrix(in,in,is,is,iorb,iorb),istate,Nchan=2) !2*[is,is]
       !
       isector  =  es_return_sector(state_list,istate)
       e_state  =  es_return_energy(state_list,istate)
       v_state  =  es_return_vector(state_list,istate)
       !
       if(MpiMaster)call build_sector(isector,sectorI)
       !
       ! c_{dw,iorb}|v> 
       jsector = getCsector(1,is,isector)
       if(jsector/=0)then
          vvinit =  apply_op_C(v_state,iorb,is,jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,1,1,istate,impGmatrix(in,in,is,is,iorb,iorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,in,is,is,iorb,iorb),istate,1,Nexc=0)
       endif
       !
       ! c^+_{dw,iorb}|v> 
       jsector = getCDGsector(1,is,isector)
       if(jsector/=0)then
          vvinit =  apply_op_CDG(v_state,iorb,is,jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,-1,2,istate,impGmatrix(in,in,is,is,iorb,iorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,in,is,is,iorb,iorb),istate,2,Nexc=0)
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    !
    if(MPIMASTER)call stop_timer()
    !
    return
  end subroutine lanc_build_gf_superc_barGdiag



  !################################################################
  !################################################################
  !################################################################
  !################################################################


  !in=1      G^{ab}_{ss} --> saved in impG(1,1,s,s,a,b)
  subroutine lanc_build_gf_superc_Gmix(ispin,iorb,jorb)
    integer      :: ispin,iorb,jorb
    integer      :: in,is
    type(sector) :: sectorI
    !
    in = 1
    is = ispin
    !
    if(ed_verbose>1)write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(iorb)//"_s"//str(is)//"_in"//str(in)
    !
    do istate=1,state_list%size
       call allocate_GFmatrix(impGmatrix(in,in,is,is,iorb,jorb),istate,Nchan=4) !2*[QdgQ,RdgR]
       !
       isector  =  es_return_sector(state_list,istate)
       e_state  =  es_return_energy(state_list,istate)
       v_state  =  es_return_vector(state_list,istate)
       !
       if(MpiMaster)call build_sector(isector,sectorI)
       !
       ! Qdg = (cdg_iorb + cdg_jorb) --> <Q.Qdg> => <C.Cdg>
       jsector = getCDGsector(1,is,isector)
       if(jsector/=0)then 
          vvinit = apply_Cops(v_state,[one,one],[ 1, 1],[iorb,jorb],[is,is],jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,1,1,istate,impGmatrix(in,in,is,is,iorb,jorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,in,is,is,iorb,jorb),istate,1,Nexc=0)
       endif
       !
       ! Q  = (c_iorb + c_jorb)       --> <Qdg.Q> => <Cdg.C>
       jsector = getCsector(1,is,isector)
       if(jsector/=0)then 
          vvinit = apply_Cops(v_state,[one,one],[-1,-1],[iorb,jorb],[is,is],jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,-1,2,istate,impGmatrix(in,in,is,is,iorb,jorb))
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,in,is,is,iorb,jorb),istate,2,Nexc=0)
       endif
       !
       !
       !
       !Rdg = (cdg_iorb -xi cdg_jorb)  --> <R.Rdg> => -xi*<Cdg.C>
       jsector = getCDGsector(1,is,isector)
       if(jsector/=0)then 
          vvinit = apply_Cops(v_state,[one,-xi],[ 1, 1],[iorb,jorb],[is,is],jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,e_state,alfa_,beta_,1,3,istate,impGmatrix(in,in,is,is,iorb,jorb))
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,in,is,is,iorb,jorb),istate,3,Nexc=0)
       endif
       !
       ! R = (c_iorb + xi c_jorb)       --> <Rdg.R> => -xi*<C.Cdg>
       jsector = getCsector(1,is,isector)
       if(jsector/=0)then 
          vvinit = apply_Cops(v_state,[one, xi],[-1,-1],[iorb,jorb],[is,is],jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,e_state,alfa_,beta_,-1,4,istate,impGmatrix(in,in,is,is,iorb,jorb))
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,in,is,is,iorb,jorb),istate,4,Nexc=0)
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    return
  end subroutine lanc_build_gf_superc_Gmix









  !in=2 \bar.G^{ab}_{ss} --> saved in impG(2,2,s,s,a,b)
  subroutine lanc_build_gf_superc_barGmix(ispin,iorb,jorb)
    integer      :: ispin,iorb,jorb
    integer      :: is,in
    type(sector) :: sectorI
    !
    in = 2
    is = 3-ispin
    !
    if(ed_verbose>1)write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(iorb)//"_s"//str(ispin)//"_in"//str(in)
    !
    !
    do istate=1,state_list%size
       call allocate_GFmatrix(impGmatrix(in,in,is,is,iorb,jorb),istate,Nchan=4) !2*[QdgQ,RdgR]
       !
       isector  =  es_return_sector(state_list,istate)
       e_state  =  es_return_energy(state_list,istate)
       v_state  =  es_return_vector(state_list,istate)
       !
       if(MpiMaster)call build_sector(isector,sectorI)
       !
       ! Qdg = (c_iorb + c_jorb)     --> <Q.Qdg> => <C.Cdg>
       jsector = getCsector(1,is,isector)
       if(jsector/=0)then 
          vvinit = apply_Cops(v_state,[one,one],[-1,-1],[iorb,jorb],[is,is],jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,1,1,istate,impGmatrix(in,in,is,is,iorb,jorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,in,is,is,iorb,jorb),istate,1,Nexc=0)
       endif
       !
       ! Q  = (cdg_iorb + cdg_jorb)  --> <Qdg.Q> => <Cdg.C>
       jsector = getCDGsector(1,is,isector)
       if(jsector/=0)then 
          vvinit = apply_Cops(v_state,[one,one],[ 1, 1],[iorb,jorb],[is,is],jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,-1,2,istate,impGmatrix(in,in,is,is,iorb,jorb))
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,in,is,is,iorb,jorb),istate,2,Nexc=0)
       endif
       !
       !
       !
       !Rdg = (c_iorb +xi c_jorb)    --> <R.Rdg> => -xi*<Cdg.C>
       jsector = getCsector(1,is,isector)
       if(jsector/=0)then 
          vvinit = apply_Cops(v_state,[one, xi],[-1,-1],[iorb,jorb],[is,is],jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,e_state,alfa_,beta_,1,3,istate,impGmatrix(in,in,is,is,iorb,jorb))
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,in,is,is,iorb,jorb),istate,3,Nexc=0)
       endif
       !
       ! R = (cdg_iorb - xi cdg_jorb) --> <Rdg.R> => -xi*<C.Cdg>
       jsector = getCDGsector(1,is,isector)
       if(jsector/=0)then 
          vvinit = apply_Cops(v_state,[one,-xi],[ 1, 1],[iorb,jorb],[is,is],jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,e_state,alfa_,beta_,-1,4,istate,impGmatrix(in,in,is,is,iorb,jorb))
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,in,is,is,iorb,jorb),istate,4,Nexc=0)
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    return
  end subroutine lanc_build_gf_superc_barGmix




  !################################################################
  !################################################################
  !################################################################
  !################################################################






  !in = 1, jn = 2:      F^{ab}_{ud} = <Cdg_aup.Cdg_bdw> --> saved in impG(1,2,s,s,a,b)
  subroutine lanc_build_gf_superc_Fmix(ispin,iorb,jorb)
    integer      :: ispin,iorb,jorb
    integer      :: is,js,in,jn
    type(sector) :: sectorI
    !
    if(ed_verbose>1)write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(iorb)//"_s"//str(ispin)//"_in"//str(in)
    !
    in = 1
    jn = 2
    !
    is = ispin 
    js = 3-ispin
    !
    do istate=1,state_list%size
       call allocate_GFmatrix(impGmatrix(in,jn,is,is,iorb,jorb),istate,Nchan=4)
       !
       isector  =  es_return_sector(state_list,istate)
       e_state  =  es_return_energy(state_list,istate)
       v_state  =  es_return_vector(state_list,istate)
       !
       if(MpiMaster)call build_sector(isector,sectorI)
       !
       !Qdg: [cdg_{up,iorb} + c_{dw,jorb}]  --> <Q.Qdg>
       isz = getsz(isector)
       if(isz<Ns)then
          jsector = getsector(isz+1,1)
          vvinit  = apply_Cops(v_state,[one,one],[ 1, -1],[iorb,jorb],[is,js],jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,1,1,istate,impGmatrix(in,jn,is,is,iorb,jorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,jn,is,is,iorb,jorb),istate,1,Nexc=0)
       endif
       !
       !Q  : [c_{up,iorb} + cdg_{dw,jorb}]   --> <Qdg.Q>
       isz = getsz(isector)
       if(isz>-Ns)then
          jsector = getsector(isz-1,1)
          vvinit  = apply_Cops(v_state,[one,one],[-1, 1],[iorb,jorb],[is,js],jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,-1,2,istate,impGmatrix(in,jn,is,is,iorb,jorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,jn,is,is,iorb,jorb),istate,2,Nexc=0)
       endif
       !
       !
       !P  : [cdg_{up,iorb} - xi*c_{dw,jorb}] --> -xi*<P.Pdg>
       isz = getsz(isector)
       if(isz<Ns)then
          jsector = getsector(isz+1,1)
          vvinit  = apply_Cops(v_state,[one,-xi],[ 1,-1],[iorb,jorb],[is,js],jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,e_state,alfa_,beta_,1,3,istate,impGmatrix(in,jn,is,is,iorb,jorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,jn,is,is,iorb,jorb),istate,3,Nexc=0)
       endif
       !
       !EVALUATE [c_{up,iorb} + xi*cdg_{dw,jorb}]|gs>  += -xi*<Pdg.P>
       isz = getsz(isector)
       if(isz>-Ns)then
          jsector = getsector(isz-1,1)
          vvinit  = apply_Cops(v_state,[one, xi],[-1, 1],[iorb,jorb],[is,js],jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,e_state,alfa_,beta_,-1,4,istate,impGmatrix(in,jn,is,is,iorb,jorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,jn,is,is,iorb,jorb),istate,4,Nexc=0)
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    !
    return
  end subroutine lanc_build_gf_superc_Fmix




  !in = 2, jn = 1: \bar.F^{ab}_{ud} = <C_adw.C_bup>     --> saved in impG(2,1,s,s,a,b)
  subroutine lanc_build_gf_superc_barFmix(ispin,iorb,jorb)
    integer      :: ispin,iorb,jorb
    integer      :: is,js,in,jn
    type(sector) :: sectorI
    !
    if(ed_verbose>1)write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(iorb)//"_s"//str(ispin)//"_in"//str(in)
    !
    in = 2
    jn = 1
    !
    is = 3-ispin 
    js = ispin
    !
    do istate=1,state_list%size
       call allocate_GFmatrix(impGmatrix(in,jn,is,is,iorb,jorb),istate,Nchan=4)
       !
       isector  =  es_return_sector(state_list,istate)
       e_state  =  es_return_energy(state_list,istate)
       v_state  =  es_return_vector(state_list,istate)
       !
       if(MpiMaster)call build_sector(isector,sectorI)
       !
       !Qdg: [c_{dw,iorb} + cdg_{up,jorb}]  --> <Q.Qdg>
       isz = getsz(isector)
       if(isz<Ns)then
          jsector = getsector(isz+1,1)
          vvinit  = apply_Cops(v_state,[one,one],[-1, 1],[iorb,jorb],[is,js],jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,1,1,istate,impGmatrix(in,jn,is,is,iorb,jorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,jn,is,is,iorb,jorb),istate,1,Nexc=0)
       endif
       !
       !Q  : [cdg_{dw,iorb} + c_{up,jorb}]   --> <Qdg.Q>
       isz = getsz(isector)
       if(isz>-Ns)then
          jsector = getsector(isz-1,1)
          vvinit  = apply_Cops(v_state,[one,one],[ 1,-1],[iorb,jorb],[is,js],jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,-1,2,istate,impGmatrix(in,jn,is,is,iorb,jorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,jn,is,is,iorb,jorb),istate,2,Nexc=0)
       endif
       !
       !
       !P  : [c_{up,iorb} - xi*cdg_{dw,jorb}] --> -xi*<P.Pdg>
       isz = getsz(isector)
       if(isz<Ns)then
          jsector = getsector(isz+1,1)
          vvinit  = apply_Cops(v_state,[one,-xi],[-1, 1],[iorb,jorb],[is,js],jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,e_state,alfa_,beta_,1,3,istate,impGmatrix(in,jn,is,is,iorb,jorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,jn,is,is,iorb,jorb),istate,3,Nexc=0)
       endif
       !
       !Pdg: [cdg_{up,iorb} + xi*c_{dw,jorb}] --> -xi*<Pdg.P>
       isz = getsz(isector)
       if(isz>-Ns)then
          jsector = getsector(isz-1,1)
          vvinit  = apply_Cops(v_state,[one, xi],[ 1,-1],[iorb,jorb],[is,js],jsector,sectorI)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,e_state,alfa_,beta_,-1,4,istate,impGmatrix(in,jn,is,is,iorb,jorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(in,jn,is,is,iorb,jorb),istate,4,Nexc=0)
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    !
    return
  end subroutine lanc_build_gf_superc_barFmix




  !################################################################
  !################################################################
  !################################################################
  !################################################################




  subroutine add_to_lanczos_gf_superc(vnorm2,Ei,As,Bs,isign,ic,istate,self)
    complex(8)                           :: vnorm2,pesoBZ,peso
    real(8)                              :: Ei,Egs,de
    integer                              :: nlanc,itype
    real(8),dimension(:)                 :: As
    real(8),dimension(size(As))          :: Bs 
    integer                              :: isign,ic,istate
    type(GFmatrix)                       :: self
    real(8),dimension(size(As),size(As)) :: Z
    real(8),dimension(size(As))          :: diag,subdiag
    integer                              :: i,j,ierr
    complex(8)                           :: iw
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(As)
    !
    pesoBZ = vnorm2/zeta_function
    if(finiteT)pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,As)
       call Bcast_MPI(MpiComm,Bs)
    endif
#endif
    !
    diag(1:Nlanc)    = As(1:Nlanc)
    subdiag(2:Nlanc) = Bs(2:Nlanc)
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    call allocate_GFmatrix(Self,istate,ic,Nlanc)
    !
    do j=1,nlanc
       de = diag(j)-Ei
       peso = pesoBZ*Z(1,j)*Z(1,j)
       !
       Self%state(istate)%channel(ic)%weight(j) = peso
       Self%state(istate)%channel(ic)%poles(j)  = isign*de
    enddo
  end subroutine add_to_lanczos_gf_superc



  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################


  function get_Gimp_superc(zeta) result(Gf)
    complex(8),dimension(:),intent(in)                                 :: zeta
    character(len=:),optional                                          :: axis
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(zeta)) :: Gf
    complex(8)                                                         :: green(size(zeta))
    complex(8)                                                         :: barG(Nspin,Norb,size(zeta))
    integer                                                            :: ispin
    integer                                                            :: ilat,jlat
    integer                                                            :: iorb,jorb
    integer                                                            :: iexc,Nexc
    integer                                                            :: ichan,Nchannel,istate,Nstates
    integer                                                            :: i,is,js,L
    complex(8)                                                         :: weight,de
    character(len=1)                                                   :: axis_
    !
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    !
    if(.not.allocated(impGmatrix))stop "ed_gf_cluster ERROR: impGmatrix not allocated!"
    !
    Gf = zero
    !
    !First get G^{ab}_{uu} and F^{ab}_{ud}
    do concurrent(ispin=1:Nspin,iorb=1:Nimp,jorb=1:Nimp,jn=1:2)
       green   = zero
       Nstates = size(impGmatrix(1,jn,ispin,ispin,iorb,jorb)%state)
       do istate=1,Nstates
          Nchannel = size(impGmatrix(1,jn,ispin,ispin,iorb,jorb)%state(istate)%channel)
          do ichan=1,Nchannel
             Nexc  = size(impGmatrix(1,jn,ispin,ispin,iorb,jorb)%state(istate)%channel(ichan)%poles)
             if(Nexc .ne. 0)then
                do iexc=1,Nexc
                   weight = impGmatrix(1,jn,ispin,ispin,iorb,jorb)%state(istate)%channel(ichan)%weight(iexc)
                   de     = impGmatrix(1,jn,ispin,ispin,iorb,jorb)%state(istate)%channel(ichan)%poles(iexc)
                   green  = green + weight/(zeta-de)
                enddo
             endif
          enddo
       enddo
       Gf(1,jn,ispin,ispin,iorb,jorb,:) = green
    enddo
    !
    !get aux \bar.G^{aa}_{dd}
    do concurrent(ispin=1:Nspin,iorb=1:Nimp)
       green   = zero
       Nstates = size(impGmatrix(2,2,ispin,ispin,iorb,iorb)%state)
       do istate=1,Nstates
          Nchannel = size(impGmatrix(2,2,ispin,ispin,iorb,iorb)%state(istate)%channel)
          do ichan=1,Nchannel
             Nexc  = size(impGmatrix(2,2,ispin,ispin,iorb,iorb)%state(istate)%channel(ichan)%poles)
             if(Nexc .ne. 0)then
                do iexc=1,Nexc
                   weight = impGmatrix(2,2,ispin,ispin,iorb,iorb)%state(istate)%channel(ichan)%weight(iexc)
                   de     = impGmatrix(2,2,ispin,ispin,iorb,iorb)%state(istate)%channel(ichan)%poles(iexc)
                   green  = green + weight/(zeta-de)
                enddo
             endif
          enddo
       enddo
       barG(ispin,iorb,:) = green
    enddo
    !
    !
    ! 
    !G^{ab}
    forall(ispin=1:Nspin,iorb=1:Nimp,jorb=1:Nimp,iorb/=jorb)
       gf(1,1,ispin,ispin,iorb,jorb,:) = 0.5d0*(gf(1,1,ispin,ispin,iorb,jorb,:) &
            - (one-xi)*gf(1,1,ispin,ispin,iorb,iorb,:) - (one-*xi)*gf(1,1,ispin,ispin,jorb,jorb,:))  
    end forall
    !
    !F^{ab}
    forall(ispin=1:Nspin,iorb=1:Nimp,jorb=1:Nimp)
       gf(1,2,ispin,ispin,iorb,jorb,:) = 0.5d0*( gf(1,2,ispin,ispin,iorb,jorb,:) &
            - (one-xi)*gf(1,1,ispin,ispin,iorb,iorb,:) - (one-*xi)*barG(ispin,jorb,:))
    end forall
    !
    !Enforce Symmetry:
    L = size(zeta)
    select case(axis_)
    case default;stop "get_Gimp_superc error: axis != [m,r]"
    case("m","M")
       do i=1,size(zeta)
          gf(2,1,:,:,:,:,i) =  conjg(gf(1,2,:,:,:,:,i)) !conjg(traspose(nn2so()))
          gf(2,2,:,:,:,:,i) = -conjg(gf(1,1,:,:,:,:,i))
       enddo
    case("r","R")
       do i=1,size(zeta)
          gf(2,1,:,:,:,:,i) =  conjg(gf(1,2,:,:,:,:,i))
          gf(2,2,:,:,:,:,i) = -conjg(gf(1,1,:,:,:,:,L-i+1))
       enddo
    end select
    !
  end function get_Gimp_superc







  function get_Sigma_superc(zeta,axis) result(Sigma)
    complex(8),dimension(:),intent(in)                                 :: zeta
    character(len=:),optional                                          :: axis
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(zeta)) :: Sigma,iG0,iG
    complex(8),dimension(Nambu*Nspin*Nimp,Nambu*Nspin*Nimp)            :: iTmp
    character(len=1)                                                   :: axis_
    !
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    !
    !
    !Get G0^-1
    iG0 = invg0_bath(zeta,dmft_bath,axis_)
    !Get G
    iG   = get_Gimp_superc(zeta,axis_)
    !
    do i=1,size(zeta)
       !Get G^-1
       iTmp  = nnn2nso_reshape(iG(:,:,:,:,:,:,i))
       call inv(iTmp)
       iG(:,:,:,:,:,:,i) = nso2nnn_reshape(iTmp)
       !
       !Get Sigma= G0^-1 - G^-1
       Sigma(1,1,:,:,:,;,i) = iG0(1,1,:,:,:,;,i) - iG(1,1,:,:,:,:,i)
       Sigma(1,2,:,:,:,;,i) = iG0(1,2,:,:,:,;,i) - iG(1,2,:,:,:,:,i)
    enddo
    !
    L = size(zeta)
    select case(axis_)
    case default;stop "get_Sigma_superc error: axis != [m,r]"
    case("m","M")
       do i=1,size(zeta)
          Sigma(2,1,:,:,:,:,i) =  conjg(Sigma(1,2,:,:,:,:,i)) !conjg(traspose(nn2so()))
          Sigma(2,2,:,:,:,:,i) = -conjg(Sigma(1,1,:,:,:,:,i))
       enddo
    case("r","R")
       do i=1,size(zeta)
          Sigma(2,1,:,:,:,:,i) =  conjg(Sigma(1,2,:,:,:,:,i))
          Sigma(2,2,:,:,:,:,i) = -conjg(Sigma(1,1,:,:,:,:,L-i+1))
       enddo
    end select
    !
  end function get_Sigma_superc




END MODULE ED_GF_SUPERC


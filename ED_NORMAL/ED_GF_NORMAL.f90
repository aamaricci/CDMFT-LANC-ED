MODULE ED_GF_NORMAL
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
  USE ED_HAMILTONIAN_NORMAL
  !
  implicit none
  private



  interface get_Gimp_normal
     module procedure :: get_Gimp_normal_array
  end interface get_Gimp_normal

  interface get_Sigma_normal
     module procedure :: get_Sigma_normal_array
  end interface get_Sigma_normal


  public :: build_gf_normal
  public :: get_Gimp_normal
  public :: get_Sigma_normal


  integer                               :: istate
  integer                               :: isector,jsector
  integer                               :: idim,idimUP,idimDW
  integer                               :: jdim,jdimUP,jdimDW
  complex(8),allocatable                :: vvinit(:)
  real(8),allocatable                   :: alfa_(:),beta_(:)
  integer                               :: r,m
  real(8)                               :: sgn,norm2,norm0
  integer                               :: Nitermax,Nlanc,vecDim
  integer                               :: in,jn
  integer                               :: inam,jnam
  integer                               :: ilat,jlat
  integer                               :: iorb,jorb
  integer                               :: ispin,jspin
  integer                               :: is,js
  integer                               :: io,jo
  integer                               :: i,j
  integer                               :: iup,idw
  integer                               :: jup,jdw  
  integer                               :: mup,mdw

  !Lanczos shared variables
  !=========================================================
  complex(8),dimension(:),allocatable   :: v_state
  real(8)                               :: e_state





contains



  !+------------------------------------------------------------------+
  !                        NORMAL
  !+------------------------------------------------------------------+
  subroutine build_gf_normal()
    !
    !Evaluates the impurity electrons Green's function :math:`G(z)` and the phonons one :math:`D(z)` using dynamical Lanczos method. The result is stored in rank-5 arrays :f:var:`impgmats`, :f:var:`impgreal` of dimensions [ |Nspin| , |Nspin| , |Norb| , |Norb| , :f:var:`Lmats` / :f:var:`Lreal` ] and rank-1 array :f:var:`impdmats`, :f:var:`impdreal`.    
    !
    !The off-diagonal components of :math:`G_{ab}` with :math:`a \neq b` are obtained using algebraic manipulation, see `j.cpc.2021.108261`_. 
    !
    !The weights and the poles obtained in this procedure are saved in a hierarchical data structure (for every state, every channel (creation or annihilation of excitations) and every degree of freedom) :f:var:`impgmatrix` of type :f:var:`gfmatrix`. 
    !
    ! .. _j.cpc.2021.108261: https://doi.org/10.1016/j.cpc.2021.108261
    !
    integer :: counter
    real(8) :: chan4
    !
    if(ed_gf_symmetric)then
       chan4=0.d0
    else
       chan4=1.d0
    endif
    !
    if(allocated(impGmatrix))deallocate(impGmatrix)
    allocate(impGmatrix(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp))
    !
    ! 
    max_exc =-huge(1d0)          !find the max excitation
    !
    if(MpiMaster)call start_timer
    !
    do ispin=1,Nspin
       do iorb=1,Nimp
          call allocate_GFmatrix(impGmatrix(1,1,ispin,ispin,iorb,iorb),Nstate=state_list%size) !2= C,Cdg
          call lanc_build_gf_normal_main(iorb,ispin)
          !
          do jorb=1,Nimp
             if(iorb==jorb)cycle
             call allocate_GFmatrix(impGmatrix(1,1,ispin,ispin,iorb,jorb),Nstate=state_list%size)!4=(Cdg_i + Cdg_j)/(CDG_i +iCDG_j)|psi>
             if(ed_gf_symmetric)then
                call lanc_build_gf_normal_mix_chan2(iorb,jorb,ispin)
             else
                call lanc_build_gf_normal_mix_chan4(iorb,jorb,ispin)
             endif
          enddo
       enddo
    enddo
    !
    if(MPIMASTER)call stop_timer()
    !
  end subroutine build_gf_normal






  !################################################################
  !################################################################
  !################################################################
  !################################################################






  subroutine lanc_build_gf_normal_main(iorb,ispin)
    integer,intent(in)          :: iorb,ispin
    type(sector)                :: sectorI
    !
    if(ed_verbose>1)write(LOGfile,*)"Get G_l"//str(iorb,3)//"_m"//str(iorb,3)//"_s"//str(ispin)
    !
    do istate=1,state_list%size
       call allocate_GFmatrix(impGmatrix(1,1,ispin,ispin,iorb,iorb),istate=istate,Nchan=2) !2= C,CDG
       !
       isector  =  es_return_sector(state_list,istate)
       e_state  =  es_return_energy(state_list,istate)
       v_state  =  es_return_vector(state_list,istate)
       !
       if(MpiMaster)call build_sector(isector,sectorI)
       !
       !ADD ONE PARTICLE:
       jsector = getCDGsector(1,ispin,isector)
       if(jsector/=0)then
          vvinit =  apply_op_CDG(v_state,iorb,ispin,jsector,sectorI)
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,e_state,alfa_,beta_,1,1,istate,impGmatrix(1,1,ispin,ispin,iorb,iorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,ispin,ispin,iorb,iorb),istate,1,Nexc=0)
       endif

       !REMOVE ONE PARTICLE:
       jsector = getCsector(1,ispin,isector)
       if(jsector/=0)then
          vvinit =  apply_op_C(v_state,iorb,ispin,jsector,sectorI)
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,e_state,alfa_,beta_,-1,2,istate,impGmatrix(1,1,ispin,ispin,iorb,iorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,ispin,ispin,iorb,iorb),istate,2,Nexc=0)
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    return
  end subroutine lanc_build_gf_normal_main







  subroutine lanc_build_gf_normal_mix_chan2(iorb,jorb,ispin)
    integer                     :: iorb,jorb,ispin
    type(sector)                :: sectorI
    !
    if(ed_verbose>1)write(LOGfile,*)"Get G_l"//str(iorb,3)//"_m"//str(jorb,3)//"_s"//str(ispin)
    !    
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(impGmatrix(1,1,ispin,ispin,iorb,jorb),istate=istate,Nchan=2) !2= add,del exc. c^+_i|psi> 
       isector  =  es_return_sector(state_list,istate)
       e_state  =  es_return_energy(state_list,istate)
       v_state  =  es_return_vector(state_list,istate)
       !
       !
       if(MpiMaster)call build_sector(isector,sectorI)
       !
       !EVALUATE (c^+_is + c^+_js)|gs> = [1,1].[C_{+1},C_{+1}].[i,j]
       jsector = getCDGsector(1,ispin,isector)
       if(jsector/=0)then
          !
          vvinit = apply_Cops(v_state,[one,one],[1,1],[iorb,jorb],[ispin,ispin],jsector,sectorI)
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,e_state,alfa_,beta_,1,1,istate,impGmatrix(1,1,ispin,ispin,iorb,jorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,ispin,ispin,iorb,jorb),istate=istate,ichan=1,Nexc=0)
       endif
       !
       !
       !EVALUATE (c_is + c_js)|gs> = [1,1].[C_{-1},C_{-1}].[i,j]
       jsector = getCsector(1,ispin,isector)
       if(jsector/=0)then
          !
          vvinit =  apply_Cops(v_state,[one,one],[-1,-1],[iorb,jorb],[ispin,ispin],jsector,sectorI)
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,e_state,alfa_,beta_,-1,2,istate,impGmatrix(1,1,ispin,ispin,iorb,jorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,ispin,ispin,iorb,jorb),istate=istate,ichan=2,Nexc=0)
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    return
  end subroutine lanc_build_gf_normal_mix_chan2






  subroutine lanc_build_gf_normal_mix_chan4(iorb,jorb,ispin)
    integer                     :: iorb,jorb,ispin
    type(sector)                :: sectorI
    !
    if(ed_verbose>1)write(LOGfile,*)"Get G_l"//str(iorb,3)//"_m"//str(jorb,3)//"_s"//str(ispin)
    !
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(impGmatrix(1,1,ispin,ispin,iorb,jorb),istate=istate,Nchan=4)
       isector  =  es_return_sector(state_list,istate)
       e_state  =  es_return_energy(state_list,istate)
       v_state  =  es_return_vector(state_list,istate)
       !
       if(MpiMaster)call build_sector(isector,sectorI)
       !
       !EVALUATE (c^+_is + c^+_js)|gs>
       jsector = getCDGsector(1,ispin,isector)
       if(jsector/=0)then
          !
          vvinit =  apply_Cops(v_state,[one,one],[1,1],[iorb,jorb],[ispin,ispin],jsector,sectorI)
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,e_state,alfa_,beta_,1,1,istate,impGmatrix(1,1,ispin,ispin,iorb,jorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,ispin,ispin,iorb,jorb),istate=istate,ichan=1,Nexc=0)
       endif
       !
       !
       !EVALUATE (c_is + c_js)|gs>
       jsector = getCsector(1,ispin,isector)
       if(jsector/=0)then
          !
          vvinit =  apply_Cops(v_state,[one,one],[-1,-1],[iorb,jorb],[ispin,ispin],jsector,sectorI)
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,e_state,alfa_,beta_,-1,2,istate,impGmatrix(1,1,ispin,ispin,iorb,jorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,ispin,ispin,iorb,jorb),istate=istate,ichan=2,Nexc=0)
       endif
       !
       !
       !EVALUATE (c^+_is + i*c^+_js)|gs>
       jsector = getCDGsector(1,ispin,isector)
       if(jsector/=0)then
          !
          vvinit =  apply_Cops(v_state,[one,xi],[1,1],[iorb,jorb],[ispin,ispin],jsector,sectorI)
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(-xi*norm2,e_state,alfa_,beta_,1,3,istate,impGmatrix(1,1,ispin,ispin,iorb,jorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,ispin,ispin,iorb,jorb),istate=istate,ichan=3,Nexc=0)
       endif
       !
       !
       !EVALUATE (c_is - i*c_js)|gs>
       jsector = getCsector(1,ispin,isector)
       if(jsector/=0)then
          !
          vvinit =  apply_Cops(v_state,[one,-xi],[-1,-1],[iorb,jorb],[ispin,ispin],jsector,sectorI)
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(-xi*norm2,e_state,alfa_,beta_,-1,4,istate,impGmatrix(1,1,ispin,ispin,iorb,jorb))
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,ispin,ispin,iorb,jorb),istate=istate,ichan=4,Nexc=0)
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    return
  end subroutine lanc_build_gf_normal_mix_chan4





  !################################################################
  !################################################################
  !################################################################
  !################################################################




  subroutine add_to_lanczos_gf_normal(vnorm2,Ei,As,Bs,isign,ichan,istate,self)
    complex(8)                           :: vnorm2,pesoBZ,peso
    real(8)                              :: Ei,Egs,de
    integer                              :: nlanc,itype
    real(8),dimension(:)                 :: As
    real(8),dimension(size(As))          :: Bs 
    integer                              :: isign,ichan,istate
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
    if((finiteT).and.(beta*(Ei-Egs) < 40))then
       pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    elseif(.not.finiteT)then
       pesoBZ = vnorm2/zeta_function
    else
       pesoBZ=0.d0
    endif
    !
    !
#ifdef _MPI
    if(MpiStatus)then
       if(MpiComm /= MPI_COMM_NULL)call Bcast_MPI(MpiComm,As)
       if(MpiComm /= MPI_COMM_NULL)call Bcast_MPI(MpiComm,Bs)
    endif
#endif
    diag             = 0.d0
    subdiag          = 0.d0
    Z                = eye(Nlanc)
    diag(1:Nlanc)    = As(1:Nlanc)
    subdiag(2:Nlanc) = Bs(2:Nlanc)
    !
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    call allocate_GFmatrix(Self,istate=istate,ichan=ichan,Nexc=Nlanc)
    !
    do j=1,nlanc
       de = diag(j)-Ei
       if (de>max_exc)max_exc=de
       peso = pesoBZ*Z(1,j)*Z(1,j)
       !
       Self%state(istate)%channel(ichan)%weight(j) = peso
       Self%state(istate)%channel(ichan)%poles(j)  = isign*de
       !
    enddo
  end subroutine add_to_lanczos_gf_normal



  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################

  function get_Gimp_normal_array(zeta,axis) result(Gf)
    complex(8),dimension(:),intent(in)                                 :: zeta
    character(len=*),optional                                          :: axis
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(zeta)) :: Gf
    complex(8),dimension(size(zeta))                                   :: green
    integer                                                            :: ispin
    integer                                                            :: ilat,jlat
    integer                                                            :: iorb,jorb
    integer                                                            :: iexc,Nexc
    integer                                                            :: ichan,Nchannel,istate,Nstates
    integer                                                            :: i,is,js
    complex(8)                                                         :: weight,de
    real(8)                                                            :: chan4
    character(len=1)                                                   :: axis_
    !
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1) !only for self-consistency, not used here
    !
    if(.not.allocated(impGmatrix))stop "ed_gf_cluster ERROR: impGmatrix not allocated!"
    !
    Gf=zero
    !
    chan4=1d0 ; if(ed_gf_symmetric)chan4=0d0
    !
    in = 1
    jn = 1
    do concurrent(ispin=1:Nspin,iorb=1:Nimp,jorb=1:Nimp)
       green   = zero
       Nstates = size(impGmatrix(in,jn,ispin,ispin,iorb,jorb)%state)
       do istate=1,Nstates
          Nchannel = size(impGmatrix(in,jn,ispin,ispin,iorb,jorb)%state(istate)%channel)
          do ichan=1,Nchannel
             Nexc  = size(impGmatrix(in,jn,ispin,ispin,iorb,jorb)%state(istate)%channel(ichan)%poles)
             if(Nexc .ne. 0)then
                do iexc=1,Nexc
                   weight = impGmatrix(in,jn,ispin,ispin,iorb,jorb)%state(istate)%channel(ichan)%weight(iexc)
                   de     = impGmatrix(in,jn,ispin,ispin,iorb,jorb)%state(istate)%channel(ichan)%poles(iexc)
                   green  = green + weight/(zeta-de)
                enddo
             endif
          enddo
       enddo
       Gf(1,1,ispin,ispin,iorb,jorb,:) = green
    enddo
    !
    do ispin=1,Nspin
       do iorb=1,Nimp
          do jorb=1,Nimp
             if(iorb==jorb)cycle
             gf(1,1,ispin,ispin,iorb,jorb,:) = 0.5d0*(gf(1,1,ispin,ispin,iorb,jorb,:) &
                  - (one-chan4*xi)*gf(1,1,ispin,ispin,iorb,iorb,:) - (one-chan4*xi)*gf(1,1,ispin,ispin,jorb,jorb,:))  
          enddo
       enddo
    enddo
    !
  end function get_Gimp_normal_array






  function get_Sigma_normal_array(zeta,axis) result(Sigma)
    complex(8),dimension(:),intent(in)                                 :: zeta
    character(len=*),optional                                          :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(zeta)) :: Sigma,invG0,invG
    complex(8),dimension(Nambu*Nspin*Nimp,Nambu*Nspin*Nimp)            :: iGzeta
    character(len=4)                                                   :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    !Get G0^-1
    invG0 = invg0_bath_function(zeta,dmft_bath,axis_)
    !
    !Get G^-1
    invG  = get_Gimp_normal_array(zeta)
    do i=1,size(zeta)
       iGzeta  = nnn2nso_reshape( invG(:,:,:,:,:,:,i) ) !rank6->rank2
       call inv(iGzeta)
       invG(:,:,:,:,:,:,i) = nso2nnn_reshape(iGzeta) !rank2->rank6
    enddo
    !
    !Get Sigma= G0^-1 - G^-1
    Sigma = invG0 - invG
    !
  end function get_Sigma_normal_array





END MODULE ED_GF_NORMAL


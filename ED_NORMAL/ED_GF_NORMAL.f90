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
     module procedure :: get_Gimp_normal_scalar
     module procedure :: get_Gimp_normal_array
  end interface get_Gimp_normal

  interface get_Sigma_normal
     module procedure :: get_Sigma_normal_scalar
     module procedure :: get_Sigma_normal_array
  end interface get_Sigma_normal


  public :: build_gf_normal
  public :: get_Gimp_normal
  public :: get_Sigma_normal


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
  complex(8),dimension(:),allocatable   :: state_cvec
  real(8)                               :: state_e



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
    counter = 0
    do ispin=1,Nspin
       do iorb=1,Nimp
          !diagonalo 
          call GFmatrix_allocate(impGmatrix(1,1,ispin,ispin,iorb,iorb),Nstate=state_list%size) !2= C,Cdg
          call lanc_build_gf_normal_main(iorb,ispin)
          !
          counter=counter+1
          if((ED_VERBOSE>1).AND.(MpiMaster))call eta(counter,Nspin*Nimp*Nimp)
          !
          !site-off-diagonal:
          do jn=1,Nambu      !==1
             do jorb=1,Norb
                if(iorb==jorb)cycle
                call GFmatrix_allocate(impGmatrix(1,1,ispin,ispin,iorb,jorb),Nstate=state_list%size)!4=(Cdg_i + Cdg_j)/(CDG_i +iCDG_j)|psi>
                if(ed_gf_symmetric)then
                   call lanc_build_gf_normal_mix_chan2(iorb,jorb,ispin)
                else
                   call lanc_build_gf_normal_mix_chan4(iorb,jorb,ispin)
                endif
                counter=counter+1
                if((ED_VERBOSE>1).and.(MpiMaster))call eta(counter,Nspin*Nimp*Nimp)
             enddo
          enddo
       enddo
    enddo
    !
    ! !nondiagonal trick
    ! do iorb=1,Norb
    !    do jorb=1,Norb
    !       if(iorb==jorb)cycle
    !       impGmats(1,1,ispin,ispin,iorb,jorb,:) = 0.5d0*(impGmats(1,1,ispin,ispin,iorb,jorb,:) &
    !            - (one-chan4*xi)*impGmats(1,1,ispin,ispin,iorb,iorb,:) - (one-chan4*xi)*impGmats(1,1,ispin,ispin,jorb,jorb,:))
    !       impGreal(1,1,ispin,ispin,iorb,jorb,:) = 0.5d0*(impGreal(1,1,ispin,ispin,iorb,jorb,:) &
    !            - (one-chan4*xi)*impGreal(isite,isite,ispin,ispin,iorb,iorb,:) - (one-chan4*xi)*impGreal(jsite,jsite,ispin,ispin,jorb,jorb,:))
    !    enddo
    ! enddo
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
    if(ed_verbose>1)write(LOGfile,*)"Get G_cluster_I"//str(iorb,3)//"_J"//str(iorb,3)
    !
    do istate=1,state_list%size
       call GFmatrix_allocate(impGmatrix(1,1,ispin,ispin,iorb,iorb),istate=istate,Nchan=2) !2= C,CDG
       !
       isector  =  es_return_sector(state_list,istate)
       state_e  =  es_return_energy(state_list,istate)
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
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>2)write(LOGfile,"(A20,I6,20I4)")&
               'From sector',isector,sectorI%Nups,sectorI%Ndws
       endif
       !
       !ADD ONE PARTICLE:
       jsector = getCDGsector(1,ispin,isector)
       if(jsector/=0)then
          vvinit =  apply_op_CDG(state_cvec,iorb,ispin,jsector,sectorI)
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,iorb,iorb,ispin,1,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,ispin,ispin,iorb,iorb),istate,1,Nexc=0)
       endif

       !REMOVE ONE PARTICLE:
       jsector = getCsector(1,ispin,isector)
       if(jsector/=0)then
          vvinit =  apply_op_C(state_cvec,iorb,ispin,jsector,sectorI)
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,iorb,iorb,ispin,2,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,ispin,ispin,iorb,iorb),istate,2,Nexc=0)
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
       if(allocated(state_cvec))deallocate(state_cvec)
       !
    enddo
    return
  end subroutine lanc_build_gf_normal_main







  subroutine lanc_build_gf_normal_mix_chan2(iorb,jorb,ispin)
    integer                     :: iorb,jorb,ispin
    type(sector)                :: sectorI
    !
    if(ed_verbose > 1)write(LOGfile,*)"Solving G_cluster_I"//str(iorb,3)//"_J"//str(jorb,3)
    !    
    do istate=1,state_list%size
       !
       call GFmatrix_allocate(impGmatrix(1,1,ispin,ispin,iorb,jorb),istate=istate,Nchan=2) !2= add,del exc. c^+_i|psi> 
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
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
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A20,I6,20I4)")&
               'From sector',isector,sectorI%Nups,sectorI%Ndws
       endif
       !
       !
       !EVALUATE (c^+_is + c^+_js)|gs> = [1,1].[C_{+1},C_{+1}].[i,j]
       jsector = getCDGsector(1,ispin,isector)
       if(jsector/=0)then
          !
          vvinit = apply_Cops(state_cvec,[one,one],[1,1],[iorb,jorb],[ispin,ispin],jsector,sectorI)
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,iorb,jorb,ispin,1,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call GFmatrix_allocate(impGmatrix(1,1,ispin,ispin,iorb,jorb),istate=istate,ichan=1,Nexc=0)
       endif
       !
       !
       !EVALUATE (c_is + c_js)|gs> = [1,1].[C_{-1},C_{-1}].[i,j]
       jsector = getCsector(1,ispin,isector)
       if(jsector/=0)then
          !
          vvinit =  apply_Cops(state_cvec,[one,one],[-1,-1],[iorb,jorb],[ispin,ispin],jsector,sectorI)
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin,2,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call GFmatrix_allocate(impGmatrix(1,1,ispin,ispin,iorb,jorb),istate=istate,ichan=2,Nexc=0)
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
       if(allocated(state_cvec))deallocate(state_cvec)
       !
    enddo
    return
  end subroutine lanc_build_gf_normal_mix_chan2






  subroutine lanc_build_gf_normal_mix_chan4(iorb,jorb,ispin)
    integer                     :: iorb,jorb,ispin
    type(sector)                :: sectorI
    !
    if(ed_verbose > 1)write(LOGfile,*)"Solving G_cluster_I"//str(iorb,3)//"_J"//str(jorb,3)
    !
    do istate=1,state_list%size
       !
       call GFmatrix_allocate(impGmatrix(1,1,ispin,ispin,iorb,jorb),istate=istate,Nchan=4)
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
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
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A20,I6,20I4)")&
               'From sector',isector,sectorI%Nups,sectorI%Ndws
       endif
       !
       !
       !EVALUATE (c^+_is + c^+_js)|gs>
       jsector = getCDGsector(1,ispin,isector)
       if(jsector/=0)then
          !
          vvinit =  apply_Cops(state_cvec,[one,one],[1,1],[iorb,jorb],[ispin,ispin],jsector,sectorI)
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,iorb,jorb,ispin,1,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call GFmatrix_allocate(impGmatrix(1,1,ispin,ispin,iorb,jorb),istate=istate,ichan=1,Nexc=0)
       endif
       !
       !
       !EVALUATE (c_is + c_js)|gs>
       jsector = getCsector(1,ispin,isector)
       if(jsector/=0)then
          !
          vvinit =  apply_Cops(state_cvec,[one,one],[-1,-1],[iorb,jorb],[ispin,ispin],jsector,sectorI)
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin,2,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call GFmatrix_allocate(impGmatrix(1,1,ispin,ispin,iorb,jorb),istate=istate,ichan=2,Nexc=0)
       endif
       !
       !
       !EVALUATE (c^+_is + i*c^+_js)|gs>
       jsector = getCDGsector(1,ispin,isector)
       if(jsector/=0)then
          !
          vvinit =  apply_Cops(state_cvec,[one,xi],[1,1],[iorb,jorb],[ispin,ispin],jsector,sectorI)
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(-xi*norm2,state_e,alfa_,beta_,1,iorb,jorb,ispin,3,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call GFmatrix_allocate(impGmatrix(1,1,ispin,ispin,iorb,jorb),istate=istate,ichan=3,Nexc=0)
       endif
       !
       !
       !EVALUATE (c_is - i*c_js)|gs>
       jsector = getCsector(ialfa,ispin,isector)
       if(jsector/=0)then
          !
          vvinit =  apply_Cops(state_cvec,[one,-xi],[-1,-1],[iorb,jorb],[ispin,ispin],jsector,sectorI)
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(-xi*norm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin,4,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call GFmatrix_allocate(impGmatrix(1,1,ispin,ispin,iorb,jorb),istate=istate,ichan=4,Nexc=0)
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
       if(allocated(state_cvec))deallocate(state_cvec)
       !
    enddo
    return
  end subroutine lanc_build_gf_normal_mix_chan4





  !################################################################
  !################################################################
  !################################################################
  !################################################################




  subroutine add_to_lanczos_gf_normal(vnorm2,Ei,alanc,blanc,isign,iorb,jorb,ispin,ichan,istate)
    complex(8)                                 :: vnorm2,pesoBZ,peso
    real(8)                                    :: Ei,Egs,de
    integer                                    :: nlanc,itype
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign,iorb,jorb,ispin,ichan,istate
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
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
       if(MpiComm /= MPI_COMM_NULL)call Bcast_MPI(MpiComm,alanc)
       if(MpiComm /= MPI_COMM_NULL)call Bcast_MPI(MpiComm,blanc)
    endif
#endif
    diag             = 0.d0
    subdiag          = 0.d0
    Z                = eye(Nlanc)
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    !
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))

    !
    call GFmatrix_allocate(impGmatrix(1,1,ispin,ispin,iorb,jorb),istate=istate,ichan=ichan,Nexc=Nlanc)
    !
    do j=1,nlanc
       de = diag(j)-Ei
       if (de>max_exc)max_exc=de
       peso = pesoBZ*Z(1,j)*Z(1,j)
       !
       impGmatrix(1,1,ispin,ispin,iorb,jorb)%state(istate)%channel(ichan)%weight(j) = peso
       impGmatrix(1,1,ispin,ispin,iorb,jorb)%state(istate)%channel(ichan)%poles(j)  = isign*de
       !
    enddo
  end subroutine add_to_lanczos_gf_normal



  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################

  function get_Gimp_normal_scalar(zeta) result(Gf)
    complex(8),intent(in)                                   :: zeta
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp) :: Gf
    complex(8)                                              :: green
    integer                                                 :: ispin
    integer                                                 :: ilat,jlat
    integer                                                 :: iorb,jorb
    integer                                                 :: iexc,Nexc
    integer                                                 :: ichan,Nchannel,istate,Nstates
    integer                                                 :: i,is,js
    complex(8)                                              :: weight,de
    real(8)                                                 :: chan4
    !
    if(.not.allocated(impGmatrix))stop "ed_gf_cluster ERROR: impGmatrix not allocated!"
    !
    Gf=zero
    !
    chan4=1d0 ; if(ed_gf_symmetric)chan4=0d0
    !
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
       Gf(1,1,ispin,ispin,iorb,jorb) = green
    enddo
    !
    do ispin=1,Nspin
       do iorb=1,Nimp
          do jorb=1,Nimp
             if(iorb==jorb)cycle
             gf(1,1,ispin,ispin,iorb,jorb) = 0.5d0*(gf(1,1,ispin,ispin,iorb,jorb) &
                  - (one-chan4*xi)*gf(1,1,ispin,ispin,iorb,iorb) - (one-chan4*xi)*gf(1,1,ispin,ispin,jorb,jorb))  
          enddo
       enddo
    enddo
    !
  end function get_Gimp_normal_scalar

  function get_Gimp_normal_array(zeta) result(Gf)
    complex(8),dimension(:),intent(in)                                 :: zeta
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(zeta)) :: Gf
    integer                                                            :: i
    do i=1,size(zeta)
       Gf(:,:,:,:,:,:,i) = get_Gimp_normal_scalar(zeta(i))
    enddo
  end function get_Gimp_normal_array







  function get_Sigma_normal_scalar(zeta) result(Sigma)
    complex(8),intent(in)                                   :: zeta
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp) :: Sigma
    complex(8),dimension(Nambu*Nspin*Nimp,Nambu*Nspin*Nimp) :: invG0,invG
    !
    if(allocated(Sigma))deallocate(Sigma)
    allocate(Sigma(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp)) ; Sigma=zero
    !
    !Get G0^-1
    invG0 = nnn2nso_reshape(invg0_bath(zeta))
    !
    !Get G^-1
    invG  = nnn2nso_reshape(get_Gimp_normal(zeta))
    call inv(invG)
    !
    !Get Sigma= G0^-1 - G^-1
    Sigma = nso2nnn_reshape(invG0 - invG)
    !
  end subroutine get_Sigma_normal

  function get_Sigma_normal_array(zeta) result(Sigma)
    complex(8),dimension(:),intent(in)                                 :: zeta
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(zeta)) :: Sigma
    integer                                                            :: i
    do i=1,size(zeta)
       Sigma(:,:,:,:,:,:,i) = get_Sigma_normal_scalar(zeta(i))
    enddo
  end function get_Sigma_normal_array





END MODULE ED_GF_NORMAL


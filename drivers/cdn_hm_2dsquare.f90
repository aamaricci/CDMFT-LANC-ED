program cdn_hm_2dsquare
  USE CDMFT_ED
  !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE MPI
  !
  implicit none
  integer                                     :: Nx,Ny,Nimp,Nlso,iloop,Nb,Nkx,Nky,Nsym
  integer                                     :: ilat,jlat,irepl,i,j,io,jo,iorb,jorb
  logical                                     :: converged
  real(8)                                     :: ts,wmixing,delta,dens_average,onsite
  real(8),allocatable,dimension(:)            :: dens_mats
  real(8),allocatable,dimension(:)            :: dens_ed
  !Bath:
  real(8),allocatable                         :: bath(:),bath_prev(:)
  !The local hybridization function:
  complex(8),allocatable                      :: Hloc(:,:)
  complex(8),allocatable,dimension(:,:,:,:,:) :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:) :: Greal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:) :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss
  !freq array:
  complex(8),allocatable                      :: Hk(:,:,:)
  complex(8),allocatable                      :: Sigma(:,:,:)
  !Density matrices:
  complex(8),allocatable,dimension(:,:)       :: reduced_density_matrix
  logical,allocatable,dimension(:)            :: orbital_mask
  !SYMMETRY BASIS for BATH:
  real(8),dimension(:,:),allocatable          :: lambdasym_vectors
  complex(8),dimension(:,:,:),allocatable     :: Hsym_basis
  character(len=64)                           :: finput,foutput
  !MPI VARIABLES (local use -> ED code has its own set of MPI variables)
  integer                                     :: comm
  integer                                     :: rank
  integer                                     :: mpi_size
  logical                                     :: master
  !
  !Init MPI: use of MPI overloaded functions in SciFor
  call init_MPI(comm,.true.)
  rank   = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  !  
  !Parse input variables
  call parse_cmd_variable(finput,"FINPUT",default='inputHM2D.conf')
  call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
  call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="hopping parameter")
  call parse_input_variable(Nx,"Nx",finput,default=2,comment="Number of cluster sites in x direction")
  call parse_input_variable(Ny,"Ny",finput,default=2,comment="Number of cluster sites in y direction")
  call parse_input_variable(Nkx,"Nkx",finput,default=100,comment="Number of kx point for BZ integration")
  call parse_input_variable(Nky,"Nky",finput,default=100,comment="Number of ky point for BZ integration")
  !
  call ed_read_input(trim(finput),comm)
  !
  !Add dmft control variables
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"Norb")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")
  call add_ctrl_var(ed_hw_band,"hwband")
  !
  !Set global variables
  if(Nlat.NE.Nx*Ny)then
     write(LOGfile,*) "                                                   "
     write(LOGfile,*) "WARNING: Nlat ≠ Nx * Ny -> Nlat will be overwritten"
     write(LOGfile,*) "                                                   "
  endif
  Nlat = Nx*Ny
  Nimp = Nlat*Norb
  Nlso = Nspin*Nimp
  !
  if(ED_VERBOSE > 0)call naming_convention()
  !  
  !Allocate Fields:
  allocate(Weiss(Nspin,Nspin,Nimp,Nimp,Lmats))
  allocate(Gmats(Nspin,Nspin,Nimp,Nimp,Lmats))
  allocate(Greal(Nspin,Nspin,Nimp,Nimp,Lreal))
  allocate(Smats(Nspin,Nspin,Nimp,Nimp,Lmats))
  allocate(Sreal(Nspin,Nspin,Nimp,Nimp,Lreal))
  allocate(Sigma(Nlso,Nlso,Lmats))


  !Build Hk and Hloc
  call generate_hk_hloc()


  !Build Hsym_basis and lambdasym_vectors
  !This trick is introduced in the case Hloc==zero ==> single site case
  Nsym = 1;if(any(abs(Hloc)/=zero))Nsym=2

  allocate(lambdasym_vectors(Nbath,Nsym))
  allocate(Hsym_basis(Nspin*Nimp,Nspin*Nimp,Nsym))

  !Replica onsite energies
  Hsym_basis(:,:,1) = zeye(Nspin*Nimp)
  !Replica hopping amplitudes
  if(Nsym==2)Hsym_basis(:,:,2) = abs(Hloc)
  !
  write(LOGfile,*) "HWBAND="//str(ed_hw_band)
  do irepl=1,Nbath
     onsite = irepl - 1 - (Nbath-1)/2d0        ![-(Nbath-1)/2:(Nbath-1)/2]
     onsite = onsite * 2*ed_hw_band/(Nbath-1)  !P-H symmetric band, -HWBAND:HWBAND
     lambdasym_vectors(irepl,1) = onsite       !Multiplies the suitable identity
  enddo
  if(mod(Nbath,2)==0)then
     lambdasym_vectors(Nbath/2,1)   = -1d-1    !Much needed small energies around
     lambdasym_vectors(Nbath/2+1,1) =  1d-1    !the fermi level. (for even Nbath)
  endif
  !
  if(Nsym==2)then
     do irepl=1,Nbath
        lambdasym_vectors(irepl,2) = 1d0+noise(0.01d0)           !Recall that TS is contained in Hloc
     enddo
  endif

  !SETUP BATH
  call ed_set_Hbath(Hsym_basis,lambdasym_vectors)
  Nb=ed_get_bath_dimension(Nsym)


  !SETUP Hloc
  call ed_set_hloc(Hloc)
  call print_matrix(Hloc)

  !SETUP SOLVER
  allocate(bath(Nb))
  call ed_init_solver(bath)


  allocate(bath_prev(Nb))
  bath_prev=zero

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath)
     call ed_get_sigma(Smats,axis='mats')


     !Retrieve ALL REDUCED DENSITY MATRICES DOWN TO THE LOCAL one
     if(dm_flag.AND.master)then
        !
        if(.not.allocated(orbital_mask))allocate(orbital_mask(Nimp))
        !
        ! All independent local RDMs (to check they are all equal)
        do ilat=1,Nlat
           do iorb=1,Norb
              io = iorb + (ilat-1)*Norb
              orbital_mask     = .false.
              orbital_mask(io) = .true.
              call ed_get_reduced_dm(reduced_density_matrix,orbital_mask,doprint=.true.)
              !Semi-analytical crosscheck of the single-orbital density matrix            
              call one_orb_benchmark(ilat,1,reduced_density_matrix)
           enddo
        enddo
        ! All independent two-site RDMs (to check NN and NNN are equal)
        do ilat=1,Nlat-1
           do iorb=1,Norb
              io = iorb + (ilat-1)*Norb
              orbital_mask = .false.
              orbital_mask(io) = .true.
              do jlat=ilat+1,Nlat
                 do jorb=1,Norb
                    jo = jorb + (jlat-1)*Norb
                    orbital_mask(io+1:Nimp) = .false.
                    orbital_mask(jo) = .true.
                    call ed_get_reduced_dm(reduced_density_matrix,orbital_mask,doprint=.true.)
                 enddo
              enddo
           enddo
        enddo
        !
     endif


     !Compute the local gfs on the imaginary axis:
     call dmft_get_gloc(Hk,Gmats,Smats,axis="m")
     call dmft_write_gf(Gmats,"Gloc",axis='matsubara',iprint=4)
     !
     !Get the Weiss field/Delta function to be fitted
     call dmft_self_consistency(Gmats,Smats,Weiss)
     call dmft_write_gf(Weiss,"Weiss",axis='mats',iprint=4)
     !
     !
     !Perform the SELF-CONSISTENCY by fitting the new bath
     call ed_chi2_fitgf(Weiss,bath)
     !
     !MIXING:
     if(iloop>1)bath = wmixing*bath + (1.d0-wmixing)*bath_prev
     bath_prev=bath
     call Bcast_MPI(comm,bath)
     !
     !
     !Check convergence (if required change chemical potential)
     converged = check_convergence(Weiss(1,1,1,1,:),dmft_error,nsuccess,nloop)
     call Bcast_MPI(comm,converged)
     !
     if(nread/=0d0)then
        allocate(dens_mats(Nimp))
        allocate(dens_ed(Nimp))
        call ed_get_dens(dens_ed)
        write(LOGfile,*)" "
        write(LOGfile,*)"Average ED-density:", sum(dens_ed)/Nimp
        !tot_density = density(up) + density(dw)
        dens_mats = zero
        do i=1,Nimp
           dens_mats(i) = dens_mats(i) + fft_get_density(Gmats(1,1,i,i,:),beta)
           dens_mats(i) = dens_mats(i) + fft_get_density(Gmats(2,2,i,i,:),beta)
        enddo
        write(LOGfile,*)" "
        write(LOGfile,*)"Average FFT-density:", sum(dens_mats)/Nimp
        !
        dens_average = sum(dens_mats)/Nimp
        call search_chemical_potential(xmu,dens_average,converged)

        call Bcast_MPI(comm,xmu)
        deallocate(dens_mats)
        deallocate(dens_ed)
     endif
     !
     !
     call end_loop
  enddo

  !Compute the local gfs on the real axis:
  call ed_get_sigma(Sreal,axis='real')
  call dmft_get_gloc(Hk,Greal,Sreal,axis='real')
  call dmft_write_gf(Greal,"Gloc",axis='real',iprint=4)

  !Compute the Kinetic Energy:
  call dmft_kinetic_energy(Hk(:,:,:),nn2so_f(Smats,Lmats))


  call ed_finalize_solver()
  call finalize_MPI()


contains


  function noise(W) result(r)
    real(8) :: W
    real(8) :: r
    r = mersenne()
    r = 2d0*r-1d0
    r = r*W
  end function noise


  !-------------------------------------------------------------------------------------------
  !PURPOSE:  Hk model for the 2d square lattice
  !-------------------------------------------------------------------------------------------
  !
  !> Hloc = H_{hop_intra_cluster}
  function hloc_model(N) result (hloc)
    integer                                     :: ix,iy,ilat,jlat,ispin,iorb,io,jo,N
    complex(8),dimension(Nspin,Nspin,Nimp,Nimp) :: hopping_matrix
    complex(8),dimension(N,N)                   :: hloc
    !
    hopping_matrix=zero
    !
    do ispin=1,Nspin
       do iorb=1,Norb
          !
          do ix=1,Nx
             do iy=1,Ny
                !
                ilat = indices2N([ix,iy])
                io   = iorb + (ilat-1)*Norb
                hopping_matrix(ispin,ispin,io,io)= zero!-mu_var
                if(ix<Nx)then !Avoid x higher outbound
                   jlat = indices2N([ix+1,iy])
                   jo   = iorb + (jlat-1)*Norb
                   hopping_matrix(ispin,ispin,io,jo)= -ts
                endif
                if(ix>1)then !Avoid x lower outbound
                   jlat = indices2N([ix-1,iy])
                   jo   = iorb + (jlat-1)*Norb
                   hopping_matrix(ispin,ispin,io,jo)= -ts
                endif
                if(iy<Ny)then !Avoid y higher outbound
                   jlat = indices2N([ix,iy+1])
                   jo   = iorb + (jlat-1)*Norb
                   hopping_matrix(ispin,ispin,io,jo)= -ts
                endif
                if(iy>1)then !Avoid y lower outbound
                   jlat = indices2N([ix,iy-1])
                   jo   = iorb + (jlat-1)*Norb
                   hopping_matrix(ispin,ispin,io,jo)= -ts
                endif
             enddo
          enddo
          !
       enddo
    enddo
    !    
    Hloc=nn2so(hopping_matrix)
    !
  end function hloc_model



  !> Hk = H_{hop_inter_cluster} + Hloc
  function hk_model(kpoint,N) result(hk)
    integer                                     :: N,ix,iy,ispin,iorb,ilat,jlat,io,jo
    real(8),dimension(:)                        :: kpoint
    real(8)                                     :: kx,ky
    complex(8),dimension(Nspin,Nspin,Nimp,Nimp) :: hopping_matrix
    complex(8),dimension(N,N)                   :: Hk
    !
    hopping_matrix=zero
    kx=kpoint(1)
    ky=kpoint(2)
    !
    do ispin=1,Nspin
       do iorb=1,Norb
          !
          do ix=1,Nx                  !Supercell bravais vector along y
             ilat = indices2N([ix,1])
             jlat = indices2N([ix,Ny])
             io   = iorb + (ilat-1)*Norb
             jo   = iorb + (jlat-1)*Norb
             hopping_matrix(ispin,ispin,io,jo)=hopping_matrix(ispin,ispin,io,jo) - ts*exp(xi*ky*Ny)
             hopping_matrix(ispin,ispin,jo,io)=hopping_matrix(ispin,ispin,jo,io) - ts*exp(-xi*ky*Ny)
          enddo
          !
          do iy=1,Ny                  !Supercell bravais vector along x
             ilat = indices2N([1,iy])
             jlat = indices2N([Nx,iy])
             io   = iorb + (ilat-1)*Norb
             jo   = iorb + (jlat-1)*Norb
             hopping_matrix(ispin,ispin,io,jo)=hopping_matrix(ispin,ispin,io,jo) - ts*exp(xi*kx*Nx)
             hopping_matrix(ispin,ispin,jo,io)=hopping_matrix(ispin,ispin,jo,io) - ts*exp(-xi*kx*Nx)
          enddo
          !
       enddo
    enddo
    !
    Hk=nn2so(hopping_matrix)+hloc_model(N)
    !
  end function hk_model






  !-------------------------------------------------------------------------------------------
  !PURPOSE:  Conventional indices for the cluster sites:
  !             y ^
  !             3 |    007 008 009
  !             2 |    004 005 006
  !             1 |    001 002 003
  !             0 |___ ___ ___ ___ _ >
  !                 0   1   2   3  x
  !-------------------------------------------------------------------------------------------
  function indices2N(indices) result(N)
    integer,dimension(2)         :: indices
    integer                      :: N
    !
    N = indices(1) + (indices(2)-1)*Nx
    !
  end function indices2N


  function N2indices(N) result(indices)
    integer,dimension(2)         :: indices
    integer                      :: N,i
    !
    indices(1)=mod(N,Nx)
    if(indices(1)==0)then
       indices(1)=Nx
       indices(2)=(N-Nx)/Nx+1
    else
       indices(2)=N/Nx+1
    endif
    !
  end function N2indices




  !-------------------------------------------------------------------------------------------
  !PURPOSE: generate Hloc and Hk
  !-------------------------------------------------------------------------------------------
  subroutine generate_hk_hloc()
    integer                                     :: ik
    real(8),dimension(Nkx*Nky,2)                :: kgrid
    real(8),dimension(2)                        :: e1,e2,bk1,bk2
    real(8)                                     :: bklen
    !
    e1 = [1d0, 0d0]
    e2 = [0d0, 1d0]
    call TB_set_ei(eix=e1,eiy=e2)
    bklen=2d0*pi
    bk1=bklen*[1d0, 0d0]
    bk2=bklen*[0d0, 1d0]
    call TB_set_bk(bkx=bk1,bky=bk2)
    !
    call TB_build_kgrid([Nkx,Nky],kgrid)
    kgrid(:,1)=kgrid(:,1)/Nx
    kgrid(:,2)=kgrid(:,2)/Ny
    !
    if(allocated(hk))deallocate(hk)
    if(allocated(hloc))deallocate(Hloc)
    !
    allocate(Hk(Nlso,Nlso,Nkx*Nky),Hloc(Nlso,Nlso))
    hk   = zero
    hloc = zero
    !
    call TB_build_model(Hk,hk_model,Nlso,kgrid)
    Hloc = hloc_model(Nlso)
    where(abs(dreal(Hloc))<1d-9)Hloc=0d0
    !
  end subroutine generate_hk_hloc







  !-------------------------------------------------------------------------------------------
  !PURPOSE: explicitate the cluster-indices convention in LOG files
  !-------------------------------------------------------------------------------------------
  subroutine naming_convention()
    integer                       :: i,j
    integer,dimension(Nx,Ny)      :: matrix
    !
    do j=1,Ny
       do i=1,Nx
          matrix(i,j)=indices2N([i,j])
       enddo
    enddo
    !
    write(LOGfile,"(A)")"The unique index of each site (on the cartesian plane) is as follows:"
    write(LOGfile,"(A)")" "
    do j=1,Ny
       write(LOGfile,"(20(I2,2x))")(matrix(i,Ny+1-j),i =1,Nx)
    enddo
    write(LOGfile,"(A)")" "
  end subroutine naming_convention




  !+---------------------------------------------------------------------------+
  !PURPOSE : check the local-dm comparing to Eq.4 in Mod.Phys.Lett.B.2013.27:05
  !+---------------------------------------------------------------------------+
  subroutine one_orb_benchmark(ilat,iorb,one_orb_dm)
    complex(8),allocatable,dimension(:,:)  :: one_orb_dm
    integer                                :: ilat,iorb
    real(8),allocatable,dimension(:,:)     :: dens,dens_up,dens_dw,docc,mag
    !
    allocate(dens(Nlat,Norb),dens_up(Nlat,Norb),dens_dw(Nlat,Norb),docc(Nlat,Norb),mag(Nlat,Norb))
    !
    call ed_get_mag(mag)
    call ed_get_dens(dens)
    call ed_get_docc(docc)
    dens_up = 0.5d0*(dens + mag)
    dens_dw = 0.5d0*(dens - mag)
    write(LOGfile,*)
    write(LOGfile,*) "LOCAL-DM BENCHMARK [Mod.Phys.Lett.B.2013.27:05]"
    write(LOGfile,*) "Semi-Analytical Estimate  |  Error"
    write(LOGfile,*) 1-dens_up(ilat,iorb)-dens_dw(ilat,iorb)+docc(ilat,iorb), "|", abs(1-dens_up(ilat,iorb)-dens_dw(ilat,iorb)+docc(ilat,iorb)-one_orb_dm(1,1))
    write(LOGfile,*) dens_up(ilat,iorb)-docc(ilat,iorb),                      "|", abs(dens_up(ilat,iorb)-docc(ilat,iorb)-one_orb_dm(2,2))
    write(LOGfile,*) dens_dw(ilat,iorb)-docc(ilat,iorb),                      "|", abs(dens_dw(ilat,iorb)-docc(ilat,iorb)-one_orb_dm(3,3))
    write(LOGfile,*) docc(ilat,iorb),                                         "|", abs(docc(ilat,iorb)-one_orb_dm(4,4))
    write(LOGfile,*)
    !
  end subroutine one_orb_benchmark



  function nn2so(In) result(Out)
    complex(8),dimension(Nspin,Nspin,Nimp,Nimp) :: In
    complex(8),dimension(Nspin*Nimp,Nspin*Nimp) :: Out
    integer :: is,js,i,j,io,jo
    do concurrent(is=1:Nspin,js=1:Nspin,i=1:Nimp,j=1:Nimp)
       io = i + (is-1)*Nimp
       jo = j + (js-1)*Nimp
       Out(io,jo) = In(is,js,i,j)
    enddo
  end function nn2so
  !
  function nn2so_f(In,L) result(Out)
    integer :: L
    complex(8),dimension(Nspin,Nspin,Nimp,Nimp,L) :: In
    complex(8),dimension(Nspin*Nimp,Nspin*Nimp,L) :: Out
    integer :: is,js,i,j,io,jo
    do concurrent(is=1:Nspin,js=1:Nspin,i=1:Nimp,j=1:Nimp)
       io = i + (is-1)*Nimp
       jo = j + (js-1)*Nimp
       Out(io,jo,:) = In(is,js,i,j,:)
    enddo
  end function nn2so_f



  function so2nn(In) result(Out)
    complex(8),dimension(Nspin*Nimp,Nspin*Nimp) :: In
    complex(8),dimension(Nspin,Nspin,Nimp,Nimp) :: Out
    integer :: is,js,i,j,io,jo
    do concurrent(is=1:Nspin,js=1:Nspin,i=1:Nimp,j=1:Nimp)
       io = i + (is-1)*Nimp
       jo = j + (js-1)*Nimp
       Out(is,js,i,j) = In(io,jo)
    enddo
  end function so2nn
  !
  function so2nn_f(In,L) result(Out)
    integer :: L
    complex(8),dimension(Nspin*Nimp,Nspin*Nimp,L) :: In
    complex(8),dimension(Nspin,Nspin,Nimp,Nimp,L) :: Out
    integer :: is,js,i,j,io,jo
    do concurrent(is=1:Nspin,js=1:Nspin,i=1:Nimp,j=1:Nimp)
       io = i + (is-1)*Nimp
       jo = j + (js-1)*Nimp
       Out(is,js,i,j,:) = In(io,jo,:)
    enddo
  end function so2nn_f

end program cdn_hm_2dsquare


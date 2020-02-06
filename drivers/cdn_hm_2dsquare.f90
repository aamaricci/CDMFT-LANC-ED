program cdn_hm_2dsquare
   USE CDMFT_ED
   !
   USE SCIFOR
   USE DMFT_TOOLS
   !
   USE MPI
   !
   implicit none
   integer                                                                :: Nx,Ny,Nlso,iloop,Nb,Nkx,Nky,iw,iii,jjj,kkk
   integer,dimension(2):: recover
   logical                                                                :: converged
   real(8)                                                                :: ts,tsp,wmixing
   !Bath:
   real(8),allocatable                                                    :: Bath(:),Bath_fitted(:),Bath_old(:)
   !The local hybridization function:
   complex(8),allocatable                                                 :: Hloc(:,:)
   complex(8),allocatable,dimension(:,:,:,:,:,:,:)                        :: Gmats,Greal,Smats,Sreal,Weiss,Weiss_old
   character(len=16)                                                      :: finput
   real(8),allocatable                                                    :: wt(:)
   complex(8),allocatable                                                 :: wm(:),wr(:)
   complex(8),allocatable                                                 :: Hk(:,:,:),Smats_lso(:,:,:)
   !custom observable example
   complex(8),dimension(:,:),allocatable                                  :: observable_matrix
   !SYMMETRIES TEST
   real(8),dimension(:),allocatable                                       :: lambdasym_vector
   complex(8),dimension(:,:,:,:,:,:,:),allocatable                        :: Hsym_basis
   !MPI VARIABLES (local use -> ED code has its own set of MPI variables)
   integer                                                                :: comm
   integer                                                                :: rank
   integer                                                                :: mpi_size
   logical                                                                :: master
   character(len=6)                                                       :: scheme
   !Init MPI: use of MPI overloaded functions in SciFor
   call init_MPI(comm,.true.)
   rank   = get_Rank_MPI(comm)
   master = get_Master_MPI(comm)

   !

   !Parse input variables
   call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
   call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
   call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="hopping parameter")
   call parse_input_variable(Nx,"Nx",finput,default=2,comment="Number of cluster sites in x direction")
   call parse_input_variable(Ny,"Ny",finput,default=2,comment="Number of cluster sites in y direction")
   call parse_input_variable(Nkx,"Nkx",finput,default=10,comment="Number of kx point for BZ integration")
   call parse_input_variable(Nky,"Nky",finput,default=10,comment="Number of ku point for BZ integration")
   call parse_input_variable(scheme,"SCHEME",finput,default="g",comment="Periodization scheme: possible g or sigma")

   !
   call ed_read_input(trim(finput),comm)
   !
   !Add dmft control variables
   !
   call add_ctrl_var(beta,"BETA")
   call add_ctrl_var(Norb,"Norb")
   call add_ctrl_var(Nspin,"Nspin")
   call add_ctrl_var(xmu,"xmu")
   call add_ctrl_var(wini,"wini")
   call add_ctrl_var(wfin,"wfin")
   call add_ctrl_var(eps,"eps")

   !set global variables
   !if (Nspin/=1.or.Norb/=1) stop "You are using too many spin-orbitals"
   !Ny=Nx:w

   !Nky=Nkx
   Nlat=Nx*Ny
   Nlso=Nlat*Nspin*Norb
   if(.not.allocated(wm))allocate(wm(Lmats))
   wm     = xi*pi/beta*real(2*arange(1,Lmats)-1,8)

   if(ED_VERBOSE > 0)call naming_convention()
   !Allocate Weiss Field:
   allocate(Weiss(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Weiss_old(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
   allocate(Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(Smats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(Smats_lso(Nlso,Nlso,Lmats))

   !Build Hk and Hloc
   call generate_hk_hloc()
   allocate(lambdasym_vector(1))
   allocate(Hsym_basis(Nlat,Nlat,Nspin,Nspin,Norb,Norb,1))
   
   !Allocate custom observable matrix (test with n_11)
   !allocate(observable_matrix(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
   !observable_matrix=zero
   !observable_matrix(1,1)=one
   !call init_custom_observables(1,Hk)
   !call add_custom_observable("test",observable_matrix)


   !Build Hbasis and lambda vector
   Hsym_basis(:,:,:,:,:,:,1)=abs(lso2nnn(Hloc))
   lambdasym_vector=[-1.d0]
   
   !setup solver
   call set_Hloc(Hsym_basis,lambdasym_vector)
   Nb=get_bath_dimension(Hsym_basis)
   allocate(bath(Nb))
   allocate(bath_old(Nb))
   Bath_old=zero
   call ed_init_solver(comm,bath)

   !DMFT loop
   iloop=0;converged=.false.
   do while(.not.converged.AND.iloop<nloop)
      iloop=iloop+1
      if(master)call start_loop(iloop,nloop,"DMFT-loop")

      !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
      call ed_solve(comm,bath) 
      call ed_get_sigma_matsubara(Smats)
      call ed_get_sigma_realaxis(Sreal)

      !Compute the local gfs:
      call dmft_gloc_matsubara(comm,Hk,Wt,Gmats,Smats)
      if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4)
      !
      !Get the Weiss field/Delta function to be fitted
      call dmft_self_consistency(comm,Gmats,Smats,Weiss,lso2nnn(Hloc),cg_scheme)
      call Bcast_MPI(comm,Weiss)
      !
      !
      !Perform the SELF-CONSISTENCY by fitting the new bath
      if(master)then
         call ed_chi2_fitgf(Weiss,bath)
         !Check convergence (if required change chemical potential)
         converged = check_convergence(Weiss(:,:,1,1,1,1,:),dmft_error,nsuccess,nloop)
      endif
      !
      !MIXING:
      if(iloop>1)Bath= wmixing*Bath + (1.d0-wmixing)*Bath_Old
      Bath_old=Bath
      !
      call Bcast_MPI(comm,bath)
      call Bcast_MPI(comm,converged)
      call Bcast_MPI(comm,xmu)
      !
      if(master)call end_loop
   enddo
   !
   !Cleanup after calculating custom observables
   call clear_custom_observables()
   !Compute the local gfs:
   call dmft_gloc_realaxis(comm,Hk,Wt,Greal,Sreal)
   if(master)call dmft_print_gf_realaxis(Greal,"Gloc",iprint=4)

   !Compute the Kinetic Energy:
   do iw=1,Lmats
      Smats_lso(:,:,iw)=nnn2lso(Smats(:,:,:,:,:,:,iw))
   enddo
   call dmft_kinetic_energy(comm,Hk(:,:,:),Wt,Smats_lso)
   
   !PERIODIZE
   !call print_periodized([Nkx,Nky],hk_model,hk_periodized,scheme)

   call finalize_MPI()


contains

   !include "auxiliary_routines.f90"

   
   !-------------------------------------------------------------------------------------------
   !PURPOSE:  Hk model for the 2d square lattice
   !-------------------------------------------------------------------------------------------
   function hloc_model(N) result (hloc)
      integer                                               :: ilat,jlat,ispin,iorb,ind1,ind2,N
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: hopping_matrix
      complex(8),dimension(N,N)                             :: hloc
      !
      hopping_matrix=zero
      !
      do ispin=1,Nspin
         do iorb=1,Norb
            do ilat=1,Nx
               do jlat=1,Ny
                  ind1=indices2N([ilat,jlat])
                  hopping_matrix(ind1,ind1,ispin,ispin,iorb,iorb)= 0.d0!-mu_var
                  if(ilat<Nx)then
                     ind2=indices2N([ilat+1,jlat])
                     hopping_matrix(ind1,ind2,ispin,ispin,iorb,iorb)= -ts
                  endif
                  if(ilat>1)then
                     ind2=indices2N([ilat-1,jlat])
                     hopping_matrix(ind1,ind2,ispin,ispin,iorb,iorb)= -ts
                  endif
                  if(jlat<Ny)then
                     ind2=indices2N([ilat,jlat+1])
                     hopping_matrix(ind1,ind2,ispin,ispin,iorb,iorb)= -ts
                  endif
                  if(jlat>1)then
                     ind2=indices2N([ilat,jlat-1])
                     hopping_matrix(ind1,ind2,ispin,ispin,iorb,iorb)= -ts
                  endif
               enddo
            enddo
         enddo
      enddo
      !
      Hloc=nnn2lso(hopping_matrix)
      !
   end function hloc_model

   function hk_model(kpoint,N) result(hk)
      integer                                               :: n,ilat,ispin,iorb,ind1,ind2
      real(8),dimension(:)                                  :: kpoint
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: hopping_matrix
      complex(8),dimension(N,N)                             :: hk
      !
      hopping_matrix=zero
      !
      do ispin=1,Nspin
         do iorb=1,Norb
            do ilat=1,Nx
               ind1=indices2N([ilat,1])
               ind2=indices2N([ilat,Ny])
               hopping_matrix(ind1,ind2,ispin,ispin,iorb,iorb)=hopping_matrix(ind1,ind2,ispin,ispin,iorb,iorb) -ts*exp(xi*kpoint(2)*Ny)
               hopping_matrix(ind2,ind1,ispin,ispin,iorb,iorb)=hopping_matrix(ind2,ind1,ispin,ispin,iorb,iorb) -ts*exp(-xi*kpoint(2)*Ny)
            enddo
            do ilat=1,Ny
               ind1=indices2N([1,ilat])
               ind2=indices2N([Nx,ilat])
               hopping_matrix(ind1,ind2,ispin,ispin,iorb,iorb)=hopping_matrix(ind1,ind2,ispin,ispin,iorb,iorb) -ts*exp(xi*kpoint(1)*Nx)
               hopping_matrix(ind2,ind1,ispin,ispin,iorb,iorb)=hopping_matrix(ind2,ind1,ispin,ispin,iorb,iorb) -ts*exp(-xi*kpoint(1)*Nx)
            enddo
         enddo
      enddo
      !
      Hk=nnn2lso(hopping_matrix)+hloc_model(N)
   end function hk_model

   function hk_periodized(kpoint,N) result(Hk)
      real(8),dimension(:)                          :: kpoint
      integer                                       :: Nlat_,Nx_,Ny_,N
      complex(8),dimension(N,N)                     :: Hk
      !
      Nlat_=Nlat
      Nx_=Nx
      Ny_=Ny
      Nlat=1
      Nx=1
      Ny=1
      !
      Hk=hk_model(kpoint,Nspin*Norb)
      !
      Nlat=Nlat_
      Nx=Nx_
      Ny=Ny_
      !
   end function hk_periodized

   !-------------------------------------------------------------------------------------------
   !PURPOSE: generate Hloc and Hk
   !-------------------------------------------------------------------------------------------

   subroutine generate_hk_hloc()
      integer                                     :: ik
      real(8),dimension(Nkx*Nky,2)                :: kgrid
      !
      call TB_build_kgrid([Nkx,Nky],kgrid)
      kgrid(:,1)=kgrid(:,1)/Nx
      kgrid(:,2)=kgrid(:,2)/Ny
      !
      if(allocated(hk))deallocate(hk)
      if(allocated(Wt))deallocate(Wt)
      if(allocated(hloc))deallocate(Hloc)
      !
      allocate(Hk(Nlso,Nlso,Nkx*Nky),Wt(Nkx*Nky),Hloc(Nlso,Nlso))
      hk=zero
      wt=zero
      hloc=zero
      !
      call TB_build_model(Hk,hk_model,Nlso,kgrid)
      Wt = 1d0/(Nkx*Nky)
      Hloc=hloc_model(Nlso)
      where(abs(dreal(Hloc))<1.d-9)Hloc=0d0
      !
   end subroutine generate_hk_hloc

   !-------------------------------------------------------------------------------------------
   !PURPOSE: auxilliary reshape functions
   !-------------------------------------------------------------------------------------------

   function lso2nnn(Hlso) result(Hnnn)
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
      integer                                               :: ilat,jlat
      integer                                               :: iorb,jorb
      integer                                               :: ispin,jspin
      integer                                               :: is,js
      Hnnn=zero
      do ilat=1,Nlat
         do jlat=1,Nlat
            do ispin=1,Nspin
               do jspin=1,Nspin
                  do iorb=1,Norb
                     do jorb=1,Norb
                        is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                        js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                        Hnnn(ilat,jlat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end function lso2nnn


   function nnn2lso(Hnnn) result(Hlso)
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
      integer                                               :: ilat,jlat
      integer                                               :: iorb,jorb
      integer                                               :: ispin,jspin
      integer                                               :: is,js
      Hlso=zero
      do ilat=1,Nlat
         do jlat=1,Nlat
            do ispin=1,Nspin
               do jspin=1,Nspin
                  do iorb=1,Norb
                     do jorb=1,Norb
                        is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                        js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                        Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end function nnn2lso

   function indices2N(indices) result(N)
      integer,dimension(2)         :: indices
      integer                      :: N,i
      !
      !
      N=Nx*(indices(2)-1)+indices(1)
   end function indices2N

   function N2indices(N) result(indices) 
      integer,dimension(2)         :: indices
      integer                      :: N,i,N_
      !
      indices(1)=mod(N,Nx)
      if(indices(1)==0)indices(1)=3
      indices(2)=N/Nx
   end function N2indices

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

end program cdn_hm_2dsquare




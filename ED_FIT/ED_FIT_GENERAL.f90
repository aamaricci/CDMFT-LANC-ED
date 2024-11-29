MODULE ED_FIT_GENERAL
  USE ED_FIT_COMMON
  USE SF_SPIN, only: pauli_tau_z,pauli_sigma_z
  USE SF_LINALG, only: kron



  implicit none
  private


  interface chi2_fitgf_general
     module procedure :: chi2_fitgf_general_normal_n4
     module procedure :: chi2_fitgf_general_superc_n4
  end interface chi2_fitgf_general

  public :: chi2_fitgf_general

  integer :: i,j,ilat,jlat,iorb,jorb,ispin,jspin,ibath,io,jo,s,l

  
contains



  !+-------------------------------------------------------------+
  !PURPOSE  : Chi^2 interface for GENERAL BATH
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_general_normal_n4(fg,bath_)
    complex(8),dimension(:,:,:,:,:)    :: fg ![Nspin][Nspin][Nimp][Nimp][Lmats]
    real(8),dimension(:),intent(inout) :: bath_
    real(8),dimension(:),allocatable   :: array_bath
    integer                            :: iter,stride,counter,Asize
    real(8)                            :: chi
    logical                            :: check
    character(len=256)                 :: suffix
    integer                            :: unit
    !
    call assert_shape(fg,[Nspin,Nspin,Nimp,Nimp,size(fg,5)],"chi2_fitgf_general_normal_n4","fg")
    !
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_general error: wrong bath dimensions"
    !
    if(cg_pow/=2.AND.cg_norm=="frobenius")then
       print *, "WARNING: CG_POW must be 2 for a meaningful definition of the Frobenius norm."
       print *, "         we'll still let you go ahead with the desired input, but please be "
       print *, "         be aware that CG_POW is not doing what you would expect for a chi^q"
    endif
    !
    call allocate_dmft_bath()
    call set_dmft_bath(bath_)
    allocate(array_bath(size(bath_)-1))
    Nlambdas   = nint(bath_(1))
    array_bath = bath_(2:)
    !
    Ldelta = Lfit ; if(Ldelta>size(fg,7))Ldelta=size(fg,7)
    !
    !
    allocate(FGmatrix(Nspin,Nspin,Nimp,Nimp,Ldelta))
    allocate(Hmask(Nspin,Nspin,Nimp,Nimp))
    allocate(Xdelta(Ldelta))
    allocate(Wdelta(Ldelta))
    !
    Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
    !
    select case(cg_weight)
    case default
       Wdelta=1d0
    case(2)
       Wdelta=1d0*arange(1,Ldelta)
    case(3)
       Wdelta=Xdelta
    end select
    !
    !
    !----------------------------------------------------------------------------------------
    Hmask=.true.
    !----------------------------------------------------------------------------------------
    !
    !
    FGmatrix = fg
    !
    !
    select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE
    case default
       if(cg_grad==0)then
          write(LOGfile,*)"  Using analytic gradient"
          select case (cg_scheme)
          case ("weiss")
             call fmin_cg(array_bath,&
                  chi2_weiss_general_normal,&
                  grad_chi2_weiss_general_normal,&
                  iter,&
                  chi, &
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cg(array_bath,&
                  chi2_delta_general_normal,&
                  grad_chi2_delta_general_normal,&
                  iter,&
                  chi, &
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_general_normal error: cg_scheme != [weiss,delta]"
          end select
       else
          write(LOGfile,*)"  Using numerical gradient"
          select case (cg_scheme)
          case ("weiss")
             call fmin_cg(array_bath,&
                  chi2_weiss_general_normal,&
                  iter,&
                  chi, &
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cg(array_bath,&
                  chi2_delta_general_normal,&
                  iter,&
                  chi, &
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_general_normal error: cg_scheme != [weiss,delta]"
          end select
       endif
       !
    case (1)
       if(cg_grad==0)then
          write(*,*) "                                                                                "
          write(*,*) "WARNING: analytic gradient not available with cg-method=1 (minimize f77 routine)"
          write(*,*) "         > we will force cg_grad=1 (so let the routine estimate the gradient)   "
          write(*,*) "                                                                                "
       endif
       select case (cg_scheme)
       case ("weiss")
          call fmin_cgminimize(array_bath,&
               chi2_weiss_general_normal,&
               iter,&
               chi, &
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               new_version=cg_minimize_ver,&
               hh_par=cg_minimize_hh,      &
               iverbose=(ed_verbose>3))
       case ("delta")
          call fmin_cgminimize(array_bath,&
               chi2_delta_general_normal,&
               iter,&
               chi, &
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               new_version=cg_minimize_ver,&
               hh_par=cg_minimize_hh,      &
               iverbose=(ed_verbose>3))
       case default
          stop "chi2_fitgf_general_normal error: cg_scheme != [weiss,delta]"
       end select
       !
    end select
    !
    write(LOGfile,"(A,ES18.9,A,I5,A)")"chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,"  <--  All Orbs, All Spins"
    !
    suffix="_ALLorb_ALLspins"//reg(ed_file_suffix)
    unit=free_unit()
    open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
    write(unit,"(ES18.9,1x,I5)") chi,iter
    close(unit)
    !
    bath_(Nbath+1:size(bath_))=array_bath
    call set_dmft_bath(bath_)               ! *** bath_ --> dmft_bath ***    (per write fit result)
    call write_dmft_bath(LOGfile)
    call save_dmft_bath()
    !
    call write_fit_result()
    !
    call get_dmft_bath(bath_)               ! ***  dmft_bath --> bath_ ***    (bath in output)
    call deallocate_dmft_bath()
    deallocate(FGmatrix,Hmask,Xdelta,Wdelta)
    deallocate(array_bath)
    !
  contains
    !
    subroutine write_fit_result()
      complex(8)        :: fgand(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,Ldelta)
      real(8)           :: w
      if(cg_scheme=='weiss')then
         fgand = g0and_bath_function(xi*Xdelta,dmft_bath)
      else
         fgand = delta_bath_function(xi*Xdelta,dmft_bath)
      endif
      !
      do ispin=1,Nspin
         do jspin=1,Nspin
            do iorb=1,Nimp
               do jorb=1,Nimp
                  suffix="_l"//reg(str(iorb))//&
                       "_m"//reg(str(jorb))//&
                       "_s"//reg(str(ispin))//&
                       "_r"//reg(str(jspin))//reg(ed_file_suffix)
                  unit=free_unit()
                  if(cg_scheme=='weiss')then
                     open(unit,file="fit_weiss"//reg(suffix)//".ed")
                  else
                     open(unit,file="fit_delta"//reg(suffix)//".ed")
                  endif
                  do i=1,Ldelta
                     w = Xdelta(i)
                     write(unit,"(5F24.15)")Xdelta(i),&
                          dimag(fg(ispin,jspin,iorb,jorb,i)),dimag(fgand(1,1,ispin,jspin,iorb,jorb,i)),&
                          dreal(fg(ispin,jspin,iorb,jorb,i)),dreal(fgand(1,1,ispin,jspin,iorb,jorb,i))
                  enddo
                  close(unit)
               enddo
            enddo
         enddo
      enddo
    end subroutine write_fit_result
    !
  end subroutine chi2_fitgf_general_normal_n4





  subroutine chi2_fitgf_general_superc_n4(fg,ff,bath_)
    complex(8),dimension(:,:,:,:,:)    :: fg ![Nspin,Nspin,Nimp,Nimp][Lmats]
    complex(8),dimension(:,:,:,:,:)    :: ff ![Nspin,Nspin,Nimp,Nimp][Lmats]
    real(8),dimension(:),intent(inout) :: bath_
    real(8),dimension(:),allocatable   :: array_bath
    integer                            :: iter,stride,counter,Asize
    real(8)                            :: chi
    logical                            :: check
    character(len=256)                 :: suffix
    integer                            :: unit
    !
    call assert_shape(fg,[Nspin,Nspin,Nimp,Nimp,size(fg,5)],"chi2_fitgf_general_superc_n4","fg")
    call assert_shape(ff,[Nspin,Nspin,Nimp,Nimp,size(fg,5)],"chi2_fitgf_general_superc_n4","ff")
    !
    !allocate(Nlambdas(Nbath))
    !
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_general_superc error: wrong bath dimensions"
    !
    if(cg_pow/=2.AND.cg_norm=="frobenius")then
       print *, "WARNING: CG_POW must be 2 for a meaningful definition of the Frobenius norm."
       print *, "         we'll still let you go ahead with the desired input, but please be "
       print *, "         be aware that CG_POW is not doing what you would expect for a chi^q"
    endif
    !
    call allocate_dmft_bath(dmft_bath)
    call set_dmft_bath(bath_,dmft_bath)
    allocate(array_bath(size(bath_)-1))
    Nlambdas  =nint(bath_(1))
    array_bath=bath_(2:)
    !
    Ldelta = Lfit ; if(Ldelta>size(fg,5))Ldelta=size(fg,5)
    !
    !
    allocate(FGmatrix(Nspin,Nspin,Nimp,Nimp,Ldelta))
    allocate(FFmatrix(Nspin,Nspin,Nimp,Nimp,Ldelta))
    allocate(Hmask(Nspin,Nspin,Nimp,Nimp))
    allocate(Xdelta(Ldelta))
    allocate(Wdelta(Ldelta))
    !
    Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
    !
    select case(cg_weight)
    case default
       Wdelta=1d0
    case(2)
       Wdelta=1d0*arange(1,Ldelta)
    case(3)
       Wdelta=Xdelta
    end select
    !
    !
    !----------------------------------------------------------------------------------------
    Hmask=.true.
    !----------------------------------------------------------------------------------------
    !
    !
    FGmatrix = fg
    FFmatrix = ff
    !
    if(cg_grad==0)then
       cg_grad=1
       write(LOGfile,*)"WARNING: CGfit Superc does not support minimization with analytical gradient (cg_grad=0)."
       write(LOGfile,*)"         To spare you some time we fall back to numerical gradient: cg_grad=1"
       call sleep(2)
    endif
    !
    select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE
    case default
       write(LOGfile,*)"  Using numerical gradient:"
       select case (cg_scheme)
       case ("weiss")
          call fmin_cg(array_bath,&
               chi2_weiss_general_superc,&
               iter,&
               chi, &
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               istop=cg_stop, &
               iverbose=(ed_verbose>3))
       case ("delta")
          call fmin_cg(array_bath,&
               chi2_delta_general_superc,&
               iter,&
               chi, &
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               istop=cg_stop, &
               iverbose=(ed_verbose>3))
       case default
          stop "chi2_fitgf_general_superc error: cg_scheme != [weiss,delta]"
       end select
       !
    case (1)
       select case (cg_scheme)
       case ("weiss")
          call fmin_cgminimize(array_bath,&
               chi2_weiss_general_superc,&
               iter,&
               chi, &
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               new_version=cg_minimize_ver,&
               hh_par=cg_minimize_hh,      &
               iverbose=(ed_verbose>3))
       case ("delta")
          call fmin_cgminimize(array_bath,&
               chi2_delta_general_superc,&
               iter,&
               chi, &
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               new_version=cg_minimize_ver,&
               hh_par=cg_minimize_hh,      &
               iverbose=(ed_verbose>3))
       case default
          stop "chi2_fitgf_general_superc error: cg_scheme != [weiss,delta]"
       end select
       !
    end select
    !
    write(LOGfile,"(A,ES18.9,A,I5,A)")"chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,"  <--  All Orbs, All Spins"
    !
    suffix="_ALLorb_ALLspins"//reg(ed_file_suffix)
    unit=free_unit()
    open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
    write(unit,"(ES18.9,1x,I5)") chi,iter
    close(unit)
    !
    bath_(Nbath+1:size(bath_))=array_bath
    call set_dmft_bath(bath_,dmft_bath) ! *** array_bath --> dmft_bath *** (per write fit result)
    call write_dmft_bath(dmft_bath,LOGfile)
    call save_dmft_bath(dmft_bath)
    !
    call write_fit_result()
    !
    call get_dmft_bath(dmft_bath,bath_)                ! ***  dmft_bath --> bath_ ***    (bath in output)
    call deallocate_dmft_bath(dmft_bath)
    deallocate(FGmatrix,FFmatrix,Hmask,Xdelta,Wdelta)
    deallocate(array_bath)
    ! deallocate(Nlambdas)
    !
  contains
    !
    subroutine write_fit_result()
      complex(8)                       :: fgand(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,Ldelta)
      real(8)                          :: w
      if(cg_scheme=='weiss')then
         fgand = g0and_bath_function(xi*Xdelta,dmft_bath)
      else
         fgand = delta_bath_function(xi*Xdelta,dmft_bath)
      endif
      !
      do ispin=1,Nspin
         do jspin=1,Nspin
            do iorb=1,Nimp
               do jorb=1,Nimp
                  suffix="_l"//reg(str(iorb))//&
                       "_m"//reg(str(jorb))//&
                       "_s"//reg(str(ispin))//&
                       "_r"//reg(str(jspin))//reg(ed_file_suffix)
                  unit=free_unit()
                  if(cg_scheme=='weiss')then
                     open(unit,file="fit_weiss"//reg(suffix)//".ed")
                  else
                     open(unit,file="fit_delta"//reg(suffix)//".ed")
                  endif
                  do i=1,Ldelta
                     w = Xdelta(i)
                     write(unit,"(9F24.15)")Xdelta(i),&
                          dimag(fg(ispin,jspin,iorb,jorb,i)),&
                          dimag(fgand(1,1,ispin,jspin,iorb,jorb,i)),&
                          dreal(fg(ispin,jspin,iorb,jorb,i)),&
                          dreal(fgand(1,1,ispin,jspin,iorb,jorb,i)),&
                          dimag(ff(ispin,jspin,iorb,jorb,i)),&
                          dimag(fgand(1,2,ispin,jspin,iorb,jorb,i)),&
                          dreal(ff(ispin,jspin,iorb,jorb,i)),&
                          dreal(fgand(1,2,ispin,jspin,iorb,jorb,i))
                  enddo
                  close(unit)
               enddo
            enddo
         enddo
      enddo
    end subroutine write_fit_result
    !
  end subroutine chi2_fitgf_general_superc_n4




  include "chi2_delta_general.f90"
  include "delta_general.f90"

  include "chi2_weiss_general.f90"
  include "weiss_general.f90"

end MODULE ED_FIT_GENERAL

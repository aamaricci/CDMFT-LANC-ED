MODULE ED_FIT_REPLICA
  USE ED_FIT_COMMON
  USE SF_SPIN
  USE SF_LINALG, only: kron,diag,inv
  !
  implicit none
  private

  interface chi2_fitgf_replica
     module procedure :: chi2_fitgf_replica_normal_n4
     module procedure :: chi2_fitgf_replica_superc_n4
  end interface chi2_fitgf_replica
  
  public :: chi2_fitgf_replica

  integer :: i,j,ilat,jlat,iorb,jorb,ispin,jspin,ibath,io,jo,s,l

contains


  !+-------------------------------------------------------------+
  !PURPOSE  : Chi^2 interface for REPLICA BATH
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_replica_normal_n4(fg,bath_)
    complex(8),dimension(:,:,:,:,:)    :: fg ![Nspin][Nspin][Nimp][Nimp][Lmats]
    real(8),dimension(:),intent(inout) :: bath_
    real(8),dimension(:),allocatable   :: array_bath
    integer                            :: iter,Asize
    real(8)                            :: chi
    logical                            :: check
    character(len=256)                 :: suffix
    integer                            :: unit
    !
    call assert_shape(fg,[Nspin,Nspin,Nimp,Nimp,size(fg,5)],"chi2_fitgf_replica_normal_n4","fg")
    !
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_replica_normal error: wrong bath dimensions"
    !
    if(cg_pow/=2.AND.cg_norm=="frobenius")then
       print *, "WARNING: CG_POW must be 2 for a meaningful definition of the Frobenius norm."
       print *, "         we'll still let you go ahead with the desired input, but please be "
       print *, "         be aware that CG_POW is not doing what you would expect for a chi^q"
    endif
    !
    call allocate_dmft_bath(dmft_bath_fit)
    call set_dmft_bath(bath_,dmft_bath_fit)
    allocate(array_bath(size(bath_)-1))
    Nlambdas   = nint(bath_(1))
    array_bath = bath_(2:)
    print*,"in:",array_bath,size(array_bath),size(bath_(2:))
    !
    Ldelta = Lfit ; if(Ldelta>size(fg,5))Ldelta=size(fg,5)
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
    !POSSIBLY TO BE USED AT SOME POINT (but beware that is meaningless for `FROBENIUS` norm)
    Hmask=.true.
    ! Aren't we sure about FG_ij = FG_ji, coming from Hbath being hermitian?
    !
    ! -> Hmask=Hbath_mask(wdiag=.false.,uplo=.true.) NO, this does more than
    !    that, putting .FALSE. where Hbath is zero, which may lead to tricky
    !    wrong fits: I have data showing this for a trimer, where at ~ 6.30d-10 
    !    dmft error, hence a very well converged point, we have a nonnegligible
    !    Im(Weiss) component for ilat=1, jlat=3, despite Hmask(1,3) = Hmask(3,1) 
    !    would result .FALSE. if computed with the Hbath_mask call above.
    !
    ! For now I'd say that we'd better dump everything inside the \chi^2, hence Hmask=.true.,
    ! but we might want to consider exploiting hermiticity at least (probably not looking at
    ! Hbath though: who guarantees zeros therein imply zeros here?). Discussion needed...
    !----------------------------------------------------------------------------------------
    !
    !
    FGmatrix = fg
    !    
    select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE
    case default
       if(cg_grad==0)then
          write(LOGfile,*)"  Using analytic gradient"
          select case (cg_scheme)
          case ("weiss")
             call fmin_cg(array_bath,&
                  chi2_weiss_replica_normal,&
                  grad_chi2_weiss_replica_normal,&
                  iter,&
                  chi, &
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cg(array_bath,&
                  chi2_delta_replica_normal,&
                  grad_chi2_delta_replica_normal,&
                  iter,&
                  chi, &
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_replica_normal error: cg_scheme != [weiss,delta]"
          end select
       else
          write(LOGfile,*)"  Using numerical gradient"
          select case (cg_scheme)
          case ("weiss")
             call fmin_cg(array_bath,&
                  chi2_weiss_replica_normal,&
                  iter,&
                  chi, &
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cg(array_bath,&
                  chi2_delta_replica_normal,&
                  iter,&
                  chi, &
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_replica_normal error: cg_scheme != [weiss,delta]"
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
          call fmin_cgminimize(array_bath,chi2_weiss_replica_normal,&
               iter,&
               chi, &
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               new_version=cg_minimize_ver,&
               hh_par=cg_minimize_hh,      &
               iverbose=(ed_verbose>3))
       case ("delta")
          call fmin_cgminimize(array_bath,chi2_delta_replica_normal,&
               iter,&
               chi, &
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               new_version=cg_minimize_ver,&
               hh_par=cg_minimize_hh,      &
               iverbose=(ed_verbose>3))
       case default
          stop "chi2_fitgf_replica_normal error: cg_scheme != [weiss,delta]"
       end select
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
    bath_(2:) = array_bath
    call set_dmft_bath(bath_,dmft_bath_fit) ! *** array_bath --> dmft_bath *** (per write fit result)
    call write_dmft_bath(dmft_bath_fit,LOGfile)
    call save_dmft_bath(dmft_bath_fit)
    !
    call write_fit_result()
    !
    call get_dmft_bath(dmft_bath_fit,bath_)                ! ***  dmft_bath --> bath_ ***    (bath in output)
    call deallocate_dmft_bath(dmft_bath_fit)
    deallocate(FGmatrix,Hmask,Xdelta,Wdelta)
    deallocate(array_bath)
    !
  contains
    !
    subroutine write_fit_result()
      real(8)    :: w
      complex(8) :: fgand(Nspin,Nspin,Nimp,Nimp,Ldelta)
      if(cg_scheme=='weiss')then
         fgand = g0and_replica_normal(array_bath)
      else
         fgand = delta_replica_normal(array_bath)
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
                          dimag(fg(ispin,jspin,iorb,jorb,i)),dimag(fgand(ispin,jspin,iorb,jorb,i)),&
                          dreal(fg(ispin,jspin,iorb,jorb,i)),dreal(fgand(ispin,jspin,iorb,jorb,i))
                  enddo
                  close(unit)
               enddo
            enddo
         enddo
      enddo
    end subroutine write_fit_result
    !
  end subroutine chi2_fitgf_replica_normal_n4





  

  subroutine chi2_fitgf_replica_superc_n4(fg,ff,bath_)
    complex(8),dimension(:,:,:,:,:)    :: fg ![Nspin,Nspin,Nimp,Nimp][Lmats]
    complex(8),dimension(:,:,:,:,:)    :: ff ![Nspin,Nspin,Nimp,Nimp][Lmats]
    real(8),dimension(:),intent(inout) :: bath_
    real(8),dimension(:),allocatable   :: array_bath
    integer                            :: iter,Asize
    real(8)                            :: chi
    logical                            :: check
    character(len=256)                 :: suffix
    integer                            :: unit
    !
    call assert_shape(fg,[Nspin,Nspin,Nimp,Nimp,size(fg,5)],"chi2_fitgf_replica_superc_n4","fg")
    call assert_shape(ff,[Nspin,Nspin,Nimp,Nimp,size(fg,5)],"chi2_fitgf_replica_superc_n4","ff")
    !
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_replica_superc error: wrong bath dimensions"
    !
    if(cg_pow/=2.AND.cg_norm=="frobenius")then
       print *, "WARNING: CG_POW must be 2 for a meaningful definition of the Frobenius norm."
       print *, "         we'll still let you go ahead with the desired input, but please be "
       print *, "         be aware that CG_POW is not doing what you would expect for a chi^q"
    endif
    !
    call allocate_dmft_bath(dmft_bath_fit)
    call set_dmft_bath(bath_,dmft_bath_fit)
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
    !POSSIBLY TO BE USED AT SOME POINT (but beware that is meaningless for `FROBENIUS` norm)
    Hmask=.true.
    ! Aren't we sure about FG_ij = FG_ji, coming from Hbath being hermitian?
    !
    ! -> Hmask=Hbath_mask(wdiag=.false.,uplo=.true.) NO, this does more than
    !    that, putting .FALSE. where Hbath is zero, which may lead to tricky
    !    wrong fits: I have data showing this for a trimer, where at ~ 6.30d-10 
    !    dmft error, hence a very well converged point, we have a nonnegligible
    !    Im(Weiss) component for ilat=1, jlat=3, despite Hmask(1,3) = Hmask(3,1) 
    !    would result .FALSE. if computed with the Hbath_mask call above.
    !
    ! For now I'd say that we'd better dump everything inside the \chi^2, hence Hmask=.true.,
    ! but we might want to consider exploiting hermiticity at least (probably not looking at
    ! Hbath though: who guarantees zeros therein imply zeros here?). Discussion needed...
    !----------------------------------------------------------------------------------------
    !
    !
    FGmatrix = fg
    FFmatrix = ff
    !
    if(cg_grad==0)then
       cg_grad=1
       write(LOGfile,*)"WARNING: REPLICA Superc does not support minimization with analytical gradient (cg_grad=0)."
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
               chi2_weiss_replica_superc,&
               iter,&
               chi, &
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               istop=cg_stop, &
               iverbose=(ed_verbose>3))
       case ("delta")
          call fmin_cg(array_bath,&
               chi2_delta_replica_superc,&
               iter,&
               chi, &
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               istop=cg_stop, &
               iverbose=(ed_verbose>3))
       case default
          stop "chi2_fitgf_replica_superc error: cg_scheme != [weiss,delta]"
       end select
       !
    case (1)
       select case (cg_scheme)
       case ("weiss")
          call fmin_cgminimize(array_bath,&
               chi2_weiss_replica_superc,&
               iter,&
               chi, &
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               new_version=cg_minimize_ver,&
               hh_par=cg_minimize_hh,      &
               iverbose=(ed_verbose>3))
       case ("delta")
          call fmin_cgminimize(array_bath,&
               chi2_delta_replica_superc,&
               iter,&
               chi, &
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               new_version=cg_minimize_ver,&
               hh_par=cg_minimize_hh,      &
               iverbose=(ed_verbose>3))
       case default
          stop "chi2_fitgf_replica_superc error: cg_scheme != [weiss,delta]"
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
    bath_(2:) = array_bath
    call set_dmft_bath(bath_,dmft_bath_fit) ! *** array_bath --> dmft_bath *** (per write fit result)
    call write_dmft_bath(dmft_bath_fit,LOGfile)
    call save_dmft_bath(dmft_bath_fit)
    !
    call write_fit_result()
    !
    call get_dmft_bath(dmft_bath_fit,bath_)                ! ***  dmft_bath --> bath_ ***    (bath in output)
    call deallocate_dmft_bath(dmft_bath_fit)
    deallocate(FGmatrix,FFmatrix,Hmask,Xdelta,Wdelta)
    deallocate(array_bath)
    !
  contains
    !
    subroutine write_fit_result()
      complex(8) :: fgand(2,Nspin,Nspin,Nimp,Nimp,Ldelta)
      integer    :: gunit,funit
      if(cg_scheme=='weiss')then
         fgand = g0and_replica_superc(array_bath)
      else
         fgand = delta_replica_superc(array_bath)
      endif
      !
      do ispin=1,Nspin
         do jspin=1,Nspin
            do iorb=1,Nimp
               do jorb=1,Nimp
                  suffix="_l"//str(iorb)//&
                       "_m"//str(jorb)//&
                       "_s"//str(ispin)//&
                       "_r"//str(jspin)//reg(ed_file_suffix)
                  gunit=free_unit()
                  funit=free_unit()
                  if(cg_scheme=='weiss')then
                     open(gunit,file="fit_weiss"//reg(suffix)//".ed")
                     open(funit,file="fit_fweiss"//reg(suffix)//".ed")
                  else
                     open(gunit,file="fit_delta"//reg(suffix)//".ed")
                     open(funit,file="fit_tdelta"//reg(suffix)//".ed")
                  endif
                  do i=1,Ldelta
                     write(gunit,"(9F24.15)")Xdelta(i),&
                          dimag(fg(ispin,jspin,iorb,jorb,i)),&
                          dimag(fgand(1,ispin,jspin,iorb,jorb,i)),&
                          dreal(fg(ispin,jspin,iorb,jorb,i)),&
                          dreal(fgand(1,ispin,jspin,iorb,jorb,i))
                     write(funit,"(9F24.15)")Xdelta(i),&
                          dimag(ff(ispin,jspin,iorb,jorb,i)),&
                          dimag(fgand(2,ispin,jspin,iorb,jorb,i)),&
                          dreal(ff(ispin,jspin,iorb,jorb,i)),&
                          dreal(fgand(2,ispin,jspin,iorb,jorb,i))
                  enddo
                  close(gunit)
                  close(funit)
               enddo
            enddo
         enddo
      enddo
    end subroutine write_fit_result
    !
  end subroutine chi2_fitgf_replica_superc_n4



  include "replica/chi2_delta_replica.f90"
  include "replica/chi2_weiss_replica.f90"
  include "replica/delta_replica.f90"
  include "replica/weiss_replica.f90"

  
end MODULE ED_FIT_REPLICA

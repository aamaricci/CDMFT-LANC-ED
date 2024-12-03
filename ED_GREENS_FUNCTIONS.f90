MODULE ED_GREENS_FUNCTIONS
  USE SCIFOR, only: reg,txtfy
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_BATH
  USE ED_GF_NORMAL
  ! USE ED_GF_SUPERC
  !
  implicit none
  private




  public :: buildGf_impurity
  public :: get_Gimp
  public :: get_Sigma

contains



  !+------------------------------------------------------------------+
  ! GF CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine buildGF_impurity()
    !
    !
    write(LOGfile,"(A)")"Get impurity Greens functions:"
    select case(ed_mode)
    case default  ;call build_gf_normal()
       ! case("superc");call build_gf_superc()
    end select
    !
    if(MPIMASTER)then
       if(ed_print_Sigma)call ed_print_impSigma()
       if(ed_print_G)call ed_print_impG()
       if(ed_print_G0)call ed_print_impG0()
       call write_szr()
    endif
    !
  end subroutine buildGF_impurity






  function get_Gimp(zeta,axis) result(self)
    complex(8),dimension(:),intent(in)                                 :: zeta
    character(len=*),optional                                          :: axis
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(zeta)) :: self
    character(len=1)                                                   :: axis_
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    select case(ed_mode)
    case default  ;self = get_gimp_normal(zeta,axis_)
       ! case("superc");self = get_gimp_superc(zeta,axis_)
    end select
  end function get_Gimp




  function get_Sigma(zeta,axis) result(Self)
    complex(8),dimension(:),intent(in)                                 :: zeta
    character(len=*),optional                                          :: axis
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(zeta)) :: self
    character(len=1)                                                   :: axis_
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    select case(ed_mode)
    case default  ;self = get_Sigma_normal(zeta,axis_)
       ! case("superc");self = get_Sigma_superc(zeta,axis_)
    end select
  end function get_Sigma








  !+------------------------------------------------------------------+
  !                         PRINT SIGMA:
  !+------------------------------------------------------------------+
  subroutine ed_print_impSigma
    character(len=64)                                             :: suffix
    integer                                                       :: i,ilat,jlat,iorb,jorb,ispin,io,jo
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,Lmats) :: Smats
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,Lreal) :: Sreal
    !
    call allocate_grids
    !
    Smats = get_Sigma(dcmplx(0d0,wm(:)))
    Sreal = get_Sigma(dcmplx(wr(:),eps))
    !
    do ispin=1,Nspin
       do ilat=1,Nlat
          do jlat=1,Nlat
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ilat-1)*Norb
                   jo = jorb + (jlat-1)*Norb
                   select case(ed_mode)
                   case default
                      suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
                      call splot("impSigma"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed",wm,Smats(1,1,ispin,ispin,io,jo,:))
                      call splot("impSigma"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,Sreal(1,1,ispin,ispin,io,jo,:))
                   case("superc")
                      suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
                      call splot("impSigma"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed",wm,Smats(1,1,ispin,ispin,io,jo,:))
                      call splot("impSelf"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed",wm,Smats(1,2,ispin,ispin,io,jo,:))
                      call splot("impSigma"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,Sreal(1,1,ispin,ispin,io,jo,:))
                      call splot("impSelf"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,Sreal(1,2,ispin,ispin,io,jo,:))
                   end select
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    call deallocate_grids
    !
  end subroutine ed_print_impSigma







  subroutine ed_print_impG
    character(len=64)                                             :: suffix
    integer                                                       :: i,ilat,jlat,iorb,jorb,ispin,io,jo
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,Lmats) :: Gmats
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,Lreal) :: Greal
    !
    call allocate_grids
    !
    Gmats = get_Gimp(dcmplx(0d0,wm(:)))
    Greal = get_Gimp(dcmplx(wr(:),eps))
    !
    do ispin=1,Nspin
       do ilat=1,Nlat
          do jlat=1,Nlat
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ilat-1)*Norb
                   jo = jorb + (jlat-1)*Norb
                   select case(ed_mode)
                   case default
                      suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
                      call splot("impG"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed",wm,Gmats(1,1,ispin,ispin,io,jo,:))
                      call splot("impG"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,Greal(1,1,ispin,ispin,io,jo,:))
                   case("superc")
                      suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
                      call splot("impG"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed",wm,Gmats(1,1,ispin,ispin,io,jo,:))
                      call splot("impF"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed",wm,Gmats(1,2,ispin,ispin,io,jo,:))
                      call splot("impG"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,Greal(1,1,ispin,ispin,io,jo,:))
                      call splot("impF"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,Greal(1,2,ispin,ispin,io,jo,:))
                   end select
                enddo
             enddo
          enddo
       enddo
    enddo
    call deallocate_grids
    !
  end subroutine ed_print_impG



  !+------------------------------------------------------------------+
  !                         PRINT G0
  !+------------------------------------------------------------------+
  subroutine ed_print_impG0
    character(len=64)                                             :: suffix
    integer                                                       :: i,ilat,jlat,iorb,jorb,ispin,io,jo
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,Lmats) :: Gmats
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,Lreal) :: Greal
    !
    call allocate_grids
    !
    Gmats = g0and_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    Greal = g0and_bath_function(dcmplx(wr(:),eps),dmft_bath)
    !
    do ispin=1,Nspin
       do ilat=1,Nlat
          do jlat=1,Nlat
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ilat-1)*Norb
                   jo = jorb + (jlat-1)*Norb
                   select case(ed_mode)
                   case default
                      suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
                      call splot("impG"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed",wm,Gmats(1,1,ispin,ispin,io,jo,:))
                      call splot("impG"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,Greal(1,1,ispin,ispin,io,jo,:))
                   case("superc")
                      suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
                      call splot("impG"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed",wm,Gmats(1,1,ispin,ispin,io,jo,:))
                      call splot("impF"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed",wm,Gmats(1,2,ispin,ispin,io,jo,:))
                      call splot("impG"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,Greal(1,1,ispin,ispin,io,jo,:))
                      call splot("impF"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,Greal(1,2,ispin,ispin,io,jo,:))
                   end select
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    call deallocate_grids
    !
  end subroutine ed_print_impG0








  !+-------------------------------------------------------------------+
  !PURPOSE  : write observables to file
  !+-------------------------------------------------------------------+
  subroutine write_szr()
    integer                                                   :: unit
    integer                                                   :: iorb,jorb,ispin,ilat
    real(8)                                                   :: wm1,wm2
    complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,2) :: Smats
    complex(8)                                                :: s1,s2
    real(8),dimension(Nspin,Nimp)                             :: zimp,simp
    !
    wm1 = pi/beta ; wm2=3d0*pi/beta
    !
    select case(ed_mode)
    case default;
       Smats = get_Sigma_normal([dcmplx(0d0,wm1),dcmplx(0d0,wm2)])
    case ("superc")
       ! Smats = get_Sigma_superc([dcmplx(0d0,wm1),dcmplx(0d0,wm2)])
    end select
    !
    !
    do ispin=1,Nspin
       do iorb=1,Nimp
          s1 = Smats(1,1,ispin,ispin,iorb,iorb,1)
          s2 = Smats(1,1,ispin,ispin,iorb,iorb,2)
          simp(ispin,iorb) = dimag(s1) - &
               wm1*( dimag(s2)-dimag(s1) )/(wm2-wm1)
          zimp(iorb,ispin) = 1.d0/( 1.d0 + abs( dimag(s1)/wm1 ))
       enddo
    enddo
    !
    open(free_unit(unit),file="zeta_info.ed")
    write(unit,"(A1,90(A10,6X))")"#",&
         ((reg(txtfy(iorb+(ispin-1)*Nimp))//"z_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Nimp),ispin=1,Nspin)
    close(unit)
    !
    open(free_unit(unit),file="sig_info.ed")
    write(unit,"(A1,90(A10,6X))")"#",&
         ((reg(txtfy(iorb+(ispin-1)*Nimp))//"sig_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Nimp),ispin=1,Nspin)
    close(unit)
    !
    open(free_unit(unit),file="zeta_all"//reg(ed_file_suffix)//".ed",position='append')
    write(unit,"(90(F15.9,1X))")&
         ((zimp(ispin,iorb),iorb=1,Nimp),ispin=1,Nspin)
    close(unit)
    open(free_unit(unit),file="zeta_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90(F15.9,1X))")&
         ((zimp(ispin,iorb),iorb=1,Nimp),ispin=1,Nspin)
    close(unit)
    !
    open(free_unit(unit),file="sig_all"//reg(ed_file_suffix)//".ed",position='append')
    write(unit,"(90(F15.9,1X))")&
         ((simp(ispin,iorb),iorb=1,Nimp),ispin=1,Nspin)
    close(unit)
    open(free_unit(unit),file="sig_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90(F15.9,1X))")&
         ((simp(ispin,iorb),iorb=1,Nimp),ispin=1,Nspin)
    close(unit)
    !
  end subroutine write_szr

end MODULE ED_GREENS_FUNCTIONS

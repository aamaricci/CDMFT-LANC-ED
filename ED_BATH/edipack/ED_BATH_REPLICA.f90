MODULE ED_BATH_REPLICA
  !In these module we implement the functions to set the matrix basis :math:`\{ \hat{O}_i \}_{i=1,\dots,N_{sym}}` and the initial variational parameters :math:`\vec{\lambda}` used to decompose each local bath hamiltonian for the  :f:var:`replica` and :f:var:`general` bath types.  
  !
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,str
  USE SF_LINALG, only: eye,inv
  USE SF_ARRAYS, only:linspace
  USE SF_MISC, only: assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  !
  USE ED_BATH_AUX
  USE ED_BATH_DIM
  USE ED_BATH_DMFT
  implicit none

  private

  interface set_Hreplica
     !
     ! This function sets the matrix basis :math:`\{ \hat{O}_i \}_{i=1,\dots,N_{sym}}` used to decompose the single bath hamiltonian :math:`h^p`. It also sets the initial values of the variational parameters :math:`\vec{\lambda}` in the :f:var:`replica` ath type.
     !
     ! Input: :f:var:`Hvec` 
     !  * rank-5, dimensions: [ |Nns| , |Nns|, |Norb| , |Norb| , |Nsym| ]
     !  * rank-3, dimensions: [ |Nnso| , |Nnso|, |Nsym| ]
     !
     ! Input: :f:var:`lambdavec`
     !  * rank-2, dimensions: [ |Nbath| , |Nsym| ]
     !  * rank-3, dimensions: [ |Nlat| , |Nbath| , |Nsym| ]
     !
     module procedure init_Hreplica_symmetries_d5
     module procedure init_Hreplica_symmetries_d3
     module procedure init_Hreplica_symmetries_lattice_d5
     module procedure init_Hreplica_symmetries_lattice_d3
     module procedure init_Hreplica_symmetries_legacy
  end interface set_Hreplica

  interface set_Hgeneral
     !
     ! This function sets the matrix basis :math:`\{ \hat{O}_i \}_{i=1,\dots,N_{sym}}` used to decompose the single bath hamiltonian :math:`h^p`. It also sets the initial values of the variational parameters :math:`\vec{\lambda}` in the :f:var:`general` ath type.
     !
     ! Input: :f:var:`Hvec` 
     !  * rank-5, dimensions: [ |Nns| , |Nns|, |Norb| , |Norb| , |Nsym| ]
     !  * rank-3, dimensions: [ |Nnso| , |Nnso|, |Nsym| ]
     !
     ! Input: :f:var:`lambdavec`
     !  * rank-2, dimensions: [ |Nbath| , |Nsym| ]
     !  * rank-3, dimensions: [ |Nlat| , |Nbath| , |Nsym| ]
     !
     module procedure init_Hgeneral_symmetries_d5
     module procedure init_Hgeneral_symmetries_d3
     module procedure init_Hgeneral_symmetries_lattice_d5
     module procedure init_Hgeneral_symmetries_lattice_d3
     module procedure init_Hgeneral_symmetries_legacy
  end interface set_Hgeneral

  !##################################################################
  !
  !     REPLICA AND GENERAL BATH ROUTINES:
  !
  !##################################################################

  public :: set_Hreplica
  public :: set_Hgeneral

  integer :: ibath

contains



  !##################################################################
  !
  !     H_REPLICA ROUTINES:
  !
  !##################################################################
  !allocate GLOBAL basis for H (used for Hbath) and vectors coefficient
  subroutine allocate_hreplica(Nsym)
    integer          :: Nsym
    integer          :: isym
    !
#ifdef _DEBUG,iorb,ispin
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG allocate_Hreplica"
#endif
    if(allocated(Hreplica_basis))deallocate(Hreplica_basis)
    if(allocated(Hreplica_lambda))deallocate(Hreplica_lambda)
    !
    allocate(Hreplica_basis(Nsym))
    allocate(Hreplica_lambda(Nbath,Nsym))
    do isym=1,Nsym
       allocate(Hreplica_basis(isym)%O(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb))
       Hreplica_basis(isym)%O=zero
       Hreplica_lambda(:,isym)=0d0
    enddo
    Hreplica_status=.true.
  end subroutine allocate_hreplica


  !deallocate GLOBAL basis for H (used for impHloc and bath) and vectors coefficient
  subroutine deallocate_hreplica()
    integer              :: isym
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG deallocate_Hreplica"
#endif
    do isym=1,size(Hreplica_basis)
       deallocate(Hreplica_basis(isym)%O)
    enddo
    deallocate(Hreplica_basis)
    deallocate(Hreplica_lambda)
    Hreplica_status=.false.
  end subroutine deallocate_hreplica



  !Initialize replica bath using \vec{H} and \vec{\lambda}
  subroutine init_Hreplica_symmetries_d5(Hvec,lambdavec)
    complex(8),dimension(:,:,:,:,:) :: Hvec      ![Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym]
    real(8),dimension(:,:)          :: lambdavec ![Nbath,Nsym]
    integer                         :: isym,Nsym
    logical                         :: bool
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_symmetries_site: from {[Hs,Lam]}_b"
#endif
    !
    if(ed_mode=="superc")Nnambu=2
    Nsym=size(lambdavec,2)
    if(size(lambdavec,1)/=Nbath)then
       write(*,*) "                                                                               "
       write(*,*) "ERROR: if you are trying to init Hreplica for inequivalent sites please note   "
       write(*,*) "       that the lambdasym array /MUST/ have [Nineq]x[Nbath]x[Nsym] shape.      "
       write(*,*) "       The legacy [Nineq]x[Nsym] is not supported anymore, for it would shadow "
       write(*,*) "       the new recommended [Nbath]x[Nsym] shape for the single impurity case.  "
       write(*,*) "                                                                               "
       stop ! This unfortunately still leaves room for nasty problems if Nbath==Nineq, but that's it...
    endif
    !
    call assert_shape(Hvec,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym],"init_Hreplica_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default
          bool = check_herm(nn2so_reshape(Hvec(:,:,:,:,isym),Nspin,Norb),Nspin*Norb)
       case("superc")
          bool = check_nambu(nn2so_reshape(Hvec(:,:,:,:,isym),Nnambu*Nspin,Norb),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hreplica_symmetries_site ERROR: not Hermitian/Nambu of replica basis O_"//str(isym)
          stop
       endif
    enddo
    !
    call allocate_hreplica(Nsym)
    !
    !
    do isym=1,Nsym
       Hreplica_lambda(:,isym)  = lambdavec(:,isym)
       Hreplica_basis(isym)%O = Hvec(:,:,:,:,isym)
    enddo
    !
    if(ed_verbose>2)then
       do ibath=1,Nbath
          write(*,*) "Hreplica #"//str(ibath)//":"
          call print_hloc(Hreplica_build(Hreplica_lambda(ibath,:)))
       enddo
    endif
    !
  end subroutine init_Hreplica_symmetries_d5

  subroutine init_Hreplica_symmetries_d3(Hvec,lambdavec)
    complex(8),dimension(:,:,:) :: Hvec       ![Nnambu*Nspin*Norb,Nnambu*Nspin*Norb,Nsym]
    real(8),dimension(:,:)      :: lambdavec ![Nbath,Nsym]
    integer                     :: isym,Nsym
    logical                     :: bool
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_symmetries_site: from {[Hs,Lam]}_b"
#endif
    if(ed_mode=="superc")Nnambu=2
    Nsym=size(lambdavec,2)
    if(size(lambdavec,1)/=Nbath)then
       write(*,*) "                                                                               "
       write(*,*) "ERROR: if you are trying to init Hreplica for inequivalent sites please note   "
       write(*,*) "       that the lambdasym array /MUST/ have [Nineq]x[Nbath]x[Nsym] shape.      "
       write(*,*) "       The legacy [Nineq]x[Nsym] is not supported anymore, for it would shadow "
       write(*,*) "       the new recommended [Nbath]x[Nsym] shape for the single impurity case.  "
       write(*,*) "                                                                               "
       stop ! This unfortunately still leaves room for nasty problems if Nbath==Nineq, but that's it...
    endif
    !
    call assert_shape(Hvec,[Nnambu*Nspin*Norb,Nnambu*Nspin*Norb,Nsym],"init_Hreplica_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default
          bool = check_herm(Hvec(:,:,isym),Nspin*Norb)
       case("superc")
          bool = check_nambu(Hvec(:,:,isym),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hreplica_symmetries_site ERROR: not Hermitian/Nambu of replica basis O_"//str(isym)
          stop
       endif
    enddo
    !
    call allocate_hreplica(Nsym)
    !
    !
    do isym=1,Nsym
       Hreplica_lambda(:,isym) = lambdavec(:,isym)
       Hreplica_basis(isym)%O  = so2nn_reshape(Hvec(:,:,isym),Nnambu*Nspin,Norb)
    enddo
    !
    if(ed_verbose>2)then
       do ibath=1,Nbath
          write(LOGfile,"(A)") "Hreplica #"//str(ibath)//":"
          call print_hloc(Hreplica_build(Hreplica_lambda(ibath,:)))
       enddo
    endif
    !
  end subroutine init_Hreplica_symmetries_d3



  subroutine init_Hreplica_symmetries_lattice_d5(Hvec,lambdavec)
    complex(8),dimension(:,:,:,:,:) :: Hvec      ![size(H),Nsym]
    real(8),dimension(:,:,:)        :: lambdavec ![Nlat,Nbath,Nsym]
    integer                         :: ilat,Nlat
    integer                         :: isym,Nsym
    logical                         :: bool
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_symmetries_lattice: from ({[Hs,Lam]}_b)_site"
#endif
    if(ed_mode=="superc")Nnambu=2
    Nlat=size(lambdavec,1)
    Nsym=size(lambdavec,3)
    !
    call assert_shape(Hvec,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym],"init_Hreplica_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default
          bool = check_herm(nn2so_reshape(Hvec(:,:,:,:,isym),Nspin,Norb),Nspin*Norb)
       case("superc")
          bool = check_nambu(nn2so_reshape(Hvec(:,:,:,:,isym),Nnambu*Nspin,Norb),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hreplica_symmetries_site ERROR: not Hermitian/Nambu of replica basis O_"//str(isym)
          stop
       endif
    enddo
    !
    if(allocated(Hreplica_lambda_ineq))deallocate(Hreplica_lambda_ineq)
    allocate(Hreplica_lambda_ineq(Nlat,Nbath,Nsym))
    call allocate_hreplica(Nsym)
    !
    do isym=1,Nsym
       Hreplica_lambda_ineq(:,:,isym)  = lambdavec(:,:,isym)
       Hreplica_basis(isym)%O          = Hvec(:,:,:,:,isym)
    enddo
    !
    if(ed_verbose>2)then
       do ilat=1,Nlat
          write(LOGfile,"(A)")"Inequivalent #"//str(ilat)//":"
          do ibath=1,Nbath
             write(LOGfile,"(A)")"> Hreplica #"//str(ibath)//":"
             call print_hloc(Hreplica_build(Hreplica_lambda_ineq(ilat,ibath,:)))
          enddo
       enddo
    endif
    !
  end subroutine init_Hreplica_symmetries_lattice_d5

  subroutine init_Hreplica_symmetries_lattice_d3(Hvec,lambdavec)
    complex(8),dimension(:,:,:) :: Hvec      ![size(H),Nsym]
    real(8),dimension(:,:,:)    :: lambdavec ![Nlat,Nbath,Nsym]
    integer                     :: ilat,Nlat
    integer                     :: isym,Nsym
    logical                     :: bool
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_symmetries_lattice: from ({[Hs,Lam]}_b)_site"
#endif
    if(ed_mode=="superc")Nnambu=2
    Nlat=size(lambdavec,1)
    Nsym=size(lambdavec,3)
    !
    call assert_shape(Hvec,[Nnambu*Nspin*Norb,Nnambu*Nspin*Norb,Nsym],"init_Hreplica_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default
          bool = check_herm(Hvec(:,:,isym),Nspin*Norb)
       case("superc")
          bool = check_nambu(Hvec(:,:,isym),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hreplica_symmetries_site ERROR: not Hermitian/Nambu of replica basis O_"//str(isym)
          stop
       endif
    enddo
    !
    if(allocated(Hreplica_lambda_ineq))deallocate(Hreplica_lambda_ineq)
    allocate(Hreplica_lambda_ineq(Nlat,Nbath,Nsym))
    call allocate_hreplica(Nsym)
    !
    do isym=1,Nsym
       Hreplica_lambda_ineq(:,:,isym)  = lambdavec(:,:,isym)
       Hreplica_basis(isym)%O          = so2nn_reshape(Hvec(:,:,isym),Nnambu*Nspin,Norb)
    enddo
    !
    if(ed_verbose>2)then
       do ilat=1,Nlat
          write(LOGfile,"(A)")"Inequivalent #"//str(ilat)//":"
          do ibath=1,Nbath
             write(LOGfile,"(A)")"> Hreplica #"//str(ibath)//":"
             call print_hloc(Hreplica_build(Hreplica_lambda_ineq(ilat,ibath,:)))
          enddo
       enddo
    endif
    !
  end subroutine init_Hreplica_symmetries_lattice_d3


  subroutine init_Hreplica_symmetries_legacy(Hvec,lambdavec)
    complex(8),dimension(:,:,:,:,:) :: Hvec      ![size(H),Nsym]
    real(8),dimension(:)            :: lambdavec ![Nsym]
    integer                         :: isym,Nsym
    logical                         :: bool
    !
    if(ed_mode=="superc")Nnambu=2
    Nsym=size(lambdavec)
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_symmetries_legacy: from {[Hs,Lam]}_b"
#endif
    call assert_shape(Hvec,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym],"init_Hreplica_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default
          bool = check_herm(nn2so_reshape(Hvec(:,:,:,:,isym),Nspin,Norb),Nspin*Norb)
       case("superc")
          bool = check_nambu(nn2so_reshape(Hvec(:,:,:,:,isym),Nnambu*Nspin,Norb),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hreplica_symmetries_site ERROR: not Hermitian/Nambu of replica basis O_"//str(isym)
          stop
       endif
    enddo
    !
    call allocate_hreplica(Nsym)
    !
    do isym=1,Nsym
       do ibath=1,Nbath
          !> BACK-COMPATIBILITY PATCH (cfr init_dmft_bath)
          Hreplica_lambda(ibath,isym) = lambdavec(isym)
       enddo
       Hreplica_basis(isym)%O = Hvec(:,:,:,:,isym)
    enddo
    !
    ! PRINT DEPRECATION MESSAGE TO LOG
    write(*,*) "                                                                               "
    write(*,*) "WARNING: Passing a single lambdasym vector to ed_set_Hreplica is /deprecated/. "
    write(*,*) "         You should instead define a different lambda for each bath component, "
    write(*,*) "         namely passing a [Nbath]x[Nsym] array instead of a [Nsym] vector.     "
    write(*,*) "         Your single lambda vector has been internally copied into the required"
    write(*,*) "         higher-rank array, so giving each replica the same set of lambdas.    "
    write(*,*) "         >>> This back-compatibility patch might be removed in a future update."
    write(*,*) "                                                                               "
    !
    if(ed_verbose>2)then
       do ibath=1,Nbath
          write(*,*) "Hreplica #"//str(ibath)//":"
          call print_hloc(Hreplica_build(Hreplica_lambda(ibath,:)))
       enddo
    endif
    !
  end subroutine init_Hreplica_symmetries_legacy






  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################








  !##################################################################
  !
  !     H_GENERAL ROUTINES:
  !
  !##################################################################
  !allocate GLOBAL basis for H (used for Hbath) and vectors coefficient
  subroutine allocate_hgeneral(Nsym)
    integer          :: Nsym
    integer          :: isym
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG allocate_Hgeneral"
#endif
    if(allocated(Hgeneral_basis))deallocate(Hgeneral_basis)
    if(allocated(Hgeneral_lambda))deallocate(Hgeneral_lambda)
    !
    allocate(Hgeneral_basis(Nsym))
    allocate(Hgeneral_lambda(Nbath,Nsym))
    do isym=1,Nsym
       allocate(Hgeneral_basis(isym)%O(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb))
       Hgeneral_basis(isym)%O=zero
       Hgeneral_lambda(:,isym)=0d0
    enddo
    Hgeneral_status=.true.
  end subroutine allocate_hgeneral


  !deallocate GLOBAL basis for H (used for impHloc and bath) and vectors coefficient
  subroutine deallocate_hgeneral()
    integer              :: isym
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG deallocate_Hgeneral"
#endif
    do isym=1,size(Hgeneral_basis)
       deallocate(Hgeneral_basis(isym)%O)
    enddo
    deallocate(Hgeneral_basis)
    deallocate(Hgeneral_lambda)
    Hgeneral_status=.false.
  end subroutine deallocate_hgeneral





  subroutine init_Hgeneral_symmetries_d5(Hvec,lambdavec)
    complex(8),dimension(:,:,:,:,:) :: Hvec      ![Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym]
    real(8),dimension(:,:)          :: lambdavec ![Nbath,Nsym]
    integer                         :: isym,Nsym
    logical                         :: bool
    !
    if(ed_mode=="superc")Nnambu=2
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hgeneral_symmetries_site: from {[Hs,Lam]}_b"
#endif
    !
    if(size(lambdavec,1)/=Nbath)then
       write(*,*) "                                                                               "
       write(*,*) "ERROR: if you are trying to init Hgeneral for inequivalent sites please note   "
       write(*,*) "       that the lambdasym array /MUST/ have [Nineq]x[Nbath]x[Nsym] shape.      "
       write(*,*) "       The legacy [Nineq]x[Nsym] is not supported anymore, for it would shadow "
       write(*,*) "       the new recommended [Nbath]x[Nsym] shape for the single impurity case.  "
       write(*,*) "                                                                               "
       stop ! This unfortunately still leaves room for nasty problems if Nbath==Nineq, but that's it...
    else
       Nsym=size(lambdavec,2)
    endif
    !
    call assert_shape(Hvec,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym],"init_Hgeneral_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default
          bool = check_herm(nn2so_reshape(Hvec(:,:,:,:,isym),Nspin,Norb),Nspin*Norb)
       case("superc")
          bool = check_nambu(nn2so_reshape(Hvec(:,:,:,:,isym),Nnambu*Nspin,Norb),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hgeneral_symmetries_site ERROR: not Hermitian/Nambu of general basis O_"//str(isym)
          stop
       endif
    enddo
    !
    call allocate_hgeneral(Nsym)
    !
    !
    do isym=1,Nsym
       Hgeneral_lambda(:,isym)  = lambdavec(:,isym)
       Hgeneral_basis(isym)%O   = Hvec(:,:,:,:,isym)
    enddo
    !
    if(ed_verbose>2)then
       do ibath=1,Nbath
          write(LOGfile,"(A)") "Hgeneral #"//str(ibath)//":"
          call print_hloc(Hgeneral_build(Hgeneral_lambda(ibath,:)))
       enddo
    endif
    !
  end subroutine init_Hgeneral_symmetries_d5


  subroutine init_Hgeneral_symmetries_d3(Hvec,lambdavec)
    complex(8),dimension(:,:,:) :: Hvec      ![Nnambu*Nspin*Norb,Nnambu*Nspin*Norb,Nsym]
    real(8),dimension(:,:)      :: lambdavec ![Nbath,Nsym]
    integer                     :: isym,Nsym
    logical                     :: bool
    !

#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hgeneral_symmetries_site: from {[Hs,Lam]}_b"
#endif
    if(ed_mode=="superc")Nnambu=2
    Nsym=size(lambdavec,2)
    if(size(lambdavec,1)/=Nbath)then
       write(*,*) "                                                                               "
       write(*,*) "ERROR: if you are trying to init Hgeneral for inequivalent sites please note   "
       write(*,*) "       that the lambdasym array /MUST/ have [Nineq]x[Nbath]x[Nsym] shape.      "
       write(*,*) "       The legacy [Nineq]x[Nsym] is not supported anymore, for it would shadow "
       write(*,*) "       the new recommended [Nbath]x[Nsym] shape for the single impurity case.  "
       write(*,*) "                                                                               "
       stop ! This unfortunately still leaves room for nasty problems if Nbath==Nineq, but that's it...
    endif
    !
    call assert_shape(Hvec,[Nnambu*Nspin*Norb,Nnambu*Nspin*Norb,Nsym],"init_Hgeneral_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default
          bool = check_herm(Hvec(:,:,isym),Nspin*Norb)
       case("superc")
          bool = check_nambu(Hvec(:,:,isym),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hgeneral_symmetries_site ERROR: not Hermitian/Nambu of general basis O_"//str(isym)
          stop
       endif
    enddo
    !
    call allocate_hgeneral(Nsym)
    !
    !
    do isym=1,Nsym
       Hgeneral_lambda(:,isym) = lambdavec(:,isym)
       Hgeneral_basis(isym)%O  = so2nn_reshape(Hvec(:,:,isym),Nnambu*Nspin,Norb)
    enddo
    !
    if(ed_verbose>2)then
       do ibath=1,Nbath
          write(LOGfile,"(A)") "Hgeneral #"//str(ibath)//":"
          call print_hloc(Hgeneral_build(Hgeneral_lambda(ibath,:)))
       enddo
    endif
    !
  end subroutine init_Hgeneral_symmetries_d3


  subroutine init_Hgeneral_symmetries_lattice_d5(Hvec,lambdavec)
    complex(8),dimension(:,:,:,:,:) :: Hvec      ![size(H),Nsym]
    real(8),dimension(:,:,:)        :: lambdavec ![Nlat,Nbath,Nsym]
    integer                         :: ilat,Nlat
    integer                         :: isym,Nsym
    logical                         :: bool
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hgeneral_symmetries_lattice: from ({[Hs,Lam]}_b)_site"
#endif
    if(ed_mode=="superc")Nnambu=2
    Nlat=size(lambdavec,1)
    Nsym=size(lambdavec,3)
    call assert_shape(Hvec,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym],"init_Hgeneral_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default
          bool = check_herm(nn2so_reshape(Hvec(:,:,:,:,isym),Nspin,Norb),Nspin*Norb)
       case("superc")
          bool = check_nambu(nn2so_reshape(Hvec(:,:,:,:,isym),Nnambu*Nspin,Norb),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hgeneral_symmetries_site ERROR: not Hermitian/Nambu of general basis O_"//str(isym)
          stop
       endif
    enddo
    !
    if(allocated(Hgeneral_lambda_ineq))deallocate(Hgeneral_lambda_ineq)
    allocate(Hgeneral_lambda_ineq(Nlat,Nbath,Nsym))
    call allocate_hgeneral(Nsym)
    !
    do isym=1,Nsym
       Hgeneral_lambda_ineq(:,:,isym)  = lambdavec(:,:,isym)
       Hgeneral_basis(isym)%O = Hvec(:,:,:,:,isym)
    enddo
    !
    if(ed_verbose>2)then
       do ilat=1,Nlat
          write(LOGfile,"(A)")"Inequivalent #"//str(ilat)//":"
          do ibath=1,Nbath
             write(LOGfile,"(A)")"> Hgeneral #"//str(ibath)//":"
             call print_hloc(Hgeneral_build(Hgeneral_lambda_ineq(ilat,ibath,:)))
          enddo
       enddo
    endif
    !
  end subroutine init_Hgeneral_symmetries_lattice_d5

  subroutine init_Hgeneral_symmetries_lattice_d3(Hvec,lambdavec)
    complex(8),dimension(:,:,:) :: Hvec      ![size(H),Nsym]
    real(8),dimension(:,:,:)    :: lambdavec ![Nlat,Nbath,Nsym]
    integer                     :: ilat,Nlat
    integer                     :: isym,Nsym
    logical                     :: bool
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hgeneral_symmetries_lattice: from ({[Hs,Lam]}_b)_site"
#endif
    if(ed_mode=="superc")Nnambu=2
    Nlat=size(lambdavec,1)
    Nsym=size(lambdavec,3)
    call assert_shape(Hvec,[Nnambu*Nspin*Norb,Nnambu*Nspin*Norb,Nsym],"init_Hgeneral_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default
          bool = check_herm(Hvec(:,:,isym),Nspin*Norb)
       case("superc")
          bool = check_nambu(Hvec(:,:,isym),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hgeneral_symmetries_site ERROR: not Hermitian/Nambu of general basis O_"//str(isym)
          stop
       endif
    enddo
    !
    if(allocated(Hgeneral_lambda_ineq))deallocate(Hgeneral_lambda_ineq)
    allocate(Hgeneral_lambda_ineq(Nlat,Nbath,Nsym))
    call allocate_hgeneral(Nsym)
    !
    do isym=1,Nsym
       Hgeneral_lambda_ineq(:,:,isym)  = lambdavec(:,:,isym)
       Hgeneral_basis(isym)%O          = so2nn_reshape(Hvec(:,:,isym),Nnambu*Nspin,Norb)
    enddo
    !
    if(ed_verbose>2)then
       do ilat=1,Nlat
          write(LOGfile,"(A)")"Inequivalent #"//str(ilat)//":"
          do ibath=1,Nbath
             write(LOGfile,"(A)")"> Hgeneral #"//str(ibath)//":"
             call print_hloc(Hgeneral_build(Hgeneral_lambda_ineq(ilat,ibath,:)))
          enddo
       enddo
    endif
    !
  end subroutine init_Hgeneral_symmetries_lattice_d3


  subroutine init_Hgeneral_symmetries_legacy(Hvec,lambdavec)
    complex(8),dimension(:,:,:,:,:) :: Hvec      ![size(H),Nsym]
    real(8),dimension(:)            :: lambdavec ![Nsym]
    integer                         :: isym,Nsym
    logical                         :: bool
    !
    if(ed_mode=="superc")Nnambu=2
    Nsym=size(lambdavec)
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hgeneral_symmetries_legacy: from {[Hs,Lam]}_b"
#endif
    call assert_shape(Hvec,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym],"init_Hgeneral_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default
          bool = check_herm(nn2so_reshape(Hvec(:,:,:,:,isym),Nspin,Norb),Nspin*Norb)
       case("superc")
          bool = check_nambu(nn2so_reshape(Hvec(:,:,:,:,isym),Nnambu*Nspin,Norb),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hgeneral_symmetries_site ERROR: not Hermitian/Nambu of general basis O_"//str(isym)
          stop
       endif
    enddo
    !
    call allocate_hgeneral(Nsym)
    !
    do isym=1,Nsym
       do ibath=1,Nbath
          !> BACK-COMPATIBILITY PATCH (cfr init_dmft_bath)
          Hgeneral_lambda(ibath,isym) = lambdavec(isym)
       enddo
       Hgeneral_basis(isym)%O = Hvec(:,:,:,:,isym)
    enddo
    !
    ! PRINT DEPRECATION MESSAGE TO LOG
    write(*,*) "                                                                               "
    write(*,*) "WARNING: Passing a single lambdasym vector to ed_set_Hgeneral is /deprecated/. "
    write(*,*) "         You should instead define a different lambda for each bath component, "
    write(*,*) "         namely passing a [Nbath]x[Nsym] array instead of a [Nsym] vector.     "
    write(*,*) "         Your single lambda vector has been internally copied into the required"
    write(*,*) "         higher-rank array, so giving each general the same set of lambdas.    "
    write(*,*) "         >>> This back-compatibility patch might be removed in a future update."
    write(*,*) "                                                                               "
    !
    if(ed_verbose>2)then
       do ibath=1,Nbath
          write(*,*) "Hgeneral #"//str(ibath)//":"
          call print_hloc(Hgeneral_build(Hgeneral_lambda(ibath,:)))
       enddo
    endif
    !
  end subroutine init_Hgeneral_symmetries_legacy




END MODULE ED_BATH_REPLICA

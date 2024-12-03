MODULE ED_OBSERVABLES
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_AUX_FUNX
  !
  USE ED_OBSERVABLES_NORMAL
  ! USE ED_OBSERVABLES_SUPERC
  !
  implicit none
  private 

  public :: observables_impurity
  public :: local_energy_impurity
  public :: density_matrix_impurity

contains


  subroutine observables_impurity()
    ! 
    ! Calculate the local observables calling the correct procedure according to the value of :f:var:`ed_mode` .
    ! Write the results on plain-text files.
    !
    ! * :code:`normal` : :f:func:`observables_normal`
    ! * :code:`superc` : :f:func:`observables_superc`
    !
    write(LOGfile,"(A)")"Get observables:"
    select case(ed_mode)
    case default  ;call observables_normal()
    ! case("superc");call observables_superc()
    end select
  end subroutine observables_impurity


  subroutine local_energy_impurity()
    ! 
    ! Calculate the local energy calling the correct procedure according to the value of :f:var:`ed_mode` .
    ! Write the results on plain-text files.
    !
    ! * :code:`normal` : :f:func:`local_energy_normal`
    ! * :code:`superc` : :f:func:`local_energy_superc`
    !
    write(LOGfile,"(A)")"Get local energy:"
    select case(ed_mode)
    case default  ;call local_energy_normal()
    ! case("superc");call local_energy_superc()
    end select
  end subroutine local_energy_impurity



  subroutine density_matrix_impurity()
    ! 
    ! Calculate the impurity reduced density matrix calling the correct procedure according to the value of :f:var:`ed_mode` .
    !
    ! * :code:`normal` : :f:func:`density_matrix_normal`
    ! * :code:`superc` : :f:func:`density_matrix_superc`
    !
    write(LOGfile,"(A)")"Get local energy:"
    select case(ed_mode)
    case default  ;call density_matrix_normal()
       ! case("superc");call density_matrix_superc()
    end select
  end subroutine density_matrix_impurity

end MODULE ED_OBSERVABLES

MODULE ED_DIAG
  USE ED_INPUT_VARS
  USE ED_DIAG_NORMAL
  USE ED_DIAG_NONSU2
  !
  implicit none
  private

  public  :: diagonalize_impurity

contains

  subroutine  diagonalize_impurity()
    !
    ! Call the correct impurity diagonalization procedure according to the value of :f:var:`ed_mode`.
    !
    ! * :f:var:`normal` : :f:var:`diagonalize_impurity_normal`
    ! * :f:var:`superc` : :f:var:`diagonalize_impurity_superc`
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG diagonalize_impurity: Start digonalization"
#endif
    !
    write(LOGfile,"(A)")"Diagonalize impurity problem:"
    select case(ed_mode)
    case default  ;call diagonalize_impurity_normal()
    case("superc");call diagonalize_impurity_superc()
    end select
#ifdef _DEBUG
    write(Logfile,"(A)")""
#endif
  end subroutine diagonalize_impurity

end MODULE ED_DIAG

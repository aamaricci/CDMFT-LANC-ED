MODULE ED_HAMILTONIAN_SUPERC_COMMON
  USE SF_MISC,      only: assert_shape
  USE SF_CONSTANTS, only: zero
  USE SF_LINALG,    only: kronecker_product,zeye
  USE SF_SP_LINALG, only: sp_lanc_tridiag
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_SECTOR
  USE ED_SETUP
  implicit none

  type(sector)       :: Hsector

  integer                                   :: Dim
  integer                                   :: i,iup,idw
  integer                                   :: m,mup,mdw
  integer                                   :: ishift,ishift_up,ishift_dw
  integer                                   :: j,ms
  integer                                   :: iorb,jorb,ispin,jspin,ibath
  integer                                   :: kp,k1,k2,k3,k4
  integer                                   :: ialfa,ibeta
  real(8)                                   :: sg1,sg2,sg3,sg4
  complex(8)                                :: htmp,htmpup,htmpdw
  logical                                   :: Jcondition
  integer                                   :: Nfoo
  complex(8),dimension(:,:,:),allocatable   :: diag_hybr ![Nspin,Norb,Nbath]
  complex(8),dimension(:,:,:,:),allocatable :: bath_diag ![Nspin,Norb/1,Nbath]
  integer,save,public                       :: iter=0
  

end MODULE ED_HAMILTONIAN_SUPERC_COMMON





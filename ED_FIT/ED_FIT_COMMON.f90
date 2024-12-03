MODULE ED_FIT_COMMON
  USE SF_CONSTANTS
  USE SF_OPTIMIZE, only:fmin_cg,fmin_cgplus,fmin_cgminimize
  USE SF_LINALG,   only:eye,zeye,inv,inv_her,operator(.x.),trace
  USE SF_IOTOOLS,  only:reg,free_unit,txtfy
  USE SF_ARRAYS,   only:arange
  USE SF_MISC,     only:assert_shape
  !
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL  
  USE ED_AUX_FUNX
  USE ED_BATH
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif

  implicit none


  integer                                     :: Ldelta
  complex(8),dimension(:,:,:,:,:),allocatable :: FGmatrix
  complex(8),dimension(:,:,:,:,:),allocatable :: FFmatrix
  logical(8),dimension(:,:,:,:),allocatable   :: Hmask
  complex(8),dimension(:,:),allocatable       :: Gdelta
  complex(8),dimension(:,:),allocatable       :: Fdelta
  real(8),dimension(:),allocatable            :: Xdelta,Wdelta
  integer                                     :: totNorb,totNspin,totNso
  integer,dimension(:),allocatable            :: getIorb,getJorb
  integer,dimension(:),allocatable            :: getIspin,getJspin
  integer,dimension(:),allocatable            :: getInambu
  integer,dimension(:),allocatable            :: getIlat,getJlat
  integer                                     :: Orb_indx,Spin_indx,Spin_mask
  type(effective_bath)                        :: chi2_bath
  integer                                     :: cg_iter_count=0
  logical                                     :: para_
  
  !location of the maximum of the chisquare over Nlso.
  integer                                       :: maxchi_loc


  !This contains the number of the lambda expansion
  !for each replica of the impurity. which are notable the same number so far.
  integer                                       :: Nlambdas
  ! integer                                     :: Nlambdas
  !
  ! THIS IS NOT NEEDED ANYMORE I GUESS, BECAUSE NLAMBDAS is FIXED NOW
  ! WE CAN ALLOCATE A RANK-2 ARRAY *dummy_lambda(Nbath,Nlambdas)*
  !
  ! !This is a dummy object which is used here to point
  ! !to the replica bath lambdas, i.e. the coefficients
  ! !of the bath item-th Hamiltonian expansion 
  ! type nsymm_vector
  !    real(8),dimension(:),allocatable         :: element          
  ! end type nsymm_vector


END MODULE ED_FIT_COMMON

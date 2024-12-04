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
  integer                                     :: maxchi_loc

  !This contains the number of the lambda expansion
  !for each replica of the impurity. which are notable the same number so far.
  integer                                     :: Nlambdas

  !This is an internal istance of the effective bath to avoid
  !overwriting the allocated internally, which is used to retrieve
  !Gimp and Sigma on the fly
  type(effective_bath)                        :: dmft_bath_fit

END MODULE ED_FIT_COMMON

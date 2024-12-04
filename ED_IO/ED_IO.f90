MODULE ED_IO
  USE SCIFOR,only: linspace,arange,str,reg,free_unit,splot,sread,assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_SETUP
  USE ED_BATH
  USE ED_BATH_FUNCTIONS
  USE ED_GREENS_FUNCTIONS
  implicit none
  private





  !Retrieve self-energy through routines:
  interface ed_get_sigma
     !| This subrotine gets from the EDIpack2 library the value of the self-energy calculated 
     ! on the Matsubara or real-frequency axis, with number of frequencies :f:var:`lmats` or :f:var:`lreal` .
     !| The self-energy is an array having the following possible dimensions:
     !
     !  * [:f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
     !  * [:f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
     !
     module procedure :: ed_get_sigma_site_n2
     module procedure :: ed_get_sigma_site_n4
     module procedure :: ed_get_sigma_site_n6
  end interface ed_get_sigma



  interface ed_get_gimp
     !This subroutine gets from the EDIpack2 library the value of the impurity Green's function calculated 
     !on the Matsubara or real-frequency axis, with number of frequencies :f:var:`lmats` or :f:var:`lreal` .
     !
     !The impurity Green's function is an array having the following possible dimensions:
     !
     !  * [:f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
     !  * [:f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
     !
     module procedure :: ed_get_gimp_site_n2
     module procedure :: ed_get_gimp_site_n4
     module procedure :: ed_get_gimp_site_n6
  end interface ed_get_gimp




  interface ed_get_g0imp
     !| This subroutine gets from the EDIpack2 library the value of the impurity non-interacting Green's function calculated 
     ! on the Matsubara or real-frequency axis, with number of frequencies :f:var:`lmats` or :f:var:`lreal` .
     !
     !The impurity non-interacting Green's function is an array having the following possible dimensions:
     ! 
     !  * [:f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]  
     !  * [:f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
     !
     !The bath is an array having the following dimension:
     !
     !  * [:f:var:`nb`] for single-impurity DMFT
     !
     !Where :f:var:`nb` is the length of the :f:var:`bath` array.
     !
     module procedure :: ed_get_g0imp_site_n2
     module procedure :: ed_get_g0imp_site_n4
     module procedure :: ed_get_g0imp_site_n6
  end interface ed_get_g0imp




  !Build Gand/Delta from a user bath
  interface ed_get_g0and
     !| This subroutine returns to the user the normal non-interacting Green's function :math:`G_0(x)` and
     ! the anomalous non-interacting Green's function :math:`F_0(x)` on a given set of frequencies. It does so
     ! by calling :f:func:`g0and_bath_function` and :f:func:`g0and_bath_function`.
     !
     !The non-interacting Green's function is an array having the following possible dimensions:
     ! 
     !  * [:f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`, :code:`size(x)`]  
     !  * [:f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :code:`size(x)`]
     !
     module procedure :: ed_get_g0and_n2
     module procedure :: ed_get_g0and_n4
     module procedure :: ed_get_g0and_n6
  end interface ed_get_g0and

  interface ed_get_delta
     !| This subroutine returns to the user the normal hybridization function :math:`\Delta(x)` and
     ! the anomalous hybridization function :math:`\Theta(x)` on a given set of frequencies. It does so
     ! by calling :f:func:`delta_bath_function` and :f:func:`fdelta_bath_function`.
     !
     !The hybridization function is an array having the following possible dimensions:
     ! 
     !  * [:f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`, :code:`size(x)`]  
     !  * [:f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :code:`size(x)`]
     !
     module procedure :: ed_get_delta_n2
     module procedure :: ed_get_delta_n4
     module procedure :: ed_get_delta_n6
  end interface ed_get_delta



  !****************************************************************************************!
  !****************************************************************************************!



  !Observables
  interface ed_get_dens
     !This subroutine gets from the EDIpack2 library the value of the charge density and passes it to the user.
     !
     !The :f:var:`self` variable can have the following dimensions:
     ! 
     !  * scalar: if :f:var:`iorb` is provided for single-impurity DMFT, density for that orbital
     !  * [:f:var:`norb`]: if no optional variable is provided for single-impurity DMFT, density for all orbitals
     !  * [:f:var:`nlat`, :f:var:`norb`]: if :f:var:`nlat` is provided for real-space DMFT, density for all impurity sites and orbitals
     !
     module procedure :: ed_get_dens_d0
     module procedure :: ed_get_dens_d1
     module procedure :: ed_get_dens_d2
  end interface ed_get_dens

  interface ed_get_mag
     !This subroutine gets from the EDIpack2 library the value of the magnetization and passes it to the user.
     !
     !The :f:var:`self` variable can have the following dimensions:
     ! 
     !  * scalar: if :f:var:`component` and :f:var:`iorb` are provided for single-impurity DMFT, given magnetization component for that orbital
     !  * [:f:var:`norb`]: for single-impurity DMFT, one magnetization component for all orbitals
     !  * [:f:var:`nlat`, :f:var:`norb`]: if :f:var:`nlat` is provided for real-space DMFT, one magnetization component for all orbitals and impurity sites
     !  * [:f:var:`nlat`, :code:`3`, :f:var:`norb`]: if :f:var:`nlat` is provided for real-space DMFT, all magnetization components for all orbitals and sites
     !
     module procedure :: ed_get_mag_d0
     module procedure :: ed_get_mag_d1
     module procedure :: ed_get_mag_d2
  end interface ed_get_mag

  interface ed_get_docc
     !This subroutine gets from the EDIpack2 library the value of the double occupation and passes it to the user.
     !
     !The :f:var:`self` variable can have the following dimensions:
     ! 
     !  * scalar: if :f:var:`iorb` is provided for single-impurity DMFT, dobule-occupation for that orbital
     !  * [:f:var:`norb`]: if no optional variable is provided for single-impurity DMFT, double-occupation for all orbitals
     !  * [:f:var:`nlat`, :f:var:`norb`]: if :f:var:`nlat` is provided for real-space DMFT, double-occupation for all impurity sites and orbitals
     !
     module procedure :: ed_get_docc_d0
     module procedure :: ed_get_docc_d1
     module procedure :: ed_get_docc_d2
  end interface ed_get_docc

  interface ed_get_phi
     !This subroutine gets from the EDIpack2 library the value of the superconducting order parameter :math:`\phi` ( :f:var:`ed_mode` = :code:`superc` ) and passes it to the user.
     !
     !The :f:var:`self` variable can have the following dimensions:
     ! 
     !  * scalar: if :f:var:`iorb` is provided for single-impurity DMFT, :math:`\phi` for that orbital
     !  * [:f:var:`norb`]: if no optional variable is provided for single-impurity DMFT, :math:`\phi` for all orbitals
     !  * [:f:var:`nlat`, :f:var:`norb`]: if :f:var:`nlat` is provided for real-space DMFT, :math:`\phi` for all impurity sites and orbitals
     !
     module procedure :: ed_get_phi_d0
     module procedure :: ed_get_phi_d1
     module procedure :: ed_get_phi_d2
  end interface ed_get_phi


  !****************************************************************************************!
  !****************************************************************************************!







  interface ed_get_epot
     !This subroutine gets from the EDIpack2 library and passes to the user the value of 
     !:f:var:`ed_epot`, the energy contribution from the interaction terms, **including** the Hartree term.
     !The returned array can have the following dimensions:
     !
     !  * scalar: for single-site DMFT
     !
     module procedure :: ed_get_epot_n0
  end interface ed_get_epot

  interface ed_get_eint
     !This subroutine gets from the EDIpack2 library and passes to the user the value of 
     !:f:var:`ed_int`, the energy contribution from the interaction terms, **excluding** the Hartree term.
     !The returned array can have the following dimensions:
     !
     !  * scalar: for single-site DMFT
     !
     module procedure :: ed_get_eint_n0
  end interface ed_get_eint

  interface ed_get_ehartree
     !This subroutine gets from the EDIpack2 library and passes to the user the value of the Hartree potential 
     !:f:var:`ed_ehartree`. The returned array can have the following dimensions:
     !
     !  * scalar: for single-site DMFT
     !
     module procedure :: ed_get_ehartree_n0
  end interface ed_get_ehartree

  interface ed_get_eknot
     !This subroutine gets from the EDIpack2 library and passes to the user the value
     !:f:var:`ed_eknot`, the kinetic term from the **local** 1-body Hamiltonian
     !The returned array can have the following dimensions:
     !
     !  * scalar: for single-site DMFT
     !
     module procedure :: ed_get_eknot_n0
  end interface ed_get_eknot


  !Get Energies
  interface ed_get_eimp
     !This subroutine gets from the EDIpack2 library and passes to the user the array [ :f:var:`ed_epot` , :f:var:`ed_eint` , :f:var:`ed_ehartree` , :f:var:`ed_eknot` ].
     !These are the expectation values various contribution to the internal energy
     !
     !  * :f:var:`ed_epot` = energy contribution from the interaction terms, **including** the Hartree term
     !  * :f:var:`ed_eint` = energy contribution from the interaction terms, **excluding** the Hartree term
     !  * :f:var:`ed_ehartree` = :math:`-\frac{U}{2} \sum_{i} \langle n_{i\uparrow} + n_{i\downarrow} \rangle 
     !    -\frac{2U^{'}-J_{H}}{2} \sum_{i < j} \langle n_{i\uparrow}+n_{i\downarrow} + n_{i\downarrow}+n_{j\downarrow} \rangle
     !    +\frac{U}{4} + \frac{2U^{'}-J_{H}}{2}` for :math:`i,j` orbitals
     !  * :f:var:`ed_eknot` = kinetic term from the **local** 1-body Hamiltonian
     !
     !The returned array can have the following dimensions:
     !
     !  * [:code:`4`]: for single-site DMFT
     !
     module procedure :: ed_get_eimp_n1
  end interface ed_get_eimp



  !****************************************************************************************!
  !****************************************************************************************!





  interface ed_get_dust
     !This subroutine gets from the EDIpack2 library and passes to the user the value of 
     !:f:var:`ed_dust` = :math:`\sum_{i < j} n_{i\uparrow}n_{j\downarrow} + n_{i\downarrow}n_{j\uparrow}` for :math:`i,j` orbitals
     !The returned array can have the following dimensions:
     !
     !  * scalar: for single-site DMFT
     !
     module procedure :: ed_get_dust_n0
  end interface ed_get_dust

  interface ed_get_dund
     !This subroutine gets from the EDIpack2 library and passes to the user the value of 
     !:f:var:`ed_dund` = :math:`\sum_{i < j} n_{i\uparrow}n_{j\uparrow}  + n_{i\downarrow}n_{j\downarrow}` for :math:`i,j` orbitals
     !The returned array can have the following dimensions:
     !
     !  * scalar: for single-site DMFT
     !
     module procedure :: ed_get_dund_n0
  end interface ed_get_dund

  interface ed_get_dse
     !This subroutine gets from the EDIpack2 library and passes to the user the value of 
     !:f:var:`ed_dse` = :math:`\sum_{i < j} c^{\dagger}_{i\uparrow}c^{\dagger}_{j\uparrow}c_{i\downarrow}c_{j\uparrow}` for :math:`i,j` orbitals
     !The returned array can have the following dimensions:
     !
     !  * scalar: for single-site DMFT
     !
     module procedure :: ed_get_dse_n0
  end interface ed_get_dse

  interface ed_get_dph
     !This subroutine gets from the EDIpack2 library and passes to the user the value of 
     !:f:var:`ed_dph` = :math:`\sum_{i < j} c^{\dagger}_{i\uparrow}c^{\dagger}_{i\downarrow}c_{j\downarrow}c_{j\uparrow}` for :math:`i,j` orbitals
     !The returned array can have the following dimensions:
     !
     !  * scalar: for single-site DMFT
     !
     module procedure :: ed_get_dph_n0
  end interface ed_get_dph

  interface ed_get_doubles
     !This subroutine gets from the EDIpack2 library and passes to the user the array [ :f:var:`ed_dust` , :f:var:`ed_dund` , :f:var:`ed_dse` , :f:var:`ed_dph` ].
     !These are the expectation values of the two-body operators associated with the density-density inter-orbital interaction (with opposite and parallel spins), 
     !spin-exchange and pair-hopping.
     !
     !  * :f:var:`ed_dust` = :math:`\sum_{i < j} n_{i\uparrow}n_{j\downarrow} + n_{i\downarrow}n_{j\uparrow}` for :math:`i,j` orbitals
     !  * :f:var:`ed_dund` = :math:`\sum_{i < j} n_{i\uparrow}n_{j\uparrow}  + n_{i\downarrow}n_{j\downarrow}` for :math:`i,j` orbitals
     !  * :f:var:`ed_dse` = :math:`\sum_{i < j} c^{\dagger}_{i\uparrow}c^{\dagger}_{j\uparrow}c_{i\downarrow}c_{j\uparrow}` for :math:`i,j` orbitals
     !  * :f:var:`ed_dph` = :math:`\sum_{i < j} c^{\dagger}_{i\uparrow}c^{\dagger}_{i\downarrow}c_{j\downarrow}c_{j\uparrow}` for :math:`i,j` orbitals
     !
     !The returned array can have the following dimensions:
     !
     !  * [:code:`4`]: for single-site DMFT
     !
     module procedure :: ed_get_doubles_n1
  end interface ed_get_doubles

  !****************************************************************************************!
  !****************************************************************************************!

  interface ed_get_cluster_dm
     module procedure :: ed_get_cluster_density_matrix_single
  end interface ed_get_cluster_dm

  !>>DA RIVEDERE<<
  interface ed_get_reduced_dm
     module procedure :: ed_get_reduced_density_matrix_single
     module procedure :: ed_get_reduced_density_matrix_single_LEGACY
  end interface ed_get_reduced_dm


  interface ed_get_sp_dm
     module procedure :: ed_get_single_particle_density_matrix_single
  end interface ed_get_sp_dm

  ! !>>DA RIVEDERE<<
  ! interface ed_gf_cluster
  !    module procedure :: ed_gf_cluster_scalar
  !    module procedure :: ed_gf_cluster_array
  ! end interface ed_gf_cluster


  ! interface ed_read_impsigma
  !    module procedure :: ed_read_impsigma_single
  ! end interface ed_read_impsigma

  interface ed_print_dm
     module procedure :: ed_print_dm_orb
     module procedure :: ed_print_dm_LEGACY
  end interface ed_print_dm


  public :: ed_get_sigma
  public :: ed_get_gimp
  public :: ed_get_g0imp
  public :: ed_get_g0and
  public :: ed_get_delta

  ! public :: ed_build_gimp
  ! public :: ed_build_sigma

  public :: ed_get_dens
  public :: ed_get_mag
  public :: ed_get_docc
  public :: ed_get_phi
  public :: ed_get_eimp
  public :: ed_get_epot
  public :: ed_get_eint 
  public :: ed_get_ehartree
  public :: ed_get_eknot
  public :: ed_get_doubles
  public :: ed_get_dust
  public :: ed_get_dund
  public :: ed_get_dse
  public :: ed_get_dph

  public :: ed_get_cluster_dm
  public :: ed_get_reduced_dm
  public :: ed_get_sp_dm
  public :: ed_print_dm



  !Frequency and time arrays:
  !=========================================================
  character(len=64)                :: suffix



  integer                                         :: ilat,jlat
  integer                                         :: iorb,jorb
  integer                                         :: ispin,jspin
  integer                                         :: is,js
  integer                                         :: io,jo
  integer                                         :: i,j
  integer                                         :: L
  complex(8),dimension(:,:,:,:,:,:,:),allocatable :: F


contains



  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity functions
  !+--------------------------------------------------------------------------+!
  include "get_sigma.f90"
  include "get_gimp.f90"
  include "get_g0imp.f90"
  include "get_gand_bath.f90"


  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the local observables
  !+--------------------------------------------------------------------------+!
  include "observables/get_dens.f90"
  include "observables/get_mag.f90"
  include "observables/get_docc.f90"
  include "observables/get_phi.f90"
  include "observables/get_energy.f90"
  include "observables/get_doubles.f90"

  !+--------------------------------------------------------------------------+!
  ! PURPOSE: RDMs
  !+--------------------------------------------------------------------------+!
  include "rdm/get_rdm.f90"
  include "rdm/get_sp_dm.f90"



  !+------------------------------------------------------------------+
  !                      PRINT DENSITY MATRICES
  !+------------------------------------------------------------------+
  subroutine ed_print_dm_orb(dm,orbital_mask,ineq)
    complex(8),dimension(:,:),intent(in)          :: dm
    logical,dimension(Nimp),intent(in)            :: orbital_mask
    integer                  ,intent(in),optional :: ineq
    integer                                       :: unit,Nsites
    character(len=64)                             :: fname,suffix
    integer                                       :: ilat,iorb,Nrdm,io
    integer,allocatable,dimension(:)              :: s1,s2,s3,s4
    !
    Nrdm = 4**count(orbital_mask)
    !
    if(size(dm,1)/=Nrdm.OR.size(dm,2)/=Nrdm)then
       stop "ERROR: reduced density matrix and orbital mask have incompatible sizes"
    endif
    !
    suffix = ""
    do ilat = 1,Nlat
       do iorb = 1,Norb
          i = iorb + (ilat-1)*Norb
          if(orbital_mask(i))then
             suffix = trim(suffix)//"_i"//reg(str(io))//"l"//reg(str(jo))
          endif
       enddo
    enddo
    if(present(ineq))then
       fname = "reduced_density_matrix"//trim(suffix)//"_ineq"//reg(str(ineq))//".dat"
    else
       fname = "reduced_density_matrix"//trim(suffix)//".dat"
    endif
    !
    unit = free_unit()
    open(unit,file=fname,action="write",position="rewind",status='unknown')
    !
    do io=1,Nrdm
       write(unit,"(*(F20.16,1X))") (dreal(dm(io,jo)),jo=1,Nrdm)
    enddo
    write(unit,*)
    !
    if(any(dimag(dm)/=0d0))then
       do io=1,Nrdm
          write(unit,"(*(F20.16,1X))") (dimag(dm(io,jo)),jo=1,Nrdm)
       enddo
       write(unit,*)
    endif
    !
    close(unit)
    !
  end subroutine ed_print_dm_orb
  !
  !
  !
  subroutine ed_print_dm_LEGACY(dm,Nrdm,ineq)
    integer                  ,intent(in)            :: Nrdm
    complex(8),dimension(:,:),intent(in)            :: dm
    integer                  ,intent(in),optional   :: ineq
    integer                                         :: unit,Nsites
    character(len=64)                               :: fname
    integer                                         :: io,jo
    !
    if(size(dm,1)/=Nrdm.OR.size(dm,2)/=Nrdm)then
       stop "ERROR: actual dm argument has incogruent size wrt explicitly passed Nrdm"
    endif
    !
    Nsites = nint( 1/Norb * log(real(Nrdm,kind=8)) / log(4d0) ) !Nrdm = 4**(Nsites*Norb)
    !
    if(present(ineq))then
       fname = "reduced_density_matrix_"//reg(str(Nsites))//"sites_ineq"//reg(str(ineq))//".dat"
    else
       fname = "reduced_density_matrix_"//reg(str(Nsites))//"sites.dat"
    endif
    !
    unit = free_unit()
    open(unit,file=fname,action="write",position="rewind",status='unknown')
    !
    do io=1,Nrdm
       write(unit,"(*(F20.16,1X))") (dreal(dm(io,jo)),jo=1,Nrdm)
    enddo
    write(unit,*)
    !
    if(any(dimag(dm)/=0d0))then
       do io=1,Nrdm
          write(unit,"(*(F20.16,1X))") (dimag(dm(io,jo)),jo=1,Nrdm)
       enddo
       write(unit,*)
    endif
    !
    close(unit)
    !
  end subroutine ed_print_dm_LEGACY








END MODULE ED_IO








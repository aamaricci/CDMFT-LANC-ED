# Cluster-DMFT Lanczos Exact Diagonalization

A Lanczos based solver for the **Cluster** Dynamical Mean-Field Theory. The code features solution of the normal (including long range order) problem as well as the superconducting state with different, user enabled, symmetries of the order parameter. 

The `ed_mode=normal` operative mode uses $N_\uparrow$:$N_\downarrow$ quantum numbers and massively parallel execution outlined here [j.cpc.2021.108261](https://doi.org/10.1016/j.cpc.2021.108261). 
The newly developed `ed_mode=superc` operative mode conserves the total spin projection $S_z=N_\uparrow-N_\downarrow$.  


The code is based on:  

* [SciFortran](https://github.com/QcmPlab/SciFortran)  

* [DMFTtools](https://github.com/QcmPlab/DMFTtools) [for the driver development only]

The code structure is as follow:  

* The set of modules compile into a top layer named `CDMFT_ED.f90`  
* The actual implementation of the DMFT equations is case by case performed in a driver program, usually placed in the directory `drivers`. 
* In the driver code the user must includes the `CDMFT_ED`, `SCIFOR` and `DMFT_TOOLS` modules and call the necessary procedures to solve the DMFT equations.

 The only two bath types implemented to date are the 'replica' and 'general' ones, in which the bath consists of multiple, noninteracting copies of the original interacting cluster, with different hybridization schemes: a unique $V$ for each replica in the former case, meaning that the i-th site in the cluster hybridizes with the j-th site in the replica, with an amplitude $V_{ij} = V\delta_{ij}$; a more general $V_{ij} = V_i \delta_{ij}$ for the latter.
 For both cases the bath can be initialized in two ways, through the unified `ed_set_Hbath` interface: one can directly pass the local hamiltonian of the cluster (going with strict semantics for 'replicas') or a symmetry-informed expansion in the form $H^\mathrm{bath} = \sum \lambda^\mathrm{sym}H^\mathrm{sym}$, where the $H^\mathrm{sym}$ terms are arrays isomorphic with $H^\mathrm{loc}$, describing a well defined symmetry of the model, and the $\lambda^\mathrm{sym}$ coefficients are contained in a `([Nineq])x[Nbath]x[Nsym]`-sized array, so to allow each replica to have a different expansion on the same unique basis. No control over the initial guess of the $V_{ij}$ values is given to the user, at the moment.   

 > An example, solving the Hubbard model on the square lattice, is contained in the file `drivers/cdn_hm_2dsquare.f90`.

## DEVELOPMENT

### MILESTONE 1

- [x] Write code for the spin-decoupled case
- [x] Code compiles and runs without errors
- [x] Check exact diagonalization comparing results with known 2x2 plaquette data (without bath)
- [x] Check single-site cluster case against single-site DMFT code
- [x] Check 1d Hubbard chain against literature

### MILESTONE 2

- [x] Include complex routines for diagonalization and bath
- [x] Test 2d BHZ model 

### MILESTONE 3

- [x] Add real-space CDMFT case for finite systems
- [ ] Test Kane-Mele model with real-space CDMFT
- [ ] Test 3d BHZ model for non spin-coupling choices of cluster
- [ ] Rewrite the code for the general spin-coupled case

### MILESTONE 4

- [x] Add cluster density matrix (CDM) computation
- [x] Test and benchmark CDM computation
- [x] Build smaller-rank reduced DMs from CDM
- [x] Check _local_ RDM against semi-analytical build (Norb=1) 
- [x] Replicate the single-site results in [Su, Dai,Tong 2013](https://doi.org/10.1142/S0217984913500346).
- [x] Replicate the (vanishing T) 2x2 results in [C. Walsh et al. 2019](https://doi.org/10.1103/PhysRevLett.122.067203).

### MILESTONE 5
- [x] rewrite the code using a better suited setup to include superconductivity, i.e. arrays of the form `[Nambu,Nambu,Nspin,Nspin,Nlat*Norb,Nlat*Norb]` where `Nambu=1,2` according to `ed_mode=normal,superc` contains the possible Nambu components of the arrays, `Nspin` is the number of independent spin channels, `Nlat` is the number of sites in  the cluster, `Norb` the number of electronic orbitals on each site. 
- [x] extend the code to include superconductive channel
- [x] split and entirely review the fit part
- [x] test `ed_mode=normal` against old code and single site results
- [ ] test `ed_mode=superc` against single site `EDIPACK2` results. 
   

--

***COPYRIGHT & LICENSING***  
The software is provided with no license, as such it is protected by copyright. The software is provided as it is and can be read and copied, in agreement with the Terms of Service of GITHUB. 
Use of the code is constrained to authors agreement.   


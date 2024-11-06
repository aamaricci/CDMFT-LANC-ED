# Code development SUPERC channel


We annotate here the strategy and the implementation of the code  support for Superconductivity. The aim of these notes is to share the workload and keep an updated list of the open/closed issues 

The aim is to update the code introducting **superc**onducting channel using `EDIpack2.0` structure as a template, there including the comments for auto-documentation.


Legend:   
`[ ] `= to be done   
`[o]` = done by to be reviewed  
`[X]` = done and checked    
## Milestone 1: code structure
[o] Update structure of the library by adopting the `EDIpack2.0` template structure (with obvious differences)

## Milestone 2: Variables, Functions and Bath  
[o] Implement changes into global and input variables `ED_INPUT_VARS` `ED_VARS_GLOBAL`: these two modules will be eventually updated again later. Special attention to the bath data structure defined here and the implementation of the `sector` type.


[o] Implement changes into `ED_AUX_FUNX` 

[ ] Import the required changes into the `ED_BATH` directory where all the bath and fit related functions are implemented. We try to follow the structure of `EDIpack` here too. Specially we separate the routines in different modules dealing with: `aux` auxiliary functions which are used in the code, `user_bath` which deals with functions for the user bath, `dmft_bath` which contains functions dealing with the internal data structure containing the bath, `bath_functions` which contains the local functions of the bath, i.e. non-interacting Anderson, its inverse and the hybridization function (these will be used later in evaluation of the GF), `fit` dealing with fit of the bath. 
**The development of this part is very import and must be done with great attention** 

## Milestone 3: symmetry sector

[ ] split the old `ED_SETUP` into `ED_SETUP` and `ED_SECTOR`. The first implement functions to setup the ED calculation, the second implements all the construction of symmetry sector through its data structure (see point above), all the functions dealing with identification of the sector information (index, quantum numbers, dimensions, etc.), the fermionic operators, etc...

## Milestone 4: Hamiltonian & Diagonalization

[ ] Using the normal case and comapring with `EDIpack2` we implement the `ED_DIAG_SUPERC` module, which performs the diagonalization sector by sector (here we assume to have all the functions that will be developed in the next point)

[ ] Implement all the Hamiltonian construction. This is a critical step. We use current structure and compare with `EDIpack2.0` concerning the main structure of the functions, the MPI split, etc.

[ ] Finally we implement the single files in direct and sparse directories, that is where the actual Hamiltonian is constructed. Comparison with `EDIpack2.0` again is essential here. 
**In a fist version we keep the s-wave structure but this part will require great attention to implement the correct symmetry of the order parameter**. 


## Milestone 5: Green's functions & Observables

[ ] Here we implement that construction of the Green's function. We adapt the algorithm from the `EDIpack2.0` code to the cluster case. **Again special attention to the symmetries should be given here**

[ ] We implement the evaluation of the impurity observables and static correlations using the superconducting symmetry sectors.


## Milestone 6: Wrap up the code, code MAIN

[ ] Implement the global modules which call the correct function based on the value of the controlling variable `ed_mode` (which must be included in `INPUT_VARS` now). 

[ ] Update the `ED_MAIN` module which wraps the ED calculations


## Milestone 7: Density Matrix

[ ] Implement the calculation of the Density Matrix. This step will require to adapt the `SPARSE_MAP` data structure to hold the information about the separation between imp and bath states used in the `DENSITY_MATRIX` calculations. 

## Milestone 8: compile and test

[ ] Update CMake compilation strategy and start testing the code. The testing will have a second checklist.  
  
  subroutine ed_get_single_particle_density_matrix_single(dm,doprint)
    complex(8),allocatable,intent(out)            :: dm(:,:)
    logical               ,intent(in) ,optional   :: doprint
    logical                                       :: doprint_
    !
    doprint_=.false.; if(present(doprint)) doprint_=doprint
    !
    if(.not.allocated(single_particle_density_matrix))then
       stop "ERROR: single_particle_density_matrix is not allocated"
    endif
    !
    !Impurity problem basis
    dm = nnn2nso_reshape(single_particle_density_matrix)
    !
    !Print to file (if requested)
    if(doprint_)then
       call ed_print_dm(dm,Nambu*Nspin*Nimp)
    endif
    !
  end subroutine ed_get_single_particle_density_matrix_single







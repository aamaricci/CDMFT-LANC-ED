subroutine ed_get_epot_(eimp)
  real(8) :: eimp(Nlat)
  eimp = ed_Epot
end subroutine ed_get_epot_

subroutine ed_get_eint_(eimp)
  real(8) :: eimp(Nlat)
  eimp = ed_Eint
end subroutine ed_get_eint_

subroutine ed_get_ehartree_(eimp)
  real(8) :: eimp(Nlat)
  eimp = ed_Ehartree
end subroutine ed_get_ehartree_

subroutine ed_get_eknot_(eimp)
  real(8) :: eimp(Nlat)
  eimp = ed_Eknot
end subroutine ed_get_eknot_




 subroutine ed_get_eimp_lattice(yii,Nineq)
   integer                      :: Nineq
   real(8),dimension(Nineq,4)    :: yii
   yii=0d0
   if(allocated(eii))then
      if(Nineq>size(eii,1)) stop "ed_get_eimp error: required N_sites > evaluated N_sites"
      yii=eii
   endif
 end subroutine ed_get_eimp_lattice

 subroutine ed_get_epot_lattice(yii,Nineq)
   integer                 :: Nineq
   real(8),dimension(Nineq) :: yii
   yii=0d0
   if(allocated(eii))then
      if(Nineq>size(eii,1)) stop "ed_get_epot error: required N_sites > evaluated N_sites"
      yii=eii(:,1)
   endif
 end subroutine ed_get_epot_lattice

 subroutine ed_get_eint_lattice(yii,Nineq)
   integer                 :: Nineq
   real(8),dimension(Nineq) :: yii
   yii=0d0
   if(allocated(eii))then
      if(Nineq>size(eii,1)) stop "ed_get_eint error: required N_sites > evaluated N_sites"
      yii=eii(:,2)
   endif
 end subroutine ed_get_eint_lattice

 subroutine ed_get_ehartree_lattice(yii,Nineq)
   integer                 :: Nineq
   real(8),dimension(Nineq) :: yii
   yii=0d0
   if(allocated(eii))then
      if(Nineq>size(eii,1)) stop "ed_get_ehartree error: required N_sites > evaluated N_sites"
      yii=eii(:,3)
   endif
 end subroutine ed_get_ehartree_lattice

 subroutine ed_get_eknot_lattice(yii,Nineq)
   integer                 :: Nineq
   real(8),dimension(Nineq) :: yii
   yii=0d0
   if(allocated(eii))then
      if(Nineq>size(eii,1)) stop "ed_get_knot error: required N_sites > evaluated N_sites"
      yii=eii(:,4)
   endif
 end subroutine ed_get_eknot_lattice

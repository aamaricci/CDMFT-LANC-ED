subroutine ed_get_phi_d0(self,iorb,ilat)
  real(8)          :: self
  integer          :: iorb
  integer,optional :: ilat
  integer          :: io
  if(.not.present(ilat))then
     if(iorb>Nimp)stop "ed_get_docc error: index > Nimp"
     io   = iorb
     self = ed_phisc(io)
  else
     if(iorb>Norb)stop "ed_get_phi error: orbital index > Norb"
     if(ilat>Nlat)stop "ed_get_phi error: lattice index > Nlat"
     io   = iorb + (ilat-1)*Norb
     self = ed_phisc(io)
  endif
end subroutine ed_get_phi_d0

subroutine ed_get_phi_d1(self)
  real(8)          :: self(Nimp)
  self = ed_phisc
end subroutine ed_get_phi_d1

subroutine ed_get_phi_d2(self)
  real(8)          :: self(Nlat,Norb)
  integer          :: ilat,iorb,io
  do concurrent(ilat=1:Nlat,iorb=1:Norb)
     io = iorb + (ilat-1)*Norb
     self(ilat,iorb) = ed_phisc(io)
  end do
end subroutine ed_get_phi_d2


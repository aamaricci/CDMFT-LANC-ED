subroutine ed_get_doubles_n1(self)
  real(8) :: self(4)
  self = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
end subroutine ed_get_doubles_n1


subroutine ed_get_dust_n0(self)
  real(8) :: self
  self = ed_Dust
end subroutine ed_get_dust_n0

subroutine ed_get_dund_n0(self)
  real(8) :: self
  self = ed_Dund
end subroutine ed_get_dund_n0

subroutine ed_get_dse_n0(self)
  real(8) :: self
  self = ed_Dse
end subroutine ed_get_dse_n0

subroutine ed_get_dph_n0(self)
  real(8) :: self
  self = ed_Dph
end subroutine ed_get_dph_n0


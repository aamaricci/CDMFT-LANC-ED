subroutine ed_get_sigma_site_n2(self,axis,type,z)
  complex(8),dimension(:,:,:),intent(inout)       :: self ! Self-energy matrix
  character(len=*),optional                       :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                       :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional                :: z
  character(len=1)                                :: axis_
  character(len=1)                                :: type_
  complex(8),dimension(:),allocatable             :: z_
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_sigma ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  !
  call assert_shape(self,[Nspin*Nimp,Nspin*Nimp,L],'ed_get_sigma','self')
  !
  if(allocated(F))deallocate(F)
  allocate(F(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,L))
  F = get_Sigma(z_)
  !
  select case(type_)
  case default; stop "ed_get_sigma ERROR: type is neither Normal, nor Anomalous"
  case ('n','N')
     self = nn2so_reshape(F(1,1,:,:,:,:,:),Nspin,Nimp,Lmats)
  case ('a','A')
     self = nn2so_reshape(F(1,2,:,:,:,:,:),Nspin,Nimp,Lmats)
  end select
  !
  deallocate(F)
  call deallocate_grids
  !
end subroutine ed_get_sigma_site_n2


subroutine ed_get_sigma_site_n4(self,axis,type,z)
  complex(8),dimension(:,:,:,:,:),intent(inout) :: self
  character(len=*),optional                     :: axis
  character(len=*),optional                     :: type
  complex(8),dimension(:),optional              :: z
  character(len=1)                              :: axis_
  character(len=1)                              :: type_
  complex(8),dimension(:),allocatable           :: z_
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_sigma ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  call assert_shape(self,[Nspin,Nspin,Nimp,Nimp,L],'ed_get_sigma','self')
  !
  if(allocated(F))deallocate(F)
  allocate(F(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,L))  
  F = get_Sigma(z_)
  !
  select case(type_)
  case default; stop "ed_get_sigma ERROR: type is neither Normal, nor Anomalous"
  case ('n','N')
     self = F(1,1,:,:,:,:,:)
  case ('a','A')
     self = F(1,2,:,:,:,:,:)
  end select
  !
  deallocate(F)
  call deallocate_grids
  !
end subroutine ed_get_sigma_site_n4


subroutine ed_get_sigma_site_n6(self,axis,z)
  complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: self
  character(len=*),optional                         :: axis
  complex(8),dimension(:),optional                  :: z
  character(len=1)                                  :: axis_
  complex(8),dimension(:),allocatable               :: z_
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_sigma ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  call assert_shape(self,[Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,L],'ed_get_sigma','self')
  !
  Self = get_Sigma(z_)
  !
  call deallocate_grids
  !
end subroutine ed_get_sigma_site_n6







subroutine ed_get_g0and_n2(x,bath_,self,axis,type)
  complex(8),dimension(:),intent(in)                              :: x !complex array of frequencies
  real(8),dimension(:)                                            :: bath_ !user-accessible bath array
  complex(8),dimension(:,:,:)                                     :: self !non-interacting Green's function
  character(len=*),optional                                       :: axis !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis 
  character(len=*),optional                                       :: type !string indicating the desired function, :code:`'n'` for normal (default), :code:`'a'` for anomalous
  !
  type(effective_bath)                                            :: dmft_bath_
  logical                                                         :: check
  character(len=1)                                                :: axis_,type_
  complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(x)) :: g0
  !
  axis_='m';if(present(axis))axis_=axis
  type_='n';if(present(type))type_=trim(type)
  check= check_bath_dimension(bath_)
  if(.not.check)stop "g0and_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  !
  g0 = g0and_bath_function(x,dmft_bath_)
  !
  call assert_shape(self,[Nspin*Nimp,Nspin*Nimp,size(x)],'ed_get_g0and','g0and')  
  select case(type_)
  case default;stop "ed_get_g0and ERROR: type is wrong: either Normal or Anomalous"
  case ('n','N');
     self = nn2so_reshape(g0(1,1,:,:,:,:,:),Nspin,Nimp,size(x))
  case ('a','A')
     self = nn2so_reshape(g0(1,2,:,:,:,:,:),Nspin,Nimp,size(x))
  end select
  !
  call deallocate_dmft_bath(dmft_bath_)
  !
end subroutine ed_get_g0and_n2


subroutine ed_get_g0and_n4(x,bath_,Self,axis,type)
  complex(8),dimension(:),intent(in)                              :: x
  real(8),dimension(:)                                            :: bath_
  complex(8),dimension(:,:,:,:,:)                                 :: Self
  character(len=*),optional                                       :: axis
  character(len=*),optional                                       :: type
  !
  type(effective_bath)                                            :: dmft_bath_
  logical                                                         :: check
  character(len=1)                                                :: axis_,type_
  complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(x)) :: g0
  !
  axis_='m';if(present(axis))axis_=axis
  type_='n';if(present(type))type_=trim(type)
  check= check_bath_dimension(bath_)
  if(.not.check)stop "g0and_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  !
  g0 = g0and_bath_function(x,dmft_bath_)
  !
  call assert_shape(self,[Nspin,Nspin,Nimp,Nimp,size(x)],'ed_get_g0and','g0and')  
  select case(type_)
  case default;stop "ed_get_g0and ERROR: type is wrong: either Normal or Anomalous"
  case ('n','N');
     Self = g0(1,1,:,:,:,:,:)
  case ('a','A')
     Self = g0(1,2,:,:,:,:,:)
  end select
  !
  call deallocate_dmft_bath(dmft_bath_)
end subroutine ed_get_g0and_n4


subroutine ed_get_g0and_n6(x,bath_,Self,axis)
  complex(8),dimension(:),intent(in)  :: x
  real(8),dimension(:)                :: bath_
  complex(8),dimension(:,:,:,:,:,:,:) :: Self
  character(len=*),optional           :: axis
  !
  type(effective_bath)                :: dmft_bath_
  logical                             :: check
  character(len=1)                    :: axis_
  !
  axis_='m';if(present(axis))axis_=axis
  check= check_bath_dimension(bath_)
  if(.not.check)stop "g0and_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  !
  call assert_shape(self,[Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(x)],'ed_get_g0and','g0and')  
  Self = g0and_bath_function(x,dmft_bath_)
  !
  call deallocate_dmft_bath(dmft_bath_)
end subroutine ed_get_g0and_n6







subroutine ed_get_delta_n2(x,bath_,self,axis,type)
  complex(8),dimension(:),intent(in)                              :: x !complex array of frequencies
  real(8),dimension(:)                                            :: bath_ !user-accessible bath array
  complex(8),dimension(:,:,:)                                     :: self !non-interacting Green's function
  character(len=*),optional                                       :: axis !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis 
  character(len=*),optional                                       :: type !string indicating the desired function, :code:`'n'` for normal (default), :code:`'a'` for anomalous
  !
  type(effective_bath)                                            :: dmft_bath_
  logical                                                         :: check
  character(len=1)                                                :: axis_,type_
  complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(x)) :: g0
  !
  axis_='m';if(present(axis))axis_=axis
  type_='n';if(present(type))type_=trim(type)
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  !
  g0 = delta_bath_function(x,dmft_bath_)
  !
  call assert_shape(self,[Nspin*Nimp,Nspin*Nimp,size(x)],'ed_get_delta','delta')  
  select case(type_)
  case default;stop "ed_get_delta ERROR: type is wrong: either Normal or Anomalous"
  case ('n','N');
     self = nn2so_reshape(g0(1,1,:,:,:,:,:),Nspin,Nimp,size(x))
  case ('a','A')
     self = nn2so_reshape(g0(1,2,:,:,:,:,:),Nspin,Nimp,size(x))
  end select
  !
  call deallocate_dmft_bath(dmft_bath_)
  !
end subroutine ed_get_delta_n2


subroutine ed_get_delta_n4(x,bath_,Self,axis,type)
  complex(8),dimension(:),intent(in)                              :: x
  real(8),dimension(:)                                            :: bath_
  complex(8),dimension(:,:,:,:,:)                                 :: Self
  character(len=*),optional                                       :: axis
  character(len=*),optional                                       :: type
  !
  type(effective_bath)                                            :: dmft_bath_
  logical                                                         :: check
  character(len=1)                                                :: axis_,type_
  complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(x)) :: g0
  !
  axis_='m';if(present(axis))axis_=axis
  type_='n';if(present(type))type_=trim(type)
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  !
  g0 = delta_bath_function(x,dmft_bath_)
  !
  call assert_shape(self,[Nspin,Nspin,Nimp,Nimp,size(x)],'ed_get_delta','delta')  
  select case(type_)
  case default;stop "ed_get_delta ERROR: type is wrong: either Normal or Anomalous"
  case ('n','N');
     Self = g0(1,1,:,:,:,:,:)
  case ('a','A')
     Self = g0(1,2,:,:,:,:,:)
  end select
  !
  call deallocate_dmft_bath(dmft_bath_)
end subroutine ed_get_delta_n4


subroutine ed_get_delta_n6(x,bath_,Self,axis)
  complex(8),dimension(:),intent(in)  :: x
  real(8),dimension(:)                :: bath_
  complex(8),dimension(:,:,:,:,:,:,:) :: Self
  character(len=*),optional           :: axis
  !
  type(effective_bath)                :: dmft_bath_
  logical                             :: check
  character(len=1)                    :: axis_
  !
  axis_='m';if(present(axis))axis_=axis
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  !
  call assert_shape(self,[Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(x)],'ed_get_delta','delta')  
  Self = delta_bath_function(x,dmft_bath_)
  !
  call deallocate_dmft_bath(dmft_bath_)
end subroutine ed_get_delta_n6





subroutine ed_get_invG0_n2(x,bath_,self,axis,type)
  complex(8),dimension(:),intent(in)                              :: x !complex array of frequencies
  real(8),dimension(:)                                            :: bath_ !user-accessible bath array
  complex(8),dimension(:,:,:)                                     :: self !non-interacting Green's function
  character(len=*),optional                                       :: axis !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis 
  character(len=*),optional                                       :: type !string indicating the desired function, :code:`'n'` for normal (default), :code:`'a'` for anomalous
  !
  type(effective_bath)                                            :: dmft_bath_
  logical                                                         :: check
  character(len=1)                                                :: axis_,type_
  complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(x)) :: g0
  !
  axis_='m';if(present(axis))axis_=axis
  type_='n';if(present(type))type_=trim(type)
  check= check_bath_dimension(bath_)
  if(.not.check)stop "invG0_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  !
  g0 = invG0_bath_function(x,dmft_bath_)
  !
  call assert_shape(self,[Nspin*Nimp,Nspin*Nimp,size(x)],'ed_get_invG0','invG0')  
  select case(type_)
  case default;stop "ed_get_invG0 ERROR: type is wrong: either Normal or Anomalous"
  case ('n','N');
     self = nn2so_reshape(g0(1,1,:,:,:,:,:),Nspin,Nimp,size(x))
  case ('a','A')
     self = nn2so_reshape(g0(1,2,:,:,:,:,:),Nspin,Nimp,size(x))
  end select
  !
  call deallocate_dmft_bath(dmft_bath_)
  !
end subroutine ed_get_invG0_n2


subroutine ed_get_invG0_n4(x,bath_,Self,axis,type)
  complex(8),dimension(:),intent(in)                              :: x
  real(8),dimension(:)                                            :: bath_
  complex(8),dimension(:,:,:,:,:)                                 :: Self
  character(len=*),optional                                       :: axis
  character(len=*),optional                                       :: type
  !
  type(effective_bath)                                            :: dmft_bath_
  logical                                                         :: check
  character(len=1)                                                :: axis_,type_
  complex(8),dimension(Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(x)) :: g0
  !
  axis_='m';if(present(axis))axis_=axis
  type_='n';if(present(type))type_=trim(type)
  check= check_bath_dimension(bath_)
  if(.not.check)stop "invG0_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  !
  g0 = invG0_bath_function(x,dmft_bath_)
  !
  call assert_shape(self,[Nspin,Nspin,Nimp,Nimp,size(x)],'ed_get_invG0','invG0')  
  select case(type_)
  case default;stop "ed_get_invG0 ERROR: type is wrong: either Normal or Anomalous"
  case ('n','N');
     Self = g0(1,1,:,:,:,:,:)
  case ('a','A')
     Self = g0(1,2,:,:,:,:,:)
  end select
  !
  call deallocate_dmft_bath(dmft_bath_)
end subroutine ed_get_invG0_n4


subroutine ed_get_invG0_n6(x,bath_,Self,axis)
  complex(8),dimension(:),intent(in)  :: x
  real(8),dimension(:)                :: bath_
  complex(8),dimension(:,:,:,:,:,:,:) :: Self
  character(len=*),optional           :: axis
  !
  type(effective_bath)                :: dmft_bath_
  logical                             :: check
  character(len=1)                    :: axis_
  !
  axis_='m';if(present(axis))axis_=axis
  check= check_bath_dimension(bath_)
  if(.not.check)stop "invG0_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  !
  call assert_shape(self,[Nambu,Nambu,Nspin,Nspin,Nimp,Nimp,size(x)],'ed_get_invG0','invG0')  
  Self = invG0_bath_function(x,dmft_bath_)
  !
  call deallocate_dmft_bath(dmft_bath_)
end subroutine ed_get_invG0_n6










program mesh_tester
	use fem_module
  use frp_generate_module
	type(struct_model) :: model
	character name_mdl*16
	real(8),allocatable :: nd_data(:,:),el_data(:,:)
	character(16) :: nd_data_name(1),el_data_name(1)
	character(32) output_file_name

  integer shift(4)

	print *, "Reading models..."; print *

  shift(1) = 2 
  shift(2) = 25
  shift(3) = 13
  shift(4) = 42
!	call read_model(model,name_mdl)
  call generate_FRP_Model(model,shift)
  call set_vf_params(model,2,4,shift)

  model%name = "gfrp_test"

	allocate(nd_data(1,model%n_nds),el_data(1,model%n_els))
	nd_data = 1d0
	el_data = 1d0
	nd_data_name(1) = "nd_data"
	el_data_name(1) = "el_data"
	output_file_name = "meshtest"

	call make_inp(nd_data,nd_data_name,el_data,el_data_name,model,output_file_name)

 	print *,"- Press ENTER to finish -"; print *
	read *

end program

program ebescg
	use fem_module
	use ebe_scg_module
	implicit none
	type(struct_model) :: model
	type(struct_bc) :: bc
	type(struct_output) :: output
	character name_mdl*16
	real(8),allocatable :: f(:),u(:)


	write(*,'(a)', advance='no') "Model: "
	read *, name_mdl; print *

	print *, "Reading models..."; print *

	call read_model(model,name_mdl)
	call set_state2d(model,0,1d0)

	allocate(u(model%n_nds*model%dim))

 	print *,"Reading boundary conditions..."; print *
	call read_cc(bc,model)

	call read_f(f,model)

	call ebe_scg(u,model,bc,f)


!!!!!!!!!!!!!!!post process
	call set_constrained_u(u,model,bc)


 	call calc_output(output,model,u)

 	call visualize_u(model,u)

 	call output_inp(model,output,"ebe_scg")


 	print *,"Calculation was completed successfully!"; print *

 	print *,"- Press ENTER to finish -"; print *

 	read *

contains

end program

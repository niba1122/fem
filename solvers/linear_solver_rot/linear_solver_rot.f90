program fem
	use rotateD_module
	implicit none
	integer i
	character name_mdl*16
	integer state2d,mat_rot
	double precision,pointer :: D_rot(:,:,:),angle(:)
	double precision,allocatable,dimension(:) :: k,u,f
	integer,allocatable,dimension(:) :: index_k_diag
	double precision thickness

	type(struct_model) :: model
	type(struct_bc) :: bc
	type(struct_output) :: output
	type(struct_data) :: data(2)

	integer iostat
	character errmsg*128

	character(32) output_file_name


	write(*,'(a)', advance='no') "Model: "
	read *, name_mdl; print *

	print *, "Reading models..."; print *

	call read_model(model,name_mdl)

	if (model%dim == 2) then
		! 厚さの読み込み
		write(*,'(a)',advance='no') "thickness: "
		read (*,*,iostat=iostat) thickness; print *
		if (iostat /= 0) call error("Invalid value")
		! 平面応力/平面ひずみ
		write(*,'(a)',advance='no') "state(0: plane stress,  >0: plane strain): "
		read (*,*,iostat=iostat) state2d; print *
		if (iostat /= 0) call error("Invalid value")
		call set_state2d(model,state2d,thickness)
	end if

	write(*,'(a)',advance='no') "which material is rotated: "
	read (*,*,iostat=iostat) mat_rot; print *
	if (iostat /= 0) call error("Invalid value")

! 	call rd_calc_angle(D_rot,model,mat_rot)
	call rd_calc_angle(angle,model,mat_rot)
	nullify(data(1)%d)
	data(2)%d(1:model%n_els,1:1,1:1) => angle

 	print *,"Reading boundary conditions..."; print *
	call read_cc(bc,model)

	call read_f(f,model)

 	print *,"Calculating K matrix..."; print *

	call init_addition_matrix_sln(k,index_k_diag,model,bc)
	call calc_addition_matrix_sln(k,model,index_k_diag,rd_calc_BDB,data)

	print *,"Solving stiffness equation..."; print *
	call set_cc(k,f,model,bc,index_k_diag)

	allocate(u(model%n_nds*2))

 	u = sl_LDL(K,f,index_k_diag)

 	call set_constrained_u(u,model,bc)

 	write(*,'(a)', advance='no') "Name of inp file: "
	read (*,*,iostat=iostat) output_file_name; print *


 	print *,"Outputting results..."; print *

 	call rd_calc_output(output,model,angle,u)
 	
 	call visualize_u(model,u)

 	call output_inp(model,output,output_file_name)

 	print *,"Calculation was completed successfully!"; print *

 	print *,"- Press ENTER to finish -"; print *

	read *

contains

end program fem
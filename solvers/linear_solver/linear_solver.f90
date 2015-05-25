program fem
	use fem_module
	implicit none
! 	external func1
	integer i
	character name_mdl*16
	integer state2d,vu,af
	double precision,allocatable,dimension(:) :: k,u,f
	integer,allocatable,dimension(:) :: index_k_diag
	double precision thickness

	type(struct_model) :: model
	type(struct_bc) :: bc
	type(struct_output) :: output

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

 	print *,"Reading boundary conditions..."; print *
	call read_cc(bc,model)




 	write(*,'(a)', advance='no') "auto equivalent nodal force? (>0: Yes, <=0: No)"
	read (*,*,iostat=iostat) af; print *
 	if (af > 0) then
		call calc_equivalent_nodal_force(f,model,(/0d0,0d0,200000d0/),0d0,0d0,1d0,-1d0)
	else
		call read_f(f,model)
	end if


 	print *,"Calculating K matrix..."; print *

	call init_addition_matrix_sln(k,index_k_diag,model,bc)

	call calc_addition_matrix_sln(k,model,index_k_diag,calc_BDB)

	print *,"Solving stiffness equation..."; print *
	call set_cc(k,f,model,bc,index_k_diag)

	allocate(u(model%n_nds*2))

 	u = sl_LDL(K,f,index_k_diag)

 	call set_constrained_u(u,model,bc)

 	write(*,'(a)', advance='no') "Name of inp file: "
	read (*,*,iostat=iostat) output_file_name; print *

 	print *,"Outputting results..."; print *

 	call calc_output(output,model,u)


 	write(*,'(a)', advance='no') "visualize u? (>0: Yes, <=0: No)"
	read (*,*,iostat=iostat) vu; print *
 	if (vu > 0) call visualize_u(model,u)

 	call output_inp(model,output,output_file_name)


 	print *,"Calculation was completed successfully!"; print *

 	print *,"- Press ENTER to finish -"; print *


!  	call model_add_d_data(model,reshape((/1d0,2d0,3d0,4d0,5d0,6d0/),(/2,3,1/)),1,"arr1")
!  	call model_add_d_data(model,reshape((/4d0,5d0,6d0,7d0,8d0,9d0/),(/1,6,1/)),1,"arr2")
!  	call model_add_d_data(model,reshape((/1d0,2d0,3d0/),(/1,3,1/)),1,"arr3")

!  	print *,model%d_data(2)%d
!  	print *,model%d_data(3)%d
!  	print *,model%d_data(4)%d

!  	call model_get_d_data(test,model,"arr2")
!  	print *, test


	read *

contains

end program fem

program fem
	use rotateD_module
	use FRPGenerate_module
	implicit none
	integer i,step
	character name_mdl*16
	integer state2d,mat_rot
	double precision,pointer :: D_rot(:,:,:),angle(:)
	double precision,allocatable,dimension(:) :: k,u,f,f_
	integer,allocatable,dimension(:) :: index_k_diag
	double precision thickness

	type(struct_model) :: model
	type(struct_bc) :: bc
	type(struct_output) :: output
	type(struct_data) :: rd_data(2)

	integer iostat
	character(32) char_step


	allocate(u(model%n_nds*2))
	allocate(f(model%n_nds*2))


!-------------------------------------------------------------------------------------------------------
!	Read model
!-------------------------------------------------------------------------------------------------------
	call generate_FRP_Model(model,bc)


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

	call rd_calc_angle(angle,model,mat_rot)


!-------------------------------------------------------------------------------------------------------
!	Read boundary conditions
!-------------------------------------------------------------------------------------------------------

	call read_cc(bc,model)

	call read_f(f_,model)


!-------------------------------------------------------------------------------------------------------
!	Initialize stiffness matrix
!-------------------------------------------------------------------------------------------------------
	call init_addition_matrix_sln(k,index_k_diag,model,bc)

!-------------------------------------------------------------------------------------------------------
!	Monte Carlo simulation
!-------------------------------------------------------------------------------------------------------

call init_genrand(1)

step = 5


nullify(rd_data(1)%d)
rd_data(2)%d(1:model%n_els,1:1,1:1) => angle

do i=1,step
	u = 0d0
	f = f_

	call calc_addition_matrix_sln(k,model,index_k_diag,rd_calc_BDB,rd_data)

	print *,"Solving stiffness equation..."; print *
	call set_cc(k,f,model,bc,index_k_diag)

 	u = sl_LDL(K,f,index_k_diag)

 	call set_constrained_u(u,model,bc)

 	print *,"Outputting results..."; print *

 	call rd_calc_output(output,model,angle,u)
 	
 	call visualize_u(model,u)


	write(char_step,'(i0)') i

 	call output_inp(model,output,trim(char_step))

 	print *,"Calculation was completed successfully!"; print *
 end do

 	print *,"- Press ENTER to finish -"; print *

	read *


contains

function gen_normal_rand()
	implicit none
	double precision gen_normal_rand,x,y,z,pi,genrand_res53

	pi = dacos(-1d0)

	do
		x = genrand_res53()
		y = genrand_res53()
		z = dsqrt(-2d0*dlog(x)) * dCOS(2*pi*y)
		if (abs(z)<=6d0) exit
	end do
!	if (z>6d0) then
!		z = 6d0
!	else if (z<-6d0) then
!		z = -6d0
!	end if

	gen_normal_rand = z
end function

end program fem
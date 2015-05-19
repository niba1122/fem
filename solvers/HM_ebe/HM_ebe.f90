program fem
	use hm_module
	use ebe_scg_module
	implicit none
	integer i,j
	character name_mdl*16
	integer state2d,mat_rot,dim
	double precision,allocatable,dimension(:) :: k
	double precision,allocatable :: BE(:,:),EB(:,:),E(:,:),B(:,:),i_mat(:,:),EBx(:,:),Ke(:,:)
	double precision,allocatable :: x(:,:)
	integer,allocatable,dimension(:) :: index_k_diag
	double precision volume
	integer auto_pbc
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

	dim = model%dim

	if (dim == 2) then
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

	write(*,'(a)',advance='no') "Use auto periodic boundary condition? (<=0: no,  >0: use): "
	read (*,*,iostat=iostat) auto_pbc; print *

	if (auto_pbc <= 0) then 
		call read_cc(bc,model)

		write(*,'(a)',advance='no') "Volume of unit cell: "	
		read (*,*,iostat=iostat) volume; print *
		if (iostat /= 0) call error("Invalid value")
	else 
		call hm_auto_periodic_bc(bc,volume,model)
		print *,volume
	end if

 	call hm_calc_BE_integral(BE,model)

	print *,"Solving stiffness equation..."; print *

   	allocate(x(model%n_nds*model%dim,model%dim*(model%dim+1)/2))
   	allocate(Ke(model%n_nds_1el*model%dim,model%n_nds_1el*model%dim))


	do i=1,dim*(dim+1)/2

		print *,"x",i

		call calc_element_integral(Ke,model,1,calc_BDB)

		print *,dsqrt(dot_product(BE(:,i),BE(:,i)))
		if (dsqrt(dot_product(BE(:,i),BE(:,i))) > 1d-10) then
			call ebe_scg(x(:,i),model,bc,BE(:,i))
		else
			print *,"EBE-SCG calculation wasn't carried out because |BE(i)| ~ 0"
			x(:,i) = 0d0
		end if
 		call set_constrained_u(x(:,i),model,bc)

 	end do

 	call hm_calc_EB_integral(EB,model)


 	call hm_calc_E_integral(E,model)
 
print *,"E:"
print *,E

  	print *,"EBx1:"
  	print *,matmul(EB,x)

print *,"DH"
print *,E-matmul(EB,x)

 	print *,"Calculation was completed successfully!"; print *

 	print *,"- Press ENTER to finish -"; print *

 	read *

contains

end program fem
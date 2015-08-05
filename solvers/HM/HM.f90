program fem
	use hm_module
	implicit none
	integer i,j
	character name_mdl*16
	integer state2d,mat_rot,dim
	double precision,allocatable,dimension(:) :: k,u
	double precision,allocatable :: BE(:,:),EB(:,:),E(:,:),B(:,:),i_mat(:,:),EBx(:,:)
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



 	print *,"Calculating K matrix..."; print *

	call init_addition_matrix_sln(k,index_k_diag,model,bc)


	print *,"Solving stiffness equation..."; print *

  allocate(x(model%n_nds*model%dim,model%dim*(model%dim+1)/2))

	do i=1,dim*(dim+1)/2
		k = 0d0
		call calc_addition_matrix_sln(k,model,index_k_diag,calc_BDB)
		call set_cc(k,BE(:,i),model,bc,index_k_diag)


do j=1,model%n_nds
if (i == 6) print *,BE((j-1)*dim+1:j*dim,6)
end do
 		x(:,i) = sl_LDL(K,BE(:,i),index_k_diag)
 		call set_constrained_u(x(:,i),model,bc)

! 		deallocate(k)
! 		deallocate(index_k_diag)
 	end do

! print *,"BE"
! do i=1,ubound(BE,1)
! 	print *, BE(i,:)
! end do


!  	print *,"x"
! do i=1,ubound(x,1)
!  	print *,x(i,:)
! end do

 	call hm_calc_EB_integral(EB,model)
! print *,"EB:"
! do i=1,ubound(EB,2)
! print *,EB(:,i)
! end do

 	call hm_calc_E_integral(E,model)
 
print *,"E:"
print *,E

  	print *,"EBx1:"
  	print *,matmul(EB,x)

!   	call hm_calc_EBx_integral(EBx,model)
! print *,"EB2x:"
! print *,EBx

print *,"DH"
print *,E-matmul(EB,x)

! print *,"DH:"

! print *,(E-matmul(EB,x))/0.035d0/0.055d0



! !  	call set_constrained_u(u,model,bc)


!  	write(*,'(a)', advance='no') "Name of inp file: "
! 	read (*,*,iostat=iostat) output_file_name; print *


!  	print *,"Outputting results..."; print *

!  	call calc_output(output,model,u)
!  	call output_inp(model,output,output_file_name)

 	print *,"Calculation was completed successfully!"; print *

 	print *,"- Press ENTER to finish -"; print *

 	read *

contains

end program fem

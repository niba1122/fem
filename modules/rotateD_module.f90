module rotateD_module
	use fem_module
	implicit none

contains

! subroutine rd_calc_angle(D_rot,model,n_mat)
! 	type(struct_model) :: model
! 	integer i,n_mat,n_els
! 	double precision,pointer,dimension(:,:) :: nodes
! 	integer,pointer,dimension(:,:) :: elements
! 	double precision,pointer,dimension(:,:,:) :: materials
! 	double precision,pointer :: D_rot(:,:,:)
! 	integer,pointer,dimension(:) :: material_nos
! 	double precision angle

! 	nodes => model%nodes
! 	elements => model%elements
! 	materials => model%materials
! 	material_nos => model%material_nos
! 	n_els = model%n_els

! 	allocate(D_rot(6,6,n_els))

! 	do i=1,n_els
! 		if (material_nos(i) == n_mat) then
!  			angle = atan2(nodes(2,elements(2,i))+nodes(2,elements(3,i))-nodes(2,elements(1,i))-nodes(2,elements(4,i)) , &
!  							&nodes(1,elements(2,i))+nodes(1,elements(3,i))-nodes(1,elements(1,i))-nodes(1,elements(4,i)) );
! 			D_rot(:,:,i) = rot_D(materials(:,:,material_nos(i)),angle)

! 		else 
! 			D_rot(:,:,i) = materials(:,:,material_nos(i))
! 		end if
! 	end do

! end subroutine

subroutine rd_calc_angle(angle,model,n_mat)
	type(struct_model) :: model
	integer i,n_mat,n_els
	double precision,pointer,dimension(:,:) :: nodes
	integer,pointer,dimension(:,:) :: elements
	double precision,pointer,dimension(:,:,:) :: materials
	double precision,pointer :: D(:,:,:)
	integer,pointer,dimension(:) :: material_nos
	double precision,pointer :: angle(:)

	nodes => model%nodes
	elements => model%elements
	materials => model%materials
	material_nos => model%material_nos
	n_els = model%n_els

	allocate(D(6,6,n_els))
	allocate(angle(n_els))

	do i=1,n_els
		if (material_nos(i) == n_mat) then
 			angle(i) = atan2(nodes(2,elements(2,i))+nodes(2,elements(3,i))-nodes(2,elements(1,i))-nodes(2,elements(4,i)) , &
 							&nodes(1,elements(2,i))+nodes(1,elements(3,i))-nodes(1,elements(1,i))-nodes(1,elements(4,i)) );
! 			D(:,:,i) = materials(:,:,material_nos(i))

		else 
! 			D(:,:,i) = materials(:,:,material_nos(i))
		end if
	end do

end subroutine

function rd_calc_BDB(model,i_el,coord,data)
	type(struct_model) :: model
	integer :: i_el
	double precision :: coord(:)
	type(struct_data),optional :: data(:)
	double precision,allocatable,dimension(:,:) :: rd_calc_BDB
	double precision,allocatable :: B(:,:),D(:,:)

	allocate(rd_calc_BDB(model%dim*model%n_nds_1el,model%dim*model%n_nds_1el))
	rd_calc_BDB = 0d0

	if (associated(data(1)%d)) then
		print *,"ok"
		call calc_D(D,model,rot_D(data(1)%d(:,:,i_el),data(2)%d(i_el,1,1)))
	else
		call calc_D(D,model,rot_D(model%materials(:,:,model%material_nos(i_el)),data(2)%d(i_el,1,1)))
	end if

	call calc_B(B,model,i_el,coord)


! 	call calc_D(D,model,i_el)

! 	print *,model%materials(:,:,model%material_nos(i_el))

	rd_calc_BDB(:,:) = matmul( matmul(transpose(B),D),B)

end function

subroutine rd_calc_output(output,model,angle,u,D_)
	type(struct_model) :: model
	type(struct_output) :: output
	integer i,l,s,t,q
	integer nn,ne
	double precision thickness
	double precision,pointer,dimension(:,:) :: nodes
	integer,pointer,dimension(:,:) :: e_n
	integer dim,m
	double precision :: u(:),angle(:)
	double precision,optional :: D_(:,:,:)

	double precision,allocatable,dimension(:) :: ue,center_coord
	double precision,allocatable,dimension(:,:) :: B
	double precision,pointer,dimension(:,:) :: eps,sig
	double precision,pointer,dimension(:) :: msig

	double precision a,detJ
	integer gn

	nn = model%n_nds
	ne = model%n_els
	dim = model%dim
	m = model%n_nds_1el

	nodes => model%nodes
	e_n => model%elements

	allocate(ue(dim*m))
	allocate(output%eps(dim*(dim+1)/2,ne))
	allocate(output%sig(6,ne))
	allocate(output%msig(ne))
	allocate(output%maxpsig(ne))
	allocate(output%maxpsigx(ne))
	allocate(output%maxpsigy(ne))
	allocate(output%u(dim*nn))

	allocate(center_coord(dim))
	center_coord = 0d0

	eps => output%eps
	sig => output%sig
	msig => output%msig

	output%u(:) = u(:)

	do i=1,ne
		ue = 0d0

		do l=1,m
			ue( (l-1)*dim+1:l*dim ) = u( (e_n(l,i)-1)*dim+1:e_n(l,i)*dim )
		end do

		call calc_B(B,model,i,center_coord)

		eps(:,i) = matmul(B,ue)

! 		call calc_D(D,model,rot_D(data(1)%d(:,:,i_el),data(2)%d(i_el,1,1)))
! 	else
! 		call calc_D(D,model,rot_D(model%materials(:,:,model%material_nos(i_el)),data(2)%d(i_el,1,1)))


		if (present(D_)) then
			call calc_sig(output,model,i,rot_D(D_(:,:,i),angle(i)))
		else 
			call calc_sig(output,model,i,rot_D(model%materials(:,:,model%material_nos(i)),angle(i)))
		end if
		deallocate(B)

	end do

end subroutine


function rot_D(D_in,theta)
	implicit none
	integer i,j,k,l,p,q,r,s
	double precision theta,pi
	double precision,dimension(:,:) :: D_in
	double precision,allocatable,dimension(:,:) :: beta,C,D2,Cr,Dr2,Cr2D,rot_D
	double precision,allocatable,dimension(:,:,:,:) :: D4,Dr4

	allocate(beta(3,3))
	allocate(C(6,6))
	allocate(D2(6,6))
	allocate(rot_D(6,6))
	allocate(D4(3,3,3,3))
	allocate(Dr4(3,3,3,3))
	allocate(Dr2(6,6))
	allocate(Cr(6,6))
	allocate(Cr2D(3,3))

	! 円周率の定義
	pi = dacos(-1d0);


	! 座標軸の回転変換テンソルの定義
	beta(1,:) = (/dcos(theta), dsin(-theta), 0d0/);
	beta(2,:) = (/-dsin(-theta), dcos(theta), 0d0/);
	beta(3,:) = (/0d0, 0d0, 1d0/);

!print *,D2
	! 弾性マトリックスDの計算(inv(C))
	D2 = D_in;

! デバッグ用出力
!open(10,file='D2.csv')
!do i=1,6
!write(10,"(f30.15, 100(',', f30.15))") D2(i,:)
!end do
!close(10);


	! Dを4階のテンソルへ変換

	D4 = 0d0;

	do i=1,3
	do j=1,3
	do k=1,3
	do l=1,3
		r = i**2+j**2
		s = k**2+l**2
		if (r==2) then
			p = 1;
		else if (r==5) then
			p = 6;
		else if (r==10) then
			p = 5;
		else if (r==8) then
			p = 2;
		else if (r==13) then
			p = 4;
		else if (r==18) then
			p = 3;
		end if
		if (s==2) then
			q = 1;
		else if (s==5) then
			q = 6;
		else if (s==10) then
			q = 5;
		else if (s==8) then
			q = 2;
		else if (s==13) then
			q = 4;
		else if (s==18) then
			q = 3;
		end if

		D4(i,j,k,l) = D2(p,q);
	end do
	end do
	end do
	end do

! デバッグ用出力
!open(10,file='D4.csv')
!	do i=1,3
!	do j=1,3
!	do k=1,3
!write(10,"(f30.15, 100(',', f30.15))", advance='no') D4(i,j,k,:)
!end do
!write(10,*)
!end do
!end do
!close(10);

	! 弾性テンソルDの座標変換

	Dr4 = 0d0;

	do i=1,3
	do j=1,3
	do k=1,3
	do l=1,3
	do p=1,3
	do q=1,3
	do r=1,3
	do s=1,3
		Dr4(i,j,k,l) = Dr4(i,j,k,l)+beta(i,p)*beta(j,q)*beta(k,r)*beta(l,s)*D4(p,q,r,s);
	end do
	end do
	end do
	end do
	end do
	end do
	end do
	end do


! デバッグ用出力
!open(10,file='Dr4.csv')
!	do i=1,3
!	do j=1,3
!	do k=1,3
!write(10,"(f30.15, 100(',', f30.15))", advance='no') Dr4(i,j,k,:)
!end do
!write(10,*)
!end do
!end do
!close(10);


	do i=1,3
	do j=1,3
	do k=1,3
	do l=1,3
		r = i**2+j**2
		s = k**2+l**2
		if (r==2) then
			p = 1;
		else if (r==5) then
			p = 6;
		else if (r==10) then
			p = 5;
		else if (r==8) then
			p = 2;
		else if (r==13) then
			p = 4;
		else if (r==18) then
			p = 3;
		end if
		if (s==2) then
			q = 1;
		else if (s==5) then
			q = 6;
		else if (s==10) then
			q = 5;
		else if (s==8) then
			q = 2;
		else if (s==13) then
			q = 4;
		else if (s==18) then
			q = 3;
		end if

		Dr2(p,q) = Dr4(i,j,k,l);
	end do
	end do
	end do
	end do

! デバッグ用出力
!open(10,file='Dr2.csv')
!do i=1,6
!write(10,"(f30.15, 100(',', f30.15))") Dr2(i,:)
!end do
!close(10);

	rot_D = Dr2;

!	Cr = inv_jordan(Dr2);

!	Er_x = 1d0/Cr(1,1);
!	Er_y = 1d0/Cr(2,2);
!	Er_z = 1d0/Cr(3,3);
!	Gr_yz = 1d0/Cr(4,4);
!	Gr_zx = 1d0/Cr(5,5);
!	Gr_xy = 1d0/Cr(6,6);
!	vr_xy = -Cr(1,2)*Er_y;
!	vr_yz = -Cr(2,3)*Er_z;
!	vr_zx = -Cr(3,1)*Er_x;


! デバッグ用出力
!open(10,file='Cr.csv')
!do i=1,6
!write(10,"(f30.15, 100(',', f30.15))") Cr(i,:)
!end do
!close(10);

!Cr

! 物性値出力
!open(10,file='material_properties.txt')
!write(10,"('Young,s Moduli(GPa)  (E11,E22,E33)  =',3('  ',f12.8))"), Er_x/1d9, Er_y/1d9, Er_z/1d9
!write(10,"('shear Moduli(GPa)    (G23,G31,G12)  =',3('  ',f12.8))"), Gr_yz/1d9, Gr_zx/1d9, Gr_xy/1d9
!write(10,"('Poisson,s Ratios     (v23,v31,v12)  =',3('  ',f12.8))"), vr_yz, vr_zx, vr_xy
!close(10)

end function

! function inv_jordan(A0)
! 	implicit none
! 	integer n,i,j,k
! 	double precision,allocatable,dimension(:,:) :: L,inv_jordan,A
! 	double precision,dimension(:,:) :: A0

! 	n = ubound(A0,1)

! 	allocate(inv_jordan(n,n))
! 	allocate(A(n,n))

! 	A = A0

! 	inv_jordan = 0d0;
! 	do i=1,n
! 		inv_jordan(i,i) = 1d0;
! 	end do

! 	do k=1,n
! 		do i=1,n
! 		if (k /= i) then
! 			do j=1,k
! 				inv_jordan(i,j) = inv_jordan(i,j) - a(i,k)*inv_jordan(k,j)/a(k,k);
! 			end do
! 			do j=k+1,n
! 				a(i,j) = a(i,j) - a(i,k)*a(k,j)/a(k,k);
! 			end do
! 			a(i,k) = 0d0;
! 		end if
! 		end do
! 	end do

! 	do i=1,n
! 	do j=1,n
! 		inv_jordan(i,j) = inv_jordan(i,j)/a(i,i)
! 	end do
! 	end do

! end function

end module
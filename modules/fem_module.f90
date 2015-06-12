module fem_module
  use config

	implicit none
!	character(1) :: slash = "\"
!	character(32) :: path_model = "..\..\models\" ! 最後にスラッシュを必ずつける
!	double precision :: pi = dacos(-1d0)

	type struct_model
		character type_el*16
		character name*16
		double precision,pointer,dimension(:,:) :: nodes
		integer,pointer,dimension(:,:) :: elements
		double precision,pointer,dimension(:,:,:) :: materials
		integer,pointer :: material_nos(:)
		integer n_nds,n_els,dim,n_nds_1el,state2d
		double precision thickness
		type(struct_data),pointer :: data(:)

	end type struct_model

	type struct_data
		double precision,pointer :: d(:,:,:)
		integer,pointer :: i(:,:,:)
	end type

	type struct_bc
!		nodes_spc,disp_spc,nodes_mpc,dir_mpc
		integer,pointer,dimension(:,:) :: nodes_spc
		double precision,pointer,dimension(:,:) :: disp_spc
! 		integer,pointer,dimension(:,:) :: nodes_mpc,dir_mpc
		integer,pointer :: mpc_master_id(:),mpc_slave_id(:),mpc_master_dir(:),mpc_slave_dir(:)
		double precision,pointer :: mpc_param(:)
		integer n_spc,n_mpc
! 		double precision,pointer,dimension(:) :: f
	end type
	type struct_output
		double precision,pointer,dimension(:,:) :: eps,sig
		double precision,pointer,dimension(:) :: u,msig,maxpsig,maxpsigx,maxpsigy
		type(struct_data),pointer :: data(:)
	end type
contains

subroutine error(msg)
	character(*) msg

	print *,"ERROR: ", msg ; print *
 	print *,"- Press ENTER to finish -"
	read *

	stop
end subroutine

subroutine read_model(model,name_mdl)
	implicit none
	character fpath*32,name_mdl*16,errmsg*128
	integer iostat
	integer i,j
	integer dim,n_nds,n_els,n_nds_1el,tmp_n,tmp_mat,n_mat
 	double precision,pointer,dimension(:,:) :: nodes
	double precision,allocatable,dimension(:) :: tmp_nd
 	integer,pointer,dimension(:,:) :: elements
 	integer,pointer,dimension(:) :: material_nos
	integer,allocatable,dimension(:) :: tmp_el
 	double precision,pointer,dimension(:,:,:) :: materials,mats_

 	type(struct_model) :: model

 	model%name = name_mdl

! double precision :: E = 200d9
! double precision :: v = 0.3d0

	fpath = trim(path_model)//trim(name_mdl)//"/"


	! 材料数の読み込み

	open(10,file=trim(fpath)//'mats.csv',status='old', iostat=iostat)
	if (iostat /= 0) call error("Cannot open 'mats.csv'")
	read (10,*,iostat=iostat) n_mat
	if (iostat /= 0) call error("Cannot read the number of materials in 'mats.csv'")


 	allocate(mats_(6,6,n_mat))
 	allocate(model%materials(6,6,n_mat))

 	materials => model%materials
 	mats_ = 0d0

	do i=1,n_mat
		read (10,*,iostat=iostat) tmp_mat
		if (iostat /= 0) then
			if (i==1) then
				call error("Cannot read the 1st ID of the material in 'mats.csv'")
			else if (i==2) then
				call error("Cannot read the 2nd ID of the material in 'mats.csv'")
			else if (i==3) then
				call error("Cannot read the 3rd ID of the material in 'mats.csv'")
			else
				write(errmsg,'(A,i0,A)') "Cannot read the ", i, "th ID of the material in 'mats.csv'"
				call error(errmsg)
			end if
		end if
		do j=1,6
			read (10,*,iostat=iostat) materials(:,j,tmp_mat)
! 			read (10,*,iostat=iostat) model%materials(:,j,tmp_mat)
			if (iostat /= 0) then
				write(errmsg,'(A,i0,A)') "Cannot read values of the material ", tmp_mat, " in 'mats.csv'"
				call error(errmsg)
			end if
		end do
! 		print *,mats_(:,:,i)
! 		materials(:,:,i) = mats_(:,:,i)
	end do

	close(10)

	! モデルの読み込み
	open(10,file=trim(fpath)//'model.csv',status='old', iostat=iostat)
	if (iostat /= 0) call error("Cannot open 'model.csv'")

	! 要素タイプの読み込み
	read (10,*,iostat=iostat) model%type_el
	if (iostat /= 0) call error("Cannot read the type of element in 'model.csv'")

	select case(trim(model%type_el))
		case("quad4")
			n_nds_1el=4
			dim = 2
		case("hexa8")
			n_nds_1el=8
			dim = 3
		case default
	end select

	model%dim = dim
	model%n_nds_1el = n_nds_1el

	! 節点数,要素数の読み込み
	read (10,*,iostat=iostat) n_nds, n_els
	if (iostat /= 0) call error("Cannot read the number of elements and nodes in 'model.csv'")


 	allocate(model%elements(n_nds_1el+1,n_els))
 	allocate(model%nodes(dim,n_nds))
!  	allocate(model%materials(n_mat,6,6))
!  	allocate(model%materials(6,6,n_els))
 	allocate(model%material_nos(n_els))
 	allocate(tmp_nd(dim))
 	allocate(tmp_el(n_nds_1el))
 	elements => model%elements
 	nodes => model%nodes
! 	materials => model%materials
	material_nos => model%material_nos

 	model%n_nds = n_nds
 	model%n_els = n_els


!  	print *,elements(:,:)

	! 節点の読み込み
 	do i=1,n_nds
		read (10,*,iostat=iostat) tmp_n,tmp_nd(:)
		if (iostat /= 0) then
			if (i==1) then
				call error("Cannot read the 1st node in 'model.csv'")
			else if (i==2) then
				call error("Cannot read the 2nd node in 'model.csv'")
			else if (i==3) then
				call error("Cannot read the 3rd node in 'model.csv'")
			else
				write(errmsg,'(a,i0,a)') "Cannot read the ", i, "th node in 'model.csv'"
				call error(errmsg)
			end if
		end if
		nodes(:,tmp_n) = tmp_nd(:)
	end do

	! 要素の読み込み
	do i=1,n_els
		read (10,*,iostat=iostat) tmp_n,tmp_mat,tmp_el(:)
		if (iostat /= 0) then
			write(errmsg,'(a,i0,a)') "Cannot read the ", i, "th elements in 'model.csv'"
			call error(errmsg)
		end if
		elements(1:n_nds_1el,tmp_n) = tmp_el(:)
! 		elements(n_nds_1el+1,tmp_n) = tmp_mat
! 		materials(:,:,i) = mats_(:,:,tmp_mat)
		material_nos(i) = tmp_mat
	end do

	close(10)


	! 各要素の材料の回転角の初期化

! 	allocate(model%angle(n_els))
! 	model%angle = 0d0

! 	materials(1,1,:) = (/1-v, v, v, 0d0, 0d0, 0d0/)
! 	materials(1,2,:) = (/v, 1-v, v, 0d0, 0d0, 0d0/)
! 	materials(1,3,:) = (/v, v, 1-v, 0d0, 0d0, 0d0/)
! 	materials(1,4,4) = (1-2*v)/2
! 	materials(1,5,5) = (1-2*v)/2
! 	materials(1,6,6) = (1-2*v)/2
! 	materials(1,:,:) = materials(1,:,:)*E / ((1+v)*(1-2*v))
! 	allocate(model%d_data(1))
! 	allocate(model%name_d_data(1))
! 	model%name_d_data(1) = "_"
! 	allocate(model%d_data(1)%d(1,1,1))
! 	model%d_data(1)%d = 0d0

end subroutine

subroutine clear_model(model)
 	type(struct_model) :: model

	model%type_el = ''
	model%name = ''
	if (associated(model%nodes)) nullify(model%nodes)
	if (associated(model%nodes)) nullify(model%elements)
	if (associated(model%materials)) nullify(model%materials)
	if (associated(model%material_nos)) nullify(model%material_nos)
	model%n_nds = 0
	model%n_els = 0
	model%dim = 0
	model%n_nds_1el = 0
	model%state2d = 0
	model%thickness = 0d0
	if (associated(model%data)) nullify(model%data)
end subroutine

subroutine clear_bc(bc)
	type(struct_bc) :: bc

	if (associated(bc%nodes_spc)) nullify(bc%nodes_spc)
	if (associated(bc%disp_spc)) nullify(bc%disp_spc)
	if (associated(bc%mpc_master_id)) nullify(bc%mpc_master_id)
	if (associated(bc%mpc_slave_id)) nullify(bc%mpc_slave_id)
	if (associated(bc%mpc_master_dir)) nullify(bc%mpc_master_dir)
	if (associated(bc%mpc_slave_dir)) nullify(bc%mpc_slave_dir)

	if (associated(bc%mpc_param)) nullify(bc%mpc_param)
	bc%n_spc = 0
	bc%n_mpc = 0
end subroutine

subroutine clear_output(output)
	type(struct_output) :: output

	if (associated(output%eps)) nullify(output%eps)
	if (associated(output%sig)) nullify(output%sig)
	if (associated(output%u)) nullify(output%u)
	if (associated(output%msig)) nullify(output%msig)
	if (associated(output%maxpsig)) nullify(output%maxpsig)
	if (associated(output%maxpsigx)) nullify(output%maxpsigx)
	if (associated(output%maxpsigy)) nullify(output%maxpsigy)
	if (associated(output%data)) nullify(output%data)
end subroutine

subroutine set_state2d(model,state2d,thickness)
	type(struct_model) :: model
	integer state2d
	double precision thickness

	model%state2d = state2d
	model%thickness = thickness

! 	if (model%dim == 2) then
! ! 		! 厚さの読み込み
! ! 		print *, "thickness:"
! ! 		read *, thickness
! ! 		! 平面応力/平面ひずみ
! ! 		print *, "state: (0: plane stress,  >0: plane strain)"
! ! 		read *, state2d
! 	else
! 		state2d = 0
! 		thickness = 0d0
! 	end if
end subroutine

subroutine init_addition_matrix_sln(k,sln,model,bc)
	implicit none
	integer i,j,s,t,p,q,r,Ki,Kj
	type(struct_model) :: model
	type(struct_bc) :: bc
	double precision,pointer,dimension(:,:) :: nodes
	integer,pointer,dimension(:,:) :: e_n

	integer m,dim
	double precision,allocatable,dimension(:) :: k
	integer,allocatable,dimension(:) :: sln,l
	integer nn,ne
	integer slm,nk

	dim = model%dim
	m = model%n_nds_1el


	nodes => model%nodes
	e_n => model%elements

	nn = model%n_nds
	ne = model%n_els


	! スカイラインマトリックスの格納の準備
	allocate(sln(nn*dim))
	allocate(l(nn*dim))

	l = 0
	do i=1,ne

		slm = nn*dim

		do s=1,m
		do t=1,dim
			if (slm> ((e_n(s,i)-1)*dim+t) ) slm = (e_n(s,i)-1)*dim+t
		end do
		end do

		do s=1,m
		do t=1,dim
			if (l((e_n(s,i)-1)*dim+t) < ((e_n(s,i)-1)*dim+t-slm+1)) then
				l((e_n(s,i)-1)*dim+t) = (e_n(s,i)-1)*dim+t-slm + 1
			end if
		end do
		end do

	end do

	nk = sum(l)


	sln = 0
	sln(1) = 1

	do i=2,nn*dim
		sln(i) = sln(i-1)+l(i)
	end do

	allocate(K(nk))
	k = 0d0

! 多点拘束の為の処理
	if (bc%n_mpc > 0) then

	do i=1,ne
		do s=1,m
		do t=1,m
			do q=1,dim
			do r=1,dim
				Ki = (e_n(s,i)-1)*dim+q
				kj = (e_n(t,i)-1)*dim+r
				if (Ki <= Kj) then
					K(sln(Kj)+Ki-Kj) = 1d0
				end if
			end do
			end do
		end do
		end do			
	end do
	do p=1,bc%n_mpc
		s = (bc%mpc_master_id(p)-1)*dim+bc%mpc_master_dir(p)
		t = (bc%mpc_slave_id(p)-1)*dim+bc%mpc_slave_dir(p)
! 		print *,s,t
		if (s>t) then
			if ( (t-l(t)+1) < (s-l(s)+1) ) then
! 				print *,l(s)
				l(s) = s - (t-l(t)+1) + 1
! 				print *,l(s)
			end if
		else if (s<t) then
			if ( (t-l(t)+1) < (s-l(s)+1) ) then
! 				print *,l(s)
				l(s) = s - (t-l(t)+1) + 1
! 				print *,l(s)
			end if
			do i=s,t-1
				if ( (sln(t)+i-t>sln(t-1)) .and. (K(sln(t)+i-t) > 0.01) .and. (i-l(i)+1 > s) ) then
					l(i) = i - s + 1
				end if
			end do
			do i=t,nn*dim
				if ( (sln(i)+t-i>sln(i-1)) .and. (K(sln(i)+t-i) > 0.01) .and. (i-l(i)+1 > s) ) then
					l(i) = i - s + 1
				end if
			end do

		end if
	end do

! Kの再定義

	nk = sum(l)
	sln(1) = 1
	do i=2,nn*dim
		sln(i) = sln(i-1)+l(i)
	end do
	deallocate(K)
	allocate(K(nk))
	k = 0d0

	end if

end subroutine

subroutine clear_addition_matrix_sln(k,sln)
	real(8),allocatable :: k(:)
	integer,allocatable :: sln(:)

	if (allocated(k)) deallocate(k)
	if (allocated(sln)) deallocate(sln)
end subroutine

! subroutine calc_element_integral(I,model,i_el,func,data)
subroutine calc_element_integral(I,model,i_el,func)
	interface
!     	function func(model,i_el,coord,data)
    	function func(model,i_el,coord)
    		import struct_model,struct_data
			type(struct_model) :: model
			integer :: i_el
			double precision :: coord(:)
! 			type(struct_data),optional :: data(:)
			double precision,allocatable,dimension(:,:) :: func
    	end function
    end interface

	type(struct_model) :: model
	integer i_el
! 	type(struct_data),optional :: data(:)
	integer k,n_gp
	double precision :: I(:,:)
	double precision,allocatable,dimension(:) :: w
	double precision,allocatable,dimension(:,:) :: coord_gp
	double precision a,detJ

	select case(trim(model%type_el))
		case("quad4")
			n_gp = 4
			allocate(coord_gp(2,4))

			! ガウス点の座標
			a = dsqrt(1d0/3d0)
			coord_gp(1,:) = (/-a, a, a, -a/)
			coord_gp(2,:) = (/-a, -a, a, a/)

			! ガウス点の重み
			w = (/1., 1., 1., 1./)

			I = 0d0

			do k=1,4
				call calc_detJ(detJ,model,i_el,coord_gp(:,k))
! 				I = I + func(model,i_el,coord_gp(:,k),data)*detJ*w(k)*model%thickness
				I = I + func(model,i_el,coord_gp(:,k))*detJ*w(k)*model%thickness
			end do

		case("hexa8")
			! ガウス点の数
			allocate(coord_gp(3,8))


			! ガウス点の座標
			a = dsqrt(1d0/3d0)
			coord_gp(1,:) = (/-a, a, a, -a, -a, a, a, -a/)
			coord_gp(2,:) = (/-a, -a, a, a, -a, -a, a, a/)
			coord_gp(3,:) = (/-a, -a, -a, -a, a, a, a, a/)

			! ガウス点の重み
			w = (/1., 1., 1., 1., 1., 1., 1., 1./)

			I = 0d0
			do k=1,8
				call calc_detJ(detJ,model,i_el,coord_gp(:,k))
! 				I = I + func(model,i_el,coord_gp(:,k),data)*detJ*w(k)
				I = I + func(model,i_el,coord_gp(:,k))*detJ*w(k)
			end do
	end select


end subroutine

subroutine calc_detJ(detJ,model,i_el,coord)
	type(struct_model) :: model
	double precision ::	coord(:)
	double precision,pointer :: nodes(:,:)
	integer,pointer :: elements(:,:)
	integer dim,n_nds_1el,i_el
	integer s,t,q

	double precision detJ
	double precision,allocatable :: J(:,:),dN_dxi(:,:)

	nodes => model%nodes
	elements => model%elements
	dim = model%dim
	n_nds_1el = model%n_nds_1el

	allocate(J(dim,dim))
	J = 0d0

	select case(trim(model%type_el))
		case("quad4")
			allocate(dN_dxi(2,4))
			dN_dxi = Nq4(coord)
		case("hexa8")
			allocate(dN_dxi(3,8))
			dN_dxi = Nh8(coord)
	end select

	! ヤコビ行列の設定
	do s=1,dim
	do t=1,dim
		do q=1,n_nds_1el
			J(s,t) = J(s,t) + dN_dxi(s,q) * nodes(t,elements(q,i_el))
		end do
	end do
	end do


	if (dim == 3) then
		! ヤコビアン
		detJ = J(1,1)*J(2,2)*J(3,3) + J(2,1)*J(3,2)*J(1,3) + J(3,1)*J(1,2)*J(2,3)
		detJ = detJ - J(1,1)*J(3,2)*J(2,3) - J(3,1)*J(2,2)*J(1,3) - J(2,1)*J(1,2)*J(3,3)

	else if (dim == 2) then
		! ヤコビアン
		detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)
	end if

end subroutine

! subroutine calc_K_sln(k,model,index_k_diag)
! 	type(struct_model) :: model
! 	double precision :: k(:)
! 	integer :: index_k_diag(:)

! 	call calc_addition_matrix_sln(k,model,index_k_diag,calc_BDB)

! end subroutine



subroutine calc_addition_matrix_sln(k,model,sln,func,data)
	integer i,s,t,q,r
	type(struct_model) :: model
	type(struct_data),optional :: data(:)
	integer n_nds,n_els,m,dim,ki,kj
	double precision,pointer,dimension(:,:) :: nodes
	integer,pointer,dimension(:,:) :: elements
	double precision,pointer,dimension(:,:,:) :: materials
	double precision,allocatable,dimension(:,:) :: ke

	double precision,dimension(:) :: k
	integer,dimension(:) :: sln

! 	double precision,pointer :: D(:,:,:)

	interface
!     	function func(model,i_el,coord,data)
    	function func(model,i_el,coord)
    		import struct_model,struct_data
			type(struct_model) :: model
			integer :: i_el
			double precision :: coord(:)
! 			type(struct_data),optional :: data(:)
			double precision,allocatable,dimension(:,:) :: func
    	end function
    end interface

	dim = model%dim
	m = model%n_nds_1el

	n_nds = model%n_nds
	n_els = model%n_els

	nodes => model%nodes
	elements => model%elements
	materials => model%materials

	do i=1,n_els
		allocate(ke(dim*m,dim*m))

! 		if (present(data)) then
! 			call calc_element_integral(Ke,model,i,func,data)
! 		else
			call calc_element_integral(Ke,model,i,func)
! 		end if

		do s=1,m
		do t=1,m
			do q=1,dim
			do r=1,dim
				Ki = (elements(s,i)-1)*dim+q
				kj = (elements(t,i)-1)*dim+r
				if (Ki <= Kj) then
					K(sln(Kj)+Ki-Kj) = K(sln(Kj)+Ki-Kj) + Ke((s-1)*dim+q,(t-1)*dim+r)
				end if

			end do
			end do
		end do
		end do			
		deallocate(ke)
	end do


end subroutine

! function calc_BDB(model,i_el,coord,data)
function calc_BDB(model,i_el,coord)
	type(struct_model) :: model
	integer :: i_el
	double precision :: coord(:)
! 	type(struct_data),optional :: data(:)
	double precision,allocatable,dimension(:,:) :: calc_BDB
	double precision,allocatable :: B(:,:),D(:,:)
	double precision :: D_(6,6)

	allocate(calc_BDB(model%dim*model%n_nds_1el,model%dim*model%n_nds_1el))
	calc_BDB = 0d0

	call calc_B(B,model,i_el,coord)

! 	if (present(data)) then
! 		if (ubound(data(1)%d,3)>1) then
! 			D_ = data(1)%d(:,:,i_el)
! 		else
! 			D_ = data(1)%d(:,:,1)
! 		end if
! 	else 
		D_ = model%materials(:,:,model%material_nos(i_el))
! 	end if

	call calc_D(D,model,D_)

	calc_BDB(:,:) = matmul( matmul(transpose(B),D),B)

	deallocate(D,B)

end function

subroutine calc_B(B,model,i_el,coord)
	type(struct_model) :: model
	integer :: i_el
	double precision :: coord(:)
	double precision,allocatable :: B(:,:),dN_dxi(:,:),Jinv(:,:),dN_dx(:,:)
	integer q

	select case(trim(model%type_el))
		case("quad4")

			allocate(dN_dxi(2,4))
			allocate(dN_dx(2,4))
			allocate(B(3,8))
			B = 0d0

			dN_dxi = Nq4(coord)

			call calc_Jinv(Jinv,model,i_el,coord,dN_dxi)

			dN_dx = matmul(Jinv, dN_dxi)

			do q=1,4
				B(:, 2*q-1) = (/ dN_dx(1,q), 0d0, dN_dx(2,q) /)
				B(:, 2*q) = (/ 0d0, dN_dx(2,q), dN_dx(1,q) /)
			end do

		case("hexa8")

			allocate(dN_dxi(3,8))
			allocate(dN_dx(3,8))
			allocate(B(6,24))
			B = 0d0

			dN_dxi = Nh8(coord)

			call calc_Jinv(Jinv,model,i_el,coord,dN_dxi)

			dN_dx = matmul(Jinv, dN_dxi)

			do q=1,8		
				B(:, 3*q-2) = (/ dN_dx(1,q), 0d0, 0d0, 0d0, dN_dx(3,q), dN_dx(2,q) /)
				B(:, 3*q-1) = (/ 0d0, dN_dx(2,q), 0d0, dN_dx(3,q), 0d0, dN_dx(1,q) /)
				B(:, 3*q) = (/ 0d0, 0d0, dN_dx(3,q), dN_dx(2,q), dN_dx(1,q), 0d0 /)
			end do
	end select

end subroutine

subroutine calc_Jinv(Jinv,model,i_el,coord,dN_dxi)
	type(struct_model) :: model
	double precision,allocatable :: Jinv(:,:)
	double precision ::	coord(:),dN_dxi(:,:)
	double precision,pointer :: nodes(:,:)
	integer,pointer :: elements(:,:)
	integer dim,n_nds_1el,i_el
	integer s,t,q

	double precision detJ
	double precision,allocatable :: J(:,:)

	nodes => model%nodes
	elements => model%elements
	dim = model%dim
	n_nds_1el = model%n_nds_1el

	allocate(Jinv(dim,dim))
	allocate(J(dim,dim))

	J = 0d0
	Jinv = 0d0

	! ヤコビ行列の設定
	do s=1,dim
	do t=1,dim
		do q=1,n_nds_1el
			J(s,t) = J(s,t) + dN_dxi(s,q) * nodes(t,elements(q,i_el))
		end do
	end do
	end do

	if (dim == 3) then
		! ヤコビアン
		detJ = J(1,1)*J(2,2)*J(3,3) + J(2,1)*J(3,2)*J(1,3) + J(3,1)*J(1,2)*J(2,3)
		detJ = detJ - J(1,1)*J(3,2)*J(2,3) - J(3,1)*J(2,2)*J(1,3) - J(2,1)*J(1,2)*J(3,3)

		! ヤコビ行列の逆行列
		Jinv(1,:) = (/J(2,2)*J(3,3)-J(2,3)*J(3,2), J(1,3)*J(3,2)-J(1,2)*J(3,3), J(1,2)*J(2,3)-J(1,3)*J(2,2)/)
		Jinv(2,:) = (/J(2,3)*J(3,1)-J(2,1)*J(3,3), J(1,1)*J(3,3)-J(1,3)*J(3,1), J(1,3)*J(2,1)-J(1,1)*J(2,3)/)
		Jinv(3,:) = (/J(2,1)*J(3,2)-J(2,2)*J(3,1), J(1,2)*J(3,1)-J(1,1)*J(3,2), J(1,1)*J(2,2)-J(1,2)*J(2,1)/)
		Jinv = Jinv/detJ

	else if (dim == 2) then
		! ヤコビアン
		detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)

		! ヤコビ行列の逆行列
		Jinv = reshape( (/J(2,2), -J(2,1), -J(1,2), J(1,1)/), (/dim,dim/) ) / detJ
	end if

end subroutine
! subroutine calc_D(D,model,)

subroutine calc_D(D,model,D_)
	type(struct_model) :: model
	integer dim,m,state2d,i
	double precision :: D_(:,:)
	integer,pointer,dimension(:,:) :: e_n
	double precision,pointer,dimension(:,:,:) :: materials
	double precision,allocatable,dimension(:,:) :: D3,C,C2,D


	dim = model%dim
	m = model%n_nds_1el

	e_n => model%elements
	materials => model%materials
	state2d = model%state2d

	allocate(D3(6,6))
	allocate(C(6,6))
	allocate(C2(3,3))

	if (dim == 2) then
		allocate(D(3,3))
		if (state2d == 0) then
! 			C = inv_jordan(materials(:,:,i))
			C = inv_jordan(D_)
			C2(1:2,1:2) = C(1:2,1:2)
			C2(3,1:2) = C(6,1:2)
			C2(1:2,3) = C(1:2,6)
			C2(3,3) = C(6,6)
			D = inv_jordan(C2)
		else if (state2d > 0) then
! 			D3 = materials(:,:,i)
			D3 = D_
			D(1:2,1:2) = D3(1:2,1:2)
			D(3,1:2) = D3(6,1:2)
			D(1:2,3) = D3(1:2,6)
			D(3,3) = D3(6,6)
		end if
	else if (dim == 3) then
		allocate(D(6,6))
! 		D(:,:) = materials(:,:,i)
		D(:,:) = D_(:,:)
	end if
end subroutine


subroutine read_f(f,model)
	type(struct_model) :: model
	integer i
	character fpath*32
	double precision,allocatable,dimension(:) :: f
	double precision,allocatable,dimension(:) :: tmp_f
	integer nn,dim,m
	integer fn,tmp_n

	integer iostat
	character errmsg*128


	dim = model%dim
	m = model%n_nds_1el
	nn = model%n_nds

	fpath = trim(path_model)//trim(model%name)//"/"

	! 等価節点力ベクトルの読み込み
	allocate(tmp_f(dim))
	allocate(f(dim*nn))

	f = 0d0
	open(10,file=trim(fpath)//'f.csv',status='old',iostat=iostat)
	if (iostat /= 0) call error("Cannot open 'f.csv'")

	read(10,*,iostat=iostat) fn

	do i=1,fn
		read (10,*,iostat=iostat) tmp_n,tmp_f(:)
		if (iostat /= 0) then
			if (i==1) then
				call error("Cannot read the 1st vector of force in 'f.csv'")
			else if (i==2) then
				call error("Cannot read the 2nd vector of force in 'f.csv'")
			else if (i==3) then
				call error("Cannot read the 3rd vector of force in 'f.csv'")
			else
				write(errmsg,'(a,i0,a)') "Cannot read the ", i, "th value of force in 'f.csv'"
				call error(errmsg)
			end if
		end if

		f((tmp_n-1)*dim+1:tmp_n*dim) = tmp_f(:)
	end do
	close(10)

end subroutine

subroutine calc_equivalent_nodal_force(f,model,f_sum,a,b,c,d)
	integer i,j,k,l
	integer :: n_els,dim,n_ptns,n_ptns_1el,n_nds_1ptn,n_nds
	integer,pointer :: elements(:,:)
	integer,allocatable :: patterns(:,:)
	real(8),pointer :: nodes(:,:)
	real(8),allocatable :: f(:),params(:),f_1ptn(:)
	real(8) :: f_sum(:),s,s_sum,a,b,c
	real(8) :: v1(3),v2(3)
	real(8),optional :: d

	type(struct_model) :: model
	logical :: flag

	dim = model%dim
	n_els = model%n_els
	n_nds = model%n_nds

	elements => model%elements
	nodes => model%nodes

	select case(dim)
	case(2) ! 2D 辺の取得

		! 要素内の辺のパターン(2D)
		n_nds_1ptn = 2
		n_ptns_1el = 4
		allocate(patterns(n_nds_1ptn,n_ptns_1el))
		allocate(params(3))
		params(1) = a
		params(2) = b
		params(3) = c
		patterns(:,1) = (/1,2/)
		patterns(:,2) = (/2,3/)
		patterns(:,3) = (/3,4/)
		patterns(:,4) = (/4,1/)

	case(3) ! 3D 面の取得

		! 要素内の辺のパターン(3D)
		n_nds_1ptn = 4
		n_ptns_1el = 6
		allocate(patterns(n_nds_1ptn,n_ptns_1el))
		allocate(params(4))
		params(1) = a
		params(2) = b
		params(3) = c
		params(4) = d
		patterns(:,1) = (/1,2,3,4/)
		patterns(:,2) = (/1,2,6,5/)
		patterns(:,3) = (/2,3,7,6/)
		patterns(:,4) = (/3,4,8,7/)
		patterns(:,5) = (/4,1,5,8/)
		patterns(:,6) = (/5,6,7,8/)

	case default
	end select

	n_ptns = 0
	s_sum = 0d0
	do i=1,n_els
	do j=1,n_ptns_1el
		flag = .true.
		do k=1,n_nds_1ptn
			! 3D: if ax+by+cz+d!=0, 2D: if ax+by+c!=0
 			if ( abs(dot_product( params(1:dim),nodes(:,elements(patterns(k,j),i)) ) + params(dim+1)) > 1d-10) then
 				flag = .false.
 			end if
 		end do

 		if (flag) then
 			if (dim == 2) then
 				! 要素の辺の長さを計算、足し込み
 				s_sum = s_sum + dsqrt( ( nodes(1,elements(patterns(2,j),i)) - nodes(1,elements(patterns(1,j),i)) )**2 + &
 								&( nodes(2,elements(patterns(2,j),i)) - nodes(2,elements(patterns(1,j),i)) )**2 )
 			else if (dim == 3) then
 				! 要素の面の面積を計算、足し込み
 				v1 = nodes(:,elements(patterns(2,j),i)) - nodes(:,elements(patterns(1,j),i))
 				v2 = nodes(:,elements(patterns(3,j),i)) - nodes(:,elements(patterns(1,j),i))
 				print *,v1,v2
 				s_sum = s_sum + 1d0/2d0*dsqrt( (dot_product(v1,v1)*dot_product(v2,v2)) - dot_product(v1,v2)**2 )
 			end if
 			n_ptns = n_ptns + 1
 		end if
	end do
	end do


	allocate(f_1ptn(dim))

! 	f_1ptn = f_sum/n_ptns

	print *,n_ptns

	allocate(f(n_nds*dim))
	f = 0d0

	do i=1,n_els
	do j=1,n_ptns_1el
		flag = .true.
		do k=1,n_nds_1ptn
			! 3D: if ax+by+cz+d!=0, 2D: if ax+by+c!=0
 			if ( abs(dot_product( params(1:dim),nodes(:,elements(patterns(k,j),i)) ) + params(dim+1)) > 1d-10) then
 				flag = .false.
 			end if
 		end do

 		if (flag) then
 			if (dim == 2) then
 				! 要素の辺の長さを計算、足し込み
 				s = dsqrt( ( nodes(1,elements(patterns(2,j),i)) - nodes(1,elements(patterns(1,j),i)) )**2 + &
 								&( nodes(2,elements(patterns(2,j),i)) - nodes(2,elements(patterns(1,j),i)) )**2 )
 			else if (dim == 3) then
 				! 要素の面の面積を計算、足し込み
 				v1 = nodes(:,elements(patterns(2,j),i)) - nodes(:,elements(patterns(1,j),i))
 				v2 = nodes(:,elements(patterns(3,j),i)) - nodes(:,elements(patterns(1,j),i))
 				s = 1d0/2d0*dsqrt( (dot_product(v1,v1)*dot_product(v2,v2)) - dot_product(v1,v2)**2 )
 			end if

 			f_1ptn(1:dim) = f_sum*s/s_sum
			do k=1,n_nds_1ptn
				l = elements(patterns(k,j),i)
				f((l-1)*dim+1:l*dim) = f((l-1)*dim+1:l*dim) + f_1ptn(1:dim)/n_nds_1ptn

 			end do
 		end if
	end do
	end do

	do i=1,n_nds
		print *,f((i-1)*dim+1:i*dim)
	end do


end subroutine

subroutine read_cc(bc,model)
	type(struct_bc) :: bc
	type(struct_model) :: model

	character fpath*32
	integer i
	integer nn,dim,m
	integer,pointer,dimension(:,:) :: cn,mpcdir,nodes_mpc
	double precision,pointer,dimension(:,:) :: cnval
	integer tmp_n,ncn,nmpc

	integer iostat
	character errmsg*128

	fpath = trim(path_model)//trim(model%name)//"/"

	dim = model%dim
	m = model%n_nds_1el
	nn = model%n_nds

	! 単点拘束条件の読み込み

	open(10,file=trim(fpath)//'cn.csv',status='old',iostat=iostat)
	if (iostat /= 0) call error("Cannot open 'cn.csv'")
	read (10,*,iostat=iostat) ncn
	if (iostat /= 0) call error("Cannot read the number of spc conditions in 'cn.csv'")

	allocate(bc%nodes_spc(dim+1,ncn))
	allocate(bc%disp_spc(dim,ncn))

	cn => bc%nodes_spc
	cnval => bc%disp_spc
	bc%n_spc = ncn

	cn = 0d0

	do i=1,ncn
		read (10,*,iostat=iostat) cn(1:(dim+1),i), cnval(1:dim,i)
		if (iostat /= 0) then
			write(errmsg,'(a,i0,a)') "Cannot read the ", i, "th spc conditions"
			call error(errmsg)
		end if
	end do
	close(10)

	! 多点拘束条件の読み込み

	open(10,file=trim(fpath)//'mpc.csv',status='old', iostat=iostat)
	if (iostat /= 0) call error("Cannot open 'mpc.csv'")
	read (10,*,iostat=iostat) nmpc
	if (iostat /= 0) call error("Cannot read the number of mpc conditions in 'cn.csv'")

	bc%n_mpc = nmpc

if (nmpc /= 0) then
! 	allocate(bc%dir_mpc(dim,nmpc))
! 	allocate(bc%nodes_mpc(2,nmpc))

! 	mpcdir => bc%dir_mpc
! 	nodes_mpc => bc%nodes_mpc

	allocate(bc%mpc_master_id(nmpc))
	allocate(bc%mpc_master_dir(nmpc))
	allocate(bc%mpc_slave_id(nmpc))
	allocate(bc%mpc_slave_dir(nmpc))
	allocate(bc%mpc_param(nmpc))

	do i=1,nmpc
! 		read (10,*,iostat=iostat) nodes_mpc(:,i),mpcdir(1:dim,i)
		read (10,*,iostat=iostat) bc%mpc_master_id(i),bc%mpc_master_dir(i),bc%mpc_slave_id(i),&
												&bc%mpc_slave_dir(i),bc%mpc_param(i)
		if (iostat /= 0) then
			write(errmsg,'(a,i0,a)') "Cannot read the ", i, "th mpc conditions"
			call error(errmsg)
		end if
	end do

end if


	close(10)

end subroutine

subroutine set_cc(k,f,model,bc,sln)
	type(struct_model) :: model
	type(struct_bc) :: bc

	integer i,l,n,q,s,t
	integer nn,dim,m
	integer,pointer,dimension(:,:) :: cn,mpcdir,nodes_mpc
	double precision,pointer,dimension(:,:) :: cnval
	integer fn,tmp_n,ncn,nmpc
	integer,allocatable,dimension(:) :: sln
	double precision,dimension(:) :: f,k

	nn = model%n_nds
	dim = model%dim
	m = model%n_nds_1el

	cn => bc%nodes_spc
	cnval => bc%disp_spc
! 	mpcdir => bc%dir_mpc
! 	nodes_mpc => bc%nodes_mpc
	nmpc = bc%n_mpc
	ncn = bc%n_spc


	! 多点拘束条件の設定

! 	print *,"add mpc"

! 	print *,bc%mpc_param

	if (nmpc /= 0) then
	do i=1,nmpc
! 	do l=1,dim
! 	if (mpcdir(l,i) /= 0) then
! 		s = (nodes_mpc(1,i)-1)*dim+l
! 		t = (nodes_mpc(2,i)-1)*dim+l

		s = (bc%mpc_master_id(i)-1)*dim+bc%mpc_master_dir(i)
		t = (bc%mpc_slave_id(i)-1)*dim+bc%mpc_slave_dir(i)

		if (s < t) then
! 			do q=1,s
! 			if (s > 1) then
! 				if ((sln(s)+q-s>sln(s-1)) .and. (sln(t)+q-t>sln(t-1))) then
! 					K(sln(s)+q-s) = K(sln(s)+q-s)+K(sln(t)+q-t)
! 				end if
! 			else
! 				if (sln(t)+q-t>sln(t-1)) then
! 					K(1) = K(1)+K(sln(t)+q-t)
! 				end if
! 			end if
! 			end do
!!!
			do q=1,s-1
				if ((sln(s)+q-s>sln(s-1)) .and. (sln(t)+q-t>sln(t-1))) then
					K(sln(s)+q-s) = K(sln(s)+q-s)+K(sln(t)+q-t)*bc%mpc_param(i)
! if (i==1) print *,K(sln(s)+q-s),q,s
				end if
			end do

			K(sln(s)) = K(sln(s))+K(sln(t))*bc%mpc_param(i)
			if (sln(t)+s-t>sln(t-1)) then
				K(sln(s)) = K(sln(s))+K(sln(t)+s-t)*bc%mpc_param(i)*2
			end if
! if (i==1) print *,K(sln(s)),s,s
!!!

			do q=(s+1),t
				if ((sln(q)+s-q>sln(q-1)) .and. (sln(t)+q-t>sln(t-1))) then
					K(sln(q)+s-q) = K(sln(q)+s-q)+K(sln(t)+q-t)*bc%mpc_param(i)
! if (i==1) print *,K(sln(q)+s-q),q,s
				end if
			end do

			do q=(t+1),nn*dim
				if ((sln(q)+s-q>sln(q-1)) .and. (sln(q)+t-q>sln(q-1))) then
					K(sln(q)+s-q) = K(sln(q)+s-q)+K(sln(q)+t-q)*bc%mpc_param(i)
! if (i==1) print *,K(sln(q)+s-q),q,s
				end if
			end do
		else if (s > t) then
			if  (sln(s)+t-s>sln(s-1)) then
				K(sln(s)) = K(sln(s))+K(sln(s)+t-s)*bc%mpc_param(i)*2
			end if
			K(sln(s)) = K(sln(s))+K(sln(t))*bc%mpc_param(i)
! if (i==1) print *,K(sln(s)),s,s

			do q=1,t
			if (t > 1) then
				if ((sln(s)+q-s>sln(s-1)) .and. (sln(t)+q-t>sln(t-1))) then
					K(sln(s)+q-s) = K(sln(s)+q-s)+K(sln(t)+q-t)*bc%mpc_param(i)
! if (i==1) print *,K(sln(s)+q-s),q,s
				end if
			else
				if (sln(s)+q-s>sln(s-1)) then
					K(sln(s)+q-s) = K(sln(s)+q-s)+K(1)*bc%mpc_param(i)
				end if
			end if
			end do

! 			do q=(t+1),s
! 				if ((sln(s)+q-s>sln(s-1)) .and. (sln(q)+t-q>sln(q-1))) then
! 					K(sln(s)+q-s) = K(sln(s)+q-s)+K(sln(q)+t-q)
! 				end if
! 			end do
!!!
			do q=(t+1),s-1
				if ((sln(s)+q-s>sln(s-1)) .and. (sln(q)+t-q>sln(q-1))) then
					K(sln(s)+q-s) = K(sln(s)+q-s)+K(sln(q)+t-q)*bc%mpc_param(i)
! if (i==1) print *,K(sln(s)+q-s),q,s
				end if
			end do

! 			K(sln(s)) = K(sln(s))+K(sln(t))*bc%mpc_param(i)

!!!

			do q=(s+1),nn*dim
				if ((sln(q)+s-q>sln(q-1)) .and. (sln(q)+t-q>sln(q-1))) then
					K(sln(q)+s-q) = K(sln(q)+s-q)+K(sln(q)+t-q)*bc%mpc_param(i)
! if (i==1) print *,K(sln(q)+s-q),q,s
				end if
			end do
		end if
		K(sln(t)) = 1d30
		f(s) = f(s) + f(t)
		f(t) = 0d0
! 	end if
! 	end do
	end do

	end if


	! 拘束条件の設定

	do i=1,ncn
	do l=1,dim
		if (cn(l+1,i) /= 0) then

			n = (cn(1,i)-1)*dim+l
			do q=1,n-1
				if ((sln(n)+q-n)>sln(n-1)) then
					f(q) = f(q) - K(sln(n)+q-n)*cnval(l,i)
				end if
			end do
			f(n) = f(n) - K(sln(n))*cnval(l,i)
			do q=n+1,nn*dim
				if ((sln(q)+n-q)>sln(q-1)) then
					f(q) = f(q) - K(sln(q)+n-q)*cnval(l,i)
				end if
			end do

			K(sln(n)) = 1d30

		end if
	end do
	end do


end subroutine


subroutine set_constrained_u(u,model,bc)
	type(struct_model) :: model
	type(struct_bc) :: bc
	integer l,i
	integer,pointer,dimension(:,:) :: nodes_spc,dir_mpc,nodes_mpc
	double precision,pointer,dimension(:,:) :: disp_spc
	double precision,dimension(:) :: u
	integer dim,n_nds_1el,ncn,nmpc,n

	dim = model%dim
	n_nds_1el = model%n_nds_1el

	ncn = bc%n_spc
	nmpc = bc%n_mpc
	nodes_spc => bc%nodes_spc
	disp_spc => bc%disp_spc
! 	nodes_mpc => bc%nodes_mpc
! 	dir_mpc => bc%dir_mpc


	! 拘束点に変位を代入し直す

	do i=1,ncn
	do l=1,dim
		n=(nodes_spc(1,i)-1)*dim+l
		if (nodes_spc(l+1,i) /= 0) then
			u(n) = disp_spc(l,i)
		end if
	end do
	end do

	do i=1,nmpc
! 	do l=1,dim
! 		if (dir_mpc(l,i) /= 0) then
! 			u( (nodes_mpc(2,i)-1)*dim+l ) = u( (nodes_mpc(1,i)-1)*dim+l )
			u( (bc%mpc_slave_id(i)-1)*dim+bc%mpc_slave_dir(i) ) = u( (bc%mpc_master_id(i)-1)*dim+bc%mpc_master_dir(i) )
! 		end if
! 	end do
	end do

end subroutine

subroutine calc_output(output,model,u,D_)
	type(struct_model) :: model
	type(struct_output) :: output
	integer i,l,s,t,q
	integer nn,ne
	double precision thickness
	double precision,pointer,dimension(:,:) :: nodes
	integer,pointer,dimension(:,:) :: e_n
	integer dim,m
	double precision :: u(:)
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

		if (present(D_)) then
			call calc_sig(output,model,i,D_(:,:,i))
		else 
			call calc_sig(output,model,i,model%materials(:,:,model%material_nos(i)))
		end if
		deallocate(B)

	end do

end subroutine

subroutine output_inp(model,output,file_name)
	type(struct_model) :: model
	type(struct_output) :: output
	character(*) file_name
	integer i
	double precision,allocatable,dimension(:,:) :: nd_data
	double precision,allocatable,dimension(:,:) :: el_data
	character(16),allocatable,dimension(:) :: nd_data_name
	character(16),allocatable,dimension(:) :: el_data_name


	if (model%dim == 2) then

	allocate(nd_data(2,model%n_nds))
	allocate(el_data(13,model%n_els))

	allocate(nd_data_name(2))
	allocate(el_data_name(13))

	nd_data_name = (/"ux","uy"/)
	el_data_name(1) = "epsx"
	el_data_name(2) = "epsy"
	el_data_name(3) = "epsxy"
	el_data_name(4) = 'sigx,'
	el_data_name(5) = 'sigy,'
	el_data_name(6) = 'sigz,'
	el_data_name(7) = 'sigyz,'
	el_data_name(8) = 'sigzx,'
	el_data_name(9) = 'sigxy,'
	el_data_name(10) = 'msig,'
	el_data_name(11) = 'maxpsig,'
	el_data_name(12) = 'mpsigx,'
	el_data_name(13) = 'mpsigy,'

	do i=1,model%n_nds
		nd_data(1,i) = output%u(i*2-1)
		nd_data(2,i) = output%u(i*2)
	end do

	el_data(1:3,:) = output%eps(:,:)
	el_data(4:9,:) = output%sig(:,:)
	el_data(10,:) = output%msig(:)
	el_data(11,:) = output%maxpsig(:)
	el_data(12,:) = output%maxpsigx(:)
	el_data(13,:) = output%maxpsigy(:)

 	else if (model%dim == 3) then

 	allocate(nd_data(3,model%n_nds))
	allocate(el_data(16,model%n_els))

	allocate(nd_data_name(3))
	allocate(el_data_name(16))

	nd_data_name = (/"ux","uy","uz"/)
	el_data_name(1) = "epsx"
	el_data_name(2) = "epsy"
	el_data_name(3) = "epsz"
	el_data_name(4) = "epsyz"
	el_data_name(5) = "epszx"
	el_data_name(6) = "epsxy"
	el_data_name(7) = 'sigx,'
	el_data_name(8) = 'sigy,'
	el_data_name(9) = 'sigz,'
	el_data_name(10) = 'sigyz,'
	el_data_name(11) = 'sigzx,'
	el_data_name(12) = 'sigxy,'
	el_data_name(13) = 'msig,'
	el_data_name(14) = 'maxpsig,'
	el_data_name(15) = 'mpsigx,'
	el_data_name(16) = 'mpsigy,'

	do i=1,model%n_nds
		nd_data(1,i) = output%u(i*3-2)
		nd_data(2,i) = output%u(i*3-1)
		nd_data(3,i) = output%u(i*3)
	end do

	el_data(1:6,:) = output%eps(:,:)
	el_data(7:12,:) = output%sig(:,:)
	el_data(13,:) = output%msig(:)
	el_data(14,:) = output%maxpsig(:)
	el_data(15,:) = output%maxpsigx(:)
	el_data(16,:) = output%maxpsigy(:)

 	end if

 	call make_inp(nd_data,nd_data_name,el_data,el_data_name,model,file_name)

end subroutine

subroutine calc_sig(output,model,i,D_)
	type(struct_model) :: model
	type(struct_output) :: output
	integer dim,m,state2d,i
	double precision :: D_(:,:)
	integer,pointer,dimension(:,:) :: e_n
	double precision,pointer,dimension(:,:,:) :: materials
	double precision,pointer,dimension(:) :: msig,maxpsig,maxpsigx,maxpsigy
	double precision,allocatable,dimension(:,:) :: D3,C,C2,D2
	double precision,pointer,dimension(:,:) :: eps,sig

	double precision psig1,psig2,the


	dim = model%dim
	m = model%n_nds_1el

	e_n => model%elements
	materials => model%materials
	state2d = model%state2d

	eps => output%eps
	sig => output%sig
	msig => output%msig
	maxpsig => output%maxpsig
	maxpsigx => output%maxpsigx
	maxpsigy => output%maxpsigy

	allocate(D2(3,3))
	allocate(D3(6,6))
	allocate(C(6,6))
	allocate(C2(3,3))

	if (dim == 2) then
		if (state2d == 0) then
! 			C = inv_jordan(rot_D(materials(:,:,i),angle(i)))
! 			C = inv_jordan(model%materials(:,:,model%material_nos(i)))
			C = inv_jordan(D_)
			C2(1:2,1:2) = C(1:2,1:2)
			C2(3,1:2) = C(6,1:2)
			C2(1:2,3) = C(1:2,6)
			C2(3,3) = C(6,6)
			D2 = inv_jordan(C2)


			sig(1,i) = D2(1,1)*eps(1,i) + D2(1,2)*eps(2,i) + D2(1,3)*eps(3,i)
			sig(2,i) = D2(2,1)*eps(1,i) + D2(2,2)*eps(2,i) + D2(2,3)*eps(3,i)
			sig(6,i) = D2(3,1)*eps(1,i) + D2(3,2)*eps(2,i) + D2(3,3)*eps(3,i)

			sig(3,i) = 0
			sig(4,i) = 0
			sig(5,i) = 0
		else if (state2d > 0) then
! 			D3 = rot_D(materials(:,:,i),angle(i))
				D3 = D_
				D2(1:2,1:2) = D3(1:2,1:2)
				D2(3,1:2) = D3(6,1:2)
				D2(1:2,3) = D3(1:2,6)
				D2(3,3) = D3(6,6)

				sig(1,i) = D3(1,1)*eps(1,i) + D3(1,2)*eps(2,i) + D3(1,6)*eps(3,i)
				sig(2,i) = D3(2,1)*eps(1,i) + D3(2,2)*eps(2,i) + D3(2,6)*eps(3,i)
				sig(3,i) = D3(3,1)*eps(1,i) + D3(3,2)*eps(2,i) + D3(3,6)*eps(3,i)
				sig(4,i) = D3(4,1)*eps(1,i) + D3(4,2)*eps(2,i) + D3(4,6)*eps(3,i)
				sig(5,i) = D3(5,1)*eps(1,i) + D3(5,2)*eps(2,i) + D3(5,6)*eps(3,i)
				sig(6,i) = D3(6,1)*eps(1,i) + D3(6,2)*eps(2,i) + D3(6,6)*eps(3,i)
		end if
	else if (dim == 3) then
		sig(:,i) = matmul(D_,eps(:,i))
	end if


	msig(i) = (sig(1,i)-sig(2,i))**2+(sig(2,i)-sig(3,i))**2+(sig(3,i)-sig(1,i))**2
	msig(i) = msig(i) + 6*(sig(4,i)**2+sig(5,i)**2+sig(6,i)**2)
	msig(i) = dsqrt(msig(i)/2)



	the = 1d0/2*datan(2d0*sig(6,i)/(sig(1,i)-sig(2,i)))
	psig1 = sig(1,i)*dcos(the)**2+sig(2,i)*dsin(the)**2+2*sig(6,i)*sin(the)*cos(the)
	psig2 = sig(1,i)*dcos(the)**2+sig(2,i)*dsin(the)**2-2*sig(6,i)*sin(the)*cos(the)

	if (psig2>psig1) then
		the = the + dacos(0d0)
		maxpsig(i) = psig2
	else
		maxpsig(i) = psig1
	end if

	maxpsigx(i) = dcos(the)
	maxpsigy(i) = dsin(the)
end subroutine


subroutine make_inp(nd_data,nd_data_name,el_data,el_data_name,model,file_name)
	type(struct_model) :: model
	integer i
	double precision,dimension(:,:) :: nd_data,el_data
	character(16),dimension(:) :: nd_data_name,el_data_name
	double precision,pointer,dimension(:,:) :: nodes
	integer,pointer,dimension(:,:) :: e_n
	integer,pointer,dimension(:) :: material_nos

	integer nn,ne,dim,n_nds_1el,n_nd_data,n_el_data
	integer,allocatable,dimension(:) :: ones_nd,ones_el
	character*16 type_el_inp
	character(*) file_name

	nn = model%n_nds
	ne = model%n_els
	dim = model%dim
	n_nds_1el = model%n_nds_1el

	n_nd_data = ubound(nd_data,1)
	n_el_data = ubound(el_data,1)

	allocate(ones_nd(n_nd_data))
	allocate(ones_el(n_el_data))
	ones_nd = 1
	ones_el = 1

	nodes => model%nodes
	e_n => model%elements
	material_nos => model%material_nos

	if (trim(model%type_el) == "quad4") then
		type_el_inp = "quad"
	else if (trim(model%type_el) == "hexa8") then
		type_el_inp = "hex"
	end if

	open(10, file=trim(path_model)//trim(model%name)//slash//trim(file_name)//'.inp')

	write (10,*) '1'
	write (10,*) 'geom'
	write (10,*) 'step1'
	write (10,*) nn, ne

	do i=1,nn
		if (dim == 2) then
			write (10,*) i, nodes(:,i), 0
		else if (dim == 3) then
			write (10,*) i, nodes(:,i)
		end if
	end do
	do i=1,ne
		write (10,*) i, material_nos(i), trim(type_el_inp), e_n(1:n_nds_1el,i)
	end do

	write (10,*) n_nd_data, n_el_data
	write (10,*) n_nd_data, ones_nd
	do i=1,n_nd_data
		write (10,*) trim(nd_data_name(i))//","
	end do
	do i=1,nn
		write (10,*) i, nd_data(:,i)
	end do

	write (10,*) n_el_data, ones_el
	do i=1,n_el_data
		write (10,*) trim(el_data_name(i))//","
	end do

	do i=1,ne
		write (10,*) i, el_data(:,i)
	end do

	close(10)


end subroutine

! subroutine visualize_u(model,u)
! 	type(struct_model) :: model
! 	real(8) :: u(:)
! 	real(8),allocatable :: umax(:),xmin(:),xmax(:),scale
! 	integer dim,n_nds,i,j
! 	real(8),pointer :: nodes(:,:)

! 	dim = model%dim
! 	n_nds = model%n_nds
! 	nodes => model%nodes

! 	allocate(umax(dim),xmin(dim),xmax(dim))

! 	umax = 0d0
! 	do i=1,n_nds
! 		do j=1,dim
! 			if (abs(u((i-1)*dim+j)) > umax(j)) umax(j) = abs(u((i-1)*dim+j))
! 		end do
! 	end do

! 	xmin = 1d30
! 	xmax = 0d0
! 	do i=1,n_nds
! 		do j=1,dim
! 			if (nodes(j,i) < xmin(j)) xmin(j) = nodes(j,i)
! 			if (nodes(j,i) > xmax(j)) xmax(j) = nodes(j,i)
! 		end do
! 	end do

! 	scale = 0d0
! 	do i=1,dim
! 		if (scale < (xmax(i)-xmin(i))/umax(i)/15d0) scale = (xmax(i)-xmin(i))/umax(i)/15d0
! 	end do

! 	do i=1,n_nds
! 		do j=1,dim
! 			nodes(j,i) = nodes(j,i) + u((i-1)*dim+j)*scale
! 		end do
! 	end do


! end subroutine

subroutine visualize_u(model,u)
	type(struct_model) :: model
	real(8) :: u(:)


	integer dim,n_nds,i,j,k,l,n_els,n_nds_1el
	real(8),pointer :: nodes(:,:)
	integer,pointer :: elements(:,:)
	real(8) :: du,dx,scale

	dim = model%dim
	n_nds = model%n_nds
	n_els = model%n_els
	n_nds_1el = model%n_nds_1el
	nodes => model%nodes
	elements => model%elements

	scale = 1d30

	do i=1,n_els
		do j=1,n_nds_1el
		do k=j+1,n_nds_1el
			du = 0d0
			dx = 0d0
			do l=1,dim
				du = du + ( u((elements(j,i)-1)*dim+l) - u((elements(k,i)-1)*dim+l) )**2
				dx = dx + ( nodes(l,elements(j,i)) - nodes(l,elements(k,i)) )**2
			end do
			du = dsqrt(du)
			dx = dsqrt(dx)
			if ((du /= 0d0) .and. (dx/du*0.9<scale)) then
				scale = dx/du*0.9
			end if
		end do
		end do
	end do

	do i=1,n_nds
		do j=1,dim
			nodes(j,i) = nodes(j,i) + u((i-1)*dim+j)*scale
		end do
	end do

end subroutine

function Nh8(xi)
	implicit none
	integer n,m,i,j
	double precision,dimension(:) :: xi
	double precision,allocatable,dimension(:,:) :: Nh8
	allocate(Nh8(3,8))

		Nh8(1,1) = -1./8.*(1-xi(2))*(1-xi(3))
		Nh8(1,2) = 1./8.*(1-xi(2))*(1-xi(3))
		Nh8(1,3) = 1./8.*(1+xi(2))*(1-xi(3))
		Nh8(1,4) = -1./8.*(1+xi(2))*(1-xi(3))
		Nh8(1,5) = -1./8.*(1-xi(2))*(1+xi(3))
		Nh8(1,6) = 1./8.*(1-xi(2))*(1+xi(3))
		Nh8(1,7) = 1./8.*(1+xi(2))*(1+xi(3))
		Nh8(1,8) = -1./8.*(1+xi(2))*(1+xi(3))
		Nh8(2,1) = -1./8.*(1-xi(3))*(1-xi(1))
		Nh8(2,2) = -1./8.*(1-xi(3))*(1+xi(1))
		Nh8(2,3) = 1./8.*(1-xi(3))*(1+xi(1))
		Nh8(2,4) = 1./8.*(1-xi(3))*(1-xi(1))
		Nh8(2,5) = -1./8.*(1+xi(3))*(1-xi(1))
		Nh8(2,6) = -1./8.*(1+xi(3))*(1+xi(1))
		Nh8(2,7) = 1./8.*(1+xi(3))*(1+xi(1))
		Nh8(2,8) = 1./8.*(1+xi(3))*(1-xi(1))
		Nh8(3,1) = -1./8.*(1-xi(1))*(1-xi(2))
		Nh8(3,2) = -1./8.*(1+xi(1))*(1-xi(2))
		Nh8(3,3) = -1./8.*(1+xi(1))*(1+xi(2))
		Nh8(3,4) = -1./8.*(1-xi(1))*(1+xi(2))
		Nh8(3,5) = 1./8.*(1-xi(1))*(1-xi(2))
		Nh8(3,6) = 1./8.*(1+xi(1))*(1-xi(2))
		Nh8(3,7) = 1./8.*(1+xi(1))*(1+xi(2))
		Nh8(3,8) = 1./8.*(1-xi(1))*(1+xi(2))

end function

function Nq4(xi)
	implicit none
	integer n,m,i,j
	double precision,dimension(:) :: xi
	double precision,allocatable,dimension(:,:) :: Nq4
	allocate(Nq4(2,4))

		Nq4(1,1) = -1./4.*(1-xi(2))
		Nq4(1,2) = 1./4.*(1-xi(2))
		Nq4(1,3) = 1./4.*(1+xi(2))
		Nq4(1,4) = -1./4.*(1+xi(2))
		Nq4(2,1) = -1./4.*(1-xi(1))
		Nq4(2,2) = -1./4.*(1+xi(1))
		Nq4(2,3) = 1./4.*(1+xi(1))
		Nq4(2,4) = 1./4.*(1-xi(1))
end function

function sl_LDL(L,b,sln)
	implicit none
	integer i,n,m,cnt,j,k,s,t
	integer,dimension(:) :: sln
	double precision lld
	double precision,dimension(:) :: b,L
	double precision,allocatable,dimension(:) :: sl_LDL,D

	n = ubound(sln,1)

	allocate(sl_LDL(n))
	allocate(D(n))

	d = 0d0
	d(1) = L(1)

	do i=2,n
		if ((i-sln(i)+sln(i-1)+1)==1) l(sln(i-1)+1) = L(sln(i-1)+1)/L(1)
	end do

	do i=2,n
	if ( (i-sln(i)+sln(i-1)+1)>2 ) then
		s = i-sln(i)+sln(i-1)+1
	else
		s = 2
	end if

	if (s<i) then
		do j=s,i-1
		if ( (i-sln(i)+sln(i-1)+1)>(j-sln(j)+sln(j-1)+1) ) then
			t = i-sln(i)+sln(i-1)+1
		else
			t = j-sln(j)+sln(j-1)+1
		end if

		lld = 0
		if (t<j) then
			do k=t,j-1
			lld = lld + l(sln(i)+k-i)*l(sln(j)+k-j)*d(k)
			end do
		end if

		l(sln(i)+j-i) = (L(sln(i)+j-i)-lld)/d(j)
		end do
	end if

	t = i-sln(i)+sln(i-1)+1
	lld = 0
	if (t<i) then
		do k=t,i-1
			lld = lld + l(sln(i)+k-i)**2*d(k)
		end do
	end if
	d(i) = L(sln(i)) - lld

	end do


	do i=1,n
		L(sln(i)) = 1d0
	end do

	sl_LDL = culc_LDL(L,D,b,sln)

end function

function culc_LDL(L,D,b,sln)
	integer i,j,n,m,k,s
	double precision lyd,lx
	integer,dimension(:) :: sln
	integer,allocatable,dimension(:) :: rowh
	double precision,dimension(:) :: b,D,L
	double precision,allocatable,dimension(:) :: culc_LDL,x,y

	n = ubound(sln,1)

	allocate(rowh(n))

	! LTの各行の列数を計算

	rowh = 1
	do j=2,n
		s = j-sln(j)+sln(j-1)+1

		do i=s,j
			if (dabs(L(sln(j)+i-j))>1d-30) rowh(i) = j-i+1
		end do
	end do

	allocate(x(n))
	allocate(y(n))
	allocate(culc_LDL(n))


	y = 0d0
	y(1) = b(1)/d(1)
	do i=2,n
		lyd = 0d0

		s = i-sln(i)+sln(i-1)+1
		do k=s,i-1
			lyd = lyd + l(sln(i)+k-i)*y(k)*d(k)
		end do
		y(i) = (b(i) - lyd)/d(i)
	end do

	x = 0d0
	x(n) = y(n)
	do i=n-1,1,-1
		lx = 0d0

		s = rowh(i)+i-1
		do k=i+1,s
			if (sln(k)-sln(k-1)>k-i) lx = lx + l(sln(k)+i-k)*x(k)
		end do
		x(i) = y(i) - lx
	end do

	culc_LDL = x

end function

! function iccg_solver(L,b,index_diag)

! end function

function inv_jordan(A0)
	implicit none
	integer n,i,j,k
	double precision,allocatable,dimension(:,:) :: L,inv_jordan,A
	double precision,dimension(:,:) :: A0

	n = ubound(A0,1)

	allocate(inv_jordan(n,n))
	allocate(A(n,n))

	A = A0

	inv_jordan = 0d0;
	do i=1,n
		inv_jordan(i,i) = 1d0;
	end do

	do k=1,n
		do i=1,n
		if (k /= i) then
			do j=1,k
				inv_jordan(i,j) = inv_jordan(i,j) - a(i,k)*inv_jordan(k,j)/a(k,k);
			end do
			do j=k+1,n
				a(i,j) = a(i,j) - a(i,k)*a(k,j)/a(k,k);
			end do
			a(i,k) = 0d0;
		end if
		end do
	end do

	do i=1,n
	do j=1,n
		inv_jordan(i,j) = inv_jordan(i,j)/a(i,i)
	end do
	end do

end function

end module fem_module

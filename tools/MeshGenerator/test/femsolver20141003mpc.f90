program fem
	implicit none
	integer dim,l,m,n,q,r,s,t,ne,i,cnt,nn,stat,fn,ncn,nk,gn,slm,Ki,Kj
	integer hole_num,step,tmp_n,tmp_mat,n_mat,auto_rot,auto_rot_num
	integer nmpc
	integer state2d
	integer hoffman_flag
	double precision a,detJ,E,v,maxu,h,rev
	double precision theta,pi
	integer,allocatable,dimension(:) :: tmp_e_n,sln,material
	integer,allocatable,dimension(:) :: nmpc1,nmpc2
	integer,allocatable,dimension(:,:) :: e_n,cn
	integer,allocatable,dimension(:,:) :: mpcdir
	double precision,allocatable,dimension(:) :: w,K,f,u,ue,msig,meps,tmp_nodes,tmp_f,tmp_d,hoffman_f
	double precision,allocatable,dimension(:,:) :: J,Jinv,B,BDB,xi_g,Nxi,Nx,nodes,Ke,eps,sig,phi,lam,cnval,D2,C,C2,d3,hoffman_prop
	double precision,allocatable,dimension(:,:,:) :: D
	character fpath*32,model*16

	double precision the,maxpsig
	double precision,dimension(2) :: psig

! 2014/10/24
!	・平面ひずみへの変更
!	・Hoffmanの破壊則の計算の追加
!

! 2014/10/03
!
!	・主方向ベクトルの追加
!

! 2014/09/24
!
!	・要素ごとに材料を定義できるように改良。
!	・要素形状による材料方向の判定機能を追加。
!


	dim = 2 ! 次元

	! 節点の数
	if (dim == 3) then
		m = 8 
	else if (dim == 2) then
		m = 4
	end if

	! モデル名
	print *, "Model:"
	read *, model
	fpath = trim(model)//"/"

	! 厚さの設定

	print *, "thickness:"
	read *, h

	! 平面応力/平面ひずみ

	print *, "state: (0: plane stress,  >0: plane strain)"
	read *, state2d

	! 回転の自動/手動

	print *, "Auto rotation of D: (>0 :auto,  <=0 :no)"

    read *, auto_rot

	! Dの回転角

	if (auto_rot <= 0) then
		print *, "rotation of D:"
		read *, theta
	end if

	! 回転を与える材料番号

	if (auto_rot > 0) then
		print *, "material No. which rotation of D is adopted:"
		read *, auto_rot_num
	end if

	! hoffman則の計算

	print *,"calculation of hoffman: (0: no,  >0:yes)"
	read *,hoffman_flag


print *,"Reading files ..."


! モデルの読み込み
	open(10,file=trim(fpath)//'model.csv',status='old')

	! 節点数,要素数の読み込み

	read (10,*) nn, ne

	allocate(e_n(m,ne))
	allocate(nodes(dim,nn))
	allocate(tmp_nodes(dim))
	allocate(tmp_e_n(m))
	allocate(material(ne))

	! 節点の読み込み
	do i=1,nn
		read (10,*) tmp_n,tmp_nodes(:)
		nodes(:,tmp_n) = tmp_nodes(:)
	end do

	! 要素の読み込み
	do i=1,ne
		read (10,*) tmp_n,tmp_mat,tmp_e_n(:)
		e_n(:,tmp_n) = tmp_e_n(:)
		material(tmp_n) = tmp_mat
	end do
	close(10)

	! 材料数の読み込み

	open(10,file=trim(fpath)//'mats.csv',status='old')
	read (10,*) n_mat

	allocate(D(n_mat,6,6))

	do q=1,n_mat
		read (10,*) tmp_mat
		do r=1,6
			read (10,*) D(tmp_mat,:,r)
		end do
	end do

	close(10)


	! 損傷に関するプロパティの読み込み

	open(10,file=trim(fpath)//'hoffman.csv',status='old')

	allocate(hoffman_prop(n_mat,9))

	do q=1,n_mat
		read (10,*) tmp_mat
		do r=1,3
			 read (10,*) hoffman_prop(tmp_mat,(r*3-2):(r*3))
		end do
	end do

print *, hoffman_prop(2,:)


	! 円周率の定義
	pi = dacos(-1d0);

	! 回転角のラジアルへの変換
	theta = theta*pi/180d0


! 変数の初期化

	allocate(xi_g(dim,m))

	allocate(w(m))
	allocate(J(dim,dim))
	allocate(Jinv(dim,dim))
	allocate(Nxi(dim,m))
	allocate(Nx(dim,m))

	allocate(B(dim*(dim+1)/2,dim*m))
	allocate(D3(6,6))
	allocate(D2(3,3))
	allocate(C(6,6))
	allocate(C2(3,3))
	allocate(BDB(dim*m,dim*m))
	allocate(Ke(dim*m,dim*m))
	allocate(tmp_f(dim))
	allocate(f(dim*nn))
	allocate(u(dim*nn))

	allocate(ue(dim*m))
	allocate(eps(dim*(dim+1)/2,ne))
!	allocate(sig(dim*(dim+1)/2,ne))
allocate(sig(6,ne))
	allocate(msig(ne))
	allocate(meps(ne))
	allocate(phi(3,ne))
	allocate(lam(3,ne))
	allocate(hoffman_f(ne))
	

!if (dim == 3) then
!	D = 0d0
!
!
!	D(1,1,:) = (/1-v, v, v, 0d0, 0d0, 0d0/)
!	D(1,2,:) = (/v, 1-v, v, 0d0, 0d0, 0d0/)
!	D(1,3,:) = (/v, v, 1-v, 0d0, 0d0, 0d0/)
!	D(1,4,4) = (1-2*v)/2
!	D(1,5,5) = (1-2*v)/2
!	D(1,6,6) = (1-2*v)/2
!	D(1,:,:) = D(1,:,:)*E / ((1+v)*(1-2*v))
!
!
!E = E*2d0
!
!	D(2,1,:) = (/1-v, v, v, 0d0, 0d0, 0d0/)
!	D(2,2,:) = (/v, 1-v, v, 0d0, 0d0, 0d0/)
!	D(2,3,:) = (/v, v, 1-v, 0d0, 0d0, 0d0/)
!	D(2,4,4) = (1-2*v)/2
!	D(2,5,5) = (1-2*v)/2
!	D(2,6,6) = (1-2*v)/2
!	D(2,:,:) = D(2,:,:)*E / ((1+v)*(1-2*v))

!end if

!else if (dim == 2) then
!	D = 0d0
!	D(1,1,:) = (/1d0, v, 0d0/)
!	D(1,2,:) = (/v, 1d0, 0d0/)
!	D(1,3,3) = (1-v)/2
!
!	D = D*E / (1-v**2)
!end if



	! ガウス点の座標の絶対値
	a = dsqrt(1d0/3d0)

if (dim == 2) then
	! ガウス点の数
	gn = 4

	! ガウス点の座標
	xi_g(1,:) = (/-a, a, a, -a/)
	xi_g(2,:) = (/-a, -a, a, a/)

	! ガウス点の重み
	w = (/1., 1., 1., 1./)

else if (dim == 3) then
	! ガウス点の数
	gn = 8

	! ガウス点の座標
	xi_g(1,:) = (/-a, a, a, -a, -a, a, a, -a/)
	xi_g(2,:) = (/-a, -a, a, a, -a, -a, a, a/)
	xi_g(3,:) = (/-a, -a, -a, -a, a, a, a, a/)

	! ガウス点の重み
	w = (/1., 1., 1., 1., 1., 1., 1., 1./)
end if



	! 等価節点力ベクトルの読み込み
	f = 0d0
	open(10,file=trim(fpath)//'f.csv',status='old')
	read(10,*) fn

	do i=1,fn
		read (10,*) tmp_n,tmp_f(:)

		f((tmp_n-1)*dim+1:tmp_n*dim) = tmp_f(:)
	end do
	close(10)


	! 単点拘束条件の読み込み

	open(10,file=trim(fpath)//'cn.csv',status='old')
	read (10,*) ncn

	allocate(cn(dim+1,ncn))
	allocate(cnval(dim,ncn))
	cn = 0d0

	do i=1,ncn
		read (10,*) cn(1:(dim+1),i), cnval(1:dim,i)
	end do
	close(10)


	! 多点拘束条件の読み込み

	open(10,file=trim(fpath)//'mpc.csv',status='old')
	read (10,*) nmpc

if (nmpc /= 0) then
	allocate(mpcdir(dim,nmpc))
	allocate(nmpc1(nmpc))
	allocate(nmpc2(nmpc))
	mpcdir = 0d0
	do i=1,nmpc
		read (10,*) nmpc1(i),nmpc2(i),mpcdir(1:dim,i)
	end do
end if

	close(10)



	! スカイラインマトリックスの格納の準備
	allocate(sln(nn*dim))
	sln = 0
	do i=1,ne

		slm = nn*dim

		do s=1,m
		do t=1,dim
			if (slm> ((e_n(s,i)-1)*dim+t) ) slm = (e_n(s,i)-1)*dim+t
		end do
		end do
		do s=1,m
		do t=1,dim
			if (sln((e_n(s,i)-1)*dim+t) < ((e_n(s,i)-1)*dim+t-slm+1)) then
				sln((e_n(s,i)-1)*dim+t) = (e_n(s,i)-1)*dim+t-slm + 1
			end if
		end do
		end do

	end do

	nk = sum(sln)

	sln(1) = 1
	do i=2,nn*dim
		sln(i) = sln(i-1)+sln(i)
	end do


	allocate(K(nk))

J=0.

K=0d0


! 各要素のループ
	do i=1,ne


	! Dの回転角の計算


	if (auto_rot > 0) then
	if (material(i) == auto_rot_num) then
		theta = atan2(nodes(2,e_n(2,i))+nodes(2,e_n(3,i))-nodes(2,e_n(1,i))-nodes(2,e_n(4,i)) , &
							&nodes(1,e_n(2,i))+nodes(1,e_n(3,i))-nodes(1,e_n(1,i))-nodes(1,e_n(4,i)) );
	else
		theta = 0d0
	end if
	end if
!	print *,theta/pi*180;


! ガウス積分のループ
	Ke=0d0

	do l=1,gn
		detJ = 0.
		J = 0.
		Jinv = 0.
		B = 0.
		BDB = 0.
		Nxi = 0.
		Nx = 0.

		! dN/dxiとdN/detaの設定

		if (dim == 2) then
			Nxi = Nq4(xi_g(:,l))
		else if (dim == 3) then
			Nxi = Nh8(xi_g(:,l))
		end if

		! ヤコビ行列の設定
		do s=1,dim
		do t=1,dim
			do q=1,m
				J(s,t) = J(s,t) + Nxi(s,q) * nodes(t,e_n(q,i))
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


		Nx = matmul(Jinv, Nxi)

		if (dim == 3) then
		do q=1,m
			B(:, dim*(q-1)+1) = (/ Nx(1,q), 0d0, 0d0, 0d0, Nx(3,q), Nx(2,q) /)
			B(:, dim*(q-1)+2) = (/ 0d0, Nx(2,q), 0d0, Nx(3,q), 0d0, Nx(1,q) /)
			B(:, dim*(q-1)+3) = (/ 0d0, 0d0, Nx(3,q), Nx(2,q), Nx(1,q), 0d0 /)

		end do
		else if (dim == 2) then
		do q=1,m
			B(:, dim*(q-1)+1) = (/ Nx(1,q), 0d0, Nx(2,q) /)
			B(:, dim*(q-1)+2) = (/ 0d0, Nx(2,q), Nx(1,q) /)
		end do
		end if


		! Dマトリックスの2D平面応力への変換

		if (dim == 2) then
			if (state2d == 0) then
				C = inv_jordan(rot_D(D(material(i),:,:),theta))
				C2(1:2,1:2) = C(1:2,1:2)
				C2(3,1:2) = C(6,1:2)
				C2(1:2,3) = C(1:2,6)
				C2(3,3) = C(6,6)
				D2 = inv_jordan(C2)
			else if (state2d > 0) then
				D3 = rot_D(D(material(i),:,:),theta)
				D2(1:2,1:2) = D3(1:2,1:2)
				D2(3,1:2) = D3(6,1:2)
				D2(1:2,3) = D3(1:2,6)
				D2(3,3) = D3(6,6)
			end if

			BDB = matmul( matmul(transpose(B),D2),B)
		else if (dim == 3) then
			BDB = matmul( matmul(transpose(B),D(material(i),:,:)),B)
		end if


		Ke = Ke + w(l)*BDB*detJ*h

	end do

	do s=1,m
	do t=1,m
		do q=1,dim
		do r=1,dim
			Ki = (e_n(s,i)-1)*dim+q
			kj = (e_n(t,i)-1)*dim+r
			if (Ki <= Kj) then
				K(sln(Kj)+Ki-Kj) = K(sln(Kj)+Ki-Kj) + Ke((s-1)*dim+q,(t-1)*dim+r)
			end if

		end do
		end do
	end do
	end do


	end do

print *, "Calculating K was completed."

	! 多点拘束条件の設定

	do i=1,nmpc
	do l=1,dim
	if (mpcdir(l,i) /= 0) then
		s = (nmpc1(i)-1)*dim+l
		t = (nmpc2(i)-1)*dim+l

		if (s < t) then
			do q=1,s
			if (s > 1) then
				if ((sln(s)+q-s>sln(s-1)) .and. (sln(t)+q-t>sln(t-1))) then
					K(sln(s)+q-s) = K(sln(s)+q-s)+K(sln(t)+q-t)
				end if
			else
				if (sln(t)+q-t>sln(t-1)) then
					K(1) = K(1)+K(sln(t)+q-t)
				end if
			end if
			end do

			do q=(s+1),t
				if ((sln(q)+s-q>sln(q-1)) .and. (sln(t)+q-t>sln(t-1))) then
					K(sln(q)+s-q) = K(sln(q)+s-q)+K(sln(t)+q-t)
				end if
			end do

			do q=(t+1),nn*dim
				if ((sln(q)+s-q>sln(q-1)) .and. (sln(q)+t-q>sln(q-1))) then
					K(sln(q)+s-q) = K(sln(q)+s-q)+K(sln(q)+t-q)
				end if
			end do
		else if (s > t) then
			do q=1,t
			if (t > 1) then
				if ((sln(s)+q-s>sln(s-1)) .and. (sln(t)+q-t>sln(t-1))) then
					K(sln(s)+q-s) = K(sln(s)+q-s)+K(sln(t)+q-t)
				end if
			else
				if (sln(s)+q-s>sln(s-1)) then
					K(sln(s)+q-s) = K(sln(s)+q-s)+K(1)
				end if
			end if
			end do

			do q=(t+1),s
				if ((sln(s)+q-s>sln(s-1)) .and. (sln(q)+t-q>sln(q-1))) then
					K(sln(s)+q-s) = K(sln(s)+q-s)+K(sln(q)+t-q)
				end if
			end do
			do q=(s+1),nn*dim
				if ((sln(q)+s-q>sln(q-1)) .and. (sln(q)+t-q>sln(q-1))) then
					K(sln(q)+s-q) = K(sln(q)+s-q)+K(sln(q)+t-q)
				end if
			end do
		end if
		K(sln(t)) = 1d30
	end if
	end do
	end do



	! 拘束条件の設定

	do i=1,ncn
	do l=1,dim
		if (cn(l+1,i) /= 0) then

			n = (cn(1,i)-1)*dim+l

			f(1) = f(1) - K(sln(n)+1-n)*cnval(l,i)
			do q=2,n-1
				if ((sln(n)+q-n)>sln(n-1)) then
					f(q) = f(q) - K(sln(n)+q-n)*cnval(l,i)
				end if
			end do
			f(n) = - K(sln(n))*cnval(l,i)
			do q=n+1,nn*dim
				if ((sln(q)+n-q)>sln(q-1)) then
					f(q) = f(q) - K(sln(q)+n-q)*cnval(l,i)
				end if
			end do

			K(sln(n)) = 1d30

		end if
	end do
	end do


print *, "Boundary conditions added."

!
! write on csv file
!	open(9, file=trim(fpath)//'K.csv')
!		write (9,"(f30.15, 1000000(',', f40.5))") K(sln(1))
!	do i=2,nn*dim
!		write (9,"(f30.15, 1000000(',', f40.5))") K(sln(i-1)+1:sln(i))
!	end do
!	close(9)
!	open(9, file=trim(fpath)//'fout.csv')
!		write (9,"(f30.15, 1000000(',', f40.5))") f(:)
!	close(9)



	! LDL分解によるuの計算

	u = 0d0

	u = sl_LDL(K,f,sln)


print *, "Displacement u was calculated."

	! 拘束点に変位を代入し直す

	do i=1,ncn
	do l=1,dim
		n=(cn(1,i)-1)*dim+l
		if (cn(l+1,i) /= 0) then
			u(n) = cnval(l,i)
		end if
	end do
	end do

	do i=1,nmpc
	do l=1,dim
		if (mpcdir(l,i) /= 0) then
			u( (nmpc2(i)-1)*dim+l ) = u( (nmpc1(i)-1)*dim+l )
		end if
	end do
	end do



! ひずみ、ミーゼス応力の計算

	if (dim == 2) then
		Nxi = Nq4((/0d0, 0d0/))
	else if (dim == 3) then
		Nxi = Nh8((/0d0, 0d0, 0d0/))
	end if

	open(10, file=trim(fpath)//'rot.csv')

	do i=1,ne

		! Dの回転角の計算

		if (auto_rot > 0) then
		if (material(i) == auto_rot_num) then
			theta = atan2(nodes(2,e_n(2,i))+nodes(2,e_n(3,i))-nodes(2,e_n(1,i))-nodes(2,e_n(4,i)) , &
								&nodes(1,e_n(2,i))+nodes(1,e_n(3,i))-nodes(1,e_n(1,i))-nodes(1,e_n(4,i)) );

		else 
			theta = 0d0
		end if
		end if

		write(10,*) i,',',theta/pi*180d0

		detJ = 0.
		J = 0.
		Jinv = 0.
		B = 0.
		Nx = 0.

		ue = 0d0
		do l=1,m
			ue( (l-1)*dim+1:l*dim ) = u( (e_n(l,i)-1)*dim+1:e_n(l,i)*dim )
		end do

		! ヤコビ行列の設定
		do s=1,dim
		do t=1,dim
			do q=1,m
				J(s,t) = J(s,t) + Nxi(s,q) * nodes(t,e_n(q,i))
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


		Nx = matmul(Jinv, Nxi)

		if (dim == 3) then
		do q=1,m
			B(:, dim*(q-1)+1) = (/ Nx(1,q), 0d0, 0d0, 0d0, Nx(3,q), Nx(2,q) /)
			B(:, dim*(q-1)+2) = (/ 0d0, Nx(2,q), 0d0, Nx(3,q), 0d0, Nx(1,q) /)
			B(:, dim*(q-1)+3) = (/ 0d0, 0d0, Nx(3,q), Nx(2,q), Nx(1,q), 0d0 /)
		end do

		else if (dim == 2) then
		do q=1,m
			B(:, dim*(q-1)+1) = (/ Nx(1,q), 0d0, Nx(2,q) /)
			B(:, dim*(q-1)+2) = (/ 0d0, Nx(2,q), Nx(1,q) /)
		end do
		end if

		eps(:,i) = matmul(B,ue)



		! Dマトリックスの2D平面応力への変換

		if (dim == 2) then
			if (state2d == 0) then
				C = inv_jordan(rot_D(D(material(i),:,:),theta))
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
				D3 = rot_D(D(material(i),:,:),theta)
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

!			sig(:,i) = matmul(D2,eps(:,i))
		else if (dim == 3) then
			sig(:,i) = matmul(D(material(i),:,:),eps(:,i))
		end if


		! ミーゼス応力
!		if (dim == 3) then
		msig(i) = (sig(1,i)-sig(2,i))**2+(sig(2,i)-sig(3,i))**2+(sig(3,i)-sig(1,i))**2
		msig(i) = msig(i) + 6*(sig(4,i)**2+sig(5,i)**2+sig(6,i)**2)
		msig(i) = dsqrt(msig(i)/2)
!		else if (dim == 2) then
!		msig(i) = sig(1,i)**2 + sig(2,i)**2 + (sig(1,i)-sig(2,i))**2
!		msig(i) = msig(i) + 6*sig(3,i)**2
!		msig(i) = dsqrt(msig(i)/2)
!		end if

		! ミーゼスひずみ
		if (dim == 3) then
		meps(i) = (eps(1,i)-eps(2,i))**2+(eps(2,i)-eps(3,i))**2+(eps(3,i)-eps(1,i))**2
		meps(i) = meps(i) + 6*(eps(4,i)**2+eps(5,i)**2+eps(6,i)**2)
		meps(i) = dsqrt(meps(i)/2)
		else if (dim == 2) then
			lam(1,i) = 1 + (eps(1,i)+eps(2,i))/2d0 + dsqrt( ((eps(1,i)+eps(2,i))/2d0)**2 - eps(1,i)*eps(2,i) + eps(3,i)**2 )
			lam(2,i) = 1 + (eps(1,i)+eps(2,i))/2d0 - dsqrt( ((eps(1,i)+eps(2,i))/2d0)**2 - eps(1,i)*eps(2,i) + eps(3,i)**2 )
			lam(3,i) = 1/lam(1,i)/lam(2,i)
			meps(i) = dsqrt(2d0/3d0* (dlog(lam(1,i))**2 + dlog(lam(2,i))**2 + dlog(lam(3,i))**2) )
		end if

		! hoffmanの破壊則の計算

		if (hoffman_flag > 0) then
			hoffman_f(i) = calc_hoffman( sig(:,i),hoffman_prop(material(i),:) );


		end if

	end do
	close(10)

	maxu = maxval(u(:))
	cnt = 0
	do
		maxu = maxu * 10
		if (maxu >= 1d0) exit
		cnt = cnt + 1
	end do
	cnt = cnt - 1

!	 変形の表示の調整

	rev = 0
	print *,"adjustment of deformation: (<0: displacement)"
	read *, rev

	open(10, file=trim(fpath)//'data.inp')

	write (10,*) '1'
	write (10,*) 'geom'
	write (10,*) 'step1'
	write (10,*) nn, ne
	do i=1,nn
		if (rev>=0) then
			write (10,*) i, nodes(:,i)+10**cnt * u((i-1)*dim+1:i*dim)*rev, 0
		else
			write (10,*) i, u((i-1)*dim+1:i*dim), 0
		end if
	end do
	do i=1,ne
		write (10,*) i, material(i), 'quad', e_n(:,i)
	end do
	write (10,*) 2, 14
	write (10,*) 2, 1, 1
	write (10,*) 'x,'
	write (10,*) 'y,'
	do i=1,nn
		write (10,*) i, 0, 0
	end do
	write (10,*) 14, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
	write (10,*) 'epsx,'
	write (10,*) 'epsy,'
	write (10,*) 'epsxy,'
	write (10,*) 'sigx,'
	write (10,*) 'sigy,'
	write (10,*) 'sigz,'
	write (10,*) 'sigyz,'
	write (10,*) 'sigzx,'
	write (10,*) 'sigxy,'
	write (10,*) 'msig,'
	write (10,*) 'maxpsig,'
	write (10,*) 'mpsigx,'
	write (10,*) 'mpsigy,'
	write (10,*) 'hoffman_f,'

print *, dacos(0d0)
	do i=1,ne


! 主応力ベクトル方向の計算

	the = 1d0/2*datan(2d0*sig(6,i)/(sig(1,i)-sig(2,i)))
	psig(1) = sig(1,i)*dcos(the)**2+sig(2,i)*dsin(the)**2+2*sig(6,i)*sin(the)*cos(the)
	psig(2) = sig(1,i)*dcos(the)**2+sig(2,i)*dsin(the)**2-2*sig(6,i)*sin(the)*cos(the)

	if (psig(2)>psig(1)) then
		the = the + dacos(0d0)
		maxpsig = psig(2)
	else
		maxpsig = psig(1)
	end if


!		write (10,*) i, eps(:,i), sig(:,i), msig(i), &
!			-( (sig(1,i)+sig(2,i)) + dsqrt((sig(1,i)-sig(2,i))**2+4*sig(3,i)) )/2  *  dsin( 1d0/2*datan(2*sig(3,i)/(sig(1,i)-sig(2,i))) ), &
!			( (sig(1,i)+sig(2,i)) + dsqrt((sig(1,i)-sig(2,i))**2+4*sig(3,i)) )/2  *  dcos( 1d0/2*datan(2*sig(3,i)/(sig(1,i)-sig(2,i))) )

		write (10,*) i, eps(:,i), sig(:,i), msig(i), maxpsig, dcos(the),dsin(the),hoffman_f(i)

	end do

	close(10)

	read *
contains

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

function sl_LDL(L0,b,sln)
	implicit none
	integer i,n,m,cnt,j,k,s,t
	integer,dimension(:) :: sln
	double precision lld
	double precision,dimension(:) :: b,L0
	double precision,allocatable,dimension(:) :: sl_LDL,D,L

	n = ubound(sln,1)

	allocate(sl_LDL(n))
	allocate(D(n))
	allocate(L(n))

	L = L0

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

	! Lの各行の列数を計算

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

function calc_hoffman(sig,strength)
	implicit none
	integer i,j,k
	double precision calc_hoffman
	double precision c1,c2,c3,c4,c5,c6,c7,c8,c9
	double precision,dimension(:) :: sig,strength  ! 強度 (x引張,y引張,z引張,x圧縮,y圧縮,z圧縮,yzせん断,zxせん断,xyせん断)

	c1 = ( 1d0/strength(2)/strength(5) + 1d0/strength(3)/strength(6) - 1d0/strength(1)/strength(4) )/2d0
	c2 = ( 1d0/strength(3)/strength(6) + 1d0/strength(1)/strength(4) - 1d0/strength(2)/strength(5) )/2d0
	c3 = ( 1d0/strength(1)/strength(4) + 1d0/strength(2)/strength(5) - 1d0/strength(3)/strength(6) )/2d0
	c4 = 1d0/strength(1) - 1d0/strength(4)
	c5 = 1d0/strength(2) - 1d0/strength(5)
	c6 = 1d0/strength(3) - 1d0/strength(6)
	c7 = 1d0/strength(7)**2
	c8 = 1d0/strength(8)**2
	c9 = 1d0/strength(9)**2

print *,c1,c2,c3,c4,c5,c6,c7,c8,c9

	calc_hoffman = c1*( (sig(2)-sig(3))**2 ) + c2*( (sig(3)-sig(1))**2 ) + c3*( (sig(1)-sig(2))**2 )&
					+ c4*sig(1) + c5*sig(2) + c6*sig(3) + c7*(sig(4)**2) + c8*(sig(5)**2) + c9*(sig(6)**2)

end function

end program fem
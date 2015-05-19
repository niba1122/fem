module ebe_scg_module
	use fem_module
contains

subroutine ebe_scg(u,model,bc,f,Ke_voxel)
	type(struct_model) :: model
	type(struct_bc) :: bc
	real(8) :: u(:),f(:),a,b
	real(8),allocatable :: s(:),Ke(:,:),r0(:),r1(:),p(:),Ku0(:),kp(:),norm_r0
	real(8),optional :: Ke_voxel(:,:)
	integer i,j,k,l,m,q,r,step
	integer n_nds,n_els,n_nds_1el,dim
	integer n_mpc,n_spc
	integer,allocatable :: mpc(:,:,:),spc(:,:)
	integer,pointer :: elements(:,:)
	integer voxel_flag

	n_nds = model%n_nds
	n_els = model%n_els
	n_nds_1el = model%n_nds_1el
	dim = model%dim
	elements => model%elements

	n_mpc = bc%n_mpc
	n_spc = bc%n_spc

	allocate(s(n_nds*dim))
	allocate(Ke(n_nds_1el*dim,n_nds_1el*dim))
	allocate(r0(n_nds*dim),r1(n_nds*dim))
	allocate(p(n_nds*dim))
	allocate(Ku0(n_nds*dim),kp(n_nds*dim))

!---------------------------------------------------------------------------------------------------
!	initialize mpc & spc conditions
!---------------------------------------------------------------------------------------------------

	!---------------------------------------------------------------------!
	! mpc(dir,slave_id,1) = master id(if 0, no mpc condition)             !
	! mpc(dir,slave_id,2) = master direction(if no mpc condition,maybe 0) !
	!---------------------------------------------------------------------!

	allocate(mpc(dim,n_nds,2))

	allocate(spc(dim,n_nds))

	mpc = 0
	spc = 0

	do i=1,n_mpc
		mpc(bc%mpc_slave_dir(i),bc%mpc_slave_id(i),1) = bc%mpc_master_id(i)
		mpc(bc%mpc_slave_dir(i),bc%mpc_slave_id(i),2) = bc%mpc_master_dir(i)
	end do

	do i=1,n_spc
		spc(1:dim,bc%nodes_spc(1,i)) = bc%nodes_spc(2:dim+1,i)
	end do

! 	do i=1,n_nds
! 		print *,mpc(1,i,:),mpc(2,i,:),mpc(3,i,:)
! 	end do
! 	do i=1,n_nds
! 		print *,spc(:,i)
! 	end do


!---------------------------------------------------------------------------------------------------
!	Calculation of daigonal parts of stiffness matrix for precondition
!---------------------------------------------------------------------------------------------------

	print *,"calculating preconditioner matrix..."
	! extract diagonal parts of stiffness matrix (stored by 's')

	s = 0d0

!$omp parallel private(Ke,i,j,k,l,m,q),default(shared)
!$omp do reduction(+:s)
	do i=1,n_els
		Ke = 0d0
		if (present(Ke_voxel)) then
				Ke = Ke_voxel
		else
			call calc_element_integral(Ke,model,i,calc_BDB)
		end if

		do j=1,n_nds_1el
		do k=1,n_nds_1el
			do l=1,dim
			do m=1,dim

			if ((j==k) .and. (l==m)) then
				s((elements(j,i)-1)*dim+l) = s((elements(j,i)-1)*dim+l) + Ke((j-1)*dim+l,(j-1)*dim+l)

! 				do q=1,bc%n_spc
! 					if ((bc%nodes_spc(l+1,q) /= 0) .and. (bc%nodes_spc(1,q) == elements(j,i))) then
! 						s((elements(j,i)-1)*dim+l) = s((elements(j,i)-1)*dim+l) + 1d30
! 					end if
! 				end do

				if (spc(l,elements(j,i))>0) then
					s((elements(j,i)-1)*dim+l) = s((elements(j,i)-1)*dim+l) + 1d30
				end if

! 				do q=1,bc%n_mpc
! 					if ((bc%mpc_slave_dir(q) == l) .and. (bc%mpc_slave_id(q) == elements(j,i))) then
! 						s((elements(j,i)-1)*dim+l) = s((elements(j,i)-1)*dim+l) + 1d30
! 						s((bc%mpc_master_id(q)-1)*dim+bc%mpc_master_dir(q))&
! 							 = s((bc%mpc_master_id(q)-1)*dim+bc%mpc_master_dir(q)) + Ke((j-1)*dim+l,(j-1)*dim+l)
! 					end if
! 				end do

				if (mpc(l,elements(j,i),1)>0) then
					s((elements(j,i)-1)*dim+l) = s((elements(j,i)-1)*dim+l) + 1d30
					s((mpc(l,elements(j,i),1)-1)*dim+mpc(l,elements(j,i),2))&
							 = s((mpc(l,elements(j,i),1)-1)*dim+mpc(l,elements(j,i),2)) + Ke((j-1)*dim+l,(j-1)*dim+l)
				end if


			else
				if ( (mpc(m,elements(k,i),1) == elements(j,i)) .and. (mpc(m,elements(k,i),2) == l) ) then
					s( (elements(j,i)-1)*dim+l ) = s( (elements(j,i)-1)*dim+l ) + Ke((j-1)*dim+l,(k-1)*dim+m)*2d0
				end if
			end if

			end do
			end do
		end do
		end do
	end do
!$omp end do
!$omp end parallel

	! calculate 's'　to the power of -1/2
	do i=1,n_nds*dim
		s(i) = 1d0/dsqrt(s(i))
	end do


!---------------------------------------------------------------------------------------------------
!	Adding spc & mpc on the right-hand side of the equation
!---------------------------------------------------------------------------------------------------

	print *,"Adding boundary condition on vector f..."

	do i=1,bc%n_mpc
		f((bc%mpc_master_id(i)-1)*dim+bc%mpc_master_dir(i)) = f((bc%mpc_master_id(i)-1)*dim+bc%mpc_master_dir(i)) &
															+ f((bc%mpc_slave_id(i)-1)*dim+bc%mpc_slave_dir(i))
		f((bc%mpc_slave_id(i)-1)*dim+bc%mpc_slave_dir(i)) = 0d0
	end do

 	do i=1,bc%n_spc
	do j=1,dim
		if (bc%nodes_spc(j+1,i) /= 0) then
			f((bc%nodes_spc(1,i)-1)*dim+j) = 0d0
		end if
	end do
	end do

 	u = 0d0

 	do i=1,bc%n_spc
	do j=1,dim
		if (bc%nodes_spc(j+1,i) /= 0) then
			u((bc%nodes_spc(1,i)-1)*dim+j) = bc%disp_spc(j,i)
		end if
	end do
	end do

 	do i=1,bc%n_mpc
	do j=1,dim
	do k=1,dim
		if ((bc%mpc_master_dir(i) == j) .and. (bc%mpc_slave_dir(i) == k)) then
			u((bc%mpc_slave_id(i)-1)*dim+k) = u((bc%mpc_master_id(i)-1)*dim+j)
		end if
	end do
	end do
	end do

	Ku0 = 0d0

!$omp parallel private(Ke,i,j,k,l,m,q),default(shared)
!$omp do reduction(+:Ku0)
	do i=1,n_els
		Ke = 0d0
		if (present(Ke_voxel)) then
			Ke = Ke_voxel
		else
			call calc_element_integral(Ke,model,i,calc_BDB)
		end if

		do j=1,n_nds_1el
		do k=1,n_nds_1el
			do l=1,dim
			do m=1,dim
				Ku0((elements(j,i)-1)*dim+l) = Ku0((elements(j,i)-1)*dim+l)&
					 + Ke((j-1)*dim+l,(k-1)*dim+m)*u((elements(k,i)-1)*dim+m)

			end do
			end do
		end do
		end do
	end do
!$omp end do
!$omp end parallel

	r0 = f - Ku0

!$omp parallel
!$omp do
	do i=1,n_nds*dim
		r0(i) = s(i)*r0(i)
	end do
!$omp end do
!$omp end parallel

! do i=1,n_nds
! 	print *,r0((i-1)*dim:i*dim)
! end do

!---------------------------------------------------------------------------------------------------
!	Calculation of CG Method
!---------------------------------------------------------------------------------------------------

	u = 0d0
	p = r0
	a = 0d0
	b = 0d0

	norm_r0 = dsqrt(dot_product(r0,r0))

	if (norm_r0 < 1d-10) then
		print *,"CAUTION: this calculation cannot be converged because |r| is very small (|r| < 1d-10)."
		print *,"|r| = ", norm_r0
	end if

	print *,"Iterative calculation start"

	do step=1,dim*n_nds
		kp = 0d0

!$omp parallel private(Ke,i,j,k,l,m,q,r),default(shared)
!$omp do reduction(+:Kp)
		do i=1,n_els
			Ke = 0d0
			if (present(Ke_voxel)) then
				Ke = Ke_voxel
			else
				call calc_element_integral(Ke,model,i,calc_BDB)
			end if

			do j=1,n_nds_1el
			do k=1,n_nds_1el

				do l=1,dim
				do m=1,dim
					Kp((elements(j,i)-1)*dim+l) = Kp((elements(j,i)-1)*dim+l)&
					 + s((elements(j,i)-1)*dim+l)*s((elements(k,i)-1)*dim+m) * Ke((j-1)*dim+l,(k-1)*dim+m)*p((elements(k,i)-1)*dim+m)

				if ((j==k).and.(l==m)) then
					if (spc(l,elements(j,i))>0) then
						Kp((elements(j,i)-1)*dim+l) = Kp((elements(j,i)-1)*dim+l) + &
							s((elements(j,i)-1)*dim+l)*s((elements(j,i)-1)*dim+l) * 1d30 * p((elements(j,i)-1)*dim+l)
					end if

					if (mpc(l,elements(j,i),1)>0) then
						Kp((elements(j,i)-1)*dim+l) = Kp((elements(j,i)-1)*dim+l) + &
							s((elements(j,i)-1)*dim+l)*s((elements(j,i)-1)*dim+l) * 1d30 * p((elements(j,i)-1)*dim+l)
					end if
				end if
					if (mpc(m,elements(k,i),1)>0) then
						Kp((elements(j,i)-1)*dim+l) = Kp((elements(j,i)-1)*dim+l)&
						 + s((elements(j,i)-1)*dim+l)*s((mpc(m,elements(k,i),1)-1)*dim+mpc(m,elements(k,i),2)) &
							 * Ke((j-1)*dim+l,(k-1)*dim+m)*p((mpc(m,elements(k,i),1)-1)*dim+mpc(m,elements(k,i),2))
					end if
					if (mpc(l,elements(j,i),1)>0) then
						Kp((mpc(l,elements(j,i),1)-1)*dim+mpc(l,elements(j,i),2)) = Kp((mpc(l,elements(j,i),1)-1)*dim+mpc(l,elements(j,i),2))&
						 + s((mpc(l,elements(j,i),1)-1)*dim+mpc(l,elements(j,i),2))*s((elements(k,i)-1)*dim+m)&
							 * Ke((j-1)*dim+l,(k-1)*dim+m)*p((elements(k,i)-1)*dim+m)
					end if

					if ( (mpc(l,elements(j,i),1)>0) .and. (mpc(m,elements(k,i),1)>0) ) then
						Kp((mpc(l,elements(j,i),1)-1)*dim+mpc(l,elements(j,i),2)) = Kp((mpc(l,elements(j,i),1)-1)*dim+mpc(l,elements(j,i),2)) + &
							s((mpc(l,elements(j,i),1)-1)*dim+mpc(l,elements(j,i),2))*s((mpc(m,elements(k,i),1)-1)*dim+mpc(m,elements(k,i),2)) &
								 * Ke((j-1)*dim+l,(k-1)*dim+m) * p((mpc(m,elements(k,i),1)-1)*dim+mpc(m,elements(k,i),2))
					end if

				end do
				end do
			end do
			end do
		end do
!$omp end do
!$omp end parallel

		a = dot_product(r0,r0)/dot_product(p,kp)
		u = u + a*p

		r1 = r0 - a*kp

! 		print *,"step = ", step, ", |r| = ", dsqrt(dot_product(r1,r1))/norm_r0
		if (dsqrt(dot_product(r1,r1))/norm_r0<1d-15) exit

		b = dot_product(r1,r1)/dot_product(r0,r0)

		p = r1 + b*p
		r0 = r1
	end do
	print *,step-1

	do i=1,n_nds*dim
		u(i) = u(i)*s(i)
	end do

end subroutine

end module
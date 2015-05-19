program ebescg
	use fem_module
	implicit none
	type(struct_model) :: model
	type(struct_bc) :: bc
	character name_mdl*16
	real(8),allocatable :: f(:),u(:)


	write(*,'(a)', advance='no') "Model: "
	read *, name_mdl; print *

	print *, "Reading models..."; print *

	call read_model(model,name_mdl)

	allocate(u(model%n_nds*model%dim))

 	print *,"Reading boundary conditions..."; print *
	call read_cc(bc,model)

	call read_f(f,model)

	call ebe_scg(u,model,bc,f)

	read *


contains

subroutine ebe_scg(u,model,bc,f)
	type(struct_model) :: model
	type(struct_bc) :: bc
	real(8) :: u(:),f(:),a,b
	real(8),allocatable :: s(:),Ke(:,:),r0(:),r1(:),p(:),Ku0(:),kp(:), Kall(:,:),norm_r0
	integer i,j,k,l,m,q,step
	integer n_nds,n_els,n_nds_1el,dim
	integer,pointer :: elements(:,:)

	n_nds = model%n_nds
	n_els = model%n_els
	n_nds_1el = model%n_nds_1el
	dim = model%dim
	elements => model%elements

	allocate(s(n_nds*dim))
	allocate(Ke(n_nds_1el*dim,n_nds_1el*dim))
	allocate(r0(n_nds*dim),r1(n_nds*dim))
	allocate(p(n_nds*dim))
	allocate(Ku0(n_nds*dim),kp(n_nds*dim))

	allocate(Kall(n_nds*dim,n_nds*dim))

	p = 0d0

!$omp parallel private(Ke)

	! extract diagonal parts of stiffness matrix (stored by 's')
	s = 0d0
	do i=1,n_els
		Ke = 0d0
		call calc_element_integral(Ke,model,i,calc_BDB)
		do j=1,n_nds_1el
			do k=1,dim

!!
				s((elements(j,i)-1)*dim+k) = s((elements(j,i)-1)*dim+k) + Ke((j-1)*dim+k,(j-1)*dim+k)
! 				s((elements(j,i)-1)*dim+k) = s((elements(j,i)-1)*dim+k) + 1d0
! print *,Ke((j-1)*dim+k,(j-1)*dim+k)

! 	do l=1,bc%n_mpc
! 		if ((bc%dir_mpc(k,l) /= 0) .and. (bc%nodes_mpc(2,l) == elements(j,i))) then
! 			s((bc%nodes_mpc(1,l)-1)*dim+k) = s((bc%nodes_mpc(1,l)-1)*dim+k) + Ke((j-1)*dim+k,(j-1)*dim+k)
! 		end if
! 	end do

	do l=1,bc%n_spc
		if ((bc%nodes_spc(k+1,l) /= 0) .and. (bc%nodes_spc(1,l) == elements(j,i))) then
! 			Ke((j-1)*dim+k,(j-1)*dim+k) = Ke((j-1)*dim+k,(j-1)*dim+k) * 1d20
			s((bc%nodes_spc(1,l)-1)*dim+k) = 1d20

		end if
	end do

	do l=1,bc%n_mpc
		if ((bc%dir_mpc(k,l) /= 0) .and. (bc%nodes_mpc(2,l) == elements(j,i))) then
! 			Ke((j-1)*dim+k,(j-1)*dim+k) = Ke((j-1)*dim+k,(j-1)*dim+k) * 1d20
			s((bc%nodes_mpc(2,l)-1)*dim+k) = 1d20
		end if
	end do

!!
			end do
		end do
	end do

print *,s

	! calculate 's'　to the power of -1/2
	do i=1,n_nds*dim
		s(i) = 1d0/dsqrt(s(i))
! 		s(i) = 1d0/s(i)
	end do

! 	s = 1d0

! print *,s

! CG method
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
		if (bc%dir_mpc(j,i) /= 0) then
			u((bc%nodes_mpc(2,i)-1)*dim+j) = u((bc%nodes_mpc(1,i)-1)*dim+j)
		end if
	end do
	end do

	Ku0 = 0d0
	do i=1,n_els
		Ke = 0d0
		call calc_element_integral(Ke,model,i,calc_BDB)
		do j=1,n_nds_1el
		do k=1,n_nds_1el
			do l=1,dim
			do m=1,dim
				Ku0((elements(j,i)-1)*dim+l) = Ku0((elements(j,i)-1)*dim+l)&
! 					 + s((elements(j,i)-1)*dim+l)*s((elements(k,i)-1)*dim+m) * Ke((j-1)*dim+l,(k-1)*dim+m)*u((elements(k,i)-1)*dim+m)
					 + Ke((j-1)*dim+l,(k-1)*dim+m)*u((elements(k,i)-1)*dim+m)
			end do
			end do
		end do
		end do
	end do


	r0 = f - Ku0

	! scaling of r0
	do i=1,n_nds*dim
		r0(i) = s(i)*r0(i)
	end do
! 	! scaling of u
! 	do i=1,n_nds*dim
! 		u(i) = u(i)/s(i)
! 	end do
u = 0d0

! 	do i=1,bc%n_mpc
! 	do j=1,dim
! 		if (bc%dir_mpc(j,i) /= 0) then
! 			r0((bc%nodes_mpc(2,i)-1)*dim+j) = 0d0
! 		end if
! 	end do
! 	end do



	do i=1,n_nds
	print *,r0((i-1)*dim+1:i*dim)
	end do

! 	r0(1) = 0d0
! 	r0(2) = 0d0
! 	r0(3) = 0d0
! 	r0(4) = 0d0
! 	r0(5) = 0d0
! 	r0(6) = 0d0
	p = r0
	a = 0d0
	b = 0d0


print *,"p = "
print *,p


! 	print *,f

	norm_r0 = dot_product(r0,r0)

! 	print *,norm_r0

	do step=1,dim*n_nds
		! calc Ke.p
		kp = 0d0
		do i=1,n_els
			Ke = 0d0
			call calc_element_integral(Ke,model,i,calc_BDB)
			do j=1,n_nds_1el
			do k=1,n_nds_1el

				do l=1,dim
				do m=1,dim


					Kp((elements(j,i)-1)*dim+l) = Kp((elements(j,i)-1)*dim+l)&
					 + s((elements(j,i)-1)*dim+l)*s((elements(k,i)-1)*dim+m) * Ke((j-1)*dim+l,(k-1)*dim+m)*p((elements(k,i)-1)*dim+m)
!!!


!!!

	do q=1,bc%n_mpc
		if ((bc%dir_mpc(l,q) /= 0) .and. (bc%nodes_mpc(2,q) == elements(k,i))) then
					Kp((elements(j,i)-1)*dim+l) = Kp((elements(j,i)-1)*dim+l)&
					 + s((elements(j,i)-1)*dim+l)*s((elements(k,i)-1)*dim+m) * Ke((j-1)*dim+l,(k-1)*dim+m)*p((bc%nodes_mpc(1,q)-1)*dim+l)
		end if
	end do


	do q=1,bc%n_mpc
		if ((bc%dir_mpc(l,q) /= 0) .and. (bc%nodes_mpc(2,q) == elements(j,i))) then
			if (bc%nodes_mpc(1,q) /= elements(k,i)) then
					Kp((bc%nodes_mpc(1,q)-1)*dim+l) = Kp((bc%nodes_mpc(1,q)-1)*dim+l)&
					 + s((elements(j,i)-1)*dim+l)*s((elements(k,i)-1)*dim+m) * Ke((j-1)*dim+l,(k-1)*dim+m)*p((elements(k,i)-1)*dim+l)
			end if
		end if
	end do


if ((j==k).and.(l==m)) then
	do q=1,bc%n_spc
		if ((bc%nodes_spc(l+1,q) /= 0) .and. (bc%nodes_spc(1,q) == elements(j,i))) then
! 			Ke((j-1)*dim+l,(j-1)*dim+l) = Ke((j-1)*dim+l,(j-1)*dim+l) * 1d20
			Kp((elements(j,i)-1)*dim+l) = 1d20 * p((elements(k,i)-1)*dim+m)
		end if
	end do

	do q=1,bc%n_mpc
		if ((bc%dir_mpc(l,q) /= 0) .and. (bc%nodes_mpc(2,q) == elements(j,i))) then
! 			Ke((j-1)*dim+l,(j-1)*dim+l) = Ke((j-1)*dim+l,(j-1)*dim+l) * 1d20
			Ke((j-1)*dim+l,(j-1)*dim+l) = 1d20 * p((elements(k,i)-1)*dim+m)
		end if
	end do
end if



				end do
				end do
			end do
			end do
		end do

if (step == 1) then
print *,"p.Kp = ",dot_product(p,kp)
print *,"r.r = ",dot_product(r0,r0)
print *,dot_product(r0,r0)/dot_product(p,kp)

end if
		a = dot_product(r0,r0)/dot_product(p,kp)
print *,"alpha = ",a
		u = u + a*p

		r1 = r0 - a*kp

		print *,"step = ", step, ", |r| = ", dsqrt(dot_product(r1,r1))/norm_r0
		if (dsqrt(dot_product(r1,r1))/norm_r0<1d-15) exit

		b = dot_product(r1,r1)/dot_product(r0,r0)

		p = r1 + b*p
		r0 = r1
	end do

	do i=1,n_nds*dim
		u(i) = u(i)*s(i)
	end do

	do i=1,n_nds
		print *,u((i-1)*dim+1:i*dim)
	end do
! print *,Ku0



end subroutine

end program
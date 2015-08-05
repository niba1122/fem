module hm_module
	use fem_module
	implicit none
contains

subroutine hm_calc_BE_integral(BE,model)
	integer i,j
	type(struct_model) :: model
	integer dim,n_nds_1el,n_nds,n_els
	integer,pointer :: elements(:,:)
	double precision,allocatable :: BE_el(:,:),BE(:,:)
	dim = model%dim
	n_nds_1el = model%n_nds_1el
	n_nds = model%n_nds
	n_els = model%n_els
	elements => model%elements

	allocate(BE(dim*n_nds,dim*(dim+1)/2))
	allocate(BE_el(dim*n_nds_1el,dim*(dim+1)/2))

	BE = 0d0
	do i=1,n_els
		BE_el = 0d0
  	call calc_element_integral(BE_el,model,i,hm_calc_BE)
  	do j=1,n_nds_1el
  		BE((elements(j,i)-1)*dim+1:elements(j,i)*dim,:) = &
  			BE((elements(j,i)-1)*dim+1:elements(j,i)*dim,:) + BE_el((j-1)*dim+1:j*dim,:)
  	end do
  end do
end subroutine

function hm_calc_BE(model,i_el,coord)
		type(struct_model) :: model
		integer i_el,dim,n_nds_1el
		double precision coord(:)
		double precision,allocatable :: hm_calc_BE(:,:),B(:,:),D(:,:)
		double precision,pointer :: D_(:,:)

		dim = model%dim
		n_nds_1el = model%n_nds_1el

		D_ => model%materials(:,:,model%material_nos(i_el))

		allocate(hm_calc_BE(dim*n_nds_1el,dim*(dim+1)/2))
		call calc_B(B,model,i_el,coord)
		call calc_D(D,model,D_)
		hm_calc_BE(:,:) = matmul(transpose(B),D)

end function


subroutine hm_calc_EB_integral(EB,model)
	integer i,j
	type(struct_model) :: model
	integer dim,n_nds_1el,n_nds,n_els
	integer,pointer :: elements(:,:)
	double precision,allocatable :: EB_el(:,:),EB(:,:)
	dim = model%dim
	n_nds_1el = model%n_nds_1el
	n_nds = model%n_nds
	n_els = model%n_els
	elements => model%elements


	allocate(EB(dim*(dim+1)/2,dim*n_nds))
	allocate(EB_el(dim*(dim+1)/2,dim*n_nds_1el))

	EB = 0d0
	do i=1,n_els

  	call calc_element_integral(EB_el,model,i,hm_calc_EB)

  	do j=1,n_nds_1el
  		EB(:,(elements(j,i)-1)*dim+1:elements(j,i)*dim) = &
  			EB(:,(elements(j,i)-1)*dim+1:elements(j,i)*dim) + EB_el(:,(j-1)*dim+1:j*dim)
  	end do
  end do
end subroutine


function hm_calc_EB(model,i_el,coord)
		type(struct_model) :: model
		integer i_el,dim,n_nds_1el
		double precision coord(:)
		double precision,allocatable :: hm_calc_EB(:,:),B(:,:),D(:,:)
		double precision,pointer :: D_(:,:)

		dim = model%dim
		n_nds_1el = model%n_nds_1el

		D_ => model%materials(:,:,model%material_nos(i_el))

		allocate(hm_calc_EB(dim*(dim+1)/2,dim*n_nds_1el))
		call calc_B(B,model,i_el,coord)
		call calc_D(D,model,D_)

		hm_calc_EB(:,:) = matmul(D,B)

end function

subroutine hm_calc_E_integral(E,model)
	integer i,j
	type(struct_model) :: model
	integer dim,n_els
	double precision,allocatable :: E_el(:,:),E(:,:)
	dim = model%dim
	n_els = model%n_els


	allocate(E(dim*(dim+1)/2,dim*(dim+1)/2))
	allocate(E_el(dim*(dim+1)/2,dim*(dim+1)/2))


	E = 0d0
	do i=1,n_els
  	call calc_element_integral(E_el,model,i,hm_calc_E)
  	E = E + E_el
  end do

end subroutine

function hm_calc_E(model,i_el,coord)
		type(struct_model) :: model
		integer i_el,dim
		double precision coord(:)
		double precision,allocatable :: hm_calc_E(:,:),D(:,:)
		double precision,pointer :: D_(:,:)

		dim = model%dim

		D_ => model%materials(:,:,model%material_nos(i_el))

		allocate(hm_calc_E(dim*(dim+1)/2,dim*(dim+1)/2))
		call calc_D(D,model,D_)

		hm_calc_E(:,:) = D

end function

subroutine hm_auto_periodic_bc(bc,volume,model) 
	type(struct_model) :: model
	type(struct_bc) :: bc
	double precision,pointer :: nodes(:,:)
	integer i,j,n_nds,dim,n_mpc
	double precision,allocatable :: min_coord(:),max_coord(:),eps(:)
	integer,allocatable :: nodes_mpc(:,:)
	double precision volume

	nodes => model%nodes 
	n_nds = model%n_nds
	dim = model%dim

	allocate(min_coord(dim),max_coord(dim),eps(dim))


	n_mpc = 0
	allocate(nodes_mpc(2,n_nds))
	select case(dim)
		! 3D
		case(3)

		! get minimum or maximum values of x,y,z
		min_coord = 1d30
		max_coord = 0d0
		do j=1,n_nds
			do i=1,dim
				if (nodes(i,j) < min_coord(i)) min_coord(i) = nodes(i,j)
				if (nodes(i,j) > max_coord(i)) max_coord(i) = nodes(i,j)
			end do
		end do

		eps = (max_coord - min_coord)*1d-9

		volume = (max_coord(1) - min_coord(1))*(max_coord(2) - min_coord(2))*(max_coord(2) - min_coord(2))

		do i=1,n_nds
			! yz面(辺上を除く)
			if (( nodes(1,i)-min_coord(1) < eps(1) ) .and. &
				&(nodes(2,i)-min_coord(2)>eps(2)) .and. (max_coord(2)-nodes(2,i)>eps(2)) .and. &
					& (nodes(3,i)-min_coord(3)>eps(3)) .and. (max_coord(3)-nodes(3,i)>eps(3))) then
				do j=1,n_nds
					if (i /= j) then
					if ((max_coord(1) - nodes(1,j) < eps(1)) &
						&.and. (abs(nodes(2,i)-nodes(2,j)) < eps(2)) .and. (abs(nodes(3,i)-nodes(3,j)) < eps(3)) )  then
						n_mpc = n_mpc + 1
						nodes_mpc(1,n_mpc) = i
						nodes_mpc(2,n_mpc) = j
					end if
					end if
				end do
			end if
			! zx面(辺上を除く)
			if (( nodes(2,i)-min_coord(2) < eps(2) ) .and. &
				&(nodes(3,i)-min_coord(3)>eps(3)) .and. (max_coord(3)-nodes(3,i)>eps(3)) .and. &
					& (nodes(1,i)-min_coord(1)>eps(1)) .and. (max_coord(1)-nodes(1,i)>eps(1))) then
				do j=1,n_nds
					if (i /= j) then
					if ((max_coord(2) - nodes(2,j) < eps(2)) &
						&.and. (abs(nodes(3,i)-nodes(3,j)) < eps(3)) .and. (abs(nodes(1,i)-nodes(1,j)) < eps(1)) )  then
						n_mpc = n_mpc + 1
						nodes_mpc(1,n_mpc) = i
						nodes_mpc(2,n_mpc) = j
					end if
					end if
				end do
			end if
			! xy面(辺上を除く)
			if (( nodes(3,i)-min_coord(3) < eps(3) ) .and. &
				&(nodes(1,i)-min_coord(1)>eps(1)) .and. (max_coord(1)-nodes(1,i)>eps(1)) .and. &
					& (nodes(2,i)-min_coord(2)>eps(2)) .and. (max_coord(2)-nodes(2,i)>eps(2))) then
				do j=1,n_nds
					if (i /= j) then
					if ((max_coord(3) - nodes(3,j) < eps(3)) &
						&.and. (abs(nodes(1,i)-nodes(1,j)) < eps(1)) .and. (abs(nodes(2,i)-nodes(2,j)) < eps(2)) )  then
						n_mpc = n_mpc + 1
						nodes_mpc(1,n_mpc) = i
						nodes_mpc(2,n_mpc) = j
					end if
					end if
				end do
			end if

			! x軸(両端を除く)
			if ( ( nodes(2,i)-min_coord(2) < eps(2) ) .and. ( nodes(3,i)-min_coord(3) < eps(3) ) .and. &
				&( nodes(1,i)-min_coord(1) > eps(1) ) .and. ( max_coord(1)-nodes(1,i) > eps(1) ) ) then
				do j=1,n_nds
					if (i /= j) then
					if ( ( abs(nodes(1,i) - nodes(1,j)) < eps(1)) .and. &
						&  ( (nodes(2,j) - min_coord(2) < eps(2)) .or. (max_coord(2) - nodes(2,j) < eps(2)) ) .and. &
							&  ( (nodes(3,j) - min_coord(3) < eps(3)) .or. (max_coord(3) - nodes(3,j) < eps(3)) ) ) then
						n_mpc = n_mpc + 1
						nodes_mpc(1,n_mpc) = i
						nodes_mpc(2,n_mpc) = j
					end if
					end if
				end do
			end if

			! y軸(両端を除く)
			if ( ( nodes(3,i)-min_coord(3) < eps(3) ) .and. ( nodes(1,i)-min_coord(1) < eps(1) ) .and. &
				&( nodes(2,i)-min_coord(2) > eps(2) ) .and. ( max_coord(2)-nodes(2,i) > eps(2) ) ) then
				do j=1,n_nds
					if (i /= j) then
					if ( ( abs(nodes(2,i) - nodes(2,j)) < eps(2)) .and. &
						&  ( (nodes(3,j) - min_coord(3) < eps(3)) .or. (max_coord(3) - nodes(3,j) < eps(3)) ) .and. &
							&  ( (nodes(1,j) - min_coord(1) < eps(1)) .or. (max_coord(1) - nodes(1,j) < eps(1)) ) ) then
						n_mpc = n_mpc + 1
						nodes_mpc(1,n_mpc) = i
						nodes_mpc(2,n_mpc) = j
					end if
					end if
				end do
			end if

			! z軸(両端を除く)
			if ( ( nodes(1,i)-min_coord(1) < eps(1) ) .and. ( nodes(2,i)-min_coord(2) < eps(2) ) .and. &
				&( nodes(3,i)-min_coord(3) > eps(3) ) .and. ( max_coord(3)-nodes(3,i) > eps(3) ) ) then
				do j=1,n_nds
					if (i /= j) then
					if ( ( abs(nodes(3,i) - nodes(3,j)) < eps(3)) .and. &
						&  ( (nodes(1,j) - min_coord(1) < eps(1)) .or. (max_coord(1) - nodes(1,j) < eps(1)) ) .and. &
							&  ( (nodes(2,j) - min_coord(2) < eps(2)) .or. (max_coord(2) - nodes(2,j) < eps(2)) ) ) then
						n_mpc = n_mpc + 1
						nodes_mpc(1,n_mpc) = i
						nodes_mpc(2,n_mpc) = j
					end if
					end if
				end do
			end if

			! 頂点
			if ( ( nodes(1,i)-min_coord(1) < eps(1) ) .and.&
				&( nodes(2,i)-min_coord(2) < eps(2) ) .and.&
					&( nodes(3,i)-min_coord(3) < eps(3) ) ) then
				do j=1,n_nds
					if (i /= j) then
					if ( ((nodes(1,j) - min_coord(1) < eps(1)) .or. (max_coord(1) - nodes(1,j) < eps(1)) ).and.&
						&((nodes(2,j) - min_coord(2) < eps(2)) .or. (max_coord(2) - nodes(2,j) < eps(2)) ).and.&
							&((nodes(3,j) - min_coord(3) < eps(3)) .or. (max_coord(3) - nodes(3,j) < eps(3)) ) ) then
						n_mpc = n_mpc + 1
						nodes_mpc(1,n_mpc) = i
						nodes_mpc(2,n_mpc) = j
					end if
					end if
				end do
			end if
		end do
	! 2D
		case(2)
	end select

! 	allocate(nodes_mpc_(2,n_mpc))
! 	nodes_mpc_ = nodes_mpc(:,1:n_mpc)
! 	deallocate(nodes_mpc)
! 	allocate(nodes_mpc(2,n_mpc))

! 	nodes_mpc(:,:) = nodes_mpc_(:,:)

! 	print *,nodes_mpc(:,:)



! 	allocate(bc%nodes_mpc(2,n_mpc))
! 	allocate(bc%dir_mpc(dim,n_mpc))

! 	bc%n_mpc = n_mpc
! 	bc%nodes_mpc(:,:) = nodes_mpc(:,1:n_mpc)
! 	bc%dir_mpc = 1

	allocate(bc%mpc_master_id(n_mpc*dim))
	allocate(bc%mpc_slave_id(n_mpc*dim))
	allocate(bc%mpc_master_dir(n_mpc*dim))
	allocate(bc%mpc_slave_dir(n_mpc*dim))
	allocate(bc%mpc_param(n_mpc*dim))

	do i=1,n_mpc
		bc%mpc_master_id((i-1)*dim+1:i*dim) = nodes_mpc(1,i)
		bc%mpc_slave_id((i-1)*dim+1:i*dim) = nodes_mpc(2,i)
		do j=1,dim
			bc%mpc_master_dir((i-1)*dim+j) = j
			bc%mpc_slave_dir((i-1)*dim+j) = j
		end do
	end do
	bc%mpc_param = 1d0
	bc%n_mpc = n_mpc*dim

	bc%n_spc = 0

end subroutine

end module hm_module

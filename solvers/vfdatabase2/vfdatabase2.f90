! warp用の物性値
program vfdatabase
  use hm_module
  implicit none
  type(struct_model) :: model
  type(struct_bc) :: bc
  type(struct_output) :: output
  integer i,n_model
  integer j,n_vf,vf_no_offset
  real(8) volume,vf,min_vf,d_vf
  real(8),allocatable :: k(:),u(:),BE(:,:),x(:,:),EB(:,:),E(:,:)
  real(8) :: DH(6,6),CH(6,6)
  integer,allocatable :: index_k_diag(:)

  min_vf = 0.38d0
  d_vf = 0.01d0
  n_vf = 40
  vf_no_offset = 10

	open(11,file=trim(path_model)//"vfdatabase2/materials.csv")
  write(11,*) n_vf

  do j=1,n_vf
    print *, "##################################################"
    print *, "# Homogenization analysis                        #"
    print *, "##################################################"

    vf = min_vf + d_vf*(j-1)
    print *, "Volume fraction = ", vf*100, " %"

    print *, "Generating unit cell ..."
    call generate_fiber_unit_cell(model,vf)
    print *, "Setting boundary conditions ..."
  	call hm_auto_periodic_bc(bc,volume,model)
    print *, "Calculating integral of BE ..."
   	call hm_calc_BE_integral(BE,model)
    print *, "Initializing k matrix ..."
  	call init_addition_matrix_sln(k,index_k_diag,model,bc)
    
    allocate(x(model%n_nds*model%dim,model%dim*(model%dim+1)/2))
  
    print *, "Solving equations ..."
  	do i=1,model%dim*(model%dim+1)/2
  		k = 0d0
  		call calc_addition_matrix_sln(k,model,index_k_diag,calc_BDB)
  		call set_cc(k,BE(:,i),model,bc,index_k_diag)
  
   		x(:,i) = sl_LDL(K,BE(:,i),index_k_diag)
   		call set_constrained_u(x(:,i),model,bc)
  
   	end do

    print *, "Calculating EB matrix ..."
   	call hm_calc_EB_integral(EB,model)
    print *, "Calculating E matrix ..."
   	call hm_calc_E_integral(E,model)

    DH = E-matmul(EB,x)
  
  CH = inv_jordan(DH)
    do i=1,6
      !print *,CH(i,:)
    end do
 ! print *,1d0/CH(3,3)

    print *, "Outputting results ..."
   	call calc_output(output,model,x(:,1))
   	call output_inp(model,output,'vfdatabase2')

    write(11,*) j+vf_no_offset-1
    do i=1,6
      write(11,*) DH(i,1),',',DH(i,2),',',DH(i,3),',',DH(i,4),',',DH(i,5),',',DH(i,6)
    end do

    print *, "Clearing variables ..."
    call clear_addition_matrix_sln(k,index_k_diag)
    call clear_model(model)
    call clear_bc(bc)
    deallocate(x)
    deallocate(EB)
    deallocate(BE)
    deallocate(E)
    print *, "Analysis finished"
  end do
  close(10)

contains

subroutine generate_fiber_unit_cell(model,vf)
  type(struct_model) :: model
  real(8) vf
  real(8),pointer :: nodes(:,:),materials(:,:,:)
  integer,pointer :: elements(:,:),material_nos(:)
  integer n_nds,n_els,n_nds_1el,dim
  integer i,j,k
  integer node1
  real(8) pi,cos225,sin225,r
  integer node_offset
  integer node_no_1_4(25)

  n_nds = 9*9*11
  n_els = 8*8*10
  n_nds_1el = 8
  dim = 3

  allocate(model%nodes(dim,n_nds))
  allocate(model%elements(n_nds_1el,n_els))
  allocate(model%material_nos(n_els))
  allocate(model%materials(6,6,2))
  nodes => model%nodes 
  elements => model%elements
  material_nos => model%material_nos
  materials => model%materials
  model%n_nds = n_nds
  model%n_els = n_els
  model%type_el = 'hexa8'
  model%n_nds_1el = 8
  model%name = 'vfdatabase2'
  model%dim = dim


  ! - Generating Elements -

  do k=1,10
    do j=1,8
      do i=1,8
        node1 = i+(j-1)*9+(k-1)*81
        elements(:,(k-1)*64+(j-1)*8+i) = (/node1,node1+1,node1+10,node1+9,node1+81,node1+82,node1+91,node1+90/)
      end do
    end do
  end do

  ! - Generating Nodes -

  pi = dacos(-1d0)
  cos225 = dcos(pi/8)
  sin225 = dsin(pi/8)
  r = dsqrt( (2d0*vf) / (16*dsin(2d0*pi/16)) )
  print *,pi,cos225,sin225,r,pi*r**2
  print *,8*(r**2)*sin225

  node_no_1_4(1:25) = (/1,2,3,4,5,10,11,12,13,14,19,20,21,22,23,28,29,30,31,32,37,38,39,40,41/)
  
  nodes(2:3,node_no_1_4(1)) = (/-0.5d0,-0.5d0/)
  nodes(2:3,node_no_1_4(2)) = (/-(0.5d0+r/dsqrt(2d0))/2,-0.5d0/)
  nodes(2:3,node_no_1_4(3)) = (/-r/dsqrt(2d0),-0.5d0/)
  nodes(2:3,node_no_1_4(4)) = (/-r*sin225,-0.5d0/)
  nodes(2:3,node_no_1_4(5)) = (/0d0,-0.5d0/)
  nodes(2:3,node_no_1_4(6)) = (/-0.5d0,-(0.5d0+r/dsqrt(2d0))/2/)
  nodes(2:3,node_no_1_4(7)) = (/-(0.5d0+r/dsqrt(2d0))/2,-(0.5d0+r/dsqrt(2d0))/2/)
  nodes(2:3,node_no_1_4(8)) = (/-r/dsqrt(2d0),-(0.5d0+r/dsqrt(2d0))/2/)
  nodes(2:3,node_no_1_4(9)) = (/-r*sin225,-(0.5d0+r*cos225)/2/)
  nodes(2:3,node_no_1_4(10)) = (/-0d0,-(0.5d0+r)/2/)
  nodes(2:3,node_no_1_4(11)) = (/-0.5d0,-r/dsqrt(2d0)/)
  nodes(2:3,node_no_1_4(12)) = (/-(0.5d0+r/dsqrt(2d0))/2,-r/dsqrt(2d0)/)
  nodes(2:3,node_no_1_4(13)) = (/-r/dsqrt(2d0),-r/dsqrt(2d0)/)
  nodes(2:3,node_no_1_4(14)) = (/-r*sin225,-r*cos225/)
  nodes(2:3,node_no_1_4(15)) = (/0d0,-r/)
  nodes(2:3,node_no_1_4(16)) = (/-0.5d0,-r*sin225/)
  nodes(2:3,node_no_1_4(17)) = (/-(0.5d0+r*cos225)/2,-r*sin225/)
  nodes(2:3,node_no_1_4(18)) = (/-r*cos225,-r*sin225/)
  nodes(2:3,node_no_1_4(19)) = (/-r*sin225,-r*sin225/)
  nodes(2:3,node_no_1_4(20)) = (/0d0,-r/2/)
  nodes(2:3,node_no_1_4(21)) = (/-0.5d0,0d0/)
  nodes(2:3,node_no_1_4(22)) = (/-(0.5d0+r)/2,0d0/)
  nodes(2:3,node_no_1_4(23)) = (/-r,0d0/)
  nodes(2:3,node_no_1_4(24)) = (/-r/2,0d0/)
  nodes(2:3,node_no_1_4(25)) = (/0d0,0d0/)

  do i=1,25
    ! 左上
    nodes(2:3,72+mod(node_no_1_4(i),9)*2-node_no_1_4(i)) = (/nodes(2,node_no_1_4(i)),-nodes(3,node_no_1_4(i))/)
    ! 右下
    nodes(2:3,(10+node_no_1_4(i)/9*18)-node_no_1_4(i)) = (/-nodes(2,node_no_1_4(i)),nodes(3,node_no_1_4(i))/)
    ! 右上
    nodes(2:3,82-node_no_1_4(i)) = -nodes(2:3,node_no_1_4(i))
  end do

  do i=1,11
    nodes(2:3,(i-1)*81+1:i*81) = nodes(2:3,1:81)
    nodes(1,(i-1)*81+1:i*81) = 0.1d0*(i-1)
  end do
  
  ! - Generating Material No. -
  do i=1,n_els
    if ( ((19<=mod(i,64)) .and. (mod(i,64)<=22)) .or. ((27<=mod(i,64)) .and. (mod(i,64)<=30)) &
      & .or. ((35<=mod(i,64)) .and. (mod(i,64)<=38)) .or. ((43<=mod(i,64)) .and. (mod(i,64)<=46)) ) then
      material_nos(i) = 1
    else
      material_nos(i) = 2
    end if
  end do

  ! - materials -

  model%materials(1,:,1) = (/0.977982954545455d+02,0.380326704545455d+02,0.380326704545455d+02,0.0d+00,0.0d+00,0.0d+00/)
  model%materials(2,:,1) = (/0.380326704545455d+02,0.977982954545455d+02,0.380326704545455d+02,0.0d+00,0.0d+00,0.0d+00/)  
  model%materials(3,:,1) = (/0.380326704545455d+02,0.380326704545455d+02,0.977982954545455d+02,0.0d+00,0.0d+00,0.0d+00/)  
  model%materials(4,:,1) = (/0.0d+00,0.0d+00,0.0d+00,0.298828125000000d+02,0.0d+00,0.0d+00/) 
  model%materials(5,:,1) = (/0.0d+00,0.0d+00,0.0d+00,0.0d+00,0.298828125000000d+02,0.0d+00/) 
  model%materials(6,:,1) = (/0.0d+00,0.0d+00,0.0d+00,0.0d+00,0.0d+00,0.298828125000000d+02/) 
  
  model%materials(1,:,2) = (/0.527692307692308d+01,0.226153846153846d+01,0.226153846153846d+01,0.0d+00,0.0d+00,0.0d+00/)
  model%materials(2,:,2) = (/0.226153846153846d+01,0.527692307692308d+01,0.226153846153846d+01,0.0d+00,0.0d+00,0.0d+00/)
  model%materials(3,:,2) = (/0.226153846153846d+01,0.226153846153846d+01,0.527692307692308d+01,0.0d+00,0.0d+00,0.0d+00/)
  model%materials(4,:,2) = (/0.0d+00,0.0d+00,0.0d+00,0.150769230769231d+01,0.0d+00,0.0d+00/)
  model%materials(5,:,2) = (/0.0d+00,0.0d+00,0.0d+00,0.0d+00,0.150769230769231d+01,0.0d+00/)
  model%materials(6,:,2) = (/0.0d+00,0.0d+00,0.0d+00,0.0d+00,0.0d+00,0.150769230769231d+01/)

  model%materials = model%materials*1d9


end subroutine


end program vfdatabase

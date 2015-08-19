!-------------------------------------------------------------------------------------------------------
!
!  Solver for 2D damage propagation analysis
!
!  Copyright (C) 2015 Nobuhito Ibaraki All Rights Reserved.
!
!-------------------------------------------------------------------------------------------------------


program fem
  use fem_module
  use frp_generate_module
  implicit none
  integer i,step,max_step,model_no
  character(8) :: i_char,step_char
  double precision,allocatable,dimension(:) :: k,u,f
  real(8),allocatable,target :: angle(:)
  real(8),allocatable,target :: damage_tensor(:,:)
  real(8) u_, du ! 損傷進展解析の変位、変位増分
  integer damage_ratio_step
!  logical damage_ratio
  real(8) damage_ratio
  real(8) f_left, f_right ! 反力の合計(左右)
  real(8) :: max_sig(6,4) ! 損傷の閾値
  real(8) :: vf !繊維束の体積含有率
  integer,allocatable,dimension(:) :: index_k_diag
  double precision thickness

  type(struct_model) :: model
  type(struct_bc) :: bc
  type(struct_output) :: output
  type(struct_data) :: data(2)

  integer shift(4)
  character errmsg*128


  character(32) output_file_name,od_data_path,od_model_name
  character(128) command

!-------------------------------------------------------------------------------------------------------
!  Config
!-------------------------------------------------------------------------------------------------------

  du = 1d-5
  model_no = 1
  max_step = 20
  od_model_name = "gfrp_damage"

!-------------------------------------------------------------------------------------------------------
!  Solver
!-------------------------------------------------------------------------------------------------------


  write(command,'(a)') "mkdir "//trim(path_model)//trim(od_model_name)
print *,"mkdir "//trim(path_model)//trim(od_model_name)
  call system(command)


  open(11, file=trim(path_model)//trim(od_model_name)//slash//'models.csv')
  write(11,*) 'du, ',du
  write(11,*) 'model_no, ',model_no
  write(11,*) 'model_no, vf, stress2L, stress2T, stress2TL, stress3T, stress3Z, stress3TZ,&
    &shift(1),shift(2),shift(3),shift(4), step of initial fraction'

  do i=1,model_no
  write(i_char, '(i0)') i

  write(*,'(a)') " ------------------------------"
  write(*,'(a,i0,a)') "  calculation of model ",i," start"
  write(*,'(a)') " ------------------------------"


  shift(1) = 0
  shift(2) = uniform_rand_integer(0,53)
  shift(3) = uniform_rand_integer(0,53)
  shift(4) = uniform_rand_integer(0,53)

  vf = 0.6d0+normal_rand_real8(0.033d0)
  print *,'vf = ',vf
  call calc_strength(max_sig,vf)

!  max_sig = 1d30
!
!! 材料2の強度
!  max_sig(1,2) = 1000.3*1d6
!  max_sig(2,2) = 34.3*1d6
!  max_sig(6,2) = 40.2*1d6
!! 材料3の強度
!  max_sig(1,3) = 34.3*1d6
!  max_sig(2,3) = 34.3*1d6
!  max_sig(6,3) = 40.2*1d6


  call generate_FRP_Model(model,shift)

  model%name = od_model_name

  call set_state2d(model,1,1d0)


  print *,'Making directory...'

  write(od_data_path,'(a,i0,a)') "model", i, slash
  write(command,'(a,a,a,a,a)') "mkdir ", trim(path_model), trim(model%name),slash, trim(od_data_path)
print *,"mkdir ", trim(path_model), trim(model%name),slash, trim(od_data_path)
  call system(command)

  ! モデルごとのデータ出力
  open(20, file=trim(path_model)//trim(model%name)//slash//trim(od_data_path)//'model'//trim(i_char)//'.csv')
  write(20,*) 'f_left,f_right'
!  damage_ratio = .false.
  damage_ratio = 0d0
  u_ = 0d0

  allocate(f(model%dim*model%n_nds))
  allocate(u(model%n_nds*model%dim))

print *,ubound(model%data(1)%i,1),ubound(model%data(2)%i,1)

!  do step=1,max_step
step = 0
  do
    step = step+1
    if (step > max_step) then
      exit
    end if

    write(step_char, '(i0)') step


    u = 0d0
    f = 0d0

!    u_ = du*step
if (step == 2) then
print *, ceiling(1 / damage_ratio)
  step = ceiling(1 / damage_ratio)
end if
u_ = du * step
  

print *,u_
    call frp_set_tensile_bc(bc,model,u_)

     print *,"Calculating K matrix..."; print *

    call init_addition_matrix_sln(k,index_k_diag,model,bc)

    call calc_addition_matrix_sln(k,model,index_k_diag,od_calc_BDB)

    print *,"Solving stiffness equation..."; print *
    call set_cc(k,f,model,bc,index_k_diag)

    u = sl_LDL(K,f,index_k_diag)

    call calc_reaction_force(f_left,f_right,model,u)
    write(20,'(d30.15,",",d30.15)') f_left, f_right

    print *,'f:',f_left,f_right

    call set_constrained_u(u,model,bc)

    print *,"Calculating stress and strain..."; print *

    call od_calc_output(output,model,u)

    print *,"Judging the damage state..."; print *

!    if (damage_ratio) then
!      damage_ratio = damage_judgment(model,output,max_sig)
!print *,'already damage_ratio'
!    else
!      damage_ratio = damage_judgment(model,output,max_sig)
!      if (damage_ratio) then
!        write(11,'(i0,7(",",d30.15),5(",",i0))') &
!          &i, vf, max_sig(1,2), max_sig(2,2), max_sig(6,2), max_sig(1,3), max_sig(2,3), max_sig(6,3),&
!            &shift(1), shift(2), shift(3), shift(4), step
!print *, 'damage_ratio!!! step',step
!      else
!        print *,'undamage_ratio'
!      end if
!    end if

    if (damage_ratio>=1) then
      damage_ratio = damage_judgment(model,output,max_sig)
print *,'already damaged'
    else
      damage_ratio = damage_judgment(model,output,max_sig)
      if (damage_ratio>=1) then
        write(11,'(i0,7(",",d30.15),5(",",i0))') &
          &i, vf, max_sig(1,2), max_sig(2,2), max_sig(6,2), max_sig(1,3), max_sig(2,3), max_sig(6,3),&
            &shift(1), shift(2), shift(3), shift(4), step
print *, 'damaged!!! step',step
      else
        print *,'undamaged'
      end if
    end if
!    call visualize_u(model,u)

    print *,"Outputting file..."; print *

    write(output_file_name,'(a,i0)') "step", step
     call od_output_inp(model,output,trim(od_data_path)//output_file_name)


    print *,"Clearing data..."; print *

    call clear_output(output)
    call clear_bc(bc)
    call clear_addition_matrix_sln(k,index_k_diag)

  end do

  close(20)
  call clear_model(model)

  deallocate(f)
  deallocate(u)

  end do

  close(11)

  print *,"Calculation was completed successfully!"; print *

  print *,"- Press ENTER to finish -"; print *

!    call model_add_d_data(model,reshape((/1d0,2d0,3d0,4d0,5d0,6d0/),(/2,3,1/)),1,"arr1")
!    call model_add_d_data(model,reshape((/4d0,5d0,6d0,7d0,8d0,9d0/),(/1,6,1/)),1,"arr2")
!    call model_add_d_data(model,reshape((/1d0,2d0,3d0/),(/1,3,1/)),1,"arr3")

!    print *,model%d_data(2)%d
!    print *,model%d_data(3)%d
!    print *,model%d_data(4)%d

!    call model_get_d_data(test,model,"arr2")
!    print *, test


  read *

contains

subroutine calc_strength(max_sig,vf)
  real(8) :: vf, max_sig(:,:)
  

  max_sig = 1d30

! 材料2の強度
  max_sig(1,2) = 1000.3d6+(vf-0.6d0)*100d6 ! I
  max_sig(2,2) = 34.3d6+(vf-0.6d0)*(-18.5d6) ! II
  max_sig(6,2) = 40.2d6+(vf-0.6d0)*(-21.8d6) ! I II
! 材料3の強度
  max_sig(2,3) = 34.3d6+(vf-0.6d0)*(-18.5d6) ! II
  max_sig(3,3) = 34.3d6+(vf-0.6d0)*(-18.5d6) ! III 
  max_sig(4,3) = 40.2d6+(vf-0.6d0)*(-21.8d6) ! II III

end subroutine

subroutine calc_reaction_force(f_left,f_right,model,u)
  type(struct_model) :: model
  real(8) :: f_left,f_right,u(:),penalty
  integer,pointer :: left_nodes(:), right_nodes(:)
  integer n_left_nodes, n_right_nodes
  integer i,dim

  n_left_nodes = ubound(model%data(1)%i,1)
  n_right_nodes = ubound(model%data(2)%i,1)
  left_nodes => model%data(1)%i(:,1,1)
  right_nodes => model%data(2)%i(:,1,1)
  dim = model%dim
  penalty = 1d30

  f_left = 0d0
  do i=1,n_left_nodes
    f_left = f_left + u(left_nodes(i)*dim-1)*penalty
!print *, u(left_nodes(i)*dim-1)*penalty
!print *, u(left_nodes(i)*dim-1)
  end do
  f_right = 0d0
  do i=1,n_right_nodes
    f_right = f_right + u(right_nodes(i)*dim-1)*penalty
!print *, u(right_nodes(i)*dim-1)*penalty
!print *, u(right_nodes(i)*dim-1)
  end do
end subroutine

function damage_judgment(model,output,max_sig)
  type(struct_model) :: model
  type(struct_output) :: output
  integer n_els,i
  real(8),pointer :: sig(:,:),angle(:),damage_tensor(:,:)
  real(8) :: max_sig(:,:),sig_rot(6)
  real(8) damage_judgment,damage_ratio ! damage_ratio: sig/maxsig 1以上で損傷発生

  n_els = model%n_els
  angle => model%data(1)%d(:,1,1)
  sig => output%sig

  damage_tensor => model%data(2)%d(:,:,1)

  do i=1,n_els
    sig_rot = rot_sig(sig(:,i),angle(i))

    if (model%material_nos(i) == 2) then ! x:L y:T z:Z
      if (max_sig(1,2) < sig_rot(1)) then
        damage_tensor(1,i) = 0.99d0
      else if (max_sig(2,2) < sig_rot(2)) then
        damage_tensor(2,i) = 0.99d0
      else if (max_sig(6,2) < sig_rot(6)) then
        damage_tensor(2,i) = 0.99d0
      end if


      damage_ratio = dabs(sig_rot(1)/max_sig(1,2))
      if (damage_ratio > damage_judgment) then
        damage_judgment = damage_ratio
      end if
      damage_ratio = dabs(sig_rot(2)/max_sig(2,2))
      if (damage_ratio > damage_judgment) then
        damage_judgment = damage_ratio
      end if
      damage_ratio = dabs(sig_rot(6)/max_sig(6,2))
      if (damage_ratio > damage_judgment) then
        damage_judgment = damage_ratio
      end if


    else if (model%material_nos(i) == 3) then ! x:T y:Z z:L
      if (max_sig(2,3) < sig_rot(1)) then
        damage_tensor(2,i) = 0.99d0
      else if (max_sig(3,3) < sig_rot(2)) then
        damage_tensor(3,i) = 0.99d0
      else if (max_sig(4,3) < sig_rot(6)) then
        damage_tensor(2,i) = 0.99d0
        damage_tensor(3,i) = 0.99d0
      end if


      damage_ratio = dabs(sig_rot(1)/max_sig(2,3))
      if (damage_ratio > damage_judgment) then
        damage_judgment = damage_ratio
      end if
      damage_ratio = dabs(sig_rot(2)/max_sig(3,3))
      if (damage_ratio > damage_judgment) then
        damage_judgment = damage_ratio
      end if
      damage_ratio = dabs(sig_rot(6)/max_sig(4,3))
      if (damage_ratio > damage_judgment) then
        damage_judgment = damage_ratio
      end if


    end if


  end do
print *,"damage_judgment = ", damage_judgment
end function


function uniform_rand_integer(min,max)
  implicit none
  integer min, max, interval
  real(8) uniform_rand_integer, genrand_res53

  interval = abs(max-min)+1

  uniform_rand_integer = int(genrand_res53()*interval)
end function


function normal_rand_real8(sd)
  implicit none
  real(8) normal_rand_real8,x,y,z,pi,genrand_res53,sd

  pi = dacos(-1d0)

  do
    x = genrand_res53()
    y = genrand_res53()
    z = dsqrt(-2d0*dlog(x)) * dCOS(2*pi*y)
    if (abs(z)<=6d0) exit
  end do
  if (z>6d0) then
    z = 6d0
  else if (z<-6d0) then
    z = -6d0
  end if

  normal_rand_real8 = z*sd
end function


subroutine frp_set_tensile_bc(bc,model,ux) ! 生成した2周期4層のFRPモデルにx方向引張の強制変位を与える
  type(struct_model) :: model
  type(struct_bc) :: bc
  real(8) ux
  integer i,n_top_nodes,n_bottom_nodes,n_left_nodes,n_right_nodes,spc_index,n_spc,dim
  integer,pointer :: left_nodes(:), right_nodes(:), bottom_nodes(:), top_nodes(:), left_center_node, right_center_node

  n_left_nodes = ubound(model%data(1)%i,1)
  n_right_nodes = ubound(model%data(2)%i,1)
  n_bottom_nodes = ubound(model%data(3)%i,1)
  n_top_nodes = ubound(model%data(4)%i,1)

  left_nodes => model%data(1)%i(:,1,1)
  right_nodes => model%data(2)%i(:,1,1)
  bottom_nodes => model%data(3)%i(:,1,1)
  top_nodes => model%data(4)%i(:,1,1)
  left_center_node => model%data(5)%i(1,1,1)
  right_center_node => model%data(6)%i(1,1,1)

  dim = model%dim

  n_spc = n_left_nodes+n_right_nodes

  bc%n_spc = n_spc


  allocate(bc%nodes_spc(dim+1,n_spc))
  allocate(bc%disp_spc(dim,n_spc))
  
  bc%nodes_spc = 0
  bc%disp_spc = 0d0

  spc_index = 1
  do i=1,n_left_nodes
    if (left_nodes(i) == left_center_node) then
      bc%nodes_spc(:,spc_index) = (/left_nodes(i),1,1/)
      bc%disp_spc(:,spc_index) = (/0d0,0d0/)
    else
      bc%nodes_spc(:,spc_index) = (/left_nodes(i),1,0/)
      bc%disp_spc(:,spc_index) = (/0d0,0d0/)
    end if
    spc_index = spc_index+1
  end do

  do i=1,n_right_nodes
    if (right_nodes(i) == right_center_node) then
      bc%nodes_spc(:,spc_index) = (/right_nodes(i),1,1/)
      bc%disp_spc(:,spc_index) = (/ux,0d0/)
    else
      bc%nodes_spc(:,spc_index) = (/right_nodes(i),1,0/)
      bc%disp_spc(:,spc_index) = (/ux,0d0/)
    end if
    spc_index=spc_index+1
  end do

  bc%n_mpc = 0

end subroutine


function od_calc_BDB(model,i_el,coord)
  type(struct_model) :: model
  integer :: i_el
  double precision :: coord(:)
  double precision,allocatable,dimension(:,:) :: od_calc_BDB
  double precision,allocatable :: B(:,:),D(:,:)
  double precision :: D_3d(6,6)


  allocate(od_calc_BDB(model%dim*model%n_nds_1el,model%dim*model%n_nds_1el))

  od_calc_BDB = 0d0

  D_3d = od_set_damage_tensor(model,i_el) ! 損傷テンソルをDマトリックスへ適用

  D_3d = rot_D(D_3d,model%data(1)%d(i_el,1,1)) ! Dマトリックスの回転

  call calc_D(D,model,D_3d)

  call calc_B(B,model,i_el,coord)


  od_calc_BDB(:,:) = matmul( matmul(transpose(B),D),B)


end function

function od_set_damage_tensor(model,i_el)
  type(struct_model) :: model
  real(8) :: dl,dt,dz
  real(8) :: od_set_damage_tensor(6,6)
  integer i_el

  dl = 1d0-model%data(2)%d(1,i_el,1)
  dt = 1d0-model%data(2)%d(2,i_el,1)
  dz = 1d0-model%data(2)%d(3,i_el,1)

  od_set_damage_tensor(1:6,1:6) = model%materials(1:6,1:6,model%material_nos(i_el)) 

  if (model%material_nos(i_el) == 2) then
    od_set_damage_tensor(1,:) = od_set_damage_tensor(1,:) * dl
    od_set_damage_tensor(2,:) = od_set_damage_tensor(2,:) * dt
    od_set_damage_tensor(3,:) = od_set_damage_tensor(3,:) * dz
    od_set_damage_tensor(:,1) = od_set_damage_tensor(:,1) * dl
    od_set_damage_tensor(:,2) = od_set_damage_tensor(:,2) * dt
    od_set_damage_tensor(:,3) = od_set_damage_tensor(:,3) * dz
    od_set_damage_tensor(4,4) = od_set_damage_tensor(4,4) * 4d0 * (dt*dz/(dt+dz))**2
    od_set_damage_tensor(5,5) = od_set_damage_tensor(5,5) * 4d0 * (dz*dl/(dz+dl))**2
    od_set_damage_tensor(6,6) = od_set_damage_tensor(6,6) * 4d0 * (dl*dt/(dl+dt))**2
  else if (model%material_nos(i_el) == 3) then
    od_set_damage_tensor(1,:) = od_set_damage_tensor(1,:) * dt
    od_set_damage_tensor(2,:) = od_set_damage_tensor(2,:) * dz
    od_set_damage_tensor(3,:) = od_set_damage_tensor(3,:) * dl
    od_set_damage_tensor(:,1) = od_set_damage_tensor(:,1) * dt
    od_set_damage_tensor(:,2) = od_set_damage_tensor(:,2) * dz
    od_set_damage_tensor(:,3) = od_set_damage_tensor(:,3) * dl
    od_set_damage_tensor(4,4) = od_set_damage_tensor(4,4) * 4d0 * (dz*dl/(dz+dl))**2
    od_set_damage_tensor(5,5) = od_set_damage_tensor(5,5) * 4d0 * (dl*dt/(dl+dt))**2
    od_set_damage_tensor(6,6) = od_set_damage_tensor(6,6) * 4d0 * (dt*dz/(dt+dz))**2
  end if

end function


function rot_sig(sig,theta)
  implicit none
  integer i,j,k,l
  double precision theta,pi
  double precision,dimension(:) :: sig,rot_sig(6)
  double precision,dimension(:,:) :: sig_tensor(6,6),sig_tensor_rot(6,6)
  double precision,allocatable,dimension(:,:) :: beta

  allocate(beta(3,3))

  ! 円周率の定義
  pi = dacos(-1d0);

  ! 座標軸の回転変換テンソルの定義
  beta(1,:) = (/dcos(theta), dsin(theta), 0d0/);
  beta(2,:) = (/-dsin(theta), dcos(theta), 0d0/);
  beta(3,:) = (/0d0, 0d0, 1d0/);

  ! 応力テンソルへ変換
  sig_tensor(1,1) = sig(1)
  sig_tensor(2,2) = sig(2)
  sig_tensor(3,3) = sig(3)
  sig_tensor(2,3) = sig(4)
  sig_tensor(3,2) = sig(4)
  sig_tensor(3,1) = sig(5)
  sig_tensor(1,3) = sig(5)
  sig_tensor(1,2) = sig(6)
  sig_tensor(2,1) = sig(6)


  ! 弾性テンソルDの座標変換

  sig_tensor_rot = 0d0;

  do i=1,3
  do j=1,3
  do k=1,3
  do l=1,3
    sig_tensor_rot(i,j) = sig_tensor_rot(i,j)+beta(i,k)*beta(j,l)*sig_tensor(i,j);
  end do
  end do
  end do
  end do

  rot_sig(1) = sig_tensor_rot(1,1)
  rot_sig(2) = sig_tensor_rot(2,2)
  rot_sig(3) = sig_tensor_rot(3,3)
  rot_sig(4) = sig_tensor_rot(2,3)
  rot_sig(5) = sig_tensor_rot(3,1)
  rot_sig(6) = sig_tensor_rot(1,2)

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

  ! 弾性マトリックスDの計算(inv(C))
  D2 = D_in;


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

  rot_D  = Dr2

!   deallocate(beta)
!   deallocate(C)
!   deallocate(D2)
!   deallocate(D4)
!   deallocate(Dr4)
!   deallocate(Cr)
!   deallocate(Cr2D)

end function


subroutine od_calc_output(output,model,u)
  type(struct_model) :: model
  type(struct_output) :: output
  integer i,l,s,t,q
  integer nn,ne
  double precision thickness
  double precision,pointer,dimension(:,:) :: nodes
  integer,pointer,dimension(:,:) :: e_n
  integer dim,m
  double precision :: u(:)

  double precision,allocatable,dimension(:) :: ue,center_coord
  double precision,allocatable,dimension(:,:) :: B
  double precision,pointer,dimension(:,:) :: eps,sig
  double precision,pointer,dimension(:) :: msig,angle

  double precision a,detJ,D_3d(6,6)
  integer gn

  nn = model%n_nds
  ne = model%n_els
  dim = model%dim
  m = model%n_nds_1el

  nodes => model%nodes
  e_n => model%elements

  angle => model%data(1)%d(:,1,1)

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


   allocate(output%data(2))
   output%data(1)%d => model%data(1)%d ! 異方性の角度
   output%data(2)%d => model%data(2)%d ! 損傷テンソル

  do i=1,ne
    ue = 0d0

    do l=1,m
      ue( (l-1)*dim+1:l*dim ) = u( (e_n(l,i)-1)*dim+1:e_n(l,i)*dim )
    end do

    call calc_B(B,model,i,center_coord)

    eps(:,i) = matmul(B,ue)

!     call calc_D(D,model,rot_D(data(1)%d(:,:,i_el),data(2)%d(i_el,1,1)))
!   else
!     call calc_D(D,model,rot_D(model%materials(:,:,model%material_nos(i_el)),data(2)%d(i_el,1,1)))

    D_3d = od_set_damage_tensor(model,i) ! 損傷テンソルをDマトリックスへ適用

    D_3d = rot_D(D_3d,model%data(1)%d(i,1,1)) ! Dマトリックスの回転

    call calc_sig(output,model,i,D_3d)

    deallocate(B)

  end do


end subroutine



subroutine od_output_inp(model,output,file_name)
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
  allocate(el_data(17,model%n_els))

  allocate(nd_data_name(2))
  allocate(el_data_name(17))

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
  el_data_name(14) = 'damage_L'
  el_data_name(15) = 'damage_T'
  el_data_name(16) = 'damage_Z'
  el_data_name(17) = 'angle'

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
  el_data(14:16,:) = output%data(2)%d(:,:,1) ! 損傷テンソル
  el_data(17,:) = output%data(1)%d(:,1,1) ! 異方性の角度

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


end program fem

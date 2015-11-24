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
  integer i,step,max_step,n_model
  integer p,q,s,t
  character(8) :: i_char,step_char
  double precision,allocatable,dimension(:) :: k,u,f
  real(8),allocatable,target :: angle(:)
  real(8),allocatable,target :: damage_tensor(:,:)
  real(8) u_, du ! 損傷進展解析の変位、変位増分
  integer damage_ratio_step
!  logical damage_ratio
  real(8) damage_ratio
  integer damaged_elem_no,damaged_elem_mat_no
  real(8) f_left, f_right ! 反力の合計(左右)
  real(8),allocatable :: max_sig(:,:) ! 損傷の閾値
  real(8) :: vf !繊維束の体積含有率

  integer n_periods_x,n_periods_y,n_areas_1weft
  real(8),allocatable :: vf_warp(:,:),vf_weft(:,:,:)
  integer,allocatable :: mat_no_warp(:,:),mat_no_weft(:,:,:)
  integer,allocatable :: focused_elements(:,:,:)

  real(8) sd_vf_meso,exp_vf_meso,min_vf_meso,max_vf_meso,&
    &sd_vf_micro,exp_vf_micro,min_vf_micro,max_vf_micro,vf_min,vf_max,vf_interval
  integer mat_no_offset_warp,mat_no_offset_weft,mat_no_offset
  real(8) sig_rot(6)

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

  intrinsic :: command_argument_count, get_command_argument
  integer length,status
  character(:),allocatable :: model_no_str

  real(8),allocatable :: tmp_params(:)

!-------------------------------------------------------------------------------------------------------
!  Config
!-------------------------------------------------------------------------------------------------------

  du = 8.89d-6
  max_step = 10
  od_model_name = "gfrp_damage"

  n_periods_x = 2
  n_periods_y = 4
  n_areas_1weft = 24


  sd_vf_meso = 0d0
  exp_vf_meso = 0.6d0
  min_vf_meso = 0.6d0
  max_vf_meso = 0.6d0
  sd_vf_micro = 0.173723d0
  exp_vf_micro = 1d0
  min_vf_micro = 1-0.37076d0
  max_vf_micro = 1+0.288101d0
  vf_interval = 0.01d0

  mat_no_offset = 5


  vf_min = min_vf_meso*min_vf_micro
  vf_max = max_vf_meso*max_vf_micro

  allocate(mat_no_warp(n_periods_x*4,n_periods_y))
  allocate(mat_no_weft(n_areas_1weft,n_periods_x*2,n_periods_y))
  allocate(vf_warp(n_periods_x*4,n_periods_y))
  allocate(vf_weft(n_areas_1weft,n_periods_x*2,n_periods_y))
  allocate(max_sig(6,mat_no_offset-1+2*( floor(vf_max/vf_interval)-floor(vf_min/vf_interval) )))
  mat_no_warp = 0
  mat_no_weft = 0
  vf_warp = 0d0
  vf_weft = 0d0
  
!-------------------------------------------------------------------------------------------------------
!  Solver
!-------------------------------------------------------------------------------------------------------


  write(command,'(a)') "mkdir "//trim(path_model)//trim(od_model_name)
print *,"mkdir "//trim(path_model)//trim(od_model_name)
  call system(command)

! 引数からモデル番号を取得
i = 1 
if (command_argument_count() > 0) then
  call get_command_argument(1, length=length, status=status)
  allocate(character(length) :: model_no_str)
  call get_command_argument(1, model_no_str, status=status)
  read(model_no_str,*) i
end if

  if (i == 1) then
    open(11, file=trim(path_model)//trim(od_model_name)//slash//'models.csv')
    write(11,*) 'du, ',du
    write(11,*) 'model_no, u when damaged, damaged element No., material No. of damaged element'
  else
    open(11, file=trim(path_model)//trim(od_model_name)//slash//'models.csv',position='append')
  end if


  write(i_char, '(i0)') i

  write(*,'(a)') " ------------------------------"
  write(*,'(a,i0,a)') "  calculation of model ",i," start"
  write(*,'(a)') " ------------------------------"

  model%name = od_model_name

shift = 0

!  shift(1) = 0
!  shift(2) = 34
!  shift(3) = 50
!  shift(4) = 22

  shift(1) = 0
  shift(2) = 17
  shift(3) = 50
  shift(4) = 27

!  shift(2) = uniform_rand_integer(0,53)
!  shift(3) = uniform_rand_integer(0,53)
!  shift(4) = uniform_rand_integer(0,53)


! 外部ファイルからパラメータ値を読み込み
  !allocate(tmp_params(n_periods_x*4 * n_periods_y + n_areas_1weft * n_periods_x*2 * n_periods_y))
  !open(22, file=trim(path_model)//trim(model%name)//slash//'params.csv')
  !read(22,*) n_model
  !read(22,'()') ! ヘッダを読み飛ばし
  !do p=1,n_model
  !  read(22,*) q, tmp_params ! qは使わない
  !  if (q==i) then
  !    exit
  !  end if
  !end do
  !close(22)
  allocate(tmp_params(n_periods_x*4 * n_periods_y + n_areas_1weft * n_periods_x*2 * n_periods_y))
  open(22, file=trim(path_model)//trim(model%name)//slash//'params'//slash//trim(i_char)//'.csv')
  read(22,*) tmp_params ! qは使わない
  close(22)

  t = 1
  do p=1,n_periods_y
    do q=1,n_periods_x*4
      vf_warp(q,p) = tmp_params(t)
      t = t + 1
    end do
  end do
  do p=1,n_periods_y
    do q=1,n_periods_x*2
      do s=1,n_areas_1weft
        vf_weft(s,q,p) = tmp_params(t)
        t = t + 1
      end do
    end do
  end do

  ! ※材料番号84と44は問題があるのでそれを用いないようにVfのリミットを+-16%に設定
  ! warpのvf
!  do q=1,n_periods_y
!    do p=1,n_periods_x*4
!      vf_warp(p,q) = normal_rand_real8(exp_vf_meso,exp_vf_meso*sd_vf_micro,0.16d0)
!    end do
!  end do
!
!  ! weftのvf
!  do q=1,n_periods_y
!    do p=1,n_periods_x*2 
!      call normal_rand_vf_micro(vf_weft(:,p,q),exp_vf_meso,exp_vf_meso*sd_vf_micro,0.16d0)
!    end do
!  end do

  vf = normal_rand_real8(0.6d0,0.033d0)
  print *,'vf = ',vf

  call calc_vf_mat_no(mat_no_warp,mat_no_weft,mat_no_offset_warp,mat_no_offset_weft,&
    &vf_warp,vf_weft,vf_min,vf_max,vf_interval,n_periods_x,n_periods_y,n_areas_1weft,mat_no_offset)
!print *,vf_min,vf_max
!print *,mat_no_warp,mat_no_weft

call calc_strength(max_sig,vf_min,vf_max,vf_interval,mat_no_offset)
  !call calc_strength(max_sig,vf)

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


  call generate_FRP_Model(model,focused_elements,n_periods_x,n_periods_y,shift)
  call read_materials(model)
  call set_vf_params(model,n_periods_x,n_periods_y,shift,mat_no_warp,mat_no_weft)
!  print *,mat_no_warp


  call set_state2d(model,1,1d0)


  print *,'Making directory...'

  write(od_data_path,'(a,i0,a)') "model", i, slash
  write(command,'(a,a,a,a,a)') "mkdir ", trim(path_model), trim(model%name),slash, trim(od_data_path)
print *,"mkdir ", trim(path_model), trim(model%name),slash, trim(od_data_path)
  call system(command)

  ! モデルごとのデータ出力
  open(20, file=trim(path_model)//trim(model%name)//slash//trim(od_data_path)//'model'//trim(i_char)//'.csv')
  write(20,*) 'f_left,f_right'
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
  if (ceiling(1 / damage_ratio) > 1) then ! ステップ1で損傷が発生していない場合
    step = ceiling(1 / damage_ratio)
  end if
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
!call visualize_u(model,u)
    print *,"Calculating stress and strain..."; print *

    call od_calc_output(output,model,u,vf_min,vf_max,vf_interval,mat_no_offset)

    print *,"Judging the damage state..."; print *

    if (step == 1) then
      ! 応力集中部の応力を一覧で取得
      open(21, file=trim(path_model)//trim(model%name)//slash//trim(od_data_path)//&
                            &'model'//trim(i_char)//'_sigTZ_stress_concentration.csv')
      do p=1,ubound(focused_elements,1)
      do q=1,ubound(focused_elements,2)
      do s=1,ubound(focused_elements,3)
!        print *,p,q,s,model%data(1)%d(focused_elements(p,q,s),1,1),model%material_nos(focused_elements(p,q,s))
        sig_rot = rot_sig(output%sig(:,focused_elements(p,q,s)), model%data(1)%d(focused_elements(p,q,s),1,1))
        write(21,*) dabs(sig_rot(6))
      end do
      end do
      end do
      close(21)
    end if

    ! 損傷判定
    if (damage_ratio>=1) then ! すでに損傷が発生している場合
      call damage_judgment(damage_ratio,damaged_elem_no,damaged_elem_mat_no,&
        &model,output,max_sig,mat_no_offset_warp,mat_no_offset_weft)
print *,'already damaged'
    else ! まだ損傷していない場合
      call damage_judgment(damage_ratio,damaged_elem_no,damaged_elem_mat_no,&
        &model,output,max_sig,mat_no_offset_warp,mat_no_offset_weft)
      if (damage_ratio>=1) then
!        write(11,'(i0,7(",",d30.15),5(",",i0))') &
!          &i, vf, max_sig(1,2), max_sig(2,2), max_sig(6,2), max_sig(1,3), max_sig(2,3), max_sig(6,3),&
!            &shift(1), shift(2), shift(3), shift(4), step
        write(11,'(i0,1(",",f30.15),2(",",i0))') &
          &i, u_/damage_ratio, damaged_elem_no, damaged_elem_mat_no
print *, 'damaged!!! step',step

      else
        print *,'undamaged'
      end if
    end if

    print *,"Outputting file..."; print *

    write(output_file_name,'(a,i0)') "step", step
     call od_output_inp(model,output,trim(od_data_path)//output_file_name,vf_min,vf_max,vf_interval,mat_no_offset,max_sig,&
                                                                                                          &focused_elements)

    print *,"Clearing data..."; print *

    call clear_output(output)
    call clear_bc(bc)
    call clear_addition_matrix_sln(k,index_k_diag)

    if (damage_ratio >= 1d0) then ! 損傷していたらループを抜ける
      exit
    end if
  end do

  close(20)
  call clear_model(model)

  deallocate(f)
  deallocate(u)
  deallocate(focused_elements)

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


!  read *

contains

subroutine calc_strength(max_sig,vf_min,vf_max,vf_interval,mat_no_offset)
  implicit none
  integer i
  integer mat_no_offset,vf_num
  real(8) vf_min,vf_max,vf_interval,vf
  real(8) :: max_sig(:,:)

  vf_num = floor(vf_max/vf_interval) - floor(vf_min/vf_interval)

  ! warp
  do i=mat_no_offset,(mat_no_offset+vf_num-1)
    vf = ceiling(vf_min/vf_interval+(i-mat_no_offset)) * vf_interval
!    max_sig(1,i) = 1000.3d6+(vf-0.6d0)*100d6 ! I
!    max_sig(2,i) = 34.3d6+(vf-0.6d0)*(-18.5d6) ! II
!    max_sig(6,i) = 40.2d6+(vf-0.6d0)*(-21.8d6) ! I II
    max_sig(1,i) = 1000.3d6/0.4d0*(1d0-vf)
    max_sig(2,i) = 272.306d6*(1d0-2*dsqrt(vf/pi)) ! II
    max_sig(6,i) = 319.146d6*(1d0-2*dsqrt(vf/pi)) ! II III
  end do

  ! weft
  do i=(mat_no_offset+vf_num),(mat_no_offset+vf_num*2-1)
    vf = ceiling(vf_min/vf_interval+(i-mat_no_offset-vf_num)) * vf_interval
!    max_sig(2,i) = 34.3d6+(vf-0.6d0)*(-18.5d6) ! II
!    max_sig(3,i) = 34.3d6+(vf-0.6d0)*(-18.5d6) ! III 
!    max_sig(4,i) = 40.2d6+(vf-0.6d0)*(-21.8d6) ! II III
    max_sig(2,i) = 272.306d6*(1d0-2*dsqrt(vf/pi)) ! II
    max_sig(3,i) = 272.306d6*(1d0-2*dsqrt(vf/pi)) ! II
    max_sig(4,i) = 319.146d6*(1d0-2*dsqrt(vf/pi)) ! II III
  end do

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

subroutine damage_judgment(damage_ratio,damaged_elem_no,damaged_elem_mat_no,&
    &model,output,max_sig,mat_no_offset_warp,mat_no_offset_weft)

  type(struct_model) :: model
  type(struct_output) :: output
  integer n_els,i,mat_no,mat_no_offset_warp,mat_no_offset_weft
  real(8),pointer :: sig(:,:),angle(:),damage_tensor(:,:)
  real(8) :: max_sig(:,:),sig_rot(6)
  real(8) damage_ratio,damage_ratio_ ! damage_ratio: sig/maxsig 1以上で損傷発生
  integer damaged_elem_no,damaged_elem_mat_no

  n_els = model%n_els
  angle => model%data(1)%d(:,1,1)
  sig => output%sig

  damage_tensor => model%data(2)%d(:,:,1)

  damage_ratio = 0d0

  do i=1,n_els
    sig_rot = rot_sig(sig(:,i),angle(i))

    mat_no = model%material_nos(i)
    ! 縦糸
    if ((mat_no_offset_warp <= mat_no) .and. (mat_no < mat_no_offset_weft)) then ! x=>L y:T z:Z
      if (max_sig(1,mat_no) < sig_rot(1)) then
        damage_tensor(1,i) = 0.99d0
      else if (max_sig(2,mat_no) < sig_rot(2)) then
        damage_tensor(2,i) = 0.99d0
      else if (max_sig(6,mat_no) < sig_rot(6)) then
        damage_tensor(2,i) = 0.99d0
      end if

      damage_ratio_ = dabs(sig_rot(1)/max_sig(1,mat_no))
      if (damage_ratio_ > damage_ratio) then
        damage_ratio = damage_ratio_
        damaged_elem_no = i
        damaged_elem_mat_no = 2
      end if
      damage_ratio_ = dabs(sig_rot(2)/max_sig(2,mat_no))
      if (damage_ratio_ > damage_ratio) then
        damage_ratio = damage_ratio_
        damaged_elem_no = i
        damaged_elem_mat_no = 2
      end if
      damage_ratio_ = dabs(sig_rot(6)/max_sig(6,mat_no))
      if (damage_ratio_ > damage_ratio) then
        damage_ratio = damage_ratio_
        damaged_elem_no = i
        damaged_elem_mat_no = 2
      end if

    ! 横糸
    else if (mat_no_offset_weft <= mat_no) then ! x:T y:Z z:L
      if (max_sig(2,mat_no) < sig_rot(1)) then
        damage_tensor(2,i) = 0.99d0
      else if (max_sig(3,mat_no) < sig_rot(2)) then
        damage_tensor(3,i) = 0.99d0
      else if (max_sig(4,mat_no) < sig_rot(6)) then
        damage_tensor(2,i) = 0.99d0
        damage_tensor(3,i) = 0.99d0
      end if


      damage_ratio_ = dabs(sig_rot(1)/max_sig(2,mat_no))
      if (damage_ratio_ > damage_ratio) then
        damage_ratio = damage_ratio_
        damaged_elem_no = i
        damaged_elem_mat_no = 3
      end if
      damage_ratio_ = dabs(sig_rot(2)/max_sig(3,mat_no))
      if (damage_ratio_ > damage_ratio) then
        damage_ratio = damage_ratio_
        damaged_elem_no = i
        damaged_elem_mat_no = 3
      end if
      damage_ratio_ = dabs(sig_rot(6)/max_sig(4,mat_no))
      if (damage_ratio_ > damage_ratio) then
        damage_ratio = damage_ratio_
        damaged_elem_no = i
        damaged_elem_mat_no = 3
      end if


    end if


  end do
print *,"damage_judgment = ", damage_ratio
end subroutine


function uniform_rand_integer(min,max)
  implicit none
  integer min, max, interval
  real(8) uniform_rand_integer, genrand_res53

  interval = abs(max-min)+1

  uniform_rand_integer = int(genrand_res53()*interval)
end function


function normal_rand_real8(ave,sd,arg_limit)
  implicit none
  real(8) normal_rand_real8,x,y,z,pi,genrand_res53,ave,sd,limit
  real(8),optional :: arg_limit

  pi = dacos(-1d0)

  if (present(arg_limit)) then
    limit = arg_limit
  else
    limit = 6d0*sd
  end if

  do
    x = genrand_res53()
    y = genrand_res53()
    z = sd*dsqrt(-2d0*dlog(x)) * dCOS(2*pi*y)
    if (abs(z)<=limit) exit
  end do

  normal_rand_real8 = z+ave
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

  D_3d = od_set_damage_tensor(model,i_el,vf_min,vf_max,vf_interval,mat_no_offset) ! 損傷テンソルをDマトリックスへ適用

  D_3d = rot_D(D_3d,model%data(1)%d(i_el,1,1)) ! Dマトリックスの回転

  call calc_D(D,model,D_3d)

  call calc_B(B,model,i_el,coord)


  od_calc_BDB(:,:) = matmul( matmul(transpose(B),D),B)


end function

!subroutine calc_strength(max_sig,vf_min,vf_max,vf_interval,mat_no_offset)
!  implicit none
!  integer i
!  integer mat_no_offset,vf_num
!  real(8) vf_min,vf_max,vf_interval,vf
!  real(8) :: max_sig(:,:)
!
!  vf_num = floor(vf_max/vf_interval) - floor(vf_min/vf_interval)
!  do i=mat_no_offset,(mat_no_offset+vf_num-1)
!    vf = ceiling(vf_min/vf_interval+(i-mat_no_offset)) * vf_interval
!    max_sig(1,i) = 1000.3d6+(vf-0.6d0)*100d6 ! I
!    max_sig(2,i) = 34.3d6+(vf-0.6d0)*(-18.5d6) ! II
!    max_sig(6,i) = 40.2d6+(vf-0.6d0)*(-21.8d6) ! I II
!  end do
!
!  ! weft
!  do i=(mat_no_offset+vf_num),(mat_no_offset+vf_num*2-1)
function od_set_damage_tensor(model,i_el,vf_min,vf_max,vf_interval,mat_no_offset)
  type(struct_model) :: model
  real(8) :: dl,dt,dz
  real(8) :: od_set_damage_tensor(6,6)
  integer mat_no_offset,vf_num,mat_no
  real(8) vf_min,vf_max,vf_interval,vf
  integer i_el

  vf_num = floor(vf_max/vf_interval) - floor(vf_min/vf_interval)

  dl = 1d0-model%data(2)%d(1,i_el,1)
  dt = 1d0-model%data(2)%d(2,i_el,1)
  dz = 1d0-model%data(2)%d(3,i_el,1)

  od_set_damage_tensor(1:6,1:6) = model%materials(1:6,1:6,model%material_nos(i_el)) 

  mat_no = model%material_nos(i_el)
  if ((mat_no_offset<=mat_no) .and. (mat_no<=(mat_no_offset+vf_num-1))) then
!  if (model%material_nos(i_el) == 2) then
    od_set_damage_tensor(1,:) = od_set_damage_tensor(1,:) * dl
    od_set_damage_tensor(2,:) = od_set_damage_tensor(2,:) * dt
    od_set_damage_tensor(3,:) = od_set_damage_tensor(3,:) * dz
    od_set_damage_tensor(:,1) = od_set_damage_tensor(:,1) * dl
    od_set_damage_tensor(:,2) = od_set_damage_tensor(:,2) * dt
    od_set_damage_tensor(:,3) = od_set_damage_tensor(:,3) * dz
    od_set_damage_tensor(4,4) = od_set_damage_tensor(4,4) * 4d0 * (dt*dz/(dt+dz))**2
    od_set_damage_tensor(5,5) = od_set_damage_tensor(5,5) * 4d0 * (dz*dl/(dz+dl))**2
    od_set_damage_tensor(6,6) = od_set_damage_tensor(6,6) * 4d0 * (dl*dt/(dl+dt))**2
!  else if (model%material_nos(i_el) == 3) then
  else if (((mat_no_offset+vf_num)<=mat_no) .and. (mat_no<=(mat_no_offset+vf_num*2-1))) then
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
    sig_tensor_rot(i,j) = sig_tensor_rot(i,j)+beta(i,k)*beta(j,l)*sig_tensor(k,l);
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


subroutine od_calc_output(output,model,u,vf_min,vf_max,vf_interval,mat_no_offset)
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

  double precision vf_min,vf_max,vf_interval
  integer mat_no_offset

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

    D_3d = od_set_damage_tensor(model,i,vf_min,vf_max,vf_interval,mat_no_offset) ! 損傷テンソルをDマトリックスへ適用

    D_3d = rot_D(D_3d,model%data(1)%d(i,1,1)) ! Dマトリックスの回転

    call calc_sig(output,model,i,D_3d)

    deallocate(B)

  end do


end subroutine



subroutine od_output_inp(model,output,file_name,vf_min,vf_max,vf_interval,mat_no_offset,max_sig,focused_elements)
  type(struct_model) :: model
  type(struct_output) :: output
  character(*) file_name
  integer i,j,k
  double precision,allocatable,dimension(:,:) :: nd_data
  double precision,allocatable,dimension(:,:) :: el_data
  character(16),allocatable,dimension(:) :: nd_data_name
  character(16),allocatable,dimension(:) :: el_data_name
  real(8) vf_min,vf_max,vf_interval
  integer vf_num,mat_no,mat_no_offset
  integer,allocatable :: tmp_material_nos(:)
  real(8) sig_rot(6)
  integer :: focused_elements(:,:,:)

  real(8) :: max_sig(:,:)


  if (model%dim == 2) then

  allocate(nd_data(2,model%n_nds))
  allocate(el_data(40,model%n_els))

  allocate(nd_data_name(2))
  allocate(el_data_name(40))

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
  el_data_name(18) = 'mat_no'
  el_data_name(19) = 'max_sig_L'
  el_data_name(20) = 'max_sig_T'
  el_data_name(21) = 'max_sig_Z'
  el_data_name(22) = 'max_sig_TZ'
  el_data_name(23) = 'max_sig_ZL'
  el_data_name(24) = 'max_sig_LT'
  el_data_name(25) = 'Vf'
  el_data_name(26) = 'sigI'
  el_data_name(27) = 'sigII'
  el_data_name(28) = 'sigIII'
  el_data_name(29) = 'sigII III'
  el_data_name(30) = 'sigIII I'
  el_data_name(31) = 'sigI II'
  el_data_name(32) = 'focused element'
  el_data_name(33) = 'damage_ratioT'
  el_data_name(34) = 'damage_ratioTZ'
  el_data_name(35) = 'abs(sigL)'
  el_data_name(36) = 'abs(sigT)'
  el_data_name(37) = 'abs(sigZ)'
  el_data_name(38) = 'abs(sigTZ)'
  el_data_name(39) = 'abs(sigZL)'
  el_data_name(40) = 'abs(sigLT)'

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
  el_data(14:16,:) = model%data(2)%d(:,:,1) ! 損傷テンソル
  el_data(17,:) = model%data(1)%d(:,1,1) ! 異方性の角度
  el_data(18,:) = model%material_nos(:) ! 物性値

  el_data(19:24,:) = 0d0

  vf_num = floor(vf_max/vf_interval) - floor(vf_min/vf_interval)
  allocate(tmp_material_nos(model%n_els))
  tmp_material_nos(:) = model%material_nos(:)
  do i=1,model%n_els
    mat_no = model%material_nos(i)
    if ((mat_no_offset<=mat_no) .and. (mat_no<=(mat_no_offset+vf_num-1))) then ! 横糸の場合
      el_data(19,i) = max_sig(1,mat_no)
      el_data(20,i) = max_sig(2,mat_no)
      el_data(21,i) = max_sig(3,mat_no)
      el_data(22,i) = max_sig(4,mat_no)
      el_data(23,i) = max_sig(5,mat_no)
      el_data(24,i) = max_sig(6,mat_no)
      el_data(25,i) = vf_interval*(ceiling(vf_min/vf_interval)+model%material_nos(i)-mat_no_offset)

      sig_rot = rot_sig(output%sig(:,i), model%data(1)%d(i,1,1))
      el_data(26:31,i) = sig_rot
      el_data(35:40,i) = dabs(sig_rot)

      if (el_data(20,i) > 0d0) then
        el_data(33,i) = dabs(el_data(27,i)/el_data(20,i))
      end if
      if (el_data(22,i) > 0d0) then
        el_data(34,i) = dabs(el_data(29,i)/el_data(22,i))
      end if

      model%material_nos(i) = 2 
    else if (((mat_no_offset+vf_num)<=mat_no) .and. (mat_no<=(mat_no_offset+vf_num*2-1))) then ! 縦糸の場合
      el_data(19,i) = max_sig(1,mat_no)
      el_data(20,i) = max_sig(2,mat_no)
      el_data(21,i) = max_sig(3,mat_no)
      el_data(22,i) = max_sig(4,mat_no)
      el_data(23,i) = max_sig(5,mat_no)
      el_data(24,i) = max_sig(6,mat_no)
      el_data(25,i) = vf_interval*(ceiling(vf_min/vf_interval)+model%material_nos(i)-mat_no_offset-vf_num)

      sig_rot = rot_sig(output%sig(:,i), model%data(1)%d(i,1,1))
      el_data(26:31,i) = (/sig_rot(3),sig_rot(1),sig_rot(2),sig_rot(6),sig_rot(4),sig_rot(5)/)
      el_data(35:40,i) = dabs((/sig_rot(3),sig_rot(1),sig_rot(2),sig_rot(6),sig_rot(4),sig_rot(5)/))

      if (el_data(20,i) > 0d0) then
        el_data(33,i) = dabs(el_data(27,i)/el_data(20,i))
      end if
      if (el_data(22,i) > 0d0) then
        el_data(34,i) = dabs(el_data(29,i)/el_data(22,i))
      end if

      model%material_nos(i) = 3
    end if
  end do

  
  el_data(32,:) = 0d0
  do i=1,ubound(focused_elements,1)
  do j=1,ubound(focused_elements,2)
  do k=1,ubound(focused_elements,3)
    el_data(32,focused_elements(i,j,k)) = 1d0
  end do
  end do
  end do


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

  call make_inp(nd_data,nd_data_name,el_data,el_data_name,model,file_name,10)

  if (model%dim == 2) then
    model%material_nos(:) = tmp_material_nos(:)
  end if

end subroutine

subroutine calc_vf_mat_no(mat_no_warp,mat_no_weft,mat_no_offset_warp,mat_no_offset_weft,&
    &vf_warp,vf_weft,vf_min,vf_max,vf_interval,n_periods_x,n_periods_y,n_areas_1weft,mat_no_offset)
  real(8) vf_min,vf_max,vf_interval
  real(8) vf_warp(:,:),vf_weft(:,:,:)
  integer i,j,k
  integer n_periods_x,n_periods_y,mat_no_offset_warp,mat_no_offset_weft,mat_no_offset,n_areas_1weft
  integer mat_no_warp(:,:), mat_no_weft(:,:,:)

  mat_no_offset_warp = mat_no_offset
  mat_no_offset_weft = mat_no_offset+floor(vf_max/vf_interval)-floor(vf_min/vf_interval)

  ! warp
  do i=1,n_periods_y
    do j=1,n_periods_x*4
      mat_no_warp(j,i) = nint(vf_warp(j,i)/vf_interval)-ceiling(vf_min/vf_interval)+mat_no_offset_warp
      if (mat_no_warp(j,i) < mat_no_offset_warp) then
        mat_no_warp(j,i) = mat_no_offset_warp
      else if (mat_no_warp(j,i) > floor(vf_max/vf_interval)-ceiling(vf_min/vf_interval)+mat_no_offset_warp) then
        mat_no_warp(j,i) = floor(vf_max/vf_interval)-ceiling(vf_min/vf_interval)+mat_no_offset_warp
      end if
    end do
  end do

  ! weft
  do i=1,n_periods_y
    do j=1,n_periods_x*2
      do k=1,n_areas_1weft
        mat_no_weft(k,j,i) = nint(vf_weft(k,j,i)/vf_interval)-ceiling(vf_min/vf_interval)+mat_no_offset_weft
        if (mat_no_weft(k,j,i) < mat_no_offset_weft) then
          mat_no_weft(k,j,i) = mat_no_offset_weft
        else if (mat_no_weft(k,j,i) > floor(vf_max/vf_interval)-ceiling(vf_min/vf_interval)+mat_no_offset_weft) then
          mat_no_weft(k,j,i) = floor(vf_max/vf_interval)-ceiling(vf_min/vf_interval)+mat_no_offset_weft
        end if
      end do
    end do
  end do

end subroutine

subroutine read_materials(model)
  type(struct_model) :: model
  real(8),pointer :: materials(:,:,:)
  integer n_mat,mat_no
  integer i,j,k

  open(20, file=trim(path_model)//trim(od_model_name)//slash//'mats.csv')
	read (20,*) n_mat

 	allocate(model%materials(6,6,n_mat))

  materials => model%materials

	do i=1,n_mat
		read (20,*) mat_no
		do j=1,6
			read (20,*) materials(:,j,mat_no)
		end do
	end do

end subroutine

subroutine normal_rand_vf_micro(vf_micro,ave,sd,limit)
  implicit none
  integer i,n_params
  real(8) ave,sd,limit,vf_micro(:)
  real(8),allocatable :: break(:)
  n_params = ubound(vf_micro,1)

  allocate(break(n_params))

  do i=1,n_params
    break(i) = normal_rand_real8(ave,sd/dsqrt(2d0),limit/2d0)
  end do

  do i=1,n_params-1
    vf_micro(i) = ave + break(i+1) - break(i)
  end do
  vf_micro(n_params) = ave + break(1) - break(n_params)

end subroutine


end program fem

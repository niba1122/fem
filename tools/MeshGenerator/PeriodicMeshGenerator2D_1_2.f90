program PeriodicMeshGenerator2D


! 2015/01/09
!  周期境界条件を出力出来るよう改良
!  
! 2014/12/31
!  ソルバーの改良に合わせて出力フォーマットを変更
! 
! 2014/12/24
! 層間移動用の境界条件を出力できるように改良
! 

	use gfrpMeshGenerator2D
	implicit none
	integer i,j,k,p,q
	double precision,allocatable,dimension(:,:) :: nodes,newNodes,nodes_,nodesOfLamina,NodesOfShiftedLamina
	integer,allocatable,dimension(:,:) :: topNodesOfLamina,bottomNodesOfLamina
	integer,allocatable,dimension(:,:) :: elements,newElements,elementsOfLamina

!	Parameters of one period
	integer,parameter :: NumOfNodesInOneElement = 4
	integer, parameter :: dim = 2

!	integer :: NumOfNodesInOnePeriod = 787
 	integer NumOfNodesInOnePeriod
	integer,parameter :: NumOfElementsInOnePeriod = 720
	integer,parameter :: NumOfNodesInOnePeriodX  = 55

	double precision,parameter :: ThicknessOfLamina = 0.6d-3
	double precision,parameter :: thicknessOfInterlaminar = 0.02d-3
	double precision :: thL,tlL,wL,thC,tlC,wC,dxC,thR,tlR,wR

	integer lastNode,lastNode_,lastElement,NumOfPeriodsX,NumOfPeriodsY,NumOfNodes,NumOfElements
	integer NumOfNodesInOnePeriodY,NumOfNodesInOneLamina,NumOfElementsInOneLamina

	integer shift

	integer numOfConstraint,NumOfElementsInOneLaminaX
	double precision ux


	NumOfPeriodsX = 2
	NumOfPeriodsY = 2

	thL = 0.11d-3
	thC = 0.11d-3
	thR = 0.11d-3
	wL = 4.04d-3
	wC = 4.04d-3
	wR = 4.04d-3
	dxC = 0.00d-3

	shift = 0



	lastNode = 0
	lastElement = 0


	allocate(nodes(1,1))
	allocate(nodesOfLamina(1,1))
	allocate(elements(NumOfNodesInOneElement+1,NumOfElementsInOnePeriod*NumOfPeriodsX*NumOfPeriodsY+&
												&(NumOfNodesInOnePeriodX-1)*NumOfPeriodsX*(NumOfPeriodsY-1)))
	allocate(newNodes(1,1))
	allocate(nodes_(1,1))
	allocate(nodesOfShiftedLamina(1,1))
	allocate(newElements(NumOfNodesInOneElement+1,NumOfElementsInOnePeriod))
	allocate(elementsOfLamina(NumOfNodesInOneElement+1,NumOfElementsInOnePeriod*NumOfPeriodsX))
	allocate(topNodesOfLamina((NumOfNodesInOnePeriodX-1)*NumOfPeriodsX+1,NumOfPeriodsY))
	allocate(bottomNodesOfLamina((NumOfNodesInOnePeriodX-1)*NumOfPeriodsX+1,NumOfPeriodsY))

	nodes = 0d0
	nodes_ = 0d0


	NumOfElementsInOneLamina = NumOfElementsInOnePeriod*NumOfPeriodsX


	do k=1,NumOfPeriodsY
		lastNode_ = lastNode

		! 1周期目の生成
		call generateGFRPMesh(newNodes,newElements,NumOfNodesInOnePeriod,thL,wL,thC,wC,dxC,thR,wR,shift)

		!　y方向の節点数のカウント
		NumOfNodesInOnePeriodY = 0
		do i=1,NumOfNodesInOnePeriod
			NumOfNodesInOnePeriodY = NumOfNodesInOnePeriodY + 1
			if (dabs(newNodes(2,i)-ThicknessOfLamina)<1d-9) exit
		end do
		! 層内の節点数の計算
		NumOfNodesInOneLamina = NumOfNodesInOnePeriod*NumOfPeriodsX-NumOfNodesInOnePeriodY*(NumOfPeriodsX-1)

		! 層内節点の配列の初期化

		deallocate(nodesOfShiftedLamina)
		allocate(nodesOfShiftedLamina(dim,NumOfNodesInOneLamina))
		deallocate(nodesOfLamina)
		allocate(nodesOfLamina(dim,NumOfNodesInOneLamina))


		nodesOfLamina(:,:) = 0d0
		nodesOfShiftedLamina(:,:) = 0d0

		! 層内節点の配列に1周期目を追加
		nodesOfLamina(:,1:NumOfNodesInOnePeriod) = newNodes(:,:)

		! 要素の配列に1周期目を追加
		elements(1:4,lastElement+1:lastElement+NumOfElementsInOnePeriod) = newElements(1:4,:) + lastNode
		elements(5,lastElement+1:lastElement+NumOfElementsInOnePeriod) = newElements(5,:)

		! 最後の節点・要素番号の変数を更新
		lastNode = lastNode + NumOfNodesInOnePeriod
		lastElement = lastElement + NumOfElementsInOnePeriod

		! 2周期目以降
		do i=2,NumOfPeriodsX
			! i周期目の生成

			call generateGFRPMesh(newNodes,newElements,NumOfNodesInOnePeriod,thL,wL,thC,wC,dxC,thR,wR,shift)

			! 層内節点の配列にi周期目を追加
			nodesOfLamina(1,lastNode-lastNode_+1:lastNode-lastNode_+NumOfNodesInOnePeriod-NumOfNodesInOnePeriodY) = &
				&newNodes(1,NumOfNodesInOnePeriodY+1:NumOfNodesInOnePeriod)+nodesOfLamina(1,lastNode-lastNode_)
			nodesOfLamina(2,lastNode-lastNode_+1:lastNode-lastNode_+NumOfNodesInOnePeriod-NumOfNodesInOnePeriodY) = &
				&newNodes(2,NumOfNodesInOnePeriodY+1:NumOfNodesInOnePeriod)


! print *,(i-1)*(NumOfNodesInOnePeriod-NumOfNodesInOnePeriodY)+NumOfNodesInOnePeriodY+1,",",&
! 				i*(NumOfNodesInOnePeriod-NumOfNodesInOnePeriodY)+NumOfNodesInOnePeriodY

			! 要素の配列にi周期目を追加
			elements(1:4,lastElement+1:lastElement+NumOfElementsInOnePeriod) = newElements(1:4,:) + lastNode - NumOfNodesInOnePeriodY
			elements(5,lastElement+1:lastElement+NumOfElementsInOnePeriod) = newElements(5,:)

			! 最後の節点・要素番号の変数を更新
			lastNode = lastNode + NumOfNodesInOnePeriod - NumOfNodesInOnePeriodY
			lastElement = lastElement + NumOfElementsInOnePeriod
		end do



!　シフト無しのコード
!
! 		! 層の上端・下端の節点番号を取得・保存
! 		p = 1
! 		q = 1
! 		do i=1,NumOfNodesInOneLamina
! 			!　下端
! 			if (dabs(nodesOfLamina(2,i)-0d0) < 1d-9) then
! 				bottomNodesOfLamina(p,k) = i
! 				p = p + 1
! 			end if
! 			! 上端
! 			if (dabs(nodesOfLamina(2,i)-ThicknessOfLamina) < 1d-9) then
! 				topNodesOfLamina(q,k) = i
! 				q = q + 1
! 			end if
! 		end do

!  		!　層の移動

! print *,"1:",NumOfNodesInOneLamina-topNodesOfLamina(shift,k)," = "&
! 							,topNodesOfLamina(shift,k)+1,":",NumOfNodesInOneLamina

! 		nodesOfShiftedLamina(1,1:NumOfNodesInOneLamina-topNodesOfLamina(shift,k)) = &
! 			nodesOfLamina(1,topNodesOfLamina(shift,k)+1:NumOfNodesInOneLamina)-nodesOfLamina(1,topNodesOfLamina(shift,k)+1)
! 		nodesOfShiftedLamina(2,1:NumOfNodesInOneLamina-topNodesOfLamina(shift,k)) = &
! 							nodesOfLamina(2,topNodesOfLamina(shift,k)+1:NumOfNodesInOneLamina)
! !-------------------------------------------------------------------------------------------------------


! print *,"1:",NumOfNodesInOneLamina-topNodesOfLamina(shift,k)
! 		nodesOfShiftedLamina(1,NumOfNodesInOneLamina-topNodesOfLamina(shift,k)+1:NumOfNodesInOneLamina)&
! 				= nodesOfLamina(1,topNodesOfLamina(1,k)+1:topNodesOfLamina(shift+1,k))&
! 					-nodesOfLamina(1,topNodesOfLamina(shift,k)+1) + nodesOfLamina(1,NumOfNodesInOneLamina)
! 		nodesOfShiftedLamina(2,NumOfNodesInOneLamina-topNodesOfLamina(shift,k)+1:NumOfNodesInOneLamina) = &
! 							nodesOfLamina(2,topNodesOfLamina(1,k)+1:topNodesOfLamina(shift+1,k))


! print *, nodesOfShiftedLamina(2,1:200)
! 		! 層の上端・下端の節点番号を取得・保存
! 		p = 1
! 		q = 1
! 		do i=1,NumOfNodesInOneLamina
! 			!　下端
! 			if (dabs(nodesOfShiftedLamina(2,i)-0d0) < 1d-9) then
! 				bottomNodesOfLamina(p,k) = i+lastNode_
! 				p = p + 1
! 			end if
! 			! 上端
! 			if (dabs(nodesOfShiftedLamina(2,i)-ThicknessOfLamina) < 1d-9) then
! 				topNodesOfLamina(q,k) = i+lastNode_
! 				q = q + 1
! 			end if
! 		end do

! 		! 節点の配列のバックアップ配列を作成
! 		deallocate(nodes_)
! 		allocate(nodes_(dim,lastNode))

! 		!　バックアップ配列に節点データを格納
! 		nodes_(:,:) = nodes(:,:)

! 		! 節点の配列の再初期化
! 		deallocate(nodes)
! 		allocate(nodes(dim,lastNode))

! 		! 節点の配列にバックアップ配列のデータを格納
! 		nodes(:,1:lastNode_) = nodes_(:,:)

! 		! 層内節点のデータを節点の配列に追加
! 		nodesOfShiftedLamina(2,:) = nodesOfShiftedLamina(2,:) + (ThicknessOfLamina+thicknessOfInterlaminar)*(k-1)
! 		nodes


		! 層の上端・下端の節点番号を取得・保存
		p = 1
		q = 1
		do i=1,NumOfNodesInOneLamina
			!　下端
			if (dabs(nodesOfLamina(2,i)-0d0) < 1d-9) then
				bottomNodesOfLamina(p,k) = i+lastNode_
				p = p + 1
			end if
			! 上端
			if (dabs(nodesOfLamina(2,i)-ThicknessOfLamina) < 1d-9) then
				topNodesOfLamina(q,k) = i+lastNode_
				q = q + 1
			end if
		end do

! 		print *,"ok1"
!  print *,topNodesOfLamina(:,k)
!  print *,"ok2"

		! 節点の配列のバックアップ配列を作成
		deallocate(nodes_)
		allocate(nodes_(dim,lastNode_))

		!　バックアップ配列に節点データを格納
		nodes_(:,:) = nodes(:,:)

		! 節点の配列の再初期化
		deallocate(nodes)
		allocate(nodes(dim,lastNode))

		! 節点の配列にバックアップ配列のデータを格納
		nodes(:,1:lastNode_) = nodes_(:,:)

		! 層内節点のデータを節点の配列に追加
		nodesOfLamina(2,:) = nodesOfLamina(2,:) + (ThicknessOfLamina+thicknessOfInterlaminar)*(k-1)
		nodes(:,lastNode_+1:lastNode) = nodesOfLamina(:,:)
! print *, bottomNodesOfLamina(:,k)
! print *, topNodesOfLamina(:,k)
shift = shift + 27
if (shift >= 54) shift = shift-54


	end do

	! 層間要素の定義
	do k=2,NumOfPeriodsY
		do i=1,(NumOfNodesInOnePeriodX-1)*NumOfPeriodsX
			elements(:,lastElement+1) = &
				&(/topNodesOfLamina(i,k-1),topNodesOfLamina(i+1,k-1),bottomNodesOfLamina(i+1,k),bottomNodesOfLamina(i,k),1/)
			lastElement = lastElement + 1
		end do
	end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	open(10,file='model.csv')
	write (10,*) "quad4"
	write (10,*) lastNode, ',',lastElement
	do i=1,lastNode
		write(10,"(i0, 3(',', e30.15))") i, nodes(:,i)
	end do

	do i=1,lastElement
		write(10,"(i0, 100(',', i0))") i, elements(5,i), elements(1:4,i)
	end do
	close(10)


	ux = 1d-5
	open(10,file='cn.csv')
	!　拘束条件を与える節点のカウント
! 	numOfConstraint=0
! 	do i=1,NumOfPeriodsY
! 		numOfConstraint = numOfConstraint + topNodesOfLamina(1,i)-bottomNodesOfLamina(1,i)+1
! 	end do
! 	numOfConstraint = numOfConstraint + (NumOfNodesInOnePeriodX-1)*NumOfPeriodsX
! 	do i=1,NumOfPeriodsY
! 		numOfConstraint = numOfConstraint + topNodesOfLamina(1,i)-bottomNodesOfLamina(1,i)+1
! 	end do
! 	numOfConstraint = numOfConstraint - 1

! 	write(10,"(i0)") numOfConstraint
! 	do i=1,NumOfPeriodsY
! 		do j=bottomNodesOfLamina(1,i),topNodesOfLamina(1,i)
! 			if (j==bottomNodesOfLamina(1,1)) then
! 				write(10,"(i0,',',i0,',',i0,',',e30.15,',',e30.15)") j,1,1,0d0,0d0
! 			else
! 				write(10,"(i0,',',i0,',',i0,',',e30.15,',',e30.15)") j,1,0,0d0,0d0
! 			end if
! 		end do
! 	end do
! 	do i=2,(NumOfNodesInOnePeriodX-1)*NumOfPeriodsX
! 		write(10,"(i0,',',i0,',',i0,',',e30.15,',',e30.15)") bottomNodesOfLamina(i,1),0,1,0d0,0d0
! 	end do
! 	NumOfElementsInOneLaminaX = (NumOfNodesInOnePeriodX-1)*NumOfPeriodsX+1
! 	do i=1,NumOfPeriodsY
! 		do j=bottomNodesOfLamina(NumOfElementsInOneLaminaX,i),topNodesOfLamina(NumOfElementsInOneLaminaX,i)
! 			if (j==bottomNodesOfLamina(NumOfElementsInOneLaminaX,1)) then
! 				write(10,"(i0,',',i0,',',i0,',',e30.15,',',e30.15)") j,1,1,ux,0d0
! 			else
! 				write(10,"(i0,',',i0,',',i0,',',e30.15,',',e30.15)") j,1,0,ux,0d0
! 			end if
! 		end do
! 	end do

	numOfConstraint=0
	do i=1,NumOfPeriodsY
		numOfConstraint = numOfConstraint + topNodesOfLamina(1,i)-bottomNodesOfLamina(1,i)+1
	end do
	do i=1,NumOfPeriodsY
		numOfConstraint = numOfConstraint + topNodesOfLamina(1,i)-bottomNodesOfLamina(1,i)+1
	end do

	write(10,"(i0)") numOfConstraint
	do i=1,NumOfPeriodsY
		do j=bottomNodesOfLamina(1,i),topNodesOfLamina(1,i)
			if (j==bottomNodesOfLamina(1,1)) then
				write(10,"(i0,',',i0,',',i0,',',e30.15,',',e30.15)") j,1,1,0d0,0d0
			else
				write(10,"(i0,',',i0,',',i0,',',e30.15,',',e30.15)") j,1,0,0d0,0d0
			end if
		end do
	end do
	NumOfElementsInOneLaminaX = (NumOfNodesInOnePeriodX-1)*NumOfPeriodsX+1
	do i=1,NumOfPeriodsY
		do j=bottomNodesOfLamina(NumOfElementsInOneLaminaX,i),topNodesOfLamina(NumOfElementsInOneLaminaX,i)
			if (j==bottomNodesOfLamina(NumOfElementsInOneLaminaX,1)) then
				write(10,"(i0,',',i0,',',i0,',',e30.15,',',e30.15)") j,1,1,ux,0d0
			else
				write(10,"(i0,',',i0,',',i0,',',e30.15,',',e30.15)") j,1,0,ux,0d0
			end if
		end do
	end do


	close(10)




	ux = 1d-5
	open(10,file='mpc.csv')
	numOfConstraint=0
	do i=1,NumOfPeriodsY
		numOfConstraint = numOfConstraint + topNodesOfLamina(1,i)-bottomNodesOfLamina(1,i)+1
	end do
	numOfConstraint = numOfConstraint + (NumOfNodesInOnePeriodX-1)*NumOfPeriodsX-1

	write(10,"(i0)") numOfConstraint
	do i=1,NumOfPeriodsY
		do j=0,topNodesOfLamina(1,i)-bottomNodesOfLamina(1,i)
			write(10,"(i0,3(',',i0))") bottomNodesOfLamina(1,i)+j,bottomNodesOfLamina((NumOfNodesInOnePeriodX-1)*NumOfPeriodsX+1,i)+j,1,1
		end do
	end do
	do i=2,(NumOfNodesInOnePeriodX-1)*NumOfPeriodsX
		write(10,"(i0,3(',',i0))") bottomNodesOfLamina(i,1),topNodesOfLamina(i,NumOfPeriodsY),1,1
	end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	close(10)



	print *,"Press Enter key to finish."
	print *,54*6+(52-1)
	read *
	stop


end program
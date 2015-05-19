program PeriodicMeshGenerator2D
	use gfrpMeshGenerator2D
	double precision,allocatable,dimension(:,:) :: nodes,newNodes,nodes_,nodesOfLamina
	integer,allocatable,dimension(:,:) :: elements,newElements

!	Parameters of one period
	integer,parameter :: NumOfNodesInOneElement = 4
	integer, parameter :: dim = 2

	integer,parameter :: NumOfNodesInOnePeriod = 787

	integer,parameter :: NumOfElementsInOnePeriod = 720

	double precision,parameter :: ThicknessOfLamina = 0.6d-3
	double precision :: thL,tlL,wL,thC,tlC,wC,dxC,thR,tlR,wR

	integer lastNode,lastNode_,lastElement,NumOfPeriodsX,NumOfPeriodsY,NumOfNodes,NumOfElements
	integer NumOfNodesInOnePeriodX,NumOfNodesInOnePeriodY




	NumOfPeriodsX = 2
	NumOfPeriodsY = 1

	thL = 0.11d-3
	thC = 0.11d-3
	thR = 0.11d-3
	wL = 4.04d-3
	wC = 4.04d-3
	wR = 4.04d-3
	dxC = 0.00d-3

	allocate(nodes(dim,NumOfNodesInOnePeriod))
	allocate(elements(NumOfNodesInOneElement+1,NumOfElementsInOnePeriod))
	allocate(newNodes(dim,NumOfNodesInOnePeriod))
	allocate(nodes_(dim,lastNode))
	allocate(newElements(NumOfNodesInOneElement+1,NumOfElementsInOnePeriod))

	nodes = 0
	nodes_ = 0



	call generateGFRPMesh(nodes,elements,thL,wL,thC,wC,dxC,thR,wR)
	lastNode = NumOfNodesInOnePeriod
	lastElement = NumOfElementsInOnePeriod

!一周期の縦の要素数をカウント
	do i=1,NumOfNodesInOnePeriod
		if (dabs(nodes(1,i))>1d-9) exit
		NumOfNodesInOnePeriodY = NumOfNodesInOnePeriodY + 1
	end do
print *,NumOfNodesInOnePeriodY


	do i=2,NumOfPeriodsX
		deallocate(nodes_)
		allocate(nodes_(dim,lastNode))

		nodes_(:,1:lastNode) = nodes(:,1:lastNode)

		call generateGFRPMesh(newNodes,newElements,thL,wL,thC,wC,dxC,thR,wR)

		lastNode_ = lastNode
		lastNode = lastNode + NumOfNodesInOnePeriod-13
		deallocate(nodes)
		allocate(nodes(dim,lastNode))

		nodes(:,1:lastNode_) = nodes_(:,:)

		nodes(1,lastNode_+1:lastNode_-13+NumOfNodesInOnePeriod) = newNodes(1,14:NumOfNodesInOnePeriod)+nodes(1,lastNode_)
		nodes(2,lastNode_+1:lastNode_-13+NumOfNodesInOnePeriod) = newNodes(2,14:NumOfNodesInOnePeriod)

	end do

	open(10,file='model.csv')

	write (10,*) lastNode, ',',lastElement
	do i=1,lastNode
		write(10,"(i0, 3(',', e30.15))") i, nodes(:,i)
	end do

	do i=1,lastElement
		write(10,"(i0, 100(',', i0))") i, elements(5,i), elements(1:4,i)
	end do
	close(10);

	print *,"Press Enter key to finish."
	read *

contains
! subroutine getNumOfNodesX

end program
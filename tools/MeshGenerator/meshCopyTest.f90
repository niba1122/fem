program test
	use gfrpMeshGenerator2D
	double precision,allocatable,dimension(:,:) :: nodes,newNodes
	integer,allocatable,dimension(:,:) :: elements,newElements

!	Parameters of one period
	integer,parameter :: NumOfNodesInOneElement = 4
	integer, parameter :: dim = 2

	integer,parameter :: NumOfNodesInOnePeriod = 787
	integer,parameter :: NumOfNodesInOnePeriodX = 55
	integer,parameter,dimension(55) :: NumOfNodesInOnePeriodY = &
		&(/13,13,13,13,13,13,13,13,13,13,13,13,13,9,&
			&9,13,13,13,13,13,13,13,13,13,13,13,13,13,&
			&13,13,13,13,13,13,13,13,13,13,13,13,9,&
			&9,13,13,13,13,13,13,13,13,13,13,13,13,13/)
	integer lastNode,lastElement,NumOfPeriodsX,NumOfPeriodsY,NumOfNodes
	integer,parameter :: NumOfElementsInOnePeriod = 720

!	parameters of morphology of GFRP
	double precision :: thL,tlL,wL,thC,tlC,wC,dxC,thR,tlR,wR

!	shift of lamina
	integer,allocatable,dimension(:) :: shiftOfLamina


!--------------------------
	NumOfPeriodsX = 2
	NumOfPeriodsY = 1
allocate(shiftOfLamina(NumOfPeriodsY))
	shiftOfLamina = (/0/)

!--------------------------
	NumOfElements = NumOfElementsInOnePeriod*NumOfPeriodsX*NumOfPeriodsY
!	NumOfNodes = ( (NumOfNodesInOnePeriodX-1)*NumOfPeriodsX + (NumOfNodesInOnePeriodY(1)-1)*NumOfPeriodsY + 1)&
! 		& + (NumOfNodesInOnePeriod-(NumOfNodesInOnePeriodX + NumOfNodesInOnePeriodY(1) - 1))&
! 					&*NumOfPeriodsX*NumOfPeriodsY

	NumOfNodes = NumOfNodesInOnePeriod*NumOfPeriodsX*NumOfPeriodsY&
		&-NumOfPeriodsX*NumOfNodesInOnePeriodX*(NumOfPeriodsY-1)&
		&-NumOfPeriodsY*NumOfNodesInOnePeriodY(1)*(NumOfPeriodsX-1)&
		&+(NumOfPeriodsX-1)*(NumOfPeriodsY-1)

print *,NumOfNodes,NumOfElements
print *,787*2-13
read *


	thL = 0.11d-3
	thC = 0.11d-3
	thR = 0.11d-3
	wL = 4.04d-3
	wC = 4.04d-3
	wR = 4.04d-3
	dxC = 0.00d-3

	allocate(nodes(dim,NumOfNodes*2-13))
	allocate(elements(NumOfNodesInOneElement+1,NumOfElements*2))
	allocate(newNodes(dim,NumOfNodes))
	allocate(newElements(NumOfNodesInOneElement+1,NumOfElements))


! 1周期目
	call generateGFRPMesh(newNodes,newElements,thL,wL,thC,wC,dxC,thR,wR)

	do i = 1,NumOfNodes
		nodes(:,i) = newNodes(:,i)
	end do
	do i = 1,NumOfElements
		elements(:,i) = newElements(:,i)
	end do
! 2周期目
	call generateGFRPMesh(newNodes,newElements,thL,wL,thC,wC,dxC,thR,wR)

	do i = 14,NumOfNodes
		nodes(1,NumOfNodes+i-13) = newNodes(1,i) + nodes(1,NumOfNodes)
		nodes(2,NumOfNodes+i-13) = newNodes(2,i)
	end do
	do i = 1,NumOfElements
		elements(1:4,NumOfElements+i) = newElements(1:4,i) + NumOfNodes-13
		elements(5,NumOfElements+i) = newElements(5,i)
	end do




	open(10,file='model.csv')

	write (10,*) NumOfNodes*2-13, ',',NumOfElements*2
	do i=1,NumOfNodes*2-13
		write(10,"(i0, 3(',', e30.15))") i, nodes(:,i)
	end do

	do i=1,NumOfElements*2
		write(10,"(i0, 100(',', i0))") i, elements(5,i), elements(1:4,i)
	end do
	close(10);

	print *,"Press Enter key to finish."
	read *

end program test
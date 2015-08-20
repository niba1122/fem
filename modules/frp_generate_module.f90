module GFRPMeshGenerator2D
contains
subroutine generateGFRPMesh(nodes,elements,NumOfNodes,thL,wL,thC,wC,dxC,thR,wR,shift)
	implicit none
	integer i,j,k,p,q
	integer lastNode
	integer lastElement
	double precision,allocatable,dimension(:,:),intent(inout) :: nodes
	integer,dimension(:,:),intent(inout) :: elements

	integer,parameter :: NumOfNodesInOneElement = 4
	integer, parameter :: dim = 2;

	integer,intent(inout) :: NumOfNodes
	integer,parameter :: NumOfElements = 704

	double precision,parameter :: widthOfModel = 8.894963d-3
	double precision,parameter :: thicknessOfInterface = 0.02d-3
	double precision,parameter :: minThicknessOfMatrix = 0.02d-3
	double precision,parameter :: ThicknessOfWarp = 0.22d-3
	double precision,parameter :: distalThicknessOfWeft = 0.04d-3
	double precision,parameter :: thicknessOfWeft = 0.22d-3
	double precision,parameter :: widthOfBlock = widthOfModel/52


	integer,dimension(53) :: bottomNodes,topNodes
	double precision,allocatable,dimension(:,:) :: shiftedNodes
  integer,allocatable,dimension(:,:) :: shiftedElements
	integer,intent(in) :: shift

	double precision :: thL,tlL,wL,thC,tlC,wC,dxC,thR,tlR,wR

	NumOfNodes = 769
	deallocate(nodes)
	allocate(nodes(dim,NumOfNodes))
	allocate(shiftedNodes(dim,NumOfNodes))
  allocate(shiftedElements(5,NumOfElements))

	tlL = thicknessOfWeft - thL
	tlC = thicknessOfWeft - thC
	tlR = thicknessOfWeft - thR

	nodes = 0d0
	elements = 0


!-------------------------------------------------------------------------------------------------------
!	Generating model
!-------------------------------------------------------------------------------------------------------

	nodes(2,2) = minThicknessOfMatrix
	nodes(2,3) = nodes(2,2) + thicknessOfInterface/2
	nodes(2,4) = nodes(2,3) + thicknessOfInterface/2
	nodes(2,5) = nodes(2,4) + tlL
	nodes(2,6) = nodes(2,5) + thL
	nodes(2,7) = nodes(2,6) + thicknessOfInterface/2
	nodes(2,8) = nodes(2,7) + thicknessOfInterface/2
	nodes(2,9) = nodes(2,8) + ThicknessOfWarp/2
	nodes(2,10) = nodes(2,9) + ThicknessOfWarp/2
	nodes(2,11) = nodes(2,10) + thicknessOfInterface/2
	nodes(2,12) = nodes(2,11) + thicknessOfInterface/2
	nodes(2,13) = nodes(2,12) + minThicknessOfMatrix

  	lastNode = 13;
 	lastElement = 0;


 	print *,"-- GFRP Mesh Generator ---------------------------"

 	print *,"paramenters:"
 	print *,"thL = ", thL
 	print *,"wL = ", wL
 	print *,"thC = ", thC
 	print *,"wC = ", wC
 	print *,"dxC= ", dxC
 	print *,"thR = ", thR
 	print *,"wR = ", wR
 	print *,""

 	do i=1,12
		call addBlock1L(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
			&minThicknessOfMatrix,thicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
 	end do
 	print *,"Block1L was generated"

	call addBlock2BR(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
		&minThicknessOfMatrix,ThicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
 	print *,"Block2BR was generated"

!  call addBlock3LC(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
!	  &minThicknessOfMatrix,thicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
!  call addBlock3LC(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
!	  &minThicknessOfMatrix,thicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
! 	print *,"Block3LC was generated"

	call addBlock2TL(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
		&minThicknessOfMatrix,ThicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
 	print *,"Block2TL was generated"

 	do i=1,24
		call addBlock1C(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
			&minThicknessOfMatrix,thicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
 	end do
 	print *,"Block1C was generated"

 	call addBlock2TR(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
		&minThicknessOfMatrix,ThicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
 	print *,"Block2TR was generated"

!  call addBlock3CR(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
!	  &minThicknessOfMatrix,thicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
!  call addBlock3CR(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
!	  &minThicknessOfMatrix,thicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
! 	print *,"Block3CR was generated"

 	call addBlock2BL(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
		&minThicknessOfMatrix,ThicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
 	print *,"Block2BL was generated"

 	do i=1,12
		call addBlock1R(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
			&minThicknessOfMatrix,thicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
 	end do
 	print *,"Block1R was generated"

 	print *,"Model was successfully generated"

 	print *,"-- End GFRP Mesh Generator ------------------------"


!-------------------------------------------------------------------------------------------------------
!	Shifting the model
!-------------------------------------------------------------------------------------------------------

if (shift > 0) then

	! 層の上端・下端の節点番号を取得・保存
	p = 1
	q = 1
	do i=1,NumOfNodes
		!　下端
		if (dabs(nodes(2,i)-0d0) < 1d-9) then
			bottomNodes(p) = i
			p = p + 1
		end if
		! 上端
		if (dabs(nodes(2,i)-nodes(2,13)) < 1d-9) then
			topNodes(q) = i
			q = q + 1
		end if
	end do

! 	shiftedNodes(1,1:NumOfNodes-topNodes(shift)) = &
! 		nodes(1,topNodes(shift)+1:NumOfNodes) - nodes(1,topNodes(shift)+1)
! 	shiftedNodes(2,1:NumOfNodes-topNodes(shift)) = nodes(2,topNodes(shift)+1:NumOfNodes)

! 	shiftedNodes(1,NumOfNodes-topNodes(shift)+1:NumOfNodes) = &
! 		nodes(1,topNodes(1)+1:topNodes(shift+1)) - nodes(1,topNodes(shift)+1) + widthOfModel
! 	shiftedNodes(2,NumOfNodes-topNodes(shift)+1:NumOfNodes) = nodes(2,topNodes(1)+1:topNodes(shift+1))

print *,"1:",NumOfNodes-topNodes(shift),",",topNodes(shift)+1,":",NumOfNodes

print *,"1:",NumOfNodes-bottomNodes(shift+1)+1,",",bottomNodes(shift+1),":",NumOfNodes
	shiftedNodes(1,1:NumOfNodes-bottomNodes(shift+1)+1) = &
		nodes(1,bottomNodes(shift+1):NumOfNodes) - nodes(1,bottomNodes(shift+1))
	shiftedNodes(2,1:NumOfNodes-bottomNodes(shift+1)+1) = nodes(2,bottomNodes(shift+1):NumOfNodes)

print *,nodes(1,topNodes(shift)+1),nodes(1,bottomNodes(shift+1))

print *,NumOfNodes-topNodes(shift)+1,":",NumOfNodes,",",topNodes(1)+1,":",topNodes(shift+1)
print *,NumOfNodes-bottomNodes(shift+1)+2,":",NumOfNodes,",",topNodes(1)+1,":",topNodes(shift+1)

	shiftedNodes(1,NumOfNodes-bottomNodes(shift+1)+2:NumOfNodes) = &
		nodes(1,topNodes(1)+1:topNodes(shift+1)) - nodes(1,bottomNodes(shift+1)) + widthOfModel
	shiftedNodes(2,NumOfNodes-bottomNodes(shift+1)+2:NumOfNodes) = nodes(2,topNodes(1)+1:topNodes(shift+1))

	NumOfNodes = 757+topNodes(shift+1)-bottomNodes(shift+1)
	deallocate(nodes)
	allocate(nodes(dim,NumOfNodes))
	print *,"numofnodes: ", NumOfNodes

	nodes(:,:) = shiftedNodes(:,1:NumOfNodes)


  j = 0
	do i=1,NumOfElements
		if (elements(1,i) < bottomNodes(shift+1)) then
			elements(1:4,i) = elements(1:4,i)+NumOfNodes-topNodes(shift+1)
      j = j+1
		else
			elements(1:4,i) = elements(1:4,i)-bottomNodes(shift+1)+1
		end if
	end do
  
  
  shiftedElements(:,(NumOfElements-j+1):NumOfElements) = elements(:,1:j)
  shiftedElements(:,1:(NumOfElements-j)) = elements(:,(j+1):NumOfElements)

  elements(:,:) = shiftedElements(:,:)

! 	do i=1,NumOfElements
! 		if (elements(1,i) < bottomNodes(shift+1)) then 
! 			elements(1:4,i) = elements(1:4,i) + topNodes(shift+1)
! 		end if
! 	end do


end if
print *,NumOfNodes, NumOfElements

!-------------------------------------------------------------------------------------------------------
!	Output of model
!-------------------------------------------------------------------------------------------------------


! 	do i=1,(12*12+1)
! 		print *, elements(:,i)
! 	end do

	! デバッグ用出力
! 	do i=1,lastNode
! 		print *,nodes(:,i)
! 	end do


! 	open(10,file='model.csv')

! 	write (10,*) lastNode, ',',lastElement
! 	do i=1,lastNode
! 		write(10,"(i0, 3(',', e30.15))") i, nodes(:,i)
! 	end do

! 	do i=1,lastElement
! 		write(10,"(i0, 100(',', i0))") i, elements(5,i), elements(1:4,i)
! 	end do
! 	close(10);

!	print *, "Press enter to finish."
!	read *

end subroutine


!-------------------------------------------------------------------------------------------------------
!	Adding Block1L(Weft is bottom of a laminar)
!-------------------------------------------------------------------------------------------------------

subroutine addBlock1L(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
	&minThicknessOfMatrix,thicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
	implicit none
	integer i,j,k
	integer,intent(inout) :: lastNode,lastElement
	integer,dimension(:,:),intent(inout) :: elements
	double precision,dimension(:,:),intent(inout) :: nodes
	double precision widthOfModel,widthOfBlock,thicknessOfInterface,minThicknessOfMatrix,ThicknessOfWarp
	double precision thicknessOfWeft,thL,tlL,wL,thC,tlC,wC,dxC,thR,tlR,wR,distalThicknessOfWeft
	double precision,parameter :: pi = dacos(-1d0)

	double precision x0,y0

	tlL = thicknessOfWeft - thL
	tlC = thicknessOfWeft - thC
	tlR = thicknessOfWeft - thR


! Generating elements
	do i=1,12
		elements(1:4,lastElement+i) = (/lastNode-13+i, lastNode+i, lastNode+1+i, lastNode-12+i/)
!		elements(:,lastElement+1) = (/lastNode-13+i, lastNode+i, 1,1/)
	end do

	! material (1:matrix,2:Warp,3:Weft,4:interface)
	elements(5,lastElement+1:lastElement+12) = (/1,4,4,3,3,4,4,2,2,4,4,1/);

	lastElement = lastElement + 12


! Generating nodes

	x0 = nodes(1,lastNode)

	nodes(1,lastNode+1) = x0+widthOfBlock
	do i=2,12
		nodes(1,lastNode+i) = nodes(1,lastNode-1)+wL/24
	end do
	nodes(1,lastNode+13) = x0+widthOfBlock

	nodes(2,lastNode+1) = 0d0;
	nodes(2,lastNode+13) = minThicknessOfMatrix*2+thicknessOfInterface*3&
							&+tlL+thL+ThicknessOfWarp

	! lower part of Weft
	nodes(2,lastNode+2) = -(tlL-distalThicknessOfWeft/2)*dcos(nodes(1,lastNode+2)*pi/wL)&
												&+tlL+minThicknessOfMatrix-distalThicknessOfWeft/2
	nodes(2,lastNode+3) = nodes(2,lastNode+2)+thicknessOfInterface/2
	nodes(2,lastNode+4) = nodes(2,lastNode+3)+thicknessOfInterface/2
	nodes(2,lastNode+5) = tlL+minThicknessOfMatrix+thicknessOfInterface

	! upper part of Weft
	nodes(2,lastNode+6) = (thL-distalThicknessOfWeft/2)*dcos(nodes(1,lastNode+6)*pi/wL)&
							&+tlL+minThicknessOfMatrix+thicknessOfInterface+distalThicknessOfWeft/2
	nodes(2,lastNode+7) = nodes(2,lastNode+6)+thicknessOfInterface/2
	nodes(2,lastNode+8) = nodes(2,lastNode+7)+thicknessOfInterface/2

	! Warp
	nodes(2,lastNode+9) = nodes(2,lastNode+8)+ThicknessOfWarp/2
	nodes(2,lastNode+10) = nodes(2,lastNode+9)+ThicknessOfWarp/2
	nodes(2,lastNode+11) = nodes(2,lastNode+10)+thicknessOfInterface/2
	nodes(2,lastNode+12) = nodes(2,lastNode+11)+thicknessOfInterface/2

	lastNode = lastNode + 13



end subroutine

!-------------------------------------------------------------------------------------------------------
!	Adding Block2BR(Weft is bottom of a laminar)
!-------------------------------------------------------------------------------------------------------

subroutine addBlock2BR(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
	&minThicknessOfMatrix,thicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
	implicit none
	integer i,j,k
	integer,intent(inout) :: lastNode,lastElement
	integer,dimension(:,:),intent(inout) :: elements
	double precision,dimension(:,:),intent(inout) :: nodes
	double precision widthOfModel,widthOfBlock,thicknessOfInterface,minThicknessOfMatrix,ThicknessOfWarp
	double precision thicknessOfWeft,thL,tlL,wL,thC,tlC,wC,dxC,thR,tlR,wR,distalThicknessOfWeft
	double precision,parameter :: pi = dacos(-1d0)

	double precision gradient
	double precision,allocatable,dimension(:,:) :: defaultNodes
	integer,allocatable,dimension(:,:) :: defaultElements


	allocate(defaultNodes(3,43))
	allocate(defaultElements(5,32))

	tlL = thicknessOfWeft - thL
	tlC = thicknessOfWeft - thC
	tlR = thicknessOfWeft - thR

! generating defaultElements
	
	defaultElements(:,1) = (/1,36,17,2,1/)
	defaultElements(:,2) = (/2,17,14,3,4/)
	defaultElements(:,3) = (/3,14,15,4,4/)
	defaultElements(:,4) = (/4,15,16,5,4/)
	defaultElements(:,5) = (/5,16,19,6,4/)
	do i=6,12
		defaultElements(1:4,i) = (/i,i+13,i+14,i+1/)
	end do
	defaultElements(5,6:12) = (/4,4,2,2,4,4,1/)

	defaultElements(:,13) = (/14,17,18,15,4/)
	defaultElements(:,14) = (/15,18,19,16,4/)

	defaultElements(:,15) = (/17,36,27,18,1/)
	do i=16,23
		defaultElements(1:4,i) = (/i+2,i+11,i+12,i+3/)
	end do
	defaultElements(5,16:23) = (/1,4,4,2,2,4,4,1/)
	defaultElements(:,24) = (/26,35,44,13,1/)

	do i=25,32
		defaultElements(1:4,i) = (/i+2,i+11,i+12,i+3/)
	end do
	defaultElements(5,25:32) = (/1,4,4,2,2,4,4,1/)

! convert to Block2BR
	do i=1,32
		Elements(1:4,lastElement+i) = defaultElements(1:4,i)+LastNode-13
		Elements(5,lastElement+i) = defaultElements(5,i)
	end do
	lastElement = lastElement+32


! x of nodes
	nodes(1,lastNode+1:lastNode+3) = nodes(1,lastNode-11)+thicknessOfInterface/2
	nodes(1,lastNode+4:lastNode+12) = nodes(1,lastNode-11)+thicknessOfInterface
nodes(1,lastNode+13) = nodes(1,lastNode-11)+(nodes(1,lastNode-12)-nodes(1,lastNode-11))/2&
&+(thicknessOfInterface*2+(widthOfBlock-thicknessOfInterface*2)/2)/3
	nodes(1,lastNode+14:lastNode+21) = nodes(1,lastNode-11)+thicknessOfInterface*2
nodes(1,lastNode+22) = nodes(1,lastNode-11)+(nodes(1,lastNode-12)-nodes(1,lastNode-11))/2&
&+(thicknessOfInterface*2+(widthOfBlock-thicknessOfInterface*2)/2)*2/3
	nodes(1,lastNode+23:lastNode+31) = nodes(1,lastNode-12) + widthOfBlock


! y of nodes
	gradient = (-thicknessOfInterface + thL - thC - distalThicknessOfWeft) / (widthOfModel/2 + dxC - wL/2 - wC/2)

	nodes(2,lastNode+1:lastNode+3) = nodes(2,lastNode-10:lastNode-8)
	nodes(2,lastNode+4) = nodes(2,lastNode-11)
	nodes(2,lastNode+5) = nodes(2,lastNode-9)
 	nodes(2,lastNode+6:lastNode+12) = nodes(2,lastNode-7:lastNode-1)+thicknessOfInterface*gradient
	nodes(2,lastNode+13) = (nodes(2,lastNode+12)+nodes(2,lastNode))/2d0
	nodes(2,lastNode+14) = nodes(2,lastNode-9)
	nodes(2,lastNode+15:lastNode+21) = nodes(2,lastNode-7:lastNode-1)+thicknessOfInterface*2d0*gradient
	nodes(2,lastNode+22) = (nodes(2,lastNode+21)+nodes(2,lastNode))/2d0
 	nodes(2,lastNode+23) = nodes(2,lastNode-12)
	nodes(2,lastNode+24:lastNode+30) = nodes(2,lastNode-7:lastNode-1)+(nodes(1,lastNode)+widthOfBlock-nodes(1,lastNode-11))*gradient
	nodes(2,lastNode+31) = nodes(2,lastNode)

	lastNode = lastNode+31

	

end subroutine

!-------------------------------------------------------------------------------------------------------
!	Adding Block3LC(Weft is bottom of a laminar)
!-------------------------------------------------------------------------------------------------------

subroutine addBlock3LC(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
	&minThicknessOfMatrix,thicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
	implicit none
	integer i,j,k
	integer,intent(inout) :: lastNode,lastElement
	integer,dimension(:,:),intent(inout) :: elements
	double precision,dimension(:,:),intent(inout) :: nodes
	double precision widthOfModel,widthOfBlock,thicknessOfInterface,minThicknessOfMatrix,ThicknessOfWarp
	double precision thicknessOfWeft,thL,tlL,wL,thC,tlC,wC,dxC,thR,tlR,wR,distalThicknessOfWeft
	double precision,parameter :: pi = dacos(-1d0)

	double precision x0,y0,gradient

	tlL = thicknessOfWeft - thL
	tlC = thicknessOfWeft - thC
	tlR = thicknessOfWeft - thR

! Generating elements
	do i=1,8
		elements(1:4,lastElement+i) = (/lastNode-9+i, lastNode+i, lastNode+1+i, lastNode-8+i/)
!		elements(:,lastElement+1) = (/lastNode-13+i, lastNode+i, 1,1/)
	end do

	! material (1:matrix,2:Warp,3:Weft,4:interface)
	elements(5,lastElement+1:lastElement+8) = (/1,4,4,2,2,4,4,1/);

	lastElement = lastElement + 8


! Generating nodes

	x0 = nodes(1,lastNode)

	nodes(1,lastNode+1:lastNode+9) = x0+widthOfBlock

	nodes(2,lastNode+1) = 0d0;
	nodes(2,lastNode+9) = minThicknessOfMatrix*2+thicknessOfInterface*3&
							&+tlL+thL+ThicknessOfWarp

	gradient = (-thicknessOfInterface + thL - thC - distalThicknessOfWeft) / (widthOfModel/2 + dxC - wL/2 - wC/2)

	nodes(2,lastNode+2:lastNode+8) =  nodes(2,lastNode-7:lastNode-1) + gradient * widthOfBlock

	lastNode = lastNode + 9

end subroutine

!-------------------------------------------------------------------------------------------------------
!	Adding Block2TL(Weft is top of a laminar)
!-------------------------------------------------------------------------------------------------------

subroutine addBlock2TL(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
	&minThicknessOfMatrix,thicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
	implicit none
	integer i,j,k
	integer,intent(inout) :: lastNode,lastElement
	integer,dimension(:,:),intent(inout) :: elements
	double precision,dimension(:,:),intent(inout) :: nodes
	double precision widthOfModel,widthOfBlock,thicknessOfInterface,minThicknessOfMatrix,ThicknessOfWarp
	double precision thicknessOfWeft,thL,tlL,wL,thC,tlC,wC,dxC,thR,tlR,wR,distalThicknessOfWeft
	double precision,parameter :: pi = dacos(-1d0)

	double precision gradient,thicknessOfModel
	double precision,allocatable,dimension(:,:) :: defaultNodes
	integer,allocatable,dimension(:,:) :: defaultElements


	allocate(defaultNodes(3,43))
	allocate(defaultElements(5,32))

	tlL = thicknessOfWeft - thL
	tlC = thicknessOfWeft - thC
	tlR = thicknessOfWeft - thR

	thicknessOfModel = 2*minThicknessOfMatrix + 3*thicknessOfInterface + thicknessOfWeft + thicknessOfWarp

! generating defaultElements
	
	defaultElements(:,1) = (/1,36,17,2,1/)
	defaultElements(:,2) = (/2,17,14,3,4/)
	defaultElements(:,3) = (/3,14,15,4,4/)
	defaultElements(:,4) = (/4,15,16,5,4/)
	defaultElements(:,5) = (/5,16,19,6,4/)
	do i=6,12
		defaultElements(1:4,i) = (/i,i+13,i+14,i+1/)
	end do
	defaultElements(5,6:12) = (/4,4,2,2,4,4,1/)

	defaultElements(:,13) = (/14,17,18,15,4/)
	defaultElements(:,14) = (/15,18,19,16,4/)

	defaultElements(:,15) = (/17,36,27,18,1/)
	do i=16,23
		defaultElements(1:4,i) = (/i+2,i+11,i+12,i+3/)
	end do
	defaultElements(5,16:23) = (/1,4,4,2,2,4,4,1/)
	defaultElements(:,24) = (/26,35,44,13,1/)

	do i=25,32
		defaultElements(1:4,i) = (/i+2,i+11,i+12,i+3/)
	end do
	defaultElements(5,25:32) = (/1,4,4,2,2,4,4,1/)

! convert to Block2TL
	do i=1,32
! 		Elements(1:4,lastElement+i) = -defaultElements(1:4,33-i)+45+LastNode-9
 		Elements(1,lastElement+i) = -defaultElements(3,33-i)+45+LastNode-9
 		Elements(2,lastElement+i) = -defaultElements(4,33-i)+45+LastNode-9
 		Elements(3,lastElement+i) = -defaultElements(1,33-i)+45+LastNode-9
 		Elements(4,lastElement+i) = -defaultElements(2,33-i)+45+LastNode-9
		Elements(5,lastElement+i) = defaultElements(5,33-i)
	end do
	lastElement = lastElement+32


! x of nodes
	nodes(1,lastNode+2:lastNode+9) = (widthOfModel/2 + dxC - wC/2) - thicknessOfInterface*2
nodes(1,lastNode+1) = nodes(1,lastNode)+(nodes(1,lastNode+2) - nodes(1,lastNode))/2+&
&(thicknessOfInterface*2+(widthOfBlock-thicknessOfInterface*2)/2)/3
	nodes(1,lastNode+11:lastnode+19) = (widthOfModel/2 + dxC - wC/2) - thicknessOfInterface
nodes(1,lastNode+10) = nodes(1,lastNode)+(nodes(1,lastNode+2) - nodes(1,lastNode))/2+&
&(thicknessOfInterface*2+(widthOfBlock-thicknessOfInterface*2)/2)*2/3
	nodes(1,lastNode+20:lastnode+22) = (widthOfModel/2 + dxC - wC/2) - thicknessOfInterface/2
 	nodes(1,lastNode+23) = nodes(1,lastNode) + widthOfBlock
	nodes(1,lastNode+24:lastnode+34) = widthOfModel/2 + dxC - wC/2
 	nodes(1,lastNode+35) = nodes(1,lastNode) + widthOfBlock

! y of nodes
	gradient = (-thicknessOfInterface + thL - thC - distalThicknessOfWeft) / (widthOfModel/2 + dxC - wL/2 - wC/2)



	nodes(2,lastNode+2:lastNode+8) = nodes(2,lastNode-7:lastNode-1) &
								& + (nodes(1,lastNode+2) - nodes(1,lastNode) )*gradient
	nodes(2,lastNode+1) = nodes(2,lastNode+2)/2d0
 	nodes(2,lastNode+9) = thicknessOfModel &
 				& - minThicknessOfMatrix - thicknessOfInterface - thC + distalThicknessOfWeft/2
 	nodes(2,lastNode+11:lastNode+17) = nodes(2,lastNode-7:lastNode-1) &
								& + ((widthOfModel/2 + dxC - wC/2 - thicknessOfInterface) - nodes(1,lastNode) )*gradient
	nodes(2,lastNode+10) = nodes(2,lastNode+11)/2d0
 	nodes(2,lastNode+18) = nodes(2,lastNode+9)
 	nodes(2,lastNode+19) = nodes(2,lastNode+9) + thicknessOfInterface
 	nodes(2,lastNode+20) = nodes(2,lastNode+9) - distalThicknessOfWeft/2
 	nodes(2,lastNode+21) = nodes(2,lastNode+9)
 	nodes(2,lastNode+22) = nodes(2,lastNode+9) + thicknessOfInterface/2
 	nodes(2,lastNode+23) = 0d0
 	nodes(2,lastNode+24:lastNode+30) = nodes(2,lastNode-7:lastNode-1) &
								& + ((widthOfModel/2 + dxC - wC/2) - nodes(1,lastNode) )*gradient
	nodes(2,lastNode+31) = nodes(2,lastNode+20)
	nodes(2,lastNode+32) = nodes(2,lastNode+21)
	nodes(2,lastNode+33) = nodes(2,lastNode+22)
	nodes(2,lastNode+34) = nodes(2,lastNode+19)
	nodes(2,lastNode+35) = thicknessOfModel

	lastNode = lastNode+35

end subroutine

!-------------------------------------------------------------------------------------------------------
!	Adding Block1C(Weft is bottom of a laminar)
!-------------------------------------------------------------------------------------------------------

subroutine addBlock1C(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
	&minThicknessOfMatrix,thicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
	implicit none
	integer i,j,k
	integer,intent(inout) :: lastNode,lastElement
	integer,dimension(:,:),intent(inout) :: elements
	double precision,dimension(:,:),intent(inout) :: nodes
	double precision widthOfModel,widthOfBlock,thicknessOfInterface,minThicknessOfMatrix,ThicknessOfWarp
	double precision thicknessOfWeft,thL,tlL,wL,thC,tlC,wC,dxC,thR,tlR,wR,distalThicknessOfWeft
	double precision,parameter :: pi = dacos(-1d0)

	double precision x0,y0,thicknessOfModel

	tlL = thicknessOfWeft - thL
	tlC = thicknessOfWeft - thC
	tlR = thicknessOfWeft - thR

	thicknessOfModel = minThicknessOfMatrix*2+thicknessOfInterface*3+thicknessOfWeft+thicknessOfWarp

! Generating elements
	do i=1,12
		elements(1:4,lastElement+i) = (/lastNode-13+i, lastNode+i, lastNode+1+i, lastNode-12+i/)
!		elements(:,lastElement+1) = (/lastNode-13+i, lastNode+i, 1,1/)
	end do

	! material (1:matrix,2:Warp,3:Weft,4:interface)
	elements(5,lastElement+1:lastElement+12) = (/1,4,4,2,2,4,4,3,3,4,4,1/);

	lastElement = lastElement + 12


! Generating nodes

	x0 = nodes(1,lastNode)
	nodes(1,lastNode+1) = x0+widthOfBlock
	do i=2,12
		nodes(1,lastNode+i) = nodes(1,lastNode-1)+wC/24
	end do
	nodes(1,lastNode+13) = x0+widthOfBlock

	nodes(2,lastNode+1) = 0d0;
	nodes(2,lastNode+13) = thicknessOfModel

 	! upper part of Weft
	nodes(2,lastNode+12) = (thC - distalThicknessOfWeft/2) * dcos(pi/wC*(nodes(1,lastNode+12)-(widthOfModel/2+dxC)))&
								& + thicknessOfModel - minThicknessOfMatrix - thC + distalThicknessOfWeft/2
	nodes(2,lastNode+11) = nodes(2,lastNode+12) - thicknessOfInterface/2
	nodes(2,lastNode+10) = nodes(2,lastNode+11) - thicknessOfInterface/2
	nodes(2,lastNode+9) = thicknessOfModel - minThicknessOfMatrix - thicknessOfInterface&
								& - thC
 	! lower part of Weft
	nodes(2,lastNode+8) = - (tlC - distalThicknessOfWeft/2) * dcos(pi/wC*(nodes(1,lastNode+8)-(widthOfModel/2+dxC)))&
					& + thicknessOfModel - minThicknessOfMatrix - thicknessOfInterface - thC - distalThicknessOfWeft/2
	nodes(2,lastNode+7) = nodes(2,lastNode+8) - thicknessOfInterface/2
	nodes(2,lastNode+6) = nodes(2,lastNode+7) - thicknessOfInterface/2

 	! Warp
 	nodes(2,lastNode+5) = nodes(2,lastNode+6) - thicknessOfWarp/2
 	nodes(2,lastNode+4) = nodes(2,lastNode+5) - thicknessOfWarp/2
	nodes(2,lastNode+3) = nodes(2,lastNode+4) - thicknessOfInterface/2
	nodes(2,lastNode+2) = nodes(2,lastNode+3) - thicknessOfInterface/2

	lastNode = lastNode + 13

end subroutine

!-------------------------------------------------------------------------------------------------------
!	Adding Block2TR(Weft is top of a laminar)
!-------------------------------------------------------------------------------------------------------

subroutine addBlock2TR(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
	&minThicknessOfMatrix,thicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
	implicit none
	integer i,j,k
	integer,intent(inout) :: lastNode,lastElement
	integer,dimension(:,:),intent(inout) :: elements
	double precision,dimension(:,:),intent(inout) :: nodes
	double precision widthOfModel,widthOfBlock,thicknessOfInterface,minThicknessOfMatrix,ThicknessOfWarp
	double precision thicknessOfWeft,thL,tlL,wL,thC,tlC,wC,dxC,thR,tlR,wR,distalThicknessOfWeft
	double precision,parameter :: pi = dacos(-1d0)

	double precision gradient,thicknessOfModel
	double precision,allocatable,dimension(:,:) :: defaultNodes
	integer,allocatable,dimension(:,:) :: defaultElements


	allocate(defaultNodes(3,43))
	allocate(defaultElements(5,32))

	tlL = thicknessOfWeft - thL
	tlC = thicknessOfWeft - thC
	tlR = thicknessOfWeft - thR

	thicknessOfModel = 2*minThicknessOfMatrix + 3*thicknessOfInterface + thicknessOfWeft + thicknessOfWarp

! generating defaultElements
	do i=1,7
		defaultElements(1:4,i) = (/i,i+16,i+17,i+1/)
	end do
	defaultElements(5,1:7) = (/1,4,4,2,2,4,4/)
	defaultElements(:,8) = (/8,24,14,9,4/)
	defaultElements(:,9) = (/9,14,15,10,4/)
	defaultElements(:,10) = (/10,15,16,11,4/)
	defaultElements(:,11) = (/11,16,26,12,4/)
	defaultElements(:,12) = (/12,26,44,13,1/)
	defaultElements(:,13) = (/14,24,25,15,4/)
	defaultElements(:,14) = (/15,25,26,16,4/)
	defaultElements(:,15) = (/1,36,27,17,1/)
	do i=16,23
		defaultElements(1:4,i) = (/i+1,i+11,i+12,i+2/)
	end do
	defaultElements(5,16:23) = (/1,4,4,2,2,4,4,1/)
	defaultElements(:,24) = (/25,35,44,26,1/)
	do i=25,32
		defaultElements(1:4,i) = (/i+2,i+11,i+12,i+3/)
	end do
	defaultElements(5,25:32) = (/1,4,4,2,2,4,4,1/)

! convert to Block2TR
	do i=1,32
		Elements(1:4,lastElement+i) = defaultElements(1:4,i)+LastNode-13
		Elements(5,lastElement+i) = defaultElements(5,i)
	end do
	lastElement = lastElement+32


! x of nodes
	nodes(1,lastNode+1) = nodes(1,lastNode-1) + thicknessOfInterface/2
	nodes(1,lastNode+2) = nodes(1,lastNode-1) + thicknessOfInterface/2
	nodes(1,lastNode+3) = nodes(1,lastNode-1) + thicknessOfInterface/2
	do i=5,13
		nodes(1,lastNode+i) = nodes(1,lastNode-1) + thicknessOfInterface
	end do
  nodes(1,lastNode+4) = nodes(1,lastNode-1)+(nodes(1,lastNode)-nodes(1,lastNode-1))/2&
  &+(thicknessOfInterface*2+(widthOfBlock-thicknessOfInterface*2)/2)/3
	do i=15,22
		nodes(1,lastNode+i) = nodes(1,lastNode-1) + thicknessOfInterface*2
	end do
  nodes(1,lastNode+14) = nodes(1,lastNode-1)+(nodes(1,lastNode)-nodes(1,lastNode-1))/2&
  &+(thicknessOfInterface*2+(widthOfBlock-thicknessOfInterface*2)/2)*2/3
	do i=23,31
		nodes(1,lastNode+i) = nodes(1,lastNode) + widthOfBlock
	end do

! y of nodes
	gradient = (thicknessOfInterface - thR + thC + distalThicknessOfWeft) / (widthOfModel/2 - dxC - wR/2 - wC/2)
	nodes(2,lastNode+1:lastNode+3) = nodes(2,lastNode-4:lastNode-2)
	nodes(2,lastNode+5:lastNode+11) = nodes(2,lastNode-11:lastNode-5) + thicknessOfInterface*gradient
	nodes(2,lastNode+4) = nodes(2,lastNode+5)/2
	nodes(2,lastNode+12) = nodes(2,lastNode+2)
	nodes(2,lastNode+13) = nodes(2,lastNode-1)
	nodes(2,lastNode+15:lastNode+21) = nodes(2,lastNode+5:lastNode+11) + thicknessOfInterface*gradient
	nodes(2,lastNode+14) = nodes(2,lastNode+15)/2
	nodes(2,lastNode+22) = nodes(2,lastNode+12)
	nodes(2,lastNode+23) = nodes(2,lastNode-12)
	nodes(2,lastNode+24:lastNode+30) = nodes(2,lastNode+15:lastNode+21)&
				& +(nodes(1,lastNode)-nodes(1,lastNode-1)+widthOfBlock-thicknessOfInterface*2)*gradient
	nodes(2,lastNode+31) = nodes(2,lastNode)

	lastNode = lastNode+31

end subroutine

!-------------------------------------------------------------------------------------------------------
!	Adding Block3CR(Weft is bottom of a laminar)
!-------------------------------------------------------------------------------------------------------

subroutine addBlock3CR(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
	&minThicknessOfMatrix,thicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
	implicit none
	integer i,j,k
	integer,intent(inout) :: lastNode,lastElement
	integer,dimension(:,:),intent(inout) :: elements
	double precision,dimension(:,:),intent(inout) :: nodes
	double precision widthOfModel,widthOfBlock,thicknessOfInterface,minThicknessOfMatrix,ThicknessOfWarp
	double precision thicknessOfWeft,thL,tlL,wL,thC,tlC,wC,dxC,thR,tlR,wR,distalThicknessOfWeft
	double precision,parameter :: pi = dacos(-1d0)

	double precision x0,y0,gradient

	tlL = thicknessOfWeft - thL
	tlC = thicknessOfWeft - thC
	tlR = thicknessOfWeft - thR

! Generating elements
	do i=1,8
		elements(1:4,lastElement+i) = (/lastNode-9+i, lastNode+i, lastNode+1+i, lastNode-8+i/)
!		elements(:,lastElement+1) = (/lastNode-13+i, lastNode+i, 1,1/)
	end do

	! material (1:matrix,2:Warp,3:Weft,4:interface)
	elements(5,lastElement+1:lastElement+8) = (/1,4,4,2,2,4,4,1/);

	lastElement = lastElement + 8


! Generating nodes

	x0 = nodes(1,lastNode)

	nodes(1,lastNode+1:lastNode+9) = x0+widthOfBlock

	nodes(2,lastNode+1) = 0d0;
	nodes(2,lastNode+9) = minThicknessOfMatrix*2+thicknessOfInterface*3&
							&+tlL+thL+ThicknessOfWarp

	gradient = (thicknessOfInterface - thR + thC + distalThicknessOfWeft) / (widthOfModel/2 - dxC - wR/2 - wC/2)

	nodes(2,lastNode+2:lastNode+8) =  nodes(2,lastNode-7:lastNode-1) + gradient * widthOfBlock

	lastNode = lastNode + 9

end subroutine

!-------------------------------------------------------------------------------------------------------
!	Adding Block2BL(Weft is top of a laminar)
!-------------------------------------------------------------------------------------------------------

subroutine addBlock2BL(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
	&minThicknessOfMatrix,thicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
	implicit none
	integer i,j,k
	integer,intent(inout) :: lastNode,lastElement
	integer,dimension(:,:),intent(inout) :: elements
	double precision,dimension(:,:),intent(inout) :: nodes
	double precision widthOfModel,widthOfBlock,thicknessOfInterface,minThicknessOfMatrix,ThicknessOfWarp
	double precision thicknessOfWeft,thL,tlL,wL,thC,tlC,wC,dxC,thR,tlR,wR,distalThicknessOfWeft
	double precision,parameter :: pi = dacos(-1d0)

	double precision gradient,thicknessOfModel
	double precision,allocatable,dimension(:,:) :: defaultNodes
	integer,allocatable,dimension(:,:) :: defaultElements


	allocate(defaultNodes(3,43))
	allocate(defaultElements(5,32))

	tlL = thicknessOfWeft - thL
	tlC = thicknessOfWeft - thC
	tlR = thicknessOfWeft - thR

	thicknessOfModel = 2*minThicknessOfMatrix + 3*thicknessOfInterface + thicknessOfWeft + thicknessOfWarp

! generating defaultElements
	do i=1,7
		defaultElements(1:4,i) = (/i,i+16,i+17,i+1/)
	end do
	defaultElements(5,1:7) = (/1,4,4,2,2,4,4/)
	defaultElements(:,8) = (/8,24,14,9,4/)
	defaultElements(:,9) = (/9,14,15,10,4/)
	defaultElements(:,10) = (/10,15,16,11,4/)
	defaultElements(:,11) = (/11,16,26,12,4/)
	defaultElements(:,12) = (/12,26,44,13,1/)
	defaultElements(:,13) = (/14,24,25,15,4/)
	defaultElements(:,14) = (/15,25,26,16,4/)
	defaultElements(:,15) = (/1,36,27,17,1/)
	do i=16,23
		defaultElements(1:4,i) = (/i+1,i+11,i+12,i+2/)
	end do
	defaultElements(5,16:23) = (/1,4,4,2,2,4,4,1/)
	defaultElements(:,24) = (/25,35,44,26,1/)
	do i=25,32
		defaultElements(1:4,i) = (/i+2,i+11,i+12,i+3/)
	end do
	defaultElements(5,25:32) = (/1,4,4,2,2,4,4,1/)

! convert to Block2BL
	do i=1,32
! 		Elements(1:4,lastElement+i) = -defaultElements(1:4,33-i)+45+LastNode-9
 		Elements(1,lastElement+i) = -defaultElements(3,33-i)+45+LastNode-9
 		Elements(2,lastElement+i) = -defaultElements(4,33-i)+45+LastNode-9
 		Elements(3,lastElement+i) = -defaultElements(1,33-i)+45+LastNode-9
 		Elements(4,lastElement+i) = -defaultElements(2,33-i)+45+LastNode-9
		Elements(5,lastElement+i) = defaultElements(5,33-i)
	end do
	lastElement = lastElement+32

! x of nodes
	nodes(1,lastNode+1:lastNode+8) = widthOfModel - wR/2 - thicknessOfInterface*2
  nodes(1,lastNode+9)=nodes(1,lastNode)+(nodes(1,lastNode+1)-nodes(1,lastNode))/2&
  &+(thicknessOfInterface*2+(widthOfBlock-thicknessOfInterface*2)/2)/3
	nodes(1,lastNode+10:lastNode+19) = widthOfModel - wR/2 - thicknessOfInterface
  nodes(1,lastNode+19)=nodes(1,lastNode)+(nodes(1,lastNode+1)-nodes(1,lastNode))/2&
  &+(thicknessOfInterface*2+(widthOfBlock-thicknessOfInterface*2)/2)*2/3
	nodes(1,lastNode+20:lastNode+22) = widthOfModel - wR/2 - thicknessOfInterface/2
	nodes(1,lastNode+23) = nodes(1,lastNode) + widthOfBlock
	nodes(1,lastNode+24:lastNode+34) = widthOfModel - wR/2
	nodes(1,lastNode+35) = nodes(1,lastNode) + widthOfBlock

! y of nodes
	gradient = (thicknessOfInterface - thR + thC + distalThicknessOfWeft) / (widthOfModel/2 - dxC - wR/2 - wC/2)

 	nodes(2,lastNode+1) = minThicknessOfMatrix + tlR + distalThicknessOfWeft/2
 	nodes(2,lastNode+1) = minThicknessOfMatrix + tlR + thicknessOfInterface-distalThicknessOfWeft/2
	nodes(2,lastNode+2:lastNode+8) = nodes(2,lastNode-7:lastNode-1) &
 				& + ((widthOfModel - wR/2) - nodes(1,lastNode-1) - thicknessOfInterface*2 )*gradient

	nodes(2,lastNode+9) = (nodes(2,lastNode) + nodes(2,lastNode+8))/2
	nodes(2,lastNode+10) = nodes(2,lastNode+1) - thicknessOfInterface
	nodes(2,lastNode+11) = nodes(2,lastNode+1)
	nodes(2,lastNode+12:lastNode+18) = nodes(2,lastNode+2:lastNode+8) &
 								& + thicknessOfInterface*gradient
 	nodes(2,lastNode+19) = (nodes(2,lastNode) + nodes(2,lastNode+18))/2
 	nodes(2,lastNode+20) = nodes(2,lastNode+1) - thicknessOfInterface/2
 	nodes(2,lastNode+21) = nodes(2,lastNode+1)
 	nodes(2,lastNode+22) = nodes(2,lastNode+1) + distalThicknessOfWeft/2
 	nodes(2,lastNode+23) = nodes(2,lastNode-8)
 	nodes(2,lastNode+24) = nodes(2,lastNode+10)
 	nodes(2,lastNode+25) = nodes(2,lastNode+20)
 	nodes(2,lastNode+26) = nodes(2,lastNode+21)
 	nodes(2,lastNode+27) = nodes(2,lastNode+22)
 	nodes(2,lastNode+28:lastNode+34) = nodes(2,lastNode+12:lastNode+18) &
 								& + thicknessOfInterface*gradient
 	nodes(2,lastNode+35) = thicknessOfModel
! 	nodes(2,lastNode+1:lastNode+3) = nodes(2,lastNode-4:lastNode-2)
! 	nodes(2,lastNode+5:lastNode+11) = nodes(2,lastNode-11:lastNode-5) + thicknessOfInterface*gradient
! 	nodes(2,lastNode+4) = nodes(2,lastNode+5)/2
! 	nodes(2,lastNode+12) = nodes(2,lastNode+2)
! 	nodes(2,lastNode+13) = nodes(2,lastNode-1)
! 	nodes(2,lastNode+15:lastNode+22) = nodes(2,lastNode+5:lastNode+12) + thicknessOfInterface*gradient
! 	nodes(2,lastNode+14) = nodes(2,lastNode+15)/2
! 	nodes(2,lastNode+23) = nodes(2,lastNode-12)
! 	nodes(2,lastNode+24:lastNode+30) = nodes(2,lastNode-11:lastNode-5)&
! 										& +(nodes(2,lastNode)-nodes(2,lastNode-1)+widthOfBlock)*gradient
! 	nodes(2,lastNode+31) = nodes(2,lastNode)

	lastNode = lastNode+35

end subroutine

!-------------------------------------------------------------------------------------------------------
!	Adding Block1R(Weft is bottom of a laminar)
!-------------------------------------------------------------------------------------------------------

subroutine addBlock1R(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
	&minThicknessOfMatrix,thicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
	implicit none
	integer i,j,k
	integer,intent(inout) :: lastNode,lastElement
	integer,dimension(:,:),intent(inout) :: elements
	double precision,dimension(:,:),intent(inout) :: nodes
	double precision widthOfModel,widthOfBlock,thicknessOfInterface,minThicknessOfMatrix,ThicknessOfWarp
	double precision thicknessOfWeft,thL,tlL,wL,thC,tlC,wC,dxC,thR,tlR,wR,distalThicknessOfWeft
	double precision,parameter :: pi = dacos(-1d0)

	double precision x0,y0

	tlL = thicknessOfWeft - thL
	tlC = thicknessOfWeft - thC
	tlR = thicknessOfWeft - thR

! Generating elements
	do i=1,12
		elements(1:4,lastElement+i) = (/lastNode-13+i, lastNode+i, lastNode+1+i, lastNode-12+i/)
!		elements(:,lastElement+1) = (/lastNode-13+i, lastNode+i, 1,1/)
	end do

	! material (1:matrix,2:Warp,3:Weft,4:interface)
	elements(5,lastElement+1:lastElement+12) = (/1,4,4,3,3,4,4,2,2,4,4,1/);

	lastElement = lastElement + 12


! Generating nodes

	x0 = nodes(1,lastNode)


	nodes(1,lastNode+1) = x0+widthOfBlock
	do i=2,12
		nodes(1,lastNode+i) = nodes(1,lastNode-1)+wR/24
	end do
	nodes(1,lastNode+13) = x0+widthOfBlock



	nodes(2,lastNode+1) = 0d0;
	nodes(2,lastNode+13) = minThicknessOfMatrix*2+thicknessOfInterface*3&
							&+thicknessOfWeft + ThicknessOfWarp

	! lower part of Weft
	nodes(2,lastNode+2) = -(tlR-distalThicknessOfWeft/2)*dcos((nodes(1,lastNode+2)-widthOfModel)*pi/wR)&
												&+tlR+minThicknessOfMatrix-distalThicknessOfWeft/2
	nodes(2,lastNode+3) = nodes(2,lastNode+2)+thicknessOfInterface/2
	nodes(2,lastNode+4) = nodes(2,lastNode+3)+thicknessOfInterface/2
	nodes(2,lastNode+5) = tlR+minThicknessOfMatrix+thicknessOfInterface

	! upper part of Weft
	nodes(2,lastNode+6) = (thR-distalThicknessOfWeft/2)*dcos((nodes(1,lastNode+6)-widthOfModel)*pi/wR)&
							&+tlR+minThicknessOfMatrix+thicknessOfInterface+distalThicknessOfWeft/2
	nodes(2,lastNode+7) = nodes(2,lastNode+6)+thicknessOfInterface/2
	nodes(2,lastNode+8) = nodes(2,lastNode+7)+thicknessOfInterface/2

	! Warp
	nodes(2,lastNode+9) = nodes(2,lastNode+8)+ThicknessOfWarp/2
	nodes(2,lastNode+10) = nodes(2,lastNode+9)+ThicknessOfWarp/2
	nodes(2,lastNode+11) = nodes(2,lastNode+10)+thicknessOfInterface/2
	nodes(2,lastNode+12) = nodes(2,lastNode+11)+thicknessOfInterface/2

	lastNode = lastNode + 13

end subroutine

end module GFRPMeshGenerator2D


module frp_generate_module

!
! made by PeriodicMeshGenerator1.3 and GFRPMeshGenerator2D1.1
!

contains

subroutine generate_FRP_Model(model,shift)
	use fem_module
	use GFRPMeshGenerator2D
	implicit none
	type(struct_model) :: model
	type(struct_bc) :: bc

	integer i,j,k,p,q
	double precision,allocatable,dimension(:,:) :: nodes,newNodes,nodes_,nodesOfLamina,NodesOfShiftedLamina
	integer,allocatable,dimension(:,:) :: topNodesOfLamina,bottomNodesOfLamina
	integer,allocatable,dimension(:,:) :: elements,newElements,elementsOfLamina

!	Parameters of one period
	integer,parameter :: NumOfNodesInOneElement = 4
	integer, parameter :: dim = 2

!	integer :: NumOfNodesInOnePeriod = 787
 	integer NumOfNodesInOnePeriod
	integer,parameter :: NumOfElementsInOnePeriod = 704
	integer,parameter :: NumOfNodesInOnePeriodX  = 53

	double precision,parameter :: ThicknessOfLamina = 0.54d-3
	double precision,parameter :: thicknessOfInterlaminar = 0.02d-3
	double precision,allocatable,dimension(:,:) :: thL,tlL,wL,thC,tlC,wC,dxC,thR,tlR,wR

	integer lastNode,lastNode_,lastElement,NumOfPeriodsX,NumOfPeriodsY,NumOfNodes,NumOfElements
	integer NumOfNodesInOnePeriodY,NumOfNodesInOneLamina,NumOfElementsInOneLamina

	integer :: shift(:)

	integer numOfConstraint,NumOfElementsInOneLaminaX

	double precision,pointer :: angle(:)


	NumOfPeriodsX = 2
	NumOfPeriodsY = 4

	allocate(thL(NumOfPeriodsX,NumOfPeriodsY))
	allocate(wL(NumOfPeriodsX,NumOfPeriodsY))
	allocate(thC(NumOfPeriodsX,NumOfPeriodsY))
	allocate(wC(NumOfPeriodsX,NumOfPeriodsY))
	allocate(dxC(NumOfPeriodsX,NumOfPeriodsY))
	allocate(thR(NumOfPeriodsX,NumOfPeriodsY))
	allocate(wR(NumOfPeriodsX,NumOfPeriodsY))

	thL = 0.113d-3
	thC = 0.113d-3
	thR = 0.113d-3
	wL = 4.206d-3
!	wC = 4.206d-3
	wC = 4.206d-3
	wR = 4.206d-3
	dxC = 0d-3

! 	thC(1,1) = 0.08d-3
! 	thC(2,1) = 0.16d-3

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
		call generateGFRPMesh(newNodes,newElements,NumOfNodesInOnePeriod,&
								&thL(1,k),wL(1,k),thC(1,k),wC(1,k),dxC(1,k),thR(1,k),wR(1,k),shift(k))

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

			call generateGFRPMesh(newNodes,newElements,NumOfNodesInOnePeriod,&
									&thL(i,k),wL(i,k),thC(i,k),wC(i,k),dxC(i,k),thR(i,k),wR(i,k),shift(k))

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
! shift = shift + 9
! if (shift >= 54) shift = shift-54


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

! 	numOfConstraint=0
! 	do i=1,NumOfPeriodsY
! 		numOfConstraint = numOfConstraint + topNodesOfLamina(1,i)-bottomNodesOfLamina(1,i)+1
! 	end do
! 	do i=1,NumOfPeriodsY
! 		numOfConstraint = numOfConstraint + topNodesOfLamina(1,i)-bottomNodesOfLamina(1,i)+1
! 	end do

! allocate(bc%nodes_spc(dim+1,numOfConstraint))
! allocate(bc%disp_spc(dim,numOfConstraint))

! k=1

! 	do i=1,NumOfPeriodsY
! 		do j=bottomNodesOfLamina(1,i),topNodesOfLamina(1,i)
! 			if (j==bottomNodesOfLamina(1,3)) then
! bc%nodes_spc(:,k) = (/j,1,1/)
! bc%disp_spc(:,k) = (/0d0,0d0/)
! 			else

! bc%nodes_spc(:,k) = (/j,1,0/)
! bc%disp_spc(:,k) = (/0d0,0d0/)
! 			end if
! k=k+1
! 		end do
! 	end do
! 	NumOfElementsInOneLaminaX = (NumOfNodesInOnePeriodX-1)*NumOfPeriodsX+1
! 	do i=1,NumOfPeriodsY
! 		do j=bottomNodesOfLamina(NumOfElementsInOneLaminaX,i),topNodesOfLamina(NumOfElementsInOneLaminaX,i)
! 			if (j==bottomNodesOfLamina(1,3)) then
! bc%nodes_spc(:,k) = (/j,1,1/)
! bc%disp_spc(:,k) = (/1d-3,0d0/)
! 			else
! bc%nodes_spc(:,k) = (/j,1,0/)
! bc%disp_spc(:,k) = (/1d-3,0d0/)
! 			end if
! k=k+1
! 		end do
! 	end do

! 	bc%n_spc = k-1


! 	bc%n_mpc = 0


! 構造体へのモデルデータの保存

	allocate(model%nodes(ubound(nodes,1),ubound(nodes,2)))
	allocate(model%elements(4,lastElement))
	allocate(model%material_nos(lastElement))
	allocate(model%materials(6,6,4))
	allocate(model%data(6))
	allocate(angle(lastElement))

!物性値


model%materials(1,:,1) = (/4038461538.46153, 1730769230.76923, 1730769230.76923, 0.00000, 0.00000, 0.00000/)
model%materials(2,:,1) = (/1730769230.76923, 4038461538.46154, 1730769230.76923, 0.00000, 0.00000, 0.00000/)
model%materials(3,:,1) = (/1730769230.76923, 1730769230.76923, 4038461538.46154, 0.00000, 0.00000, 0.00000/)
model%materials(4,:,1) = (/0.00000, 0.00000, 0.00000, 1153846153.84615, 0.00000, 0.00000/)
model%materials(5,:,1) = (/0.00000, 0.00000, 0.00000, 0.00000, 1153846153.84615, 0.00000/)
model%materials(6,:,1) = (/0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 1153846153.84615/)

model%materials(1,:,2) = (/41240000000.00000, 3843000000.00000, 3844000000.00000, 30390.00000, -28.14000, -48.82000/)
model%materials(2,:,2) = (/3843000000.00000, 10050000000.00000, 2763000000.00000, -2425000.00000, -8.02100, -134.10000/)
model%materials(3,:,2) = (/3844000000.00000, 2763000000.00000, 10050000000.00000, 2527000.00000, -85.76000, -28.65000/)
model%materials(4,:,2) = (/30390.00000, -2425000.00000, 2527000.00000, 2326000000.00000, -0.00011, 24.37000/)
model%materials(5,:,2) = (/-28.14000, -8.02100, -85.76000, -0.00011, 3183000000.00000, 36.71000/)
model%materials(6,:,2) = (/-48.82000, -134.10000, -28.65000, 24.37000, 36.71000, 3182000000.00000/)

model%materials(1,:,3) = (/10050000000.00000, 2763000000.00000, 3843000000.00000, 8.02100, 134.10000, -2425000.00000/)
model%materials(2,:,3) = (/2763000000.00000, 10050000000.00000, 3844000000.00000, 85.76000, 28.65000, 2527000.00000/)
model%materials(3,:,3) = (/3843000000.00000, 3844000000.00000, 41240000000.00000, 28.14000, 48.82000, 30390.00000/) 
model%materials(4,:,3) = (/8.02100, 85.76000, 28.14000, 3183000000.00000, 36.71000, 0.00011/)
model%materials(5,:,3) = (/134.10000, 28.65000, 48.82000, 36.71000, 3182000000.00000, -24.37000/)
model%materials(6,:,3) = (/-2425000.00000, 2527000.00000, 30390.00000, 0.00011, -24.37000, 2326000000.00000/)

model%materials(1,:,4) = (/4038461538.46153, 1730769230.76923, 1730769230.76923, 0.00000, 0.00000, 0.00000/)
model%materials(2,:,4) = (/1730769230.76923, 4038461538.46154, 1730769230.76923, 0.00000, 0.00000, 0.00000/)
model%materials(3,:,4) = (/1730769230.76923, 1730769230.76923, 4038461538.46154, 0.00000, 0.00000, 0.00000/)
model%materials(4,:,4) = (/0.00000, 0.00000, 0.00000, 1153846153.84615, 0.00000, 0.00000/)
model%materials(5,:,4) = (/0.00000, 0.00000, 0.00000, 0.00000, 1153846153.84615, 0.00000/)
model%materials(6,:,4) = (/0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 1153846153.84615/)


	model%nodes(:,:) = nodes(:,:)
	model%elements(1:4,:) = elements(1:4,:)
	model%material_nos(:) = elements(5,:)
	model%n_els = lastElement
	model%n_nds = lastNode
	model%dim = dim

	model%type_el = 'quad4'
	model%n_nds_1el = 4

	! 損傷テンソルの初期化
	allocate(model%data(2)%d(3,model%n_els,1))
	model%data(2)%d = 0d0


	! 異方性の回転
	angle = 0d0
	do i=1,lastElement
		if (elements(5,i) == 2) then
 			angle(i) = atan2(nodes(2,elements(2,i))+nodes(2,elements(3,i))-nodes(2,elements(1,i))-nodes(2,elements(4,i)) , &
 							&nodes(1,elements(2,i))+nodes(1,elements(3,i))-nodes(1,elements(1,i))-nodes(1,elements(4,i)) );
! 			D(:,:,i) = materials(:,:,material_nos(i))

		else 
! 			D(:,:,i) = materials(:,:,material_nos(i))
		end if
	end do


	model%data(1)%d(1:model%n_els,1:1,1:1) => angle(:) ! 異方性の角度をモデルに保存



! 上下左右端のノードを保存


! 左端

k=0
	do i=1,NumOfPeriodsY
		do j=bottomNodesOfLamina(1,i),topNodesOfLamina(1,i)
			k=k+1
		end do
	end do

	allocate(model%data(1)%i(k,1,1))


k=1
	do i=1,NumOfPeriodsY
		do j=bottomNodesOfLamina(1,i),topNodesOfLamina(1,i)
			model%data(1)%i(k,1,1) = j
			k=k+1
		end do
	end do

! 左中央
	allocate(model%data(5)%i(1,1,1))
	model%data(5)%i(1,1,1) = topNodesOfLamina(1,NumOfPeriodsY/2)

! 右端
k=0
	NumOfElementsInOneLaminaX = (NumOfNodesInOnePeriodX-1)*NumOfPeriodsX+1
	do i=1,NumOfPeriodsY
		do j=bottomNodesOfLamina(NumOfElementsInOneLaminaX,i),topNodesOfLamina(NumOfElementsInOneLaminaX,i)
			k=k+1
		end do
	end do

	allocate(model%data(2)%i(k,1,1))
K=1
	do i=1,NumOfPeriodsY
		do j=bottomNodesOfLamina(NumOfElementsInOneLaminaX,i),topNodesOfLamina(NumOfElementsInOneLaminaX,i)
			model%data(2)%i(k,1,1) = j
			k=k+1
		end do
	end do



! 右中央
	allocate(model%data(6)%i(1,1,1))
	model%data(6)%i(1,1,1) = topNodesOfLamina(NumOfElementsInOneLaminaX,NumOfPeriodsY/2)


! 下端
	allocate(model%data(3)%i(NumOfElementsInOneLaminaX,1,1))

	do i=1,NumOfElementsInOneLaminaX
		model%data(3)%i(i,1,1) = bottomNodesOfLamina(i,1)
	end do


! 上端
	allocate(model%data(4)%i(NumOfElementsInOneLaminaX,1,1))

	do i=1,NumOfElementsInOneLaminaX
		model%data(4)%i(i,1,1) = topNodesOfLamina(i,NumOfPeriodsY)
	end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	print *,"Press Enter key to finish."

! 	type struct_model
! 		character type_el*16
! 		character name*16
! 		double precision,pointer,dimension(:,:) :: nodes
! 		integer,pointer,dimension(:,:) :: elements
! 		double precision,pointer,dimension(:,:,:) :: materials
! 		integer,pointer :: material_nos(:)
! 		integer n_nds,n_els,dim,n_nds_1el,state2d
! 		double precision thickness
! 		character(16),pointer :: name_d_data(:)




end subroutine

! Vfのパラメータ値をセット
!subroutine set_vf_params(model,n_periods_x,n_periods_y,shift,mat_no_warp,mat_no_weft)
subroutine set_vf_params(model,n_periods_x,n_periods_y,shift)
  use fem_module
  implicit none
  type(struct_model) :: model
!  integer :: shift(:),mat_no_warp(:,:),mat_no_weft(:,:)
integer :: shift(:)
integer mat_no_warp(4,4),mat_no_weft(4,4)
  integer n_els_1period,n_els_1period_x_weft,n_els_1period_x_warp,n_periods_x,n_periods_y,n_els_1lamina,param_no
  integer block_no,cnt_warp,cnt_weft,warp_no,mat_no_1weft,vfmicro_no
  integer n_els
  integer,allocatable :: offset_warp(:),offset_weft(:)
  integer shift2cnm_warp(52),shift2cnm_weft(52)
  integer i,j
  integer,pointer :: material_nos(:)
mat_no_warp = 0

shift2cnm_warp(1:14) = (/0,1,2,3,4,5,6,7,8,9,10,11,12,15/)*2
shift2cnm_warp(15:27) = (/18,19,20,21,22,23,24,25,26,27,28,29,30/)*2
shift2cnm_warp(28:52) = shift2cnm_warp(2:26)+60

shift2cnm_weft(1:14) = (/0,1,2,3,4,5,6,7,8,9,10,11,12,12/)*2
shift2cnm_weft(15:27) = (/12,13,14,15,16,17,18,19,20,21,22,23,24/)*2
shift2cnm_weft(28:52) = shift2cnm_weft(2:26)+48

  n_els_1period_x_warp = 120
  n_els_1period_x_weft = 96

  n_els_1period = 704
  n_els_1lamina = n_els_1period*n_periods_x
  warp_no = 2
  !weft_no = 3

  n_els = model%n_els

  material_nos => model%material_nos


  allocate(offset_warp(size(shift)))
  do i=1,size(shift)
    offset_warp(i) = n_els_1period_x_warp/2 - mod((n_els_1period_x_warp/4+shift2cnm_warp(shift(i)+1)),n_els_1period_x_warp/2)
  end do
  allocate(offset_weft(size(shift)))
  do i=1,size(shift)
    offset_weft(i) = n_els_1period_x_weft/2 - mod((n_els_1period_x_weft/4+shift2cnm_weft(shift(i)+1)),n_els_1period_x_weft/2)
  end do


  do i=1,n_periods_y
    cnt_warp = 1
    cnt_weft = 1
    block_no = 0
    do j=1,n_els_1lamina
      if (material_nos((i-1)*n_els_1lamina+j) == 2) then
        param_no = ((cnt_warp-offset_warp(i) + n_els_1period_x_warp/2-1)/(n_els_1period_x_warp/2))+1
        if (param_no > 2*n_periods_x) then
          param_no = 1
        end if
        vfmicro_no = mod(cnt_warp,2)

        material_nos((i-1)*n_els_1lamina+j) = vfmicro_no+param_no*10+100
        cnt_warp = cnt_warp + 1
      else if (material_nos((i-1)*n_els_1lamina+j) == 3) then
        param_no = ((cnt_weft-offset_weft(i) + n_els_1period_x_weft/2-1)/(n_els_1period_x_weft/2))+1
        if (param_no > 2*n_periods_x) then
          param_no = 1
        end if

        mat_no_1weft = mod( (cnt_weft + n_els_1period_x_weft/2 - offset_weft(i)),n_els_1period_x_weft/2)
        if (mat_no_1weft == 0) then
          mat_no_1weft = n_els_1period_x_weft/2
        end if

        vfmicro_no = ((mat_no_1weft-1)/4+1)*2
        if ((mod(mat_no_1weft,4) == 1) .or. (mod(mat_no_1weft,4) == 3)) then
          vfmicro_no = vfmicro_no - 1
        end if
        
        !material_nos((i-1)*n_els_1lamina+j) = 20+vfmicro_no
        material_nos((i-1)*n_els_1lamina+j) = 20+param_no
        cnt_weft = cnt_weft + 1
      end if 
    end do 
  end do

end subroutine


! program test
! 	use GFRPMeshGenerator2D
! 	double precision,allocatable,dimension(:,:) :: nodes
! 	integer,allocatable,dimension(:,:) :: elements
! 	integer,parameter :: NumOfNodesInOneElement = 4
! 	integer, parameter :: dim = 2;
! 	integer,parameter :: NumOfNodes = 787
! 	integer,parameter :: NumOfElements = 720
! 	double precision :: thL,tlL,wL,thC,tlC,wC,dxC,thR,tlR,wR

! 	thL = 0.11d-3
! 	thC = 0.11d-3
! 	thR = 0.11d-3
! 	wL = 4.04d-3
! 	wC = 4.04d-3
! 	wR = 4.04d-3
! 	dxC = 0.00d-3

! 	allocate(nodes(dim,NumOfNodes))
! 	allocate(elements(NumOfNodesInOneElement+1,NumOfElements))

! 	call generateGFRPMesh(nodes,elements,thL,wL,thC,wC,dxC,thR,wR)

! end program test

end module



! program test
! 	use GFRPMeshGenerator2D
! 	double precision,allocatable,dimension(:,:) :: nodes
! 	integer,allocatable,dimension(:,:) :: elements
! 	integer,parameter :: NumOfNodesInOneElement = 4
! 	integer, parameter :: dim = 2;
! 	integer,parameter :: NumOfNodes = 787
! 	integer,parameter :: NumOfElements = 720
! 	double precision :: thL,tlL,wL,thC,tlC,wC,dxC,thR,tlR,wR

! 	thL = 0.11d-3
! 	thC = 0.11d-3
! 	thR = 0.11d-3
! 	wL = 4.04d-3
! 	wC = 4.04d-3
! 	wR = 4.04d-3
! 	dxC = 0.00d-3

! 	allocate(nodes(dim,NumOfNodes))
! 	allocate(elements(NumOfNodesInOneElement+1,NumOfElements))

! 	call generateGFRPMesh(nodes,elements,thL,wL,thC,wC,dxC,thR,wR)

! end program test

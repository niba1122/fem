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
! 	integer,parameter :: NumOfNodes = 787
	integer,parameter :: NumOfElements = 720

	double precision,parameter :: widthOfModel = 9.12d-3
	double precision,parameter :: thicknessOfInterface = 0.04d-3
	double precision,parameter :: minThicknessOfMatrix = 0.02d-3
	double precision,parameter :: ThicknessOfWarp = 0.22d-3
	double precision,parameter :: distalThicknessOfWeft = 0.04d-3
	double precision,parameter :: thicknessOfWeft = 0.22d-3
	double precision,parameter :: widthOfBlock = widthOfModel/54


	integer,dimension(55) :: bottomNodes,topNodes
	double precision,allocatable,dimension(:,:) :: shiftedNodes
	integer,intent(in) :: shift

	double precision :: thL,tlL,wL,thC,tlC,wC,dxC,thR,tlR,wR

	NumOfNodes = 787
	deallocate(nodes)
	allocate(nodes(dim,NumOfNodes))
	allocate(shiftedNodes(dim,NumOfNodes))

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

	call addBlock3LC(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
		&minThicknessOfMatrix,ThicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
 	print *,"Block3LC was generated"

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

 	call addBlock3CR(nodes,elements,lastNode,lastElement,widthOfModel,widthOfBlock,thicknessOfInterface,&
		&minThicknessOfMatrix,ThicknessOfWarp,thicknessOfWeft,thL,wL,thC,wC,dxC,thR,wR,distalThicknessOfWeft)
 	print *,"Block3CR was generated"

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

	NumOfNodes = 775+topNodes(shift+1)-bottomNodes(shift+1)
	deallocate(nodes)
	allocate(nodes(dim,NumOfNodes))
	print *,"numofnodes: ", NumOfNodes

	nodes(:,:) = shiftedNodes(:,1:NumOfNodes)


	do i=1,NumOfElements
		if (elements(1,i) < bottomNodes(shift+1)) then
			elements(1:4,i) = elements(1:4,i)+NumOfNodes-topNodes(shift+1)
		else
			elements(1:4,i) = elements(1:4,i)-bottomNodes(shift+1)+1
		end if
	end do

! 	do i=1,NumOfElements
! 		if (elements(1,i) < bottomNodes(shift+1)) then 
! 			elements(1:4,i) = elements(1:4,i) + topNodes(shift+1)
! 		end if
! 	end do

end if

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
	nodes(2,lastNode+2) = -(tlL-distalThicknessOfWeft/2)*dcos((x0+widthOfBlock)*pi/wL)&
												&+tlL+minThicknessOfMatrix-distalThicknessOfWeft/2
	nodes(2,lastNode+3) = nodes(2,lastNode+2)+thicknessOfInterface/2
	nodes(2,lastNode+4) = nodes(2,lastNode+3)+thicknessOfInterface/2
	nodes(2,lastNode+5) = tlL+minThicknessOfMatrix+thicknessOfInterface

	! upper part of Weft
	nodes(2,lastNode+6) = (thL-distalThicknessOfWeft/2)*dcos((x0+widthOfBlock)*pi/wL)&
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
	nodes(1,lastNode+4:lastNode+13) = nodes(1,lastNode-11)+thicknessOfInterface
	nodes(1,lastNode+14:lastNode+22) = nodes(1,lastNode-11)+thicknessOfInterface*2
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
	nodes(1,lastNode+1:lastNode+9) = (widthOfModel/2 + dxC - wC/2) - thicknessOfInterface*2
	nodes(1,lastNode+10:lastnode+19) = (widthOfModel/2 + dxC - wC/2) - thicknessOfInterface
	nodes(1,lastNode+20:lastnode+22) = (widthOfModel/2 + dxC - wC/2) - thicknessOfInterface/2
 	nodes(1,lastNode+23) = nodes(1,lastNode) + widthOfBlock
	nodes(1,lastNode+24:lastnode+34) = widthOfModel/2 + dxC - wC/2
 	nodes(1,lastNode+35) = nodes(1,lastNode) + widthOfBlock

! y of nodes
	gradient = (-thicknessOfInterface + thL - thC - distalThicknessOfWeft) / (widthOfModel/2 + dxC - wL/2 - wC/2)



	nodes(2,lastNode+2:lastNode+8) = nodes(2,lastNode-7:lastNode-1) &
								& + ((widthOfModel/2 + dxC - wC/2 - thicknessOfInterface*2) - nodes(1,lastNode) )*gradient
	nodes(2,lastNode+1) = nodes(2,lastNode+2)/2d0
 	nodes(2,lastNode+9) = thicknessOfModel &
 				& - minThicknessOfMatrix - thicknessOfInterface - thC + distalThicknessOfWeft/2
 	nodes(2,lastNode+11:lastNode+17) = nodes(2,lastNode-7:lastNode-1) &
								& + ((widthOfModel/2 + dxC - wC/2 - thicknessOfInterface) - nodes(1,lastNode) )*gradient
	nodes(2,lastNode+10) = nodes(2,lastNode+11)/2d0
 	nodes(2,lastNode+18) = nodes(2,lastNode+9)
 	nodes(2,lastNode+19) = nodes(2,lastNode+9) + thicknessOfInterface
 	nodes(2,lastNode+20) = nodes(2,lastNode+9) - thicknessOfInterface/2
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
	nodes(2,lastNode+12) = (thC - distalThicknessOfWeft/2) * dcos(pi/wC*(x0+widthOfBlock-(widthOfModel/2+dxC)))&
								& + thicknessOfModel - minThicknessOfMatrix - thC + distalThicknessOfWeft/2
	nodes(2,lastNode+11) = nodes(2,lastNode+12) - thicknessOfInterface/2
	nodes(2,lastNode+10) = nodes(2,lastNode+11) - thicknessOfInterface/2
	nodes(2,lastNode+9) = thicknessOfModel - minThicknessOfMatrix - thicknessOfInterface&
								& - thC
 	! lower part of Weft
	nodes(2,lastNode+8) = - (tlC - distalThicknessOfWeft/2) * dcos(pi/wC*(x0+widthOfBlock-(widthOfModel/2+dxC)))&
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
	do i=4,13
		nodes(1,lastNode+i) = nodes(1,lastNode-1) + thicknessOfInterface
	end do
	do i=14,22
		nodes(1,lastNode+i) = nodes(1,lastNode-1) + thicknessOfInterface*2
	end do
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
	nodes(1,lastNode+1:lastNode+9) = widthOfModel - wR/2 - thicknessOfInterface*2
	nodes(1,lastNode+10:lastNode+19) = widthOfModel - wR/2 - thicknessOfInterface
	nodes(1,lastNode+20:lastNode+22) = widthOfModel - wR/2 - thicknessOfInterface/2
	nodes(1,lastNode+23) = nodes(1,lastNode) + widthOfBlock
	nodes(1,lastNode+24:lastNode+34) = widthOfModel - wR/2
	nodes(1,lastNode+35) = nodes(1,lastNode) + widthOfBlock

! y of nodes
	gradient = (thicknessOfInterface - thR + thC + distalThicknessOfWeft) / (widthOfModel/2 - dxC - wR/2 - wC/2)

 	nodes(2,lastNode+1) = minThicknessOfMatrix + tlR + distalThicknessOfWeft/2
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
 	nodes(2,lastNode+22) = nodes(2,lastNode+1) + thicknessOfInterface/2
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
	nodes(2,lastNode+2) = -(tlR-distalThicknessOfWeft/2)*dcos((x0+widthOfBlock-widthOfModel)*pi/wR)&
												&+tlR+minThicknessOfMatrix-distalThicknessOfWeft/2
	nodes(2,lastNode+3) = nodes(2,lastNode+2)+thicknessOfInterface/2
	nodes(2,lastNode+4) = nodes(2,lastNode+3)+thicknessOfInterface/2
	nodes(2,lastNode+5) = tlR+minThicknessOfMatrix+thicknessOfInterface

	! upper part of Weft
	nodes(2,lastNode+6) = (thR-distalThicknessOfWeft/2)*dcos((x0+widthOfBlock-widthOfModel)*pi/wR)&
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
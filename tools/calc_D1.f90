program calc_d1
	integer i
	double precision E,v
	double precision :: D(6,6),D1(6,6),C1(6,6)

	print *,"Young's modulus?"
	read (*,*,iostat=iostat) E; print *
	if (iostat /= 0) then
		print *,"Error: Bad value."
		print *,""
		print *,"Press Enter to finish."
		read *
		stop
	end if

	print *,"Poisson's ratio?"
	read (*,*,iostat=iostat) v; print *
	if (iostat /= 0) then
		print *,"Error: Bad value."
		print *,""
		print *,"- Press Enter to finish -"
		read *
		stop
	end if

	D = 0d0
	D(1,:) = (/1-v, v, v, 0d0, 0d0, 0d0/)
	D(2,:) = (/v, 1-v, v, 0d0, 0d0, 0d0/)
	D(3,:) = (/v, v, 1-v, 0d0, 0d0, 0d0/)
	D(4,4) = (1-2*v)/2
	D(5,5) = (1-2*v)/2
	D(6,6) = (1-2*v)/2
	D(:,:) = D(:,:)*E / ((1+v)*(1-2*v))

	C1 = 0d0
	C1(1,1:3) = (/0d0,-1d0,-1d0/)
	C1(2,1:3) = (/-1d0,0d0,-1d0/)
	C1(3,1:3) = (/-1d0,-1d0,0d0/)
	C1(4,4) = 2d0
	C1(5,5) = 2d0
	C1(6,6) = 2d0
	C1 = C1*v/E

	D1 = -matmul(matmul(D,C1),D)

	open(10,file='D.csv')
		do i=1,6
			write(10,'(f30.15,5(",",f30.15))') D(i,:)
		end do
		write(10,'(a)') ""
		do i=1,6
			write(10,'(f30.15,5(",",f30.15))') D1(i,:)
		end do

		D = D - 0.1d0*D1
		write(10,'(a)') ""
		do i=1,6
			write(10,'(f30.15,5(",",f30.15))') D(i,:)
		end do
	close(10)

	print *,"D matrix was made successfully!"
	print *,""
	print *,"- Press Enter to finish -"

	read *

end program
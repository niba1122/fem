program D_matrix_maker
	integer i
	double precision E,v
	double precision :: D(6,6)

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

	open(10,file='D.csv')
		do i=1,6
			write(10,'(E21.15e2,6(",",E21.15e2))') D(i,:)
!			write(10,'(E9.3e2,6(",",E9.3e2))') D(i,:)
		end do
	close(10)

	print *,"D matrix was made successfully!"
	print *,""
	print *,"- Press Enter to finish -"

	read *

end program

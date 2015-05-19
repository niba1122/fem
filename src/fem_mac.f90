program fem
  implicit none
  integer file_num, keyword_pos, iostat, i
  integer length, status
  character(256) :: command
  character(32) :: tmp
  character(64) :: file_path
  character(:),allocatable :: solver_name
  intrinsic :: command_argument_count, get_command_argument

  if (command_argument_count() > 0) then
    call get_command_argument(1, length=length, status=status)

    allocate(character(length) :: solver_name)
    call get_command_argument(1, solver_name, status=status)

    open(10,file="../solvers/"//trim(solver_name)//"/sources.dat",status='old', iostat=iostat)
    if (iostat /= 0) call error("Cannot read solvers/"//solver_name//"/sources.dat")
    read (10,*,iostat=iostat) file_num
    if (iostat /= 0) call error("sources.dat is badly formated")

    command = "gfortran -o ../solvers/"//solver_name//"/"//solver_name//".out" 
    do i=1,file_num
      read (10,*,iostat=iostat) tmp
      if (iostat /= 0) call error("sources.dat is badly formated")

      keyword_pos = index(tmp, "module:") 
      if (keyword_pos>0) then
        file_path = "../modules/"//trim(tmp(keyword_pos+7:32))//".f90"
      else
        file_path = "../solvers/"//solver_name//"/"//tmp
      end if

      command = trim(command)//" "//trim(file_path)

    end do

    print *,trim(command)
    call system(command)
    call system("rm *.mod")
    close(10)
  end if

contains

subroutine error(msg)
	character(*) msg

	print *,"ERROR: ", msg ; print *
 	print *,"- Press ENTER to finish -"
	read *

	stop
end subroutine

end program

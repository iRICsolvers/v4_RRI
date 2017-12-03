module iric
  implicit none
  private
  character(len=64):: cgns_name
  integer:: cgns_f

  public iric_cgns_open, iric_cgns_close
  public iric_read_input_condition

contains
  subroutine iric_cgns_open(cname, ierr)
    include 'cgnslib_f.h'
    implicit none

    character(len=*),intent(in):: cname
    integer,intent(out):: ierr

    cgns_name = cname
    call cg_open_f(cgns_name, CG_MODE_MODIFY, cgns_f, ierr)
    call cg_iric_init_f(cgns_f, ierr)

  end subroutine

  subroutine iric_cgns_close()
    implicit none
    integer:: ierr

    call cg_close_f(cgns_f, ierr)
  end subroutine

  subroutine iric_read_input_condition()
    call iric_read_grid()
  end subroutine

  subroutine iric_read_grid()
    use globals

    integer:: isize, jsize, ierr
    double precision, dimension(:,:), allocatable:: grid_x, grid_y

    call cg_iric_gotogridcoord2d_f(isize, jsize, ierr)
    if (ierr /= 0) stop "CGNS grid read error"
    allocate(grid_x(isize, jsize), grid_y(isize, jsize))
    call cg_iric_getgridcoord2d_f(grid_x, grid_y, ierr)

    xllcorner = grid_x(1, 1)
    yllcorner = grid_y(1, 1)
    cellsize = grid_x(2, 1) - grid_x(1, 1)

    ! ONLY FOR DEBUGGING
    print *, "ISIZE, JSIZE = ", isize, jsize
    print *, "XLLCORNER, YLLCORNER = ", xllcorner, yllcorner
    print *, "CELLSIZE = ", cellsize

    deallocate(grid_x, grid_y)
  end subroutine
end module

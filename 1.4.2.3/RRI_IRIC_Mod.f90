module iric
  implicit none
  private
  character(len=64):: cgns_name
  integer:: cgns_f

  public iric_cgns_open, iric_cgns_close
  public iric_read_input_condition

contains
  subroutine iric_cgns_open()
    implicit none
    include 'cgnslib_f.h'

    integer:: ierr

    cgns_name = "Case1.cgn"
    call cg_open_f(cgns_name, CG_MODE_MODIFY, cgns_f, ierr)
    if (ierr /= 0) stop "cg_open_f failed"
    call cg_iric_init_f(cgns_f, ierr)
    if (ierr /= 0) stop "cg_iric_init_f failed"

  end subroutine

  subroutine iric_cgns_close()
    implicit none
    integer:: ierr

    call cg_close_f(cgns_f, ierr)
  end subroutine

  subroutine iric_read_input_condition()
    call iric_read_grid()
  end subroutine

  subroutine iric_read_cell_attr_int(name, v)
    use globals
    implicit none

    character(len=*),intent(in):: name
    integer, dimension(:,:), allocatable, intent(out):: v

    integer, dimension(:,:), allocatable:: tmpv
    integer:: i, j, ierr

    allocate(tmpv(nx, ny))
    call cg_iric_read_grid_integer_cell_f(cgns_f, name, tmpv, ierr)

    do i = 1, nx
      do j = 1, ny
        v(j, i) = tmpv(i, j)
      end do
    end do
  end subroutine

  subroutine iric_read_cell_attr_real(name, v)
    use globals
    implicit none

    character(len=*),intent(in):: name
    double precision, dimension(:,:), allocatable, intent(out):: v

    double precision, dimension(:,:), allocatable:: tmpv
    integer:: i, j, ierr

    allocate(tmpv(nx, ny))
    call cg_iric_read_grid_real_cell_f(cgns_f, name, tmpv, ierr)

    do i = 1, nx
      do j = 1, ny
        v(j, i) = tmpv(i, j)
      end do
    end do
  end subroutine

  subroutine iric_write_result_int(name, v)
    use globals
    implicit none

    character(len=*),intent(in):: name
    integer, dimension(:,:), allocatable, intent(in):: v

    integer, dimension(:,:), allocatable:: tmpv
    integer:: i, j, ierr

    allocate(tmpv(nx, ny))

    do i = 1, nx
      do j = 1, ny
        tmpv(i, j) = v(j, i)
      end do
    end do

    call cg_iric_write_sol_cell_integer_f(cgns_f, name, tmpv, ierr)
  end subroutine

  subroutine iric_write_result_real(name, v)
    use globals
    implicit none

    character(len=*),intent(in):: name
    double precision, dimension(:,:), allocatable, intent(in):: v

    double precision, dimension(:,:), allocatable:: tmpv
    integer:: i, j, ierr

    allocate(tmpv(nx, ny))

    do i = 1, nx
      do j = 1, ny
        tmpv(i, j) = v(j, i)
      end do
    end do

    call cg_iric_write_sol_cell_integer_f(cgns_f, name, tmpv, ierr)
  end subroutine

  subroutine iric_read_grid()
    use globals

    integer:: isize, jsize, ierr
    double precision, dimension(:,:), allocatable:: grid_x, grid_y

    call cg_iric_gotogridcoord2d_f(isize, jsize, ierr)
    if (ierr /= 0) stop "CGNS grid read error"
    allocate(grid_x(isize, jsize), grid_y(isize, jsize))
    call cg_iric_getgridcoord2d_f(grid_x, grid_y, ierr)

    ! iRIC から読み込んだ格子はメートル単位の座標系
    utm = 1
    xllcorner = grid_x(1, 1)
    yllcorner = grid_y(1, 1)
    cellsize = grid_x(2, 1) - grid_x(1, 1)
    nx = isize - 1
    ny = jsize - 1

    allocate(zs(ny, nx), zb(ny, nx), zb_riv(ny, nx), domain(ny, nx))
    allocate(riv(ny, nx), acc(ny, nx))
    allocate(dir(ny, nx))
    allocate(land(ny, nx))

    call iric_read_cell_attr_real("Elevation", zs)
    call iric_read_cell_attr_int("Acc", acc)
    call iric_read_cell_attr_int("Dir", dir)
    land = 1
    if (land_switch.eq.1) then
      call iric_read_cell_attr_int("Land", land)
    endif

    ! ONLY FOR DEBUGGING
    print *, "ISIZE, JSIZE = ", isize, jsize
    print *, "XLLCORNER, YLLCORNER = ", xllcorner, yllcorner
    print *, "CELLSIZE = ", cellsize

    deallocate(grid_x, grid_y)
  end subroutine

  subroutine iric_write_qu(qs_ave)
    use globals
    implicit none

    double precision, dimension(:,:,:), allocatable, intent(in):: qs_ave
    double precision, dimension(:,:), allocatable:: v
    integer:: i, j

    allocate(v(ny, nx))
    do i = 1, ny
      do j = 1, nx
        v(i, j) = ((qs_ave(1, i, j) + (qs_ave(3, i, j) - qs_ave(4, i, j)) / 2.d0) * area)
      end do
    end do
    call iric_write_result_real('qu', v)
    deallocate(v)
  end subroutine

  subroutine iric_write_qv(qs_ave)
    use globals
    implicit none

    double precision, dimension(:,:,:), allocatable, intent(in):: qs_ave
    double precision, dimension(:,:), allocatable:: v
    integer:: i, j

    allocate(v(ny, nx))
    do i = 1, ny
      do j = 1, nx
        v(i, j) = ((qs_ave(2, i, j) + (qs_ave(3, i, j) + qs_ave(4, i, j)) / 2.d0) * area)
      end do
    end do
    call iric_write_result_real('qv', v)
    deallocate(v)
  end subroutine

  subroutine iric_write_gu(qg_ave)
    use globals
    implicit none

    double precision, dimension(:,:,:), allocatable, intent(in):: qg_ave
    double precision, dimension(:,:), allocatable:: v
    integer:: i, j

    allocate(v(ny, nx))
    do i = 1, ny
      do j = 1, nx
        v(i, j) = ((qg_ave(1, i, j) + (qg_ave(3, i, j) - qg_ave(4, i, j)) / 2.d0) * area)
      end do
    end do
    call iric_write_result_real('gu', v)
    deallocate(v)
  end subroutine

  subroutine iric_write_gv(qg_ave)
    use globals
    implicit none

    double precision, dimension(:,:,:), allocatable, intent(in):: qg_ave
    double precision, dimension(:,:), allocatable:: v
    integer:: i, j

    allocate(v(ny, nx))
    do i = 1, ny
      do j = 1, nx
        v(i, j) = ((qg_ave(2, i, j) + (qg_ave(3, i, j) + qg_ave(4, i, j)) / 2.d0) * area)
      end do
    end do
    call iric_write_result_real('gu', v)
    deallocate(v)
  end subroutine

  subroutine iric_cgns_output_result( &
    hs,hr,hg,ar_ave,qs_ave,qg_ave)
    use globals

    double precision, dimension(:,:), allocatable, intent(in):: &
      hs, hr, hg, ar_ave
    double precision, dimension(:,:,:), allocatable, intent(in):: &
      qs_ave, qg_ave
    integer:: ierr

    call cg_iric_write_sol_time_f(cgns_f, time, ierr)

    if (outswitch_hs /= 0) then
      call iric_write_result_real('hs', hs)
    end if
    if (outswitch_hr /= 0) then
      call iric_write_result_real('hr', hr)
    end if
    if (outswitch_hg /= 0) then
      call iric_write_result_real('hg', hg)
    end if
    if (outswitch_qr /= 0) then
      call iric_write_result_real('qr', ar_ave)
    end if
    if (outswitch_qu /= 0) then
      call iric_write_qu(qs_ave)
    end if
    if (outswitch_qv /= 0) then
      call iric_write_qv(qs_ave)
    end if
    if (outswitch_gu /= 0) then
      call iric_write_gu(qg_ave)
    end if
    if (outswitch_gv /= 0) then
      call iric_write_gv(qg_ave)
    end if
    if (outswitch_gampt_ff /= 0) then
      call iric_write_result_real('gampt_ff', gampt_ff)
    end if
  end subroutine

end module

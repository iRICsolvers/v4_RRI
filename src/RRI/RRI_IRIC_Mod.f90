module RRI_iric

    use iric

    implicit none
    character(len=64):: cgns_name
    integer:: cgns_f, icount

    public iric_cgns_open, iric_cgns_close
    public iric_read_input_condition
    public iric_cgns_output_result
    public iric_write_cell_real, iric_write_cell_integer
    public iric_read_cell_attr_int, iric_read_cell_attr_real

contains
    subroutine iric_cgns_open()
        implicit none

        integer:: ierr

        !引数取得
        !icount = nargs()
        icount = iargc()
        !if (icount ==  2) then
        if (icount == 1) then
            !call getarg(1, cgns_name, ierr)
            call getarg(1, cgns_name)
        else
            write (*, "(a)") "You should specify an argument."
            stop
        end if

        call cg_iric_open(cgns_name, IRIC_MODE_MODIFY, cgns_f, ierr)
        if (ierr /= 0) stop "cg_iric_open failed"

        ! guiにcgnsファイルを読込みであることを知らせるファイルを生成
        call iric_initoption(IRIC_OPTION_CANCEL, ierr)

    end subroutine

    subroutine iric_cgns_close()
        implicit none
        integer:: ierr

        call cg_iric_close(cgns_f, ierr)
    end subroutine

    subroutine iric_read_input_condition()
        call iric_read_grid()
    end subroutine

    subroutine iric_read_cell_attr_int(name, v)
        use globals
        implicit none

        character(len=*), intent(in):: name
        integer, intent(out):: v(ny, nx)

        integer, dimension(:, :), allocatable:: tmpv
        integer:: i, j, jj, ierr

        allocate (tmpv(nx, ny))
        call cg_iric_read_grid_integer_cell(cgns_f, name, tmpv, ierr)

        do i = 1, nx
            do j = 1, ny
                jj = ny - j + 1
                v(j, i) = tmpv(i, jj)
            end do
        end do
    end subroutine

    subroutine iric_read_cell_attr_real(name, v)
        use globals
        implicit none

        character(len=*), intent(in):: name
        double precision, intent(out):: v(ny, nx)

        double precision, dimension(:, :), allocatable:: tmpv
        integer:: i, j, jj, ierr

        allocate (tmpv(nx, ny))
        call cg_iric_read_grid_real_cell(cgns_f, name, tmpv, ierr)

        do i = 1, nx
            do j = 1, ny
                jj = ny - j + 1
                v(j, i) = tmpv(i, jj)
            end do
        end do
    end subroutine

    subroutine iric_write_result_real(name, v)
        use globals
        implicit none

        character(len=*), intent(in):: name
        double precision, dimension(:, :), allocatable, intent(in):: v

        double precision, dimension(:, :), allocatable:: tmpv
        integer:: i, j, ierr

        allocate (tmpv(nx, ny))

        do i = 1, nx
            do j = 1, ny
                tmpv(i, j) = v(ny - j + 1, i)
                !tmpv(i, j) = v(j, i)
                !tmpv(i, j) = v(i, j)
            end do
        end do

        call cg_iric_write_sol_cell_real(cgns_f, name, tmpv, ierr)
    end subroutine

    subroutine iric_write_result_integer(name, v)
        use globals
        implicit none

        character(len=*), intent(in):: name
        integer, dimension(:, :), allocatable, intent(in):: v

        integer, dimension(:, :), allocatable:: tmpv
        integer:: i, j, ierr

        allocate (tmpv(nx, ny))

        do i = 1, nx
            do j = 1, ny
                tmpv(i, j) = v(ny - j + 1, i)
                !tmpv(i, j) = v(j, i)
                !tmpv(i, j) = v(i, j)
            end do
        end do

        call cg_iric_write_sol_cell_integer(cgns_f, name, tmpv, ierr)
    end subroutine

    subroutine iric_write_cell_real(name, v)
        use globals
        implicit none

        character(len=*), intent(in):: name
        double precision, dimension(:, :), allocatable, intent(in):: v

        double precision, dimension(:, :), allocatable:: tmpv
        integer:: i, j, ierr

        allocate (tmpv(nx, ny))

        do i = 1, nx
            do j = 1, ny
                tmpv(i, j) = v(ny - j + 1, i)
                !tmpv(i, j) = v(j, i)
                !tmpv(i, j) = v(i, j)
            end do
        end do

        call cg_iric_write_grid_real_cell(cgns_f, name, tmpv, ierr)
    end subroutine

    subroutine iric_write_cell_integer(name, v)
        use globals
        implicit none

        character(len=*), intent(in):: name
        integer, dimension(:, :), allocatable, intent(in):: v

        integer, dimension(:, :), allocatable:: tmpv
        integer:: i, j, ierr

        allocate (tmpv(nx, ny))

        do i = 1, nx
            do j = 1, ny
                tmpv(i, j) = v(ny - j + 1, i)
                !tmpv(i, j) = v(j, i)
                !tmpv(i, j) = v(i, j)
            end do
        end do

        call cg_iric_write_grid_integer_cell(cgns_f, name, tmpv, ierr)
    end subroutine

    subroutine iric_read_grid()
        use globals

        integer:: isize, jsize, ierr
        double precision, dimension(:, :), allocatable:: grid_x, grid_y

        call cg_iRIC_Read_Grid2d_Str_Size(cgns_f, isize, jsize, ierr)
        print *, "isize, jsize", isize, jsize
        if (ierr /= 0) stop "CGNS grid read error"
        allocate (grid_x(isize, jsize), grid_y(isize, jsize))
        print *, "grid_x, grid_y allocated"
        call cg_iRIC_Read_Grid2d_Coords(cgns_f, grid_x, grid_y, ierr)
        print *, "cg_iRIC_Read_Grid2d_Coords called"

        ! iRIC から読み込んだ格子はメートル単位の座標系
        ! utm = 1
        xllcorner = grid_x(1, 1)
        yllcorner = grid_y(1, 1)
        cellsize = grid_x(2, 1) - grid_x(1, 1)
        nx = isize - 1
        ny = jsize - 1

        allocate (zs(ny, nx), zb(ny, nx), zb_riv(ny, nx), domain(ny, nx))
        allocate (riv(ny, nx), acc(ny, nx))
        allocate (dir(ny, nx))
        allocate (land(ny, nx))

        print *, "iric_read_cell_attr_real Elevation"
        call iric_read_cell_attr_real("Elevation", zs)
        print *, "iric_read_cell_attr_int Acc"
        call iric_read_cell_attr_int("Acc", acc)
        print *, "iric_read_cell_attr_int Dir"
        call iric_read_cell_attr_int("Dir", dir)
        land = 1
        if (land_switch .eq. 1) then
            print *, "iric_read_cell_attr_int Land"
            call iric_read_cell_attr_int("Land", land)
        end if
        print *, "iric_read_cell_attr_real all ok"

        ! ONLY FOR DEBUGGING
        print *, "ISIZE, JSIZE = ", isize, jsize
        print *, "XLLCORNER, YLLCORNER = ", xllcorner, yllcorner
        print *, "CELLSIZE = ", cellsize

        deallocate (grid_x, grid_y)
    end subroutine

    subroutine iric_write_qu(qs_ave)
        use globals
        implicit none

        double precision, dimension(:, :, :), allocatable, intent(in):: qs_ave
        double precision, dimension(:, :), allocatable:: v
        integer:: i, j

        allocate (v(ny, nx))
        do i = 1, ny
            do j = 1, nx
                v(i, j) = ((qs_ave(1, i, j) + (qs_ave(3, i, j) - qs_ave(4, i, j))/2.d0)*area)
            end do
        end do
        call iric_write_result_real('qu', v)
        deallocate (v)
    end subroutine

    subroutine iric_write_qv(qs_ave)
        use globals
        implicit none

        double precision, dimension(:, :, :), allocatable, intent(in):: qs_ave
        double precision, dimension(:, :), allocatable:: v
        integer:: i, j

        allocate (v(ny, nx))
        do i = 1, ny
            do j = 1, nx
                v(i, j) = ((qs_ave(2, i, j) + (qs_ave(3, i, j) + qs_ave(4, i, j))/2.d0)*area)
            end do
        end do
        call iric_write_result_real('qv', v)
        deallocate (v)
    end subroutine

    subroutine iric_write_gu(qg_ave)
        use globals
        implicit none

        double precision, dimension(:, :, :), allocatable, intent(in):: qg_ave
        double precision, dimension(:, :), allocatable:: v
        integer:: i, j

        allocate (v(ny, nx))
        do i = 1, ny
            do j = 1, nx
                v(i, j) = ((qg_ave(1, i, j) + (qg_ave(3, i, j) - qg_ave(4, i, j))/2.d0)*area)
            end do
        end do
        call iric_write_result_real('gu', v)
        deallocate (v)
    end subroutine

    subroutine iric_write_gv(qg_ave)
        use globals
        implicit none

        double precision, dimension(:, :, :), allocatable, intent(in):: qg_ave
        double precision, dimension(:, :), allocatable:: v
        integer:: i, j

        allocate (v(ny, nx))
        do i = 1, ny
            do j = 1, nx
                v(i, j) = ((qg_ave(2, i, j) + (qg_ave(3, i, j) + qg_ave(4, i, j))/2.d0)*area)
            end do
        end do
        call iric_write_result_real('gv', v)
        deallocate (v)
    end subroutine

    subroutine iric_cgns_output_result( &
        sum_qp_t, qp_t, hs, hr, hg, qr_ave, qs_ave, qg_ave)
        use globals

        double precision, dimension(:, :), allocatable, intent(in):: &
            sum_qp_t, qp_t, hs, hr, hg, qr_ave
        double precision, dimension(:, :, :), allocatable, intent(in):: &
            qs_ave, qg_ave

        double precision, dimension(:, :), allocatable :: rain_rate, rain_vol
        double precision, dimension(:, :), allocatable :: tmpv
        integer :: i, j
        integer:: ierr

        call cg_iric_write_sol_start(cgns_f, ierr)
        call cg_iric_write_sol_time(cgns_f, time, ierr)

        !node valueが出力されていないとグラフが表示できないための仮出力
        allocate (tmpv(1:ny + 1, 1:nx + 1))
        tmpv = 0.0d0
        call cg_iRIC_Write_Sol_Node_Real(cgns_f, "dummy", tmpv, ierr)
        deallocate (tmpv)

        allocate (rain_rate(1:ny, 1:nx), rain_vol(1:ny, 1:nx))
        do i = 1, ny
            do j = 1, nx
                rain_rate(i, j) = qp_t(i, j)*3600.d0*1000.d0
                rain_vol(i, j) = sum_qp_t(i, j)*1000.d0
            end do
        end do
        call iric_write_result_real('total_qp_t[mm]', rain_vol)
        call iric_write_result_real('qp_t[mm_h]', rain_rate)
        call iric_write_result_real('hs[m]', hs)
        call iric_write_result_real('hr[m]', hr)
        call iric_write_result_real('qr[m3_s]', qr_ave)
        call iric_write_qu(qs_ave)
        call iric_write_qv(qs_ave)
        call iric_write_result_real('hg[m]', hg)
        call iric_write_gu(qg_ave)
        call iric_write_gv(qg_ave)
        call iric_write_result_real('gampt_ff', gampt_ff)
        deallocate (rain_rate)
        deallocate (rain_vol)

        !if (outswitch_hg /= 0) then
        !  call iric_write_result_real('hg', hg)
        !end if
        !if (outswitch_gu /= 0) then
        !  call iric_write_gu(qg_ave)
        !end if
        !if (outswitch_gv /= 0) then
        !  call iric_write_gv(qg_ave)
        !end if
        !if (outswitch_gampt_ff /= 0) then
        !  call iric_write_result_real('gampt_ff', gampt_ff)
        !end if

        call cg_iric_write_sol_end(cgns_f, ierr)

    end subroutine

end module

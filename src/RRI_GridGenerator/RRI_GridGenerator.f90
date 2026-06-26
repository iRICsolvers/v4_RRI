! demAdjust2.f90
!
! coded by T.Sayama on April 28, 2010
! v2.1 Jan 10, 2013

! DEM Adjustment Program
!
program RRI_GridGenerator

    use RRI_iric
    use iric
    implicit none

    character(len=1024) :: infile_dem, infile_dir, infile_acc

    integer acc_thresh
    real(8) increment, carve, lift
    parameter(acc_thresh=0)    ! adjustment is done for cells with acc > acc_thresh
    parameter(lift=500.d0)     ! [m]
    parameter(carve=5.d0)      ! [m]
    parameter(increment=0.01d0) ! [m]

    integer utm

    integer nx, ny, nx_file, ny_file
    real(8) xllcorner, yllcorner, cellsize, nodata
    real(8) xll_file, yll_file, cellsize_file, nodata_file
    real(8) dem_lowest

    real(8) d1, d2, d3, d4
    real(8) x1, x2, y1, y2, dx, dy
    real(8) length, dis

    integer i, j, k, l, ii, jj, iii, jjj, ios

    real(8), dimension(:, :), allocatable :: dem, adem
    integer, dimension(:, :), allocatable :: dir, acc, riv
    integer, dimension(:, :), allocatable :: downstream_i, downstream_j
    real(8), dimension(:, :), allocatable :: downstream_distance
    integer, dimension(:, :), allocatable :: landuse, gsd_s, gsd_c
    real(8), dimension(:, :), allocatable :: width, depth, height
    integer, dimension(:, :), allocatable :: upstream
    real(8), dimension(:), allocatable :: total_length
    integer, dimension(:), allocatable :: upstream_i, upstream_j, upstream_i_temp, upstream_j_temp
    integer, dimension(:), allocatable :: upstream_order, upstream_order_work
    integer switch, s1_i, s1_j, s2_i, s2_j, numup
    integer s1_down_i, s1_down_j, adjustment_count, rise_count, fall_count

    integer riv_thresh
    real(8) width_param_c, width_param_s
    real(8) depth_param_c, depth_param_s
    integer height_limit_param
    real(8) height_param

    character*256 ctemp
    integer :: ier, icount
    real(8), allocatable, save :: gxx(:, :), gyy(:, :)
    !--------------------------------------------------

    ! Open the grid generation CGNS file and read conditions.
    call iric_cgns_open()
    call cg_iric_read_string(cgns_f, "demfile", infile_dem, ier)
    if (ier /= 0) error stop "Failed to read demfile"
    call cg_iric_read_string(cgns_f, "accfile", infile_acc, ier)
    if (ier /= 0) error stop "Failed to read accfile"
    call cg_iric_read_string(cgns_f, "dirfile", infile_dir, ier)
    if (ier /= 0) error stop "Failed to read dirfile"
    call cg_iric_read_integer(cgns_f, "utm", utm, ier)
    if (ier /= 0) error stop "Failed to read coordinate system"

    open (10, file=trim(infile_dem), status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Failed to open DEM file"
    open (20, file=trim(infile_dir), status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Failed to open DIR file"
    open (30, file=trim(infile_acc), status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Failed to open ACC file"

    ! STEP 1 : Reading File
    read (10, *) ctemp, nx
    read (10, *) ctemp, ny
    read (10, *) ctemp, xllcorner
    read (10, *) ctemp, yllcorner
    read (10, *) ctemp, cellsize
    read (10, *) ctemp, nodata

    allocate (dem(ny, nx), dir(ny, nx), acc(ny, nx), adem(ny, nx), upstream(ny, nx))
    allocate (downstream_i(ny, nx), downstream_j(ny, nx), downstream_distance(ny, nx))

    ! Start generate grid shape for CGNS output
    allocate (gxx(nx + 1, ny + 1), gyy(nx + 1, ny + 1))
    gxx(1, 1) = xllcorner; gyy(1, 1) = yllcorner
    do j = 2, ny + 1
        gxx(1, j) = gxx(1, 1)
        gyy(1, j) = gyy(1, j - 1) + cellsize
    end do

    do i = 2, nx + 1
        gxx(i, 1) = gxx(i - 1, 1) + cellsize
        gyy(i, 1) = gyy(1, 1)
    end do

    do i = 2, nx + 1
        do j = 2, ny + 1
            gxx(i, j) = gxx(i - 1, j) + cellsize
            gyy(i, j) = gyy(i, j - 1) + cellsize
        end do
    end do
    ! End generate grid shape for CGNS output

    rewind (10)

    read (10, '(a30)') ctemp

    read (10, '(a30)') ctemp

    read (10, '(a30)') ctemp

    read (10, '(a30)') ctemp

    read (10, '(a30)') ctemp

    read (10, '(a30)') ctemp

    do i = 1, ny
        read (10, *) (dem(i, j), j=1, nx)
    end do

    read (20, *) ctemp, nx_file
    read (20, *) ctemp, ny_file
    read (20, *) ctemp, xll_file
    read (20, *) ctemp, yll_file
    read (20, *) ctemp, cellsize_file
    read (20, *) ctemp, nodata_file
    call check_grid_header("DIR", nx, ny, xllcorner, yllcorner, cellsize, &
                           nx_file, ny_file, xll_file, yll_file, cellsize_file)

    do i = 1, ny
        read (20, *) (dir(i, j), j=1, nx)
    end do

    read (30, *) ctemp, nx_file
    read (30, *) ctemp, ny_file
    read (30, *) ctemp, xll_file
    read (30, *) ctemp, yll_file
    read (30, *) ctemp, cellsize_file
    read (30, *) ctemp, nodata_file
    call check_grid_header("ACC", nx, ny, xllcorner, yllcorner, cellsize, &
                           nx_file, ny_file, xll_file, yll_file, cellsize_file)

    do i = 1, ny
        read (30, *) (acc(i, j), j=1, nx)
    end do

    write (*, *) "Done STEP 1"

    ! STEP 2 : dx, dy

    ! d1 : South
    x1 = xllcorner
    y1 = yllcorner
    x2 = xllcorner + nx*cellsize
    y2 = yllcorner
    if (utm .eq. 0) call hubeny_sub(x1, y1, x2, y2, d1)

    ! d2 : North
    x1 = xllcorner
    y1 = yllcorner + ny*cellsize
    x2 = xllcorner + nx*cellsize
    y2 = yllcorner + ny*cellsize
    if (utm .eq. 0) call hubeny_sub(x1, y1, x2, y2, d2)

    ! d3 : West
    x1 = xllcorner
    y1 = yllcorner
    x2 = xllcorner
    y2 = yllcorner + ny*cellsize
    if (utm .eq. 0) call hubeny_sub(x1, y1, x2, y2, d3)

    ! d4 : East
    x1 = xllcorner + nx*cellsize
    y1 = yllcorner
    x2 = xllcorner + nx*cellsize
    y2 = yllcorner + ny*cellsize
    if (utm .eq. 0) call hubeny_sub(x1, y1, x2, y2, d4)

    if (utm .eq. 1) then
        dx = cellsize
        dy = cellsize
    else
        dx = (d1 + d2)/2.d0/real(nx)
        dy = (d3 + d4)/2.d0/real(ny)
    end if
    write (*, *) "dx [m] : ", dx, "dy [m] : ", dy

    length = (dx + dy)/2.d0

    write (*, *) "Done STEP 2"

    ! STEP 3 : Set "dir = 0"
    do i = 1, ny
        do j = 1, nx

            if (dem(i, j) .lt. -100.d0) cycle

            if (dir(i, j) .le. -1) then
                write (*, *) "dir(", i, ",", j, ") is detected."
                write (*, *) "dem(", i, ",", j, ") is set to be nodata"
                dem(i, j) = nodata
                write (*, *) "dir(", i, ",", j, ") is set to be zero"
                dir(i, j) = 0
            end if

            call down(dir, nx, ny, i, j, length, ii, jj, dis)
            if (ii .eq. 0 .or. jj .eq. 0 .or. ii .gt. ny .or. jj .gt. nx) then
                dir(i, j) = 0
                cycle
            end if
            if (dem(ii, jj) .lt. -100.d0) dir(i, j) = 0

        end do
    end do
    write (*, *) "Done STEP 3"

    ! STEP 4 : Slope Calculation (Steepest Slope Among the Eight Directions)
    ! The slope array was never used by the DEM adjustment or grid output.
    ! Precompute downstream references here instead of recalculating them
    ! repeatedly along every upstream path.
    do i = 1, ny
        do j = 1, nx
            call down(dir, nx, ny, i, j, length, downstream_i(i, j), &
                      downstream_j(i, j), downstream_distance(i, j))
        end do
    end do

    write (*, *) "Done STEP 4"

    ! STEP 5 : Find The Most Upstream Cell
    upstream(:, :) = 1
    do i = 1, ny
        do j = 1, nx
            if (dem(i, j) .lt. -100.d0 .or. dir(i, j) .eq. 0 .or. acc(i, j) .lt. acc_thresh) then
                upstream(i, j) = 0
                cycle
            end if
            ii = downstream_i(i, j)
            jj = downstream_j(i, j)
            upstream(ii, jj) = 0
        end do
    end do
    numup = sum(upstream(:, :))

    allocate (upstream_i(numup), upstream_j(numup))
    allocate (upstream_i_temp(numup), upstream_j_temp(numup))
    allocate (upstream_order(numup), upstream_order_work(numup))

    numup = sum(upstream(:, :))
    k = 0
    do i = 1, ny
        do j = 1, nx
            if (upstream(i, j) .eq. 1) then
                k = k + 1
                upstream_i(k) = i
                upstream_j(k) = j
            end if
        end do
    end do

    write (*, *) "Done STEP 5"

    ! STEP 6 : Calc Total Length from The Most Upstream Cell
    allocate (total_length(numup))

    total_length(:) = 0.d0
    do k = 1, numup
        i = upstream_i(k)
        j = upstream_j(k)
        total_length(k) = 0.d0
        do
            ii = downstream_i(i, j)
            jj = downstream_j(i, j)
            total_length(k) = total_length(k) + downstream_distance(i, j)
            if (dir(ii, jj) .eq. 0) exit
            i = ii
            j = jj
        end do
    end do

    write (*, *) "Done STEP 6"

    ! STEP 7 : Decide the Order from the Longest Total Length

    upstream_i_temp(:) = upstream_i(:)
    upstream_j_temp(:) = upstream_j(:)
    call sort_upstream_order(total_length, upstream_order, upstream_order_work, numup)
    do k = 1, numup
        upstream_i(k) = upstream_i_temp(upstream_order(k))
        upstream_j(k) = upstream_j_temp(upstream_order(k))
    end do

    write (*, *) "Done STEP 7"

    ! STEP 8 : Adjust DEM
    adem = dem

    ! adem > 0
    where (adem .gt. -50.d0 .and. adem .le. 0.d0) adem = 0.d0

    ! lifting
    do k = 1, numup
        !write(*,*) "l: ", k, "/", numup
1110    continue
        i = upstream_i(k)
        j = upstream_j(k)
        do
            ii = downstream_i(i, j)
            jj = downstream_j(i, j)
            iii = downstream_i(ii, jj)
            jjj = downstream_j(ii, jj)
            if (dir(ii, jj) .eq. 0) exit
            if ((adem(i, j) - adem(ii, jj)) .gt. lift .and. (adem(iii, jjj) - adem(ii, jj)) .gt. lift) then
                adem(ii, jj) = adem(i, j)
                goto 1110
            end if
            i = ii
            j = jj
        end do
    end do

    ! carving
    do k = 1, numup
        !write(*,*) "c: ", k, "/", numup
        i = upstream_i(k)
        j = upstream_j(k)
        do
            ii = downstream_i(i, j)
            jj = downstream_j(i, j)
            if (adem(i, j) .lt. adem(ii, jj) - carve) then
                adem(ii, jj) = adem(i, j)
            end if
            if (dir(ii, jj) .eq. 0) exit
            i = ii
            j = jj
        end do
    end do

    ! lifting and carving
    do k = 1, numup
        !write(*,*) "cf: ", k, "/", numup
1112    continue
        i = upstream_i(k)
        j = upstream_j(k)
        switch = 0
        do
            ii = downstream_i(i, j)
            jj = downstream_j(i, j)

            if (switch .eq. 0 .and. adem(i, j) .lt. adem(ii, jj)) then
                s1_i = i
                s1_j = j
                s1_down_i = ii
                s1_down_j = jj
                switch = 1
            end if

            if (switch .eq. 1 .and. adem(i, j) .gt. adem(ii, jj)) then
                s2_i = i
                s2_j = j
                switch = 2
            end if

            if (switch .eq. 2) then
                rise_count = ceiling((adem(s1_down_i, s1_down_j) - adem(s1_i, s1_j))/increment)
                fall_count = ceiling((adem(s2_i, s2_j) - adem(ii, jj))/increment)
                adjustment_count = max(1, min(rise_count, fall_count))
                adem(s1_i, s1_j) = adem(s1_i, s1_j) + adjustment_count*increment
                adem(s2_i, s2_j) = adem(s2_i, s2_j) - adjustment_count*increment
                goto 1112
            end if

            if (dir(ii, jj) .eq. 0) exit
            i = ii
            j = jj
        end do
    end do

    write (*, *) "Done STEP 8"

    !! Sample Output
    !open(1000, file = "longest_line.txt")
    !i = upstream_i(1)
    !j = upstream_j(1)
    !do
    ! call down(dir, nx, ny, i, j, length, ii, jj, dis)
    ! write(1000, *) dem(i, j), adem(i, j)
    ! if( dir(ii, jj) .eq. 0 ) exit
    ! i = ii
    ! j = jj
    !enddo
    !close(1000)

    ! STEP 9 : Make River width, depth, leavy height
    ! river widhth, depth, leavy height, river length, river area ratio
    allocate (riv(ny, nx), width(ny, nx), depth(ny, nx), height(ny, nx))
    allocate (landuse(ny, nx), gsd_s(ny, nx), gsd_c(ny, nx))
    width = 0.d0
    depth = 0.d0
    height = 0.d0
    landuse = 1
    gsd_s = 1
    gsd_c = 1

    call cg_iric_read_integer(cgns_f, "riv_thresh", riv_thresh, ier)
    if (ier /= 0) error stop "Failed to read riv_thresh"
    call cg_iric_read_real(cgns_f, "width_param_c", width_param_c, ier)
    if (ier /= 0) error stop "Failed to read width_param_c"
    call cg_iric_read_real(cgns_f, "width_param_s", width_param_s, ier)
    if (ier /= 0) error stop "Failed to read width_param_s"
    call cg_iric_read_real(cgns_f, "depth_param_c", depth_param_c, ier)
    if (ier /= 0) error stop "Failed to read depth_param_c"
    call cg_iric_read_real(cgns_f, "depth_param_s", depth_param_s, ier)
    if (ier /= 0) error stop "Failed to read depth_param_s"
    call cg_iric_read_real(cgns_f, "height_param", height_param, ier)
    if (ier /= 0) error stop "Failed to read height_param"
    call cg_iric_read_integer(cgns_f, "height_limit_param", height_limit_param, ier)
    if (ier /= 0) error stop "Failed to read height_limit_param"

    if (riv_thresh < 0) error stop "riv_thresh must be zero or greater"
    if (width_param_c < 0.d0 .or. depth_param_c < 0.d0) &
        error stop "River width and depth coefficients must be zero or greater"

    riv = 0 ! slope cell
    if (riv_thresh .gt. 0) then
        where (acc .gt. riv_thresh) riv = 1 ! river cell
    end if

    where (riv .eq. 1) width = width_param_c*(acc*dx*dy*1.d-6)**width_param_s
    where (riv .eq. 1) depth = depth_param_c*(acc*dx*dy*1.d-6)**depth_param_s
    where (riv .eq. 1 .and. acc .gt. height_limit_param) height = height_param

    write (*, *) "Done STEP 9"

    ! STEP 10 : scaleUP

    write (*, *) "Done STEP 10"

    ! Recreate the structured grid before writing cell attributes.
    call cg_iRIC_Write_Grid2d_Coords(cgns_f, nx + 1, ny + 1, gxx, gyy, ier)
    if (ier /= 0) error stop "Failed to create the calculation grid"
    call iric_write_cell_real("elevation_c", nx, ny, adem)
    call iric_write_cell_integer("dir_c", nx, ny, dir)
    call iric_write_cell_integer("acc_c", nx, ny, acc)
    call iric_write_cell_real("width_c", nx, ny, width)
    call iric_write_cell_real("depth_c", nx, ny, depth)
    call iric_write_cell_real("height_c", nx, ny, height)
    call iric_write_cell_integer("landuse_c", nx, ny, landuse)
    call iric_write_cell_integer("GSD_s", nx, ny, gsd_s)
    call iric_write_cell_integer("GSD_c", nx, ny, gsd_c)
    call iric_cgns_close()

end program RRI_GridGenerator

subroutine sort_upstream_order(total_length, order, work, count)
    implicit none

    integer, intent(in) :: count
    real(8), intent(in) :: total_length(count)
    integer, intent(out) :: order(count)
    integer, intent(inout) :: work(count)
    integer :: width, left, middle, right, left_pos, right_pos, output_pos
    integer :: left_index, right_index, k
    logical :: take_left

    do k = 1, count
        order(k) = k
    end do

    width = 1
    do while (width < count)
        left = 1
        do while (left <= count)
            middle = min(left + width, count + 1)
            right = min(left + 2*width, count + 1)
            left_pos = left
            right_pos = middle
            output_pos = left

            do while (left_pos < middle .and. right_pos < right)
                left_index = order(left_pos)
                right_index = order(right_pos)
                take_left = total_length(left_index) > total_length(right_index)
                if (total_length(left_index) == total_length(right_index)) &
                    take_left = left_index > right_index

                if (take_left) then
                    work(output_pos) = left_index
                    left_pos = left_pos + 1
                else
                    work(output_pos) = right_index
                    right_pos = right_pos + 1
                end if
                output_pos = output_pos + 1
            end do

            do while (left_pos < middle)
                work(output_pos) = order(left_pos)
                left_pos = left_pos + 1
                output_pos = output_pos + 1
            end do
            do while (right_pos < right)
                work(output_pos) = order(right_pos)
                right_pos = right_pos + 1
                output_pos = output_pos + 1
            end do
            left = left + 2*width
        end do
        order(:) = work(:)
        width = width*2
    end do
end subroutine sort_upstream_order

subroutine check_grid_header(label, nx, ny, xll, yll, cellsize, &
                             nx_file, ny_file, xll_file, yll_file, cellsize_file)
    implicit none
    character(len=*), intent(in) :: label
    integer, intent(in) :: nx, ny, nx_file, ny_file
    real(8), intent(in) :: xll, yll, cellsize, xll_file, yll_file, cellsize_file
    real(8) :: coordinate_tolerance, cellsize_tolerance

    ! ASCII raster headers can differ slightly after GIS export and decimal conversion.
    coordinate_tolerance = max(1.d-8, abs(cellsize)*1.d-6)
    cellsize_tolerance = max(1.d-12, abs(cellsize)*1.d-6)
    if (nx_file /= nx .or. ny_file /= ny) then
        write(*, *) trim(label), " dimensions do not match DEM:", nx_file, ny_file, nx, ny
        error stop "Input raster dimensions do not match"
    end if
    if (abs(xll_file - xll) > coordinate_tolerance .or. &
        abs(yll_file - yll) > coordinate_tolerance .or. &
        abs(cellsize_file - cellsize) > cellsize_tolerance) then
        write(*, *) trim(label), " grid geometry does not match DEM"
        error stop "Input raster geometry does not match"
    end if
end subroutine check_grid_header

! Search Downstream
subroutine down(dir, nx, ny, i, j, length, ii, jj, dis)

    integer nx, ny, i, j, ii, jj
    real(8) length, dis
    integer dir(ny, nx)

    ! right
    if (dir(i, j) .eq. 1) then
        ii = i
        jj = j + 1
        dis = length
        ! right down
    elseif (dir(i, j) .eq. 2) then
        ii = i + 1
        jj = j + 1
        dis = length*sqrt(2.d0)
        ! down
    elseif (dir(i, j) .eq. 4) then
        ii = i + 1
        jj = j
        dis = length
        ! left down
    elseif (dir(i, j) .eq. 8) then
        ii = i + 1
        jj = j - 1
        dis = length*sqrt(2.d0)
        ! left
    elseif (dir(i, j) .eq. 16) then
        ii = i
        jj = j - 1
        dis = length
        ! left up
    elseif (dir(i, j) .eq. 32) then
        ii = i - 1
        jj = j - 1
        dis = length*sqrt(2.d0)
        ! up
    elseif (dir(i, j) .eq. 64) then
        ii = i - 1
        jj = j
        dis = length
        ! right up
    elseif (dir(i, j) .eq. 128) then
        ii = i - 1
        jj = j + 1
        dis = length*sqrt(2.d0)
        ! zero
    else
        ii = i
        jj = j
        dis = length
    end if

end subroutine

! Hubeny_sub.f90
subroutine hubeny_sub(x1_deg, y1_deg, x2_deg, y2_deg, d)
    implicit none

    real(8) x1_deg, y1_deg, x2_deg, y2_deg
    real(8) x1, y1, x2, y2
    real(8) pi, dx, dy, mu, a, b, e, W, N, M, d

    !read(*,*) x1_deg, y1_deg, x2_deg, y2_deg ! in degree

    pi = 3.1415926535897d0

    x1 = x1_deg*pi/180.d0
    y1 = y1_deg*pi/180.d0
    x2 = x2_deg*pi/180.d0
    y2 = y2_deg*pi/180.d0

    dy = y1 - y2
    dx = x1 - x2
    mu = (y1 + y2)/2.

    a = 6378137.000d0 ! Semi-Major Axis
    b = 6356752.314d0 ! Semi-Minor Axis

    e = sqrt((a**2.d0 - b**2.d0)/(a**2.d0))

    W = sqrt(1.-e**2.d0*(sin(mu))**2.d0)

    N = a/W

    M = a*(1.-e**2.d0)/W**(3.d0)

    d = sqrt((dy*M)**2.d0 + (dx*N*cos(mu))**2.d0)

    !write(*,*) d ! in m

end subroutine

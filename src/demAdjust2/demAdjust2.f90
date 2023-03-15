! demAdjust2.f90
!
! coded by T.Sayama on April 28, 2010
! v2.1 Jan 10, 2013

! DEM Adjustment Program
!
program demAdjust2

    use RRI_iric
    use iric
    implicit none

    character*256 infile, infile_dem, infile_dir, infile_acc
    character*256 outfile_adem, outfile_adir, outfile_slo

    parameter(infile="demAdjust2.txt")

    integer acc_thresh
    real(8) increment, carve, lift
    parameter(acc_thresh=0)    ! adjustment is done for cells with acc > acc_thresh
    parameter(lift=500.d0)     ! [m]
    parameter(carve=5.d0)      ! [m]
    parameter(increment=0.01d0) ! [m]

    integer utm
    parameter(utm=0) ! 1: UTM, 0 : LatLon

    integer nx, ny
    real(8) xllcorner, yllcorner, cellsize, nodata
    real(8) dem_lowest

    real(8) d1, d2, d3, d4
    real(8) x1, x2, y1, y2, dx, dy
    real(8) length, dis, dif, slope

    integer i, j, k, l, ii, jj, iii, jjj, ios

    real(8), dimension(:, :), allocatable :: dem, adem, slo
    integer, dimension(:, :), allocatable :: dir, acc, adir, riv
    real(8), dimension(:, :), allocatable :: width, depth, height
    integer, dimension(:, :), allocatable :: upstream
    real(8), dimension(:), allocatable :: total_length
    integer, dimension(:), allocatable :: upstream_i, upstream_j, upstream_i_temp, upstream_j_temp
    integer switch, s1_i, s1_j, s2_i, s2_j, numup, kk, longest_kk
    real(8) longest_length

    integer riv_thresh
    real(8) width_param_c, width_param_s
    real(8) depth_param_c, depth_param_s
    integer height_limit_param
    real(8) height_param

    character*256 ctemp
    character*20 ctemp2

    integer :: ier, icount
    real(8), allocatable, save :: gxx(:, :), gyy(:, :)
    !--------------------------------------------------

    ! Start opne CGNS and Read file name
    call iric_cgns_open()
    call cg_iric_read_string(cgns_f, "demfile", infile_dem, ier)
    call cg_iric_read_string(cgns_f, "accfile", infile_acc, ier)
    call cg_iric_read_string(cgns_f, "dirfile", infile_dir, ier)
    ! End opne CGNS and Read file name

    !Output file name
    outfile_adem = "./topo/adem.txt"
    outfile_adir = "./topo/adir.txt"

    !! STEP 0 : Open Files
    !open(1, file = infile, status = "old")
    !read(1, '(a)') infile_dem
    !read(1, '(a)') infile_dir
    !read(1, '(a)') infile_acc
    !read(1, '(a)') outfile_adem
    !read(1, '(a)') outfile_adir
    !!read(1, '(a)') outfile_slo

    open (10, file=trim(infile_dem), status="old")
    open (20, file=trim(infile_dir), status="old")
    open (30, file=trim(infile_acc), status="old")
    open (40, file=trim(outfile_adem))
    open (50, file=trim(outfile_adir))
    !open(60, file = outfile_slo)

    ! STEP 1 : Reading File
    read (10, *) ctemp, nx
    read (10, *) ctemp, ny
    read (10, *) ctemp, xllcorner
    read (10, *) ctemp, yllcorner
    read (10, *) ctemp, cellsize
    read (10, *) ctemp, nodata

    allocate (dem(ny, nx), dir(ny, nx), acc(ny, nx), adem(ny, nx), adir(ny, nx), slo(ny, nx), upstream(ny, nx))

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
    write (40, '(a30)') ctemp
    write (50, '(a30)') ctemp
    !write(60, '(a30)') ctemp

    read (10, '(a30)') ctemp
    write (40, '(a30)') ctemp
    write (50, '(a30)') ctemp
    !write(60, '(a30)') ctemp

    read (10, '(a30)') ctemp
    write (40, '(a30)') ctemp
    write (50, '(a30)') ctemp
    !write(60, '(a30)') ctemp

    read (10, '(a30)') ctemp
    write (40, '(a30)') ctemp
    write (50, '(a30)') ctemp
    !write(60, '(a30)') ctemp

    read (10, '(a30)') ctemp
    write (40, '(a30)') ctemp
    write (50, '(a30)') ctemp
    !write(60, '(a30)') ctemp

    read (10, '(a30)') ctemp
    write (40, '(a30)') ctemp
    write (50, '(a30)') ctemp
    !write(60, '(a30)') ctemp

    do i = 1, ny
        read (10, *) (dem(i, j), j=1, nx)
    end do

    read (20, *) ctemp, nx
    read (20, *) ctemp, ny
    read (20, *) ctemp, xllcorner
    read (20, *) ctemp, yllcorner
    read (20, *) ctemp, cellsize
    read (20, *) ctemp, nodata

    do i = 1, ny
        read (20, *) (dir(i, j), j=1, nx)
    end do

    read (30, *) ctemp, nx
    read (30, *) ctemp, ny
    read (30, *) ctemp, xllcorner
    read (30, *) ctemp, yllcorner
    read (30, *) ctemp, cellsize
    read (30, *) ctemp, nodata

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
    do i = 1, ny
        write (50, '(<nx>i6)') (dir(i, j), j=1, nx)
    end do

    write (*, *) "Done STEP 3"

    ! STEP 4 : Slope Calculation (Steepest Slope Among the Eight Directions)
    ! (Not necessary for demAdjust2 algorithm)
    slo(:, :) = nodata
    do i = 1, ny
        do j = 1, nx
            if (dem(i, j) .lt. -100.d0) cycle ! outside the domain
            if (dir(i, j) .eq. 0) cycle ! outlet point
            slope = 0.d0
            do ii = i - 1, i + 1
                do jj = j - 1, j + 1
                    if (ii .eq. 0 .or. jj .eq. 0) cycle
                    if (ii .gt. ny .or. jj .gt. nx) cycle
                    if (dem(ii, jj) .lt. -100.) cycle
                    dis = length*sqrt(2.d0)
                    if (ii .eq. 0 .or. jj .eq. 0) dis = length
                    slope = max(slope, abs(dem(i, j) - dem(ii, jj))/dis)
                end do
            end do
            slo(i, j) = slope
        end do
    end do
    !do i = 1, ny
    ! write(60, '(<nx>f10.3)') (slo(i, j), j = 1, nx)
    !enddo

    write (*, *) "Done STEP 4"

    ! STEP 5 : Find The Most Upstream Cell
    upstream(:, :) = 1
    do i = 1, ny
        do j = 1, nx
            if (dem(i, j) .lt. -100.d0 .or. dir(i, j) .eq. 0 .or. acc(i, j) .lt. acc_thresh) then
                upstream(i, j) = 0
                cycle
            end if
            call down(dir, nx, ny, i, j, length, ii, jj, dis)
            upstream(ii, jj) = 0
        end do
    end do
    numup = sum(upstream(:, :))

    allocate (upstream_i(numup), upstream_j(numup))
    allocate (upstream_i_temp(numup), upstream_j_temp(numup))

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
            call down(dir, nx, ny, i, j, length, ii, jj, dis)
            total_length(k) = total_length(k) + dis
            if (dir(ii, jj) .eq. 0) exit
            i = ii
            j = jj
        end do
    end do

    write (*, *) "Done STEP 6"

    ! STEP 7 : Decide the Order from the Longest Total Length

    upstream_i_temp(:) = upstream_i(:)
    upstream_j_temp(:) = upstream_j(:)
    do k = 1, numup
        !write(*,*) "S7", k, "/", numup
        longest_length = 0.d0
        do kk = 1, numup
            if (total_length(kk) .ge. longest_length) then
                longest_length = total_length(kk)
                longest_kk = kk
            end if
        end do
        upstream_i(k) = upstream_i_temp(longest_kk)
        upstream_j(k) = upstream_j_temp(longest_kk)
        total_length(longest_kk) = -1*total_length(longest_kk)
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
            call down(dir, nx, ny, i, j, length, ii, jj, dis)
            call down(dir, nx, ny, ii, jj, length, iii, jjj, dis)
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
1111    continue
        i = upstream_i(k)
        j = upstream_j(k)
        do
            call down(dir, nx, ny, i, j, length, ii, jj, dis)
            if (adem(i, j) .lt. adem(ii, jj) - carve) then
                adem(ii, jj) = adem(i, j)
                goto 1111
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
            call down(dir, nx, ny, i, j, length, ii, jj, dis)

            if (switch .eq. 0 .and. adem(i, j) .lt. adem(ii, jj)) then
                s1_i = i
                s1_j = j
                switch = 1
            end if

            if (switch .eq. 1 .and. adem(i, j) .gt. adem(ii, jj)) then
                s2_i = i
                s2_j = j
                switch = 2
            end if

            if (switch .eq. 2) then
                adem(s1_i, s1_j) = adem(s1_i, s1_j) + increment
                adem(s2_i, s2_j) = adem(s2_i, s2_j) - increment
                goto 1112
            end if

            if (dir(ii, jj) .eq. 0) exit
            i = ii
            j = jj
        end do
    end do

    do i = 1, ny
        write (40, '(10000f13.5)') (adem(i, j), j=1, nx)
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
    width = 0.d0
    depth = 0.d0
    height = 0.d0

    call cg_iric_read_integer(cgns_f, "riv_thresh", riv_thresh, ier)
    call cg_iric_read_real(cgns_f, "width_param_c", width_param_c, ier)
    call cg_iric_read_real(cgns_f, "width_param_s", width_param_s, ier)
    call cg_iric_read_real(cgns_f, "depth_param_c", depth_param_c, ier)
    call cg_iric_read_real(cgns_f, "depth_param_s", depth_param_s, ier)
    call cg_iric_read_real(cgns_f, "height_param", height_param, ier)
    call cg_iric_read_integer(cgns_f, "height_limit_param", height_limit_param, ier)

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

    !CGNSファイルに出力
    call cg_iRIC_Write_Grid2d_Coords(cgns_f, nx + 1, ny + 1, gxx, gyy, ier)
    call iric_write_cell_real("elevation_c", nx, ny, adem)
    call iric_write_cell_integer("dir_c", nx, ny, dir)
    call iric_write_cell_integer("acc_c", nx, ny, acc)
    call iric_write_cell_real("width_c", nx, ny, width)
    call iric_write_cell_real("depth_c", nx, ny, depth)
    call iric_write_cell_real("height_c", nx, ny, height)
    call iric_cgns_close()

end program demAdjust2

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
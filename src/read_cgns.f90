!****************************************************************************************************
!
!   Subroutine : Read Calculation Conditions
!
!****************************************************************************************************
subroutine read_cgns(fid)
    use iric
    use global_variables
    
    integer, intent(in) :: fid
    integer :: ier
    
    integer :: n, count
    integer :: ind_size, size
    integer, dimension(:), allocatable :: indices
    integer :: ii, jj
    
    integer :: esize
    integer, dimension(:,:), allocatable :: edges
    
    integer :: nsize
    real(8), dimension(:), allocatable :: xtmp, ytmp
    
    ! variables for time condition
    call cg_iric_read_real(fid, "dt", dt, ier)
    call cg_iric_read_real(fid, "stime", stime, ier)
    call cg_iric_read_real(fid, "etime", etime, ier)
    call cg_iric_read_real(fid, "toutput", toutput, ier)  
    
    istime = int(stime/dt)
    ietime = int(etime/dt)
    ioutput = int(toutput/dt)
    
    ! grid conditions
    call cg_iric_read_grid2d_str_size(fid, isize, jsize, ier)
    if(ier /= 0)then
        write(*,*) "Error: no grid data."
        stop
    end if
    
    call alloc_variables(isize, jsize)
    call cg_iric_read_grid2d_coords(fid, xx, yy, ier)
    call cg_iric_read_grid_real_cell(fid, "elevation", zc, ier)
    call cg_iric_read_grid_real_cell(fid, "roughness", sn, ier)
    call cg_iric_read_grid_integer_cell(fid, "obstacle", iobst, ier)
    
    ! boundary conditions
    call cg_iric_read_bc_count(fid, "inflow", count)
    do n = 1, count
        
        call cg_iric_read_bc_indicessize(fid, "inflow", n, size, ier)
        
        !indices = point or cell
        ! i,j for a point or a cell
        allocate(indices(size*2))
        call cg_iric_read_bc_indices(fid, "inflow", n, indices, ier)
        
        
        ! position = edge
        esize = size / 2
        allocate(edges(2, esize))
        do m = 1, esize
            ii = 4 * (m - 1) + 1 ; jj = ii + 1
            edges(1, m) = indices(ii)    !i-index
            edges(2, m) = indices(jj)    !j-index
            
            write(*,*) edges(1, m), edges(2, m)
        end do
        
        
        
        
        !read discharge function the bc id =n
        call cg_iric_read_bc_functionalsize(fid, "inflow", n, "discharge", nsize, ier)    
        allocate(xtmp(1:nsize), ytmp(1:nsize))
        call cg_iric_read_bc_functional(fid, "inflow", n, "discharge", xtmp, ytmp, ier)
    end do

end subroutine
!****************************************************************************************************
!
!   Subroutine : Write Calculation Results
!
!****************************************************************************************************
subroutine write_cgns(fid)
    use iric
    use global_variables
    integer, intent(in) :: fid
    integer :: ier
    
    ! write time
    call cg_iric_write_sol_time(fid, ctime, ier)
    
    ! write grid
    
    ! write result such as discharge and waterlevel
    !call cg_iric_write_sol_baseiterative_real(fid, "QUP[m3_s]", qup, ier)
    !call cg_iric_write_sol_baseiterative_real(fid, "HDW[m]", hdw, ier)
    
    ! write node result such as elevation, depth, waterlevel
    call cg_iric_write_sol_cell_real(fid, "Elevation[m]", zz(1:isize-1, 1:jsize-1), ier)
    
    ! write cell result such as elevation, depth, waterlevel
    call cg_iric_write_sol_cell_real(fid, "Elevation_C[m]", zc(1:isize-1, 1:jsize-1), ier)
    call cg_iric_write_sol_cell_real(fid, "Depth[m]", depth(1:isize-1, 1:jsize-1), ier)
    call cg_iric_write_sol_cell_real(fid, "WaterLevel[m]", wl(1:isize-1, 1:jsize-1), ier)
    
    ! write result such as velocity iface
    !call cg_iric_write_sol_iface_real(fid, "Vx", vx(1:isize, 1:jsize-1), ier)
    
    ! write result such as velocity jface
    !call cg_iric_write_sol_jface_real(fid, "Vy", vy(1:isize-1, 1:jsize), ier)
    
end subroutine
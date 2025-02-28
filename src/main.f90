!****************************************************************************************************
!
!   Sample program for the iRIC version 4
!
!****************************************************************************************************
    program samaple
    use iric
    use global_variables

    implicit none
    integer, parameter :: strmax = 1024
    
    integer :: fid, ier
    character(len=strmax) :: cgnsfile
    
    integer :: icount
    
!======================================================================
    write(*,'(a)') "---------- START Sample Solver ---------"
    
    ! get argument
    if(iargc() == 1) then
        call getarg(1, cgnsfile)
        
    else
        write(*,"(a)") "Error: no argument."
        write(*,"(a)") "usage: "
        write(*,"(a)") "***.exe Case1.cgn"
        stop
    end if
    
    ! open cgns file
    call cg_iric_open(cgnsfile, IRIC_MODE_MODIFY, fid, ier)
    call cg_iric_clear_sol(fid, ier)
    call iric_initoption(IRIC_OPTION_CANCEL, ier)
    call iric_initoption(IRIC_OPTION_DIVIDESOLUTIONS, ier)
    
    
    ! set default
    call set_default()
    
    ! read calculation conditions
    call read_cgns(fid)
    
    ! set initial
    call set_initial()
    
    ! start time loop
    do icount = istime, ietime
        
        ! set time
        ctime = icount * dt
        
        ! the update button is clicked or not
        call cg_iric_check_update(fid, ier)
        
        ! the cancle button is clicked or not
        call iric_check_cancel(ier)
        if (ier == 1) then
            write(*,*) ""
            write(*,*) "solver stops because the cancel button was clicked."
            exit
        end if
        
        ! output result
        if(mod(icount, ioutput) == 0) then
            write(*,"(a7, f20.3, a)") "time = ", ctime, "[sec]"
            call cg_iric_write_sol_start(fid, ier)
            call write_cgns(fid)
            call cg_iric_write_sol_end(fid, ier)
        end if
        
        ! calculation main code for one step        
        call calc_main()
        
    end do
    
    ! close cgns file
    call cg_iric_close(fid, ier)
    write(*,*) ""
    write(*,'(a)') "---------- END: Sample Solver ---------"
    end program


program check_cgns

    use iric

    implicit none

    character(len=1024) :: fname
    integer :: ier, icount, idf
    integer :: run_type
    
    ! get argument
    icount = iargc()
    if (icount == 1) then
        call getarg(1, fname)
    else
        write(*,"(a)") "You should specify an argument."
        stop
    endif

    ! file open
    call cg_iric_open(fname, IRIC_MODE_READ, idf, ier)
    if (ier /= 0) stop "cg_iric_open failed"
    
    ! get run type
    call cg_iric_read_integer(idf, "run_type", run_type, ier)
    if (ier /= 0) stop "Failed to read run_type"
    
    open(1, file="runtype.txt", status="replace", action="write")
    write(1, "(i0)") run_type
    close(1)
    
    call cg_iric_close(idf, ier)    
    stop

end program
    
    
    
    

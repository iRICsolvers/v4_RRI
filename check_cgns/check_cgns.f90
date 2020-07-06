
    program check_cgns
    implicit none
    include 'cgnslib_f.h'
    include 'iriclib_f.h'

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
    call cg_open_f(fname, CG_MODE_READ, idf, ier)
    if (ier /= 0) stop "cg_open_f failed"
    call cg_iric_init_f(idf, ier)

    
    ! get run type
    call cg_iric_read_integer_f("run_type", run_type, ier)
    
    open(1,file ="runtype.txt")
    write(1,*) run_type
    close(1)
    
    call cg_close_f(idf, ier)    
    stop
    end program
    
    
    
    
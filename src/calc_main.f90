!****************************************************************************************************
!
!   Subroutine : Calculation Main
!
!****************************************************************************************************
subroutine calc_main()
    use iric
    use global_variables
    
    real(8) :: tmp
    
    
    ! update 
    do j = 1, jsize-1
        do i = 1, isize-1
            
            call RANDOM_NUMBER(tmp)
            zc(i,j) = zc(i,j) + 0.01 * tmp
            
            call RANDOM_NUMBER(tmp)
            wl(i,j) = wl(i,j) + 0.01 * tmp
            depth(i,j) = wl(i,j) - zc(i,j)
            
        end do
    end do
    
    
end subroutine
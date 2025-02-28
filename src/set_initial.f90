!****************************************************************************************************
!
!   Subroutine : Set Initial Conditions
!
!****************************************************************************************************
subroutine set_initial()
    use iric
    use global_variables
    
    real(8) :: zsum
    
    ! cross-sectional averaged elevation and set water level and depth
    do i = 1, isize-1
        zsum = 0.0d0
        do j = 1, jsize-1
            zsum = zsum + zc(i,j)
        end do
        zave(i) = dble(zsum / (jsize-1))
        
        !write(*,*) i, zave(i)
        do j = 1, jsize-1
            wl(i,j) = zave(i) + 2.0d0
            depth(i,j) = wl(i,j) - zc(i,j)
            if(depth(i,j) < 0.0d0)then
                depth(i,j) = 0.0d0    
            end if
        end do
    end do    
    
end subroutine
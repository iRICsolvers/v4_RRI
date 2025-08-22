! RRI_Dam.f90  !this dam file is totally modified for RSR model

! reading dam control file
subroutine dam_read
    use globals
    use dam_mod
    use RRI_iric
    use iric

    implicit none

    integer :: i
    integer :: size, ier
    integer, dimension(2) :: tmp_idx
    character*256 :: outfile

    !allocate (damflg(riv_count), dam_qin(riv_count))
    !damflg(:) = 0

    !if (dam_switch .eq. 1) then

        !open(99, file=Damfile, status="old" )
        !read(99,*) dam_num
        allocate (dam_name(dam_num), dam_kind(dam_num) &
                  , dam_ix(dam_num), dam_iy(dam_num) &
                  , dam_vol(dam_num), dam_vol_temp(dam_num) &
                  , dam_volmax(dam_num), dam_state(dam_num) &
                  , dam_qout(dam_num), dam_loc(dam_num) &
                  , dam_floodq(dam_num))
        allocate (dam_w_vol(dam_num), dam_qin(dam_num), dam_reserv_area(dam_num),reserv_elev(dam_num),dam_discharge_switch(dam_num),outswitch_dam(dam_num),outfile_dam(dam_num)) !RSR model 20240724
        dam_vol(:) = 0.d0
        dam_state(:) = 0
        dam_qout(:) = 0.d0  !RSR model 20240724
        dam_qin(:) = 0.d0   !RSR model 20240724
       
        reserv_elev(:) = -9999 !for dam
        dam_discharge_switch(:) =0   !for dam
        !!need to add function to read dam_reserv_area and reserv_elev  !RSR model 20240724

        do i = 1, dam_num
            !read(99,*) dam_name(i), dam_iy(i), dam_ix(i), dam_volmax(i), dam_floodq(i)

            call cg_iric_read_bc_indicessize(cgns_f, "dam", i, size, ier)
            if (size /= 1) then
                write (*, *) "Error: A boundary condition for dam can be specified on only one cell for one."
                call iric_cgns_close()
                stop
            end if
            call cg_iric_read_bc_indices(cgns_f, "dam", i, tmp_idx, ier)

            dam_iy(i) = ny - tmp_idx(2) + 1
            dam_ix(i) = tmp_idx(1)
            write(*,'(a,i)') "for check03, i =", i
            call cg_iric_read_bc_string(cgns_f,"dam",i,"name", dam_name(i), ier)
            call cg_iric_read_bc_real(cgns_f, "dam", i, "dam_volmax", dam_volmax(i), ier)
            call cg_iric_read_bc_real(cgns_f, "dam", i, "dam_floodq", dam_floodq(i), ier)
            call cg_iric_read_bc_real(cgns_F,"dam",i, "dam_reserv_area",dam_reserv_area(i), ier)
            call cg_iric_read_bc_real(cgns_f,"dam",i, "reserv_elev",reserv_elev(i), ier)
            if(reserv_elev(i)>-9999)then
            write(*,*)trim(dam_name(i))," elevation has been modified as ", reserv_elev(i)
            endif
            call cg_iric_read_bc_integer(cgns_f, "dam",i, "dam_discharge_switch",dam_discharge_switch(i),ier)
            outfile_dam(i) = ''
            call cg_iric_read_bc_integer(cgns_f,"dam", i, "outswitch_dam", outswitch_dam(i), ier)
            if (outswitch_dam(i) == 1)then
             call cg_iric_read_bc_string(cgns_f, "dam",i,"outfile_dam", outfile, ier)
               outfile_dam(i) =  trim(outfile)// trim(dam_name(i)) // ".txt"
            endif

        end do
        !close(99)
   ! end if
end subroutine dam_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----read dam discharge operation file; added by Qin 2021/5/24
!for dam
subroutine dam_discharge_operation(i)
 use globals
 use dam_mod
 use sediment_mod!, only: xtmp, ytmp !for dam
 use RRI_iric
 use iric
 implicit none
    
 integer :: i,t,t_count,tt,ier
 real(8) tdiff
 !do i= 1,dam_num
    call cg_iric_read_bc_functionalsize(cgns_f, 'dam', i, 'dam_outflow', t_count, ier) 
    if(t_count==0) then
        write(*,*) "Please set the hourly dam ourflow discharge data"
        stop
    endif
    write(*,'(a,2i)') "for check04, i =", i, t_count
    allocate(xtmp(t_count),ytmp(t_count))
    call cg_iric_read_bc_functionalwithname(cgns_f,'dam', i, 'dam_outflow','dam_discharge_period', xtmp,ier)
    max_dam_discharge_t = xtmp(t_count)
    if(i==1) allocate(  dam_discharge(max_dam_discharge_t,dam_num), dam_discharge_period(max_dam_discharge_t))

    if(xtmp(1)<1)then
        write(*,*)"The given outflow start time of ", trim(dam_name(i))," should be >= 1 hour"
        stop
    endif    
    if(xtmp(t_count)/=lasth) then
        write(*,*)"The outflow duration of ",trim(dam_name(i))," should be same with the end of computation duration"
        stop
    endif

    call cg_iric_read_bc_functionalwithname(cgns_f,'dam', i,'dam_outflow','dam_discharge', ytmp,ier)

    do t = 1, t_count

     if(t==1)then
        tt=t
        if(xtmp(t)>1)then
         do
          dam_discharge(tt,i) = ytmp(t)
          if(i==1) dam_discharge_period(tt) = tt*3600.
          tt= tt+1
          if(tt > xtmp(t)) exit
         enddo
        else
          dam_discharge(tt,i) = ytmp(t)
          if(i==1)  dam_discharge_period(tt) = xtmp(t)*3600.
        endif
     else
        tdiff = xtmp(t)-xtmp(t-1)
        tt= xtmp(t-1)+1
        if(tdiff>1.)then
         do
          dam_discharge(tt,i) = ytmp(t)
          if(i==1) dam_discharge_period(tt) = tt*3600.
          tt= tt+1
          if(tt > xtmp(t)) exit
         enddo
        else
          dam_discharge(tt,i) = ytmp(t)
          if(i==1) dam_discharge_period(tt) = xtmp(t)*3600.
        endif          

     endif

    enddo
 !enddo    
    DEALLOCATE(xtmp, STAT = ier)
    DEALLOCATE(ytmp, STAT = ier)

 !    write(*,'(a,i)') "End time of dam outflow discharge(hour) : ",max_dam_discharge_t
end subroutine dam_discharge_operation   

subroutine dam_prepare(qr_idx)
    
use globals
use dam_mod!, only :dam_qin
implicit none
        
integer :: k, kk
real(8) :: qr_idx(riv_count)
        
dam_qin(:) = 0.d0
do k = 1, riv_count
    kk = down_riv_idx(k)
    if(damflg(kk) == 0) cycle
    dam_qin(damflg(kk)) = dam_qin(damflg(kk)) + qr_idx(k)
enddo

    
end subroutine dam_prepare
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!for dam
subroutine dam_operation(i)

!use globals, only :ddt,area
use globals!, only :time, ddt ! v1.4
use dam_mod
use sediment_mod, only: 
!use sediment_mod!, only:sed_switch,dam_sedi_totalV
implicit none

integer  t
real(8)  qdiff
integer i    
dam_qout(i) = 0.d0
if (dam_discharge_switch(i) == 1) then
 
   t = int((time+ddt)/3600)+1
   !----
   if((time+ddt)>dam_discharge_period(t))then
    write(*,*) "error: time of dam outflow < present comupation time"
    write(*,'("Time of outflow (s): ",f15.5)') dam_discharge_period(t)
    write(*,'("Present time (s): ",f15.5)') time+ddt
   stop
   elseif ((time+ddt) .lt. dam_discharge_period(t-1))then
    write(*,*) "error: present comupation time <  previous time section of dam outflow data"
    write(*,'("Time of outflow (s): ",f15.5)') dam_discharge_period(t)
    write(*,'("Time of previous outflow  time section (s): ",f15.5)')dam_discharge_period(t-1)
    write(*,'("Present time (s): ",f15.5)') time+ddt
   endif
      if(t.ge.max_dam_discharge_t) t = max_dam_discharge_t
  if(dam_state(i)==1 )then
   !if( dam_discharge(ttemp,i) > dam_qin(i)) then
   if( dam_discharge(t,i) .ge. 0.d0) then ! Attention: follow the observed dam discharge even dam is full 
      dam_qout(i) = dam_discharge(t,i)
   elseif(dam_discharge(t,i)<-88.d0 .and. dam_floodq(i)>dam_qin(i)) then
      dam_qout(i) =  dam_floodq(i)
   else
      dam_qout(i) = dam_qin(i)
   endif   
  else ! dam_state(i) == 0
      if(dam_discharge(t,i) .ge. 0.d0)then    
         dam_qout(i) = dam_discharge(t,i)          
      elseif(dam_discharge(t, i) == -88.d0)then
            !dam_w_vol(i) = 0.d0
            dam_qout(i) = dam_qin(i)  
      else  
         if ( dam_qin(i) .lt. dam_floodq(i) ) then ! v1.4
            dam_qout(i) = dam_qin(i)   
         else        
            dam_qout(i) = dam_floodq(i) ! v1.4
         endif
      endif
  endif
else
  if(dam_state(i)==1)then
      if(dam_qin(i) .lt. dam_floodq(i))then
         dam_qout(i) = dam_floodq(i)
      else
         dam_qout(i) = dam_qin(i)
      endif
  else 
      if ( dam_qin(i) .lt. dam_floodq(i) ) then ! v1.4
         dam_qout(i) = dam_qin(i)
      else        
         dam_qout(i) = dam_floodq(i) ! v1.4
      endif
  endif
endif

   qdiff = dam_qin(i) - dam_qout(i) ! v1.4
!if(qdiff.lt.0.d0) qdiff = 0.d0
  dam_vol_temp(i) = dam_vol_temp(i) + qdiff * ddt
!write(*,*) dam_qin(i), dam_qout(i), dam_discharge(ttemp,i)

end subroutine dam_operation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dam_checkstate(t)

use globals!,only :riv_count,area
use dam_mod
use sediment_mod!, only:sed_switch,dam_sedi_totalV
implicit none
    
integer :: i,t
    
do i=1, dam_num
 dam_vol_temp(i) = dam_vol_temp(i) / 6.d0
 dam_w_vol(i) = dam_w_vol(i) + dam_vol_temp(i) 
 if(dam_w_vol(i).lt.0.) dam_w_vol(i) = 0.d0 !modified 2021/6/9
 if(sed_switch.ge.1)then !modified by Qin 
 if(t.gt.10)then
 dam_vol(i) = dam_w_vol(i)+ dam_sedi_totalV(i)
 else
 dam_vol(i) = dam_w_vol(i)
 endif   
 else
 !dam_vol(i) = dam_vol(i) + dam_vol_temp(i)
 dam_vol(i) = dam_w_vol(i)
 endif 
 if (dam_vol(i) .gt. dam_volmax(i))then
 dam_state(i) = 1
 else
 dam_state(i) = 0
 endif
 
end do
    
end subroutine dam_checkstate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dam_write

use globals
use dam_mod
use sediment_mod
implicit none
integer i, k,m
!character*256 dam_outfile(100)
!---modified by Qin
!write(1001,'(100f20.5)') (dam_vol(i), dam_qin(dam_loc(i)) * area, dam_qout(i) * area, i = 1, dam_num) 
if(sed_switch > 0) then
do i = 1, dam_num
 !  k= dam_loc(i)
   write(900+i,'(100f20.5)')time, dam_vol(i), dam_sedi_totalV(i),dam_sedi_total_b(i)/0.6d0,dam_sedi_total_s(i)/0.6d0,dam_qin(i), dam_qout(i),dam_outflow_qb(i),dam_outflow_qs(i),dam_outflow_cs(i),(dam_sedi_qsi(i,m), m = 1,Np) ! v1.4
!   close(1001)
enddo    
else
do i = 1, dam_num
write(900+i,'(100f20.5)') time, dam_vol(i), dam_qin(i), dam_qout(i)
enddo
endif
end subroutine dam_write
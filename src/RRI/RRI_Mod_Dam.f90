module dam_mod                                   ! this file is totally modified for the RSR model

integer dam_switch
!integer dam_discharge_switch                     ! added by Qin 2021/5/24
character*256 damfile                            ! dam control file
character*256 dam_dischargefile                  ! dam discharge operation file ; added by Qin 2021/5/24
    
integer :: dam_num                               ! number of dam
character*256, allocatable :: dam_name(:)        ! dam name

integer, allocatable :: dam_kind(:)              ! dam id number
integer, allocatable :: dam_ix(:), dam_iy(:)     ! dam location (x, y)
integer, allocatable :: dam_loc(:)               ! dam location (k)
integer, allocatable :: dam_state(:)             ! dam status (full:1, not full:0)
integer, allocatable :: damflg(:)                ! dam exist at each grid-cell (exist:>=1, not exist:0)
integer, allocatable :: dam_discharge_switch(:) !for dam

real(8), allocatable :: dam_qin(:)               ! dam inflow
real(8), allocatable :: dam_vol(:)               ! storage volume
real(8), allocatable :: dam_w_vol(:)             ! storaged water volume added by Qin
real(8), allocatable :: dam_vol_temp(:)          ! storage volume (temporary)
real(8), allocatable :: dam_volmax(:)            ! maximum storage [m3]
real(8), allocatable :: dam_qout(:)              ! dam outflow
real(8), allocatable :: dam_floodq(:)            ! flood discharge[m3/s]
real(8), allocatable :: dam_reserv_area(:)       ! The area of dam reservoir [m2];  added by Qin
real(8), allocatable :: reserv_elev(:)           ! The bed elevation of the reservior [m]; added by Qin 2021/6/13 
real(8), allocatable :: dam_discharge(:,:)       ! Operated dam discharge   [m3/s]; addedby Qin 2021/5/14
real, allocatable :: dam_discharge_period(:)  ! Dam discharge time [s]; added by Qin 2021/5/14
integer:: max_dam_discharge_t                     !added by Qin 2021/5/14
character*256, allocatable:: outfile_dam(:) !for dam
integer, allocatable:: outswitch_dam(:) !for dam

end module dam_mod

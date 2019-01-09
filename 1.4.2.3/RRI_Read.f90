subroutine RRI_Read

use globals
use dam_mod
use tecout_mod

implicit none
include '..\include\cgnslib_f.h'

integer i
character*256 format_version

integer :: ier, cgns_f, icount
character(len=64):: cgns_name

!CGNSファイルを開く
!cgns_name = "input.cgn"

!引数取得
icount = nargs()
if (icount ==  2) then
    call getarg(1, cgns_name, ier)
else
    write(*,"(a)") "You should specify an argument."
    stop
endif

call cg_open_f(cgns_name, CG_MODE_READ, cgns_f, ier)
if (ier /= 0) stop "cg_open_f failed"
call cg_iric_init_f(cgns_f, ier)
!if (ier /= 0) stop "cg_iric_init_f failed"


call cg_iric_read_string_f("rrifile", rrifile, ier)

write(*,*) index(rrifile,"\", back=.true.) 
rri_dir = rrifile(1:index(rrifile,"\", back=.true.) )

write(*,"(a,a)") "rri_dir = ", rri_dir

!open(1, file = "..\RRI_Input.txt", status = 'old')
open(1, file = trim(rrifile), status = 'old')

read(1,'(a)') format_version
write(*,'("format_version : ", a)') trim(adjustl(format_version))

if( format_version .ne. "RRI_Input_Format_Ver1_4_2" ) stop "This RRI model requires RRI_Input_Format_Ver1_4_2"

read(1,*)
write(*,*)

read(1,'(a)') rainfile
read(1,'(a)') demfile
read(1,'(a)') accfile
read(1,'(a)') dirfile

rainfile = trim(rri_dir)//rainfile(3:len(rainfile))
demfile = trim(rri_dir)//demfile(3:len(demfile))
accfile = trim(rri_dir)//accfile(3:len(accfile))
dirfile = trim(rri_dir)//dirfile(3:len(dirfile))

write(*,'("rainfile : ", a)') trim(adjustl(rainfile))
write(*,'("demfile : ", a)') trim(adjustl(demfile))
write(*,'("accfile : ", a)') trim(adjustl(accfile))
write(*,'("dirfile : ", a)') trim(adjustl(dirfile))

read(1,*)
write(*,*)

read(1,*) utm
read(1,*) eight_dir
read(1,*) lasth
read(1,*) dt
read(1,*) dt_riv
read(1,*) outnum

!----------　ここから（未済） ----------
!降雨データの左上コーナ座標値とセルサイズ
!雨の読込みをiRICで行うと不要になるはず　　→　要確認
read(1,*) xllcorner_rain
read(1,*) yllcorner_rain
read(1,*) cellsize_rain_x, cellsize_rain_y

!----------　ここまで（未済） ----------

write(*,'("utm : ", i5)') utm
write(*,'("eight_dir : ", i5)') eight_dir
write(*,'("lasth : ", i8)') lasth
write(*,'("dt : ", i12)') dt
write(*,'("dt_riv : ", i8)') dt_riv
write(*,'("outnum : ", i8)') outnum
write(*,'("xllcorner_rain : ", f15.5)') xllcorner_rain
write(*,'("yllcorner_rain : ", f15.5)') yllcorner_rain
write(*,'("cellsize_rain_x : ", f15.5, "  cellsize_rain_y : ", f15.5)') cellsize_rain_x, cellsize_rain_y

read(1,*)
write(*,*)

read(1,*) ns_river

!----------　ここから（未済） ----------
!cell属性のcomplex型で処理する
read(1,*) num_of_landuse

allocate( dif(num_of_landuse) )
allocate( ns_slope(num_of_landuse), soildepth(num_of_landuse) )
allocate( gammaa(num_of_landuse) )

read(1,*) (dif(i), i = 1, num_of_landuse)
read(1,*) (ns_slope(i), i = 1, num_of_landuse)
read(1,*) (soildepth(i), i = 1, num_of_landuse)
read(1,*) (gammaa(i), i = 1, num_of_landuse)

write(*,'("ns_river : ", f12.3)') ns_river
write(*,'("num_of_landuse : ", i5)') num_of_landuse
write(*,'("dif : ", 100i5)') (dif(i), i = 1, num_of_landuse)
write(*,'("ns_slope : ", 100f12.3)') (ns_slope(i), i = 1, num_of_landuse)
write(*,'("soildepth : ", 100f12.3)') (soildepth(i), i = 1, num_of_landuse)
write(*,'("gammaa : ", 100f12.3)') (gammaa(i), i = 1, num_of_landuse)


read(1,*)
write(*,*)

allocate( ksv(num_of_landuse), faif(num_of_landuse) )

read(1,*) (ksv(i), i = 1, num_of_landuse)
read(1,*) (faif(i), i = 1, num_of_landuse)

write(*,'("ksv : ", 100e12.3)') (ksv(i), i = 1, num_of_landuse)
write(*,'("faif : ", 100f12.3)') (faif(i), i = 1, num_of_landuse)

read(1,*)
write(*,*)

allocate( ka(num_of_landuse), gammam(num_of_landuse), beta(num_of_landuse) )

read(1,*) (ka(i), i = 1, num_of_landuse)
read(1,*) (gammam(i), i = 1, num_of_landuse)
read(1,*) (beta(i), i = 1, num_of_landuse)

write(*,'("ka : ", 100e12.3)') (ka(i), i = 1, num_of_landuse)
write(*,'("gammam : ", 100f12.3)') (gammam(i), i = 1, num_of_landuse)
write(*,'("beta : ", 100f12.3)') (beta(i), i = 1, num_of_landuse)

do i = 1, num_of_landuse
 if( gammam(i) .gt. gammaa(i) ) stop "gammag must be smaller than gammaa"
enddo

read(1,*)
write(*,*)

allocate( ksg(num_of_landuse), gammag(num_of_landuse), kg0(num_of_landuse), fpg(num_of_landuse), rgl(num_of_landuse) )

read(1,*) (ksg(i), i = 1, num_of_landuse)
read(1,*) (gammag(i), i = 1, num_of_landuse)
read(1,*) (kg0(i), i = 1, num_of_landuse)
read(1,*) (fpg(i), i = 1, num_of_landuse)
read(1,*) (rgl(i), i = 1, num_of_landuse)

write(*,'("ksg : ", 100e12.3)') (ksg(i), i = 1, num_of_landuse)
write(*,'("gammag : ", 100f12.3)') (gammag(i), i = 1, num_of_landuse)
write(*,'("kg0 : ", 100e12.3)') (kg0(i), i = 1, num_of_landuse)
write(*,'("fpg : ", 100f12.3)') (fpg(i), i = 1, num_of_landuse)
write(*,'("rgl : ", 100e12.3)') (rgl(i), i = 1, num_of_landuse)

!----------　ここまで（未済） ----------


read(1,*)
write(*,*)

read(1,*) riv_thresh
read(1,*) width_param_c
read(1,*) width_param_s
read(1,*) depth_param_c
read(1,*) depth_param_s
read(1,*) height_param
read(1,*) height_limit_param

read(1,*)
write(*,*)

read(1,*) rivfile_switch

!------- ここから（ファイルで与える？要件等） ------
read(1,'(a)') widthfile
read(1,'(a)') depthfile
read(1,'(a)') heightfile

widthfile = trim(rri_dir)//widthfile(3:len(widthfile))
depthfile = trim(rri_dir)//depthfile(3:len(depthfile))
heightfile = trim(rri_dir)//heightfile(3:len(heightfile))

!------- ここまで（ファイルで与える？要件等） ------


if(rivfile_switch.eq.0) then
 write(*,'("riv_thresh : ", i7)') riv_thresh
 write(*,'("width_param_c : ", f12.2)') width_param_c
 write(*,'("width_param_s : ", f12.2)') width_param_s
 write(*,'("depth_param_c : ", f12.2)') depth_param_c
 write(*,'("depth_param_s : ", f12.2)') depth_param_s
 write(*,'("height_param : ", f12.2)') height_param
 write(*,'("height_limit_param : ", i10)') height_limit_param
else
 write(*,'("widthfile : ", a)') trim(adjustl(widthfile))
 write(*,'("depthfile : ", a)') trim(adjustl(depthfile))
 write(*,'("heightfile : ", a)') trim(adjustl(heightfile))
endif

read(1,*)
write(*,*)

read(1,*) init_slo_switch, init_riv_switch, init_gw_switch, init_gampt_ff_switch

!------- ここから（ファイルで与える？要件等） ------
read(1,"(a)") initfile_slo
read(1,'(a)') initfile_riv
read(1,'(a)') initfile_gw
read(1,'(a)') initfile_gampt_ff

initfile_slo = trim(rri_dir)//initfile_slo(3:len(initfile_slo))
initfile_riv = trim(rri_dir)//initfile_riv(3:len(initfile_riv))
initfile_gw = trim(rri_dir)//initfile_gw(3:len(initfile_gw))
initfile_gampt_ff = trim(rri_dir)//initfile_gampt_ff(3:len(initfile_gampt_ff))

!------- ここまで（ファイルで与える？要件等） ------

if(init_slo_switch.ne.0) write(*,'("initfile_slo : ", a)') trim(adjustl(initfile_slo))
if(init_riv_switch.ne.0) write(*,'("initfile_riv : ", a)') trim(adjustl(initfile_riv))
if(init_gw_switch.ne.0) write(*,'("initfile_gw : ", a)') trim(adjustl(initfile_gw))
if(init_gampt_ff_switch.ne.0) write(*,'("initfile_gampt_ff : ", a)') trim(adjustl(initfile_gampt_ff))

read(1,*)
write(*,*)

read(1,*) bound_slo_wlev_switch, bound_riv_wlev_switch

!------- ここから（ファイルで与える？要件等） ------
read(1,'(a)') boundfile_slo_wlev
read(1,'(a)') boundfile_riv_wlev

boundfile_slo_wlev = trim(rri_dir)//boundfile_slo_wlev(3:len(boundfile_slo_wlev))
boundfile_riv_wlev = trim(rri_dir)//boundfile_riv_wlev(3:len(boundfile_riv_wlev))

!------- ここまで（ファイルで与える？要件等） ------

if(bound_slo_wlev_switch.ne.0) write(*,'("boundfile_slo_wlev : ", a)') trim(adjustl(boundfile_slo_wlev))
if(bound_riv_wlev_switch.ne.0) write(*,'("boundfile_riv_wlev : ", a)') trim(adjustl(boundfile_riv_wlev))

read(1,*)
write(*,*)

read(1,*) bound_slo_disc_switch, bound_riv_disc_switch

!------- ここから（ファイルで与える？要件等） ------
read(1,'(a)') boundfile_slo_disc
read(1,'(a)') boundfile_riv_disc

boundfile_slo_disc = trim(rri_dir)//boundfile_slo_disc(3:len(boundfile_slo_disc))
boundfile_riv_disc = trim(rri_dir)//boundfile_riv_disc(3:len(boundfile_riv_disc))

!------- ここまで（ファイルで与える？要件等） ------

if(bound_slo_disc_switch.ne.0) write(*,'("boundfile_slo_disc : ", a)') trim(adjustl(boundfile_slo_disc))
if(bound_riv_disc_switch.ne.0) write(*,'("boundfile_riv_disc : ", a)') trim(adjustl(boundfile_riv_disc))

read(1,*)
write(*,*)

read(1,*) land_switch


!------- ここから（ファイルで与える？要件等） ------
read(1,'(a)') landfile

landfile = trim(rri_dir)//landfile(3:len(landfile))
!------- ここまで（ファイルで与える？要件等） ------
if(land_switch.eq.1) write(*,'("landfile : ", a)') trim(adjustl(landfile))

read(1,*)
write(*,*)

read(1,*) dam_switch
!call cg_iric_read_integer_f("dam_switch", dam_switch, ier)

!------- ここから（ファイルで与える？要件等） ------
read(1,'(a)') damfile

damfile = trim(rri_dir)//damfile(3:len(damfile))
!------- ここまで（ファイルで与える？要件等） ------
if(dam_switch.eq.1) write(*,'("damfile : ", a)') trim(adjustl(damfile))

read(1,*)
write(*,*)

read(1,*) div_switch
!call cg_iric_read_integer_f("div_switch", div_switch, ier)

!------- ここから（ファイルで与える？要件等） ------
read(1,'(a)') divfile
divfile = trim(rri_dir)//divfile(3:len(divfile))
!------- ここまで（ファイルで与える？要件等） ------
if(div_switch.eq.1) write(*,'("divfile : ", a)') trim(adjustl(divfile))

read(1,*)
write(*,*)

read(1,*) evp_switch
read(1,'(a)') evpfile
read(1,*) xllcorner_evp
read(1,*) yllcorner_evp
read(1,*) cellsize_evp_x, cellsize_evp_y

evpfile = trim(rri_dir)//evpfile(3:len(evpfile))


if( evp_switch .ne. 0 ) then
 write(*,'("evpfile : ", a)') trim(adjustl(evpfile))
 write(*,'("xllcorner_evp : ", f15.5)') xllcorner_evp
 write(*,'("yllcorner_evp : ", f15.5)') yllcorner_evp
 write(*,'("cellsize_evp_x : ", f15.5, " cellsize_evp_y : ", f15.5)') cellsize_evp_x, cellsize_evp_y
endif

read(1,*)
write(*,*)

read(1,*) sec_length_switch

!------- ここから（ファイルで与える？要件等） ------
read(1,'(a)') sec_length_file

sec_length_file = trim(rri_dir)//sec_length_file(3:len(sec_length_file))

!------- ここまで（ファイルで与える？要件等） ------

if(sec_length_switch.eq.1) write(*,'("sec_length : ", a)') trim(adjustl(sec_length_file))

read(1,*)
write(*,*)

read(1,*) sec_switch

!------- ここから（ファイルで与える？要件等） ------
read(1,'(a)') sec_map_file
read(1,'(a)') sec_file

sec_map_file = trim(rri_dir)//sec_map_file(3:len(sec_map_file))
sec_file = trim(rri_dir)//sec_file(3:len(sec_file))

!------- ここまで（ファイルで与える？要件等） ------

if(sec_switch.eq.1) write(*,'("sec_map_file : ", a)') trim(adjustl(sec_map_file))
if(sec_switch.eq.1) write(*,'("sec_file : ", a)') trim(adjustl(sec_file))

read(1,*)
write(*,*)

!read(1,*) emb_switch
!read(1,'(a)') embrfile
!read(1,'(a)') embbfile
!if(emb_switch.eq.1) write(*,'("embrfile : ", a)') trim(adjustl(embrfile))
!if(emb_switch.eq.1) write(*,'("embbfile : ", a)') trim(adjustl(embbfile))

read(1,*) outswitch_hs, outswitch_hr, outswitch_hg, outswitch_qr, outswitch_qu, outswitch_qv, &
          outswitch_gu, outswitch_gv, outswitch_gampt_ff, outswitch_storage

read(1,'(a)') outfile_hs
read(1,'(a)') outfile_hr
read(1,'(a)') outfile_hg
read(1,'(a)') outfile_qr
read(1,'(a)') outfile_qu
read(1,'(a)') outfile_qv
read(1,'(a)') outfile_gu
read(1,'(a)') outfile_gv
read(1,'(a)') outfile_gampt_ff
read(1,'(a)') outfile_storage



outfile_hs = trim(rri_dir)//outfile_hs(3:len(outfile_hs))
outfile_hr = trim(rri_dir)//outfile_hr(3:len(outfile_hr))
outfile_hg = trim(rri_dir)//outfile_hg(3:len(outfile_hg))
outfile_qr = trim(rri_dir)//outfile_qr(3:len(outfile_qr))
outfile_qu = trim(rri_dir)//outfile_qu(3:len(outfile_qu))
outfile_qv = trim(rri_dir)//outfile_qv(3:len(outfile_qv))
outfile_gu = trim(rri_dir)//outfile_gu(3:len(outfile_gu))
outfile_gv = trim(rri_dir)//outfile_gv(3:len(outfile_gv))
outfile_gampt_ff = trim(rri_dir)//outfile_gampt_ff(3:len(outfile_gampt_ff))
outfile_storage = trim(rri_dir)//outfile_storage(3:len(outfile_storage))


if(outswitch_hs .ne. 0) write(*,'("outfile_hs : ", a)') trim(adjustl(outfile_hs))
if(outswitch_hr .ne. 0) write(*,'("outfile_hr : ", a)') trim(adjustl(outfile_hr))
if(outswitch_hg .ne. 0) write(*,'("outfile_hg : ", a)') trim(adjustl(outfile_hg))
if(outswitch_qr .ne. 0) write(*,'("outfile_qr : ", a)') trim(adjustl(outfile_qr))
if(outswitch_qu .ne. 0) write(*,'("outfile_qu : ", a)') trim(adjustl(outfile_qu))
if(outswitch_qv .ne. 0) write(*,'("outfile_qv : ", a)') trim(adjustl(outfile_qv))
if(outswitch_gu .ne. 0) write(*,'("outfile_gu : ", a)') trim(adjustl(outfile_gu))
if(outswitch_gv .ne. 0) write(*,'("outfile_gv : ", a)') trim(adjustl(outfile_gv))
if(outswitch_gampt_ff .ne. 0) write(*,'("outfile_gampt_ff : ", a)') trim(adjustl(outfile_gampt_ff))
if(outswitch_storage .ne. 0) write(*,'("outfile_storage : ", a)') trim(adjustl(outfile_storage))

read(1,*)
write(*,*)

read(1,*) hydro_switch
read(1,'(a)') location_file
location_file = trim(rri_dir)//location_file(3:len(location_file))

if(hydro_switch .eq. 1) write(*,'("location_file : ", a)') trim(adjustl(location_file))

write(*,*)

close(1)
call cg_close_f(cgns_f, ier)

! Parameter Check
do i = 1, num_of_landuse
 if( ksv(i) .gt. 0.d0 .and. ka(i) .gt. 0.d0 ) &
  stop "Error: both ksv and ka are non-zero."
 if( gammam(i) .gt. gammaa(i) ) &
  stop "Error: gammam must be smaller than gammaa."
enddo

! Set da, dm and infilt_limit
allocate( da(num_of_landuse), dm(num_of_landuse), infilt_limit(num_of_landuse))
da(:) = 0.d0
dm(:) = 0.d0
infilt_limit(:) = 0.d0
do i = 1, num_of_landuse
 if( soildepth(i) .gt. 0.d0 .and. ksv(i) .gt. 0.d0 ) infilt_limit(i) = soildepth(i) * gammaa(i)
 if( soildepth(i) .gt. 0.d0 .and. ka(i) .gt. 0.d0 ) da(i)= soildepth(i) * gammaa(i)
 if( soildepth(i) .gt. 0.d0 .and. ka(i) .gt. 0.d0 .and. gammam(i) .gt. 0.d0 ) &
  dm(i) = soildepth(i) * gammam(i)
enddo

! if ksg(i) = 0.d0 -> no gw calculation
gw_switch = 0
do i = 1, num_of_landuse
 if( ksg(i) .gt. 0.d0 ) then
  gw_switch = 1
 else
  gammag(i) = 0.d0
  kg0(i) = 0.d0
  fpg(i) = 0.d0
  rgl(i) = 0.d0
 endif
enddo

end subroutine RRI_Read

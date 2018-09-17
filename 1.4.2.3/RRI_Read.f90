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
cgns_name = "Case1.cgn"

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

open(1, file = "RRI_Input.txt", status = 'old')

read(1,'(a)') format_version
write(*,'("format_version : ", a)') trim(adjustl(format_version))

if( format_version .ne. "RRI_Input_Format_Ver1_4_2" ) stop "This RRI model requires RRI_Input_Format_Ver1_4_2"

read(1,*)
write(*,*)

read(1,'(a)') rainfile
read(1,'(a)') demfile
read(1,'(a)') accfile
read(1,'(a)') dirfile

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


call cg_iric_read_integer_f("utm", utm, ier)
call cg_iric_read_integer_f("eight_dir", eight_dir, ier)
call cg_iric_read_integer_f("lasth", lasth, ier)
call cg_iric_read_integer_f("dt", dt, ier)
call cg_iric_read_integer_f("dt_riv", dt_riv, ier)
call cg_iric_read_integer_f("outnum", outnum, ier)

!----------　ここから（未済） ----------
!降雨データの左上コーナ座標値とセルサイズ
!雨の読込みをiRICで行うと不要になるはず　　→　要確認
read(1,*) xllcorner_rain
read(1,*) yllcorner_rain
read(1,*) cellsize_rain_x, cellsize_rain_y

call cg_iric_read_real_f("xllcorner_rain", xllcorner_rain, ier)
call cg_iric_read_real_f("yllcorner_rain", yllcorner_rain, ier)
call cg_iric_read_real_f("cellsize_rain_x", cellsize_rain_x, ier)
call cg_iric_read_real_f("cellsize_rain_y", cellsize_rain_y, ier)
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
call cg_iric_read_real_f("ns_river", ns_river, ier)


!----------　ここから（未済） ----------
!cell属性のcomplex型で処理する
read(1,*) num_of_landuse
call cg_iric_read_integer_f("num_of_landuse", num_of_landuse, ier)

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

call cg_iric_read_integer_f("riv_thresh", riv_thresh, ier)
call cg_iric_read_real_f("width_param_c", width_param_c, ier)
call cg_iric_read_real_f("width_param_s", width_param_s, ier)
call cg_iric_read_real_f("depth_param_c", depth_param_c, ier)
call cg_iric_read_real_f("depth_param_s", depth_param_s, ier)
call cg_iric_read_real_f("height_param", height_param, ier)
call cg_iric_read_real_f("height_limit_param", height_limit_param, ier)


read(1,*)
write(*,*)

read(1,*) rivfile_switch
call cg_iric_read_integer_f("rivfile_switch", rivfile_switch, ier)

!------- ここから（ファイルで与える？要件等） ------
read(1,'(a)') widthfile
read(1,'(a)') depthfile
read(1,'(a)') heightfile
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
call cg_iric_read_integer_f("init_slo_switch", init_slo_switch, ier)
call cg_iric_read_integer_f("init_riv_switch", init_riv_switch, ier)
call cg_iric_read_integer_f("init_gw_switch", init_gw_switch, ier)
call cg_iric_read_integer_f("init_gampt_ff_switch", init_gampt_ff_switch, ier)

!------- ここから（ファイルで与える？要件等） ------
read(1,"(a)") initfile_slo
read(1,'(a)') initfile_riv
read(1,'(a)') initfile_gw
read(1,'(a)') initfile_gampt_ff
!------- ここまで（ファイルで与える？要件等） ------

if(init_slo_switch.ne.0) write(*,'("initfile_slo : ", a)') trim(adjustl(initfile_slo))
if(init_riv_switch.ne.0) write(*,'("initfile_riv : ", a)') trim(adjustl(initfile_riv))
if(init_gw_switch.ne.0) write(*,'("initfile_gw : ", a)') trim(adjustl(initfile_gw))
if(init_gampt_ff_switch.ne.0) write(*,'("initfile_gampt_ff : ", a)') trim(adjustl(initfile_gampt_ff))

read(1,*)
write(*,*)

read(1,*) bound_slo_wlev_switch, bound_riv_wlev_switch
call cg_iric_read_integer_f("bound_slo_wlev_switch", bound_slo_wlev_switch, ier)
call cg_iric_read_integer_f("bound_riv_wlev_switch", bound_riv_wlev_switch, ier)

!------- ここから（ファイルで与える？要件等） ------
read(1,'(a)') boundfile_slo_wlev
read(1,'(a)') boundfile_riv_wlev
!------- ここまで（ファイルで与える？要件等） ------

if(bound_slo_wlev_switch.ne.0) write(*,'("boundfile_slo_wlev : ", a)') trim(adjustl(boundfile_slo_wlev))
if(bound_riv_wlev_switch.ne.0) write(*,'("boundfile_riv_wlev : ", a)') trim(adjustl(boundfile_riv_wlev))

read(1,*)
write(*,*)

read(1,*) bound_slo_disc_switch, bound_riv_disc_switch
call cg_iric_read_integer_f("bound_slo_disc_switch", bound_slo_disc_switch, ier)
call cg_iric_read_integer_f("bound_riv_disc_switch", bound_riv_disc_switch, ier)

!------- ここから（ファイルで与える？要件等） ------
read(1,'(a)') boundfile_slo_disc
read(1,'(a)') boundfile_riv_disc
!------- ここまで（ファイルで与える？要件等） ------

if(bound_slo_disc_switch.ne.0) write(*,'("boundfile_slo_disc : ", a)') trim(adjustl(boundfile_slo_disc))
if(bound_riv_disc_switch.ne.0) write(*,'("boundfile_riv_disc : ", a)') trim(adjustl(boundfile_riv_disc))

read(1,*)
write(*,*)

read(1,*) land_switch
call cg_iric_read_integer_f("land_switch", land_switch, ier)

!------- ここから（ファイルで与える？要件等） ------
read(1,'(a)') landfile
!------- ここまで（ファイルで与える？要件等） ------
if(land_switch.eq.1) write(*,'("landfile : ", a)') trim(adjustl(landfile))

read(1,*)
write(*,*)

read(1,*) dam_switch
call cg_iric_read_integer_f("dam_switch", dam_switch, ier)

!------- ここから（ファイルで与える？要件等） ------
read(1,'(a)') damfile
!------- ここまで（ファイルで与える？要件等） ------
if(dam_switch.eq.1) write(*,'("damfile : ", a)') trim(adjustl(damfile))

read(1,*)
write(*,*)

read(1,*) div_switch
call cg_iric_read_integer_f("div_switch", div_switch, ier)

!------- ここから（ファイルで与える？要件等） ------
read(1,'(a)') divfile
!------- ここまで（ファイルで与える？要件等） ------
if(div_switch.eq.1) write(*,'("divfile : ", a)') trim(adjustl(divfile))

read(1,*)
write(*,*)

read(1,*) evp_switch
read(1,'(a)') evpfile
read(1,*) xllcorner_evp
read(1,*) yllcorner_evp
read(1,*) cellsize_evp_x, cellsize_evp_y

call cg_iric_read_integer_f("evp_switch", evp_switch, ier)
call cg_iric_read_string_f("evpfile", evpfile, ier)
call cg_iric_read_real_f("xllcorner_evp", xllcorner_evp, ier)
call cg_iric_read_real_f("yllcorner_evp", yllcorner_evp, ier)
call cg_iric_read_real_f("cellsize_evp_x", cellsize_evp_x, ier)
call cg_iric_read_real_f("cellsize_evp_y", cellsize_evp_y, ier)


if( evp_switch .ne. 0 ) then
 write(*,'("evpfile : ", a)') trim(adjustl(evpfile))
 write(*,'("xllcorner_evp : ", f15.5)') xllcorner_evp
 write(*,'("yllcorner_evp : ", f15.5)') yllcorner_evp
 write(*,'("cellsize_evp_x : ", f15.5, " cellsize_evp_y : ", f15.5)') cellsize_evp_x, cellsize_evp_y
endif

read(1,*)
write(*,*)

read(1,*) sec_length_switch
call cg_iric_read_integer_f("sec_length_switch", sec_length_switch, ier)

!------- ここから（ファイルで与える？要件等） ------
read(1,'(a)') sec_length_file
!------- ここまで（ファイルで与える？要件等） ------

if(sec_length_switch.eq.1) write(*,'("sec_length : ", a)') trim(adjustl(sec_length_file))

read(1,*)
write(*,*)

read(1,*) sec_switch
call cg_iric_read_integer_f("sec_switch", sec_switch, ier)

!------- ここから（ファイルで与える？要件等） ------
read(1,'(a)') sec_map_file
read(1,'(a)') sec_file
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

call cg_iric_read_integer_f("outswitch_hs", outswitch_hs, ier)
call cg_iric_read_integer_f("outswitch_hr", outswitch_hr, ier)
call cg_iric_read_integer_f("outswitch_hg", outswitch_hg, ier)
call cg_iric_read_integer_f("outswitch_qr", outswitch_qr, ier)
call cg_iric_read_integer_f("outswitch_qu", outswitch_qu, ier)
call cg_iric_read_integer_f("outswitch_qv", outswitch_qv, ier)
call cg_iric_read_integer_f("outswitch_gu", outswitch_gu, ier)
call cg_iric_read_integer_f("outswitch_gv", outswitch_gv, ier)
call cg_iric_read_integer_f("outswitch_gampt_ff", outswitch_gampt_ff, ier)
call cg_iric_read_integer_f("outswitch_storage", outswitch_storage, ier)

outfile_hs =""
call cg_iric_read_string_f("outfile_hs_folder", outfile_hs, ier)
outfile_hs = trim(adjustl(outfile_hs))//"/hs_"

outfile_hr =""
call cg_iric_read_string_f("outfile_hr_folder", outfile_hr, ier)
outfile_hr = trim(adjustl(outfile_hr))//"/hr_"

outfile_hg = ""
call cg_iric_read_string_f("outfile_hg_folder", outfile_hg, ier)
outfile_hg = trim(adjustl(outfile_hg))//"/hg_"

outfile_qr=""
call cg_iric_read_string_f("outfile_qr_folder", outfile_qr, ier)
outfile_qr = trim(adjustl(outfile_qr))//"/qr_"

outfile_qu=""
call cg_iric_read_string_f("outfile_qu_folder", outfile_qu, ier)
outfile_qu = trim(adjustl(outfile_qu))//"/qu_"

outfile_qv=""
call cg_iric_read_string_f("outfile_qv_folder", outfile_qv, ier)
outfile_qv = trim(adjustl(outfile_qv))//"/qv_"

outfile_gu=""
call cg_iric_read_string_f("outfile_gu_folder", outfile_gu, ier)
outfile_gu = trim(adjustl(outfile_gu))//"/gu_"

outfile_gv=""
call cg_iric_read_string_f("outfile_gv_folder", outfile_gv, ier)
outfile_gv = trim(adjustl(outfile_gv))//"/gv_"

outfile_gampt_ff=""
call cg_iric_read_string_f("outfile_gampt_ff_folder", outfile_gampt_ff, ier)
outfile_gampt_ff = trim(adjustl(outfile_gampt_ff))//"/gampt_ff_"

outfile_storage=""
call cg_iric_read_string_f("outfile_storage_folder", outfile_storage, ier)
outfile_storage = trim(adjustl(outfile_storage))//"/storage.dat"


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

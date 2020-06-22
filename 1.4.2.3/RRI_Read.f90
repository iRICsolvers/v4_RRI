subroutine RRI_Read

use globals
use dam_mod
use tecout_mod

implicit none
include 'cgnslib_f.h'

integer i
character*256 format_version

integer :: ier, cgns_f, icount
character(len=64):: cgns_name

!--------------------------------------------------
!CGNSファイルを開く
!--------------------------------------------------
icount = iargc()
if (icount ==  1) then
    call getarg(1, cgns_name)
else
    write(*,"(a)") "You should specify an argument."
    stop
endif
call cg_open_f(cgns_name, CG_MODE_READ, cgns_f, ier)
if (ier /= 0) stop "cg_open_f failed"
call cg_iric_init_f(cgns_f, ier)

!--------------------------------------------------
!RRI バージョン情報
!--------------------------------------------------
!read(1,'(a)') format_version
call cg_iric_read_string_f("rri_ver", format_version, ier)

write(*,'("format_version : ", a)') trim(adjustl(format_version))
if( format_version .ne. "Ver1_4_2 for iRIC" ) stop "This RRI model requires RRI_Input_Format_Ver1_4_2"
write(*,*)

!--------------------------------------------------
!外部ファイル　→　格子属性として設定
!--------------------------------------------------
!read(1,'(a)') rainfile
!read(1,'(a)') demfile
!read(1,'(a)') accfile
!read(1,'(a)') dirfile
!
!!rainfile = trim(rri_dir)//rainfile(3:len(rainfile))
!!demfile = trim(rri_dir)//demfile(3:len(demfile))
!!accfile = trim(rri_dir)//accfile(3:len(accfile))
!!dirfile = trim(rri_dir)//dirfile(3:len(dirfile))
!
!write(*,'("rainfile : ", a)') trim(adjustl(rainfile))
!write(*,'("demfile : ", a)') trim(adjustl(demfile))
!write(*,'("accfile : ", a)') trim(adjustl(accfile))
!write(*,'("dirfile : ", a)') trim(adjustl(dirfile))
!
!read(1,*)
!write(*,*)

!--------------------------------------------------
!基本オプション
!--------------------------------------------------
!read(1,*) utm
call cg_iric_read_integer_f("utm", utm, ier)

!read(1,*) eight_dir
call cg_iric_read_integer_f("eight_dir", eight_dir, ier)

write(*,'("utm : ", i5)') utm
write(*,'("eight_dir : ", i5)') eight_dir
write(*,*)

!--------------------------------------------------
!時間条件
!--------------------------------------------------
!read(1,*) lasth
call cg_iric_read_integer_f("lasth", lasth, ier)

!read(1,*) dt
call cg_iric_read_integer_f("dt", dt, ier)

!read(1,*) dt_riv
call cg_iric_read_integer_f("dt_riv", dt_riv, ier)

!read(1,*) outnum
call cg_iric_read_integer_f("outnum", outnum, ier)

write(*,'("lasth : ", i8)') lasth
write(*,'("dt : ", i12)') dt
write(*,'("dt_riv : ", i8)') dt_riv
write(*,'("outnum : ", i8)') outnum
write(*,*)

!--------------------------------------------------
!降雨条件　格子属性として設定
!--------------------------------------------------
!read(1,*) xllcorner_rain
!read(1,*) yllcorner_rain
!read(1,*) cellsize_rain_x, cellsize_rain_y
!write(*,'("xllcorner_rain : ", f15.5)') xllcorner_rain
!write(*,'("yllcorner_rain : ", f15.5)') yllcorner_rain
!write(*,'("cellsize_rain_x : ", f15.5, "  cellsize_rain_y : ", f15.5)') cellsize_rain_x, cellsize_rain_y
write(*,"(a)") "Rain data can be set as the grid attribute.　"
write(*,*)

!--------------------------------------------------
!slopeパラメータ　→　格子セル属性　複合条件として設定
!--------------------------------------------------
!read(1,*) num_of_landuse
call cg_iric_read_complex_count_f("landuse_c", num_of_landuse, ier)
!
allocate( dif(num_of_landuse) )
allocate( ns_slope(num_of_landuse), soildepth(num_of_landuse) )
allocate( gammaa(num_of_landuse) )
!
do i = 1, num_of_landuse
	!read(1,*) (dif(i), i = 1, num_of_landuse)
	call cg_iric_read_complex_integer_f("landuse_c", i, "dif", dif(i), ier)

	!read(1,*) (ns_slope(i), i = 1, num_of_landuse)
	call cg_iric_read_complex_real_f("landuse_c", i, "ns_slope", ns_slope(i), ier)

	!read(1,*) (soildepth(i), i = 1, num_of_landuse)
	call cg_iric_read_complex_real_f("landuse_c", i, "soildepth", soildepth(i), ier)

	!read(1,*) (gammaa(i), i = 1, num_of_landuse)
	call cg_iric_read_complex_real_f("landuse_c", i, "gammaa", gammaa(i), ier)
end do
!

write(*,'("num_of_landuse : ", i5)') num_of_landuse
write(*,'("dif : ", 100i5)') (dif(i), i = 1, num_of_landuse)
write(*,'("ns_slope : ", 100f12.3)') (ns_slope(i), i = 1, num_of_landuse)
write(*,'("soildepth : ", 100f12.3)') (soildepth(i), i = 1, num_of_landuse)
write(*,'("gammaa : ", 100f12.3)') (gammaa(i), i = 1, num_of_landuse)
!
!
!read(1,*)
write(*,*)
!
allocate( ksv(num_of_landuse), faif(num_of_landuse) )
!
do i = 1, num_of_landuse
	!read(1,*) (ksv(i), i = 1, num_of_landuse)
	call cg_iric_read_complex_real_f("landuse_c", i, "ksv", ksv(i), ier)

	!read(1,*) (faif(i), i = 1, num_of_landuse)
	call cg_iric_read_complex_real_f("landuse_c", i, "faif", faif(i), ier)
end do
!
write(*,'("ksv : ", 100e12.3)') (ksv(i), i = 1, num_of_landuse)
write(*,'("faif : ", 100f12.3)') (faif(i), i = 1, num_of_landuse)
!
!read(1,*)
write(*,*)
!
allocate( ka(num_of_landuse), gammam(num_of_landuse), beta(num_of_landuse) )
!

do i = 1, num_of_landuse
	!read(1,*) (ka(i), i = 1, num_of_landuse)
	call cg_iric_read_complex_real_f("landuse_c", i, "ka", ka(i), ier)

	!read(1,*) (gammam(i), i = 1, num_of_landuse)
	call cg_iric_read_complex_real_f("landuse_c", i, "gammam", gammam(i), ier)

	!read(1,*) (beta(i), i = 1, num_of_landuse)
	call cg_iric_read_complex_real_f("landuse_c", i, "beta", beta(i), ier)
end do

write(*,'("ka : ", 100e12.3)') (ka(i), i = 1, num_of_landuse)
write(*,'("gammam : ", 100f12.3)') (gammam(i), i = 1, num_of_landuse)
write(*,'("beta : ", 100f12.3)') (beta(i), i = 1, num_of_landuse)

do i = 1, num_of_landuse
 if( gammam(i) .gt. gammaa(i) ) stop "gammag must be smaller than gammaa"
enddo
!
!read(1,*)
write(*,*)
!
allocate( ksg(num_of_landuse), gammag(num_of_landuse), kg0(num_of_landuse), fpg(num_of_landuse), rgl(num_of_landuse) )
!
do i = 1, num_of_landuse
	!read(1,*) (ksg(i), i = 1, num_of_landuse)
	call cg_iric_read_complex_real_f("landuse_c", i, "ksg", ksg(i), ier)
	
	!read(1,*) (gammag(i), i = 1, num_of_landuse)
	call cg_iric_read_complex_real_f("landuse_c", i, "gammag", gammag(i), ier)
	
	!read(1,*) (kg0(i), i = 1, num_of_landuse)
	call cg_iric_read_complex_real_f("landuse_c", i, "kg0", kg0(i), ier)
	
	!read(1,*) (fpg(i), i = 1, num_of_landuse)
	call cg_iric_read_complex_real_f("landuse_c", i, "fpg", fpg(i), ier)
	
	!read(1,*) (rgl(i), i = 1, num_of_landuse)
	call cg_iric_read_complex_real_f("landuse_c", i, "rgl", rgl(i), ier)
end do
!
write(*,'("ksg : ", 100e12.3)') (ksg(i), i = 1, num_of_landuse)
write(*,'("gammag : ", 100f12.3)') (gammag(i), i = 1, num_of_landuse)
write(*,'("kg0 : ", 100e12.3)') (kg0(i), i = 1, num_of_landuse)
write(*,'("fpg : ", 100f12.3)') (fpg(i), i = 1, num_of_landuse)
write(*,'("rgl : ", 100e12.3)') (rgl(i), i = 1, num_of_landuse)
!read(1,*)
write(*,*)


!--------------------------------------------------
!河道パラメータ
!--------------------------------------------------
!read(1,*) ns_river
call cg_iric_read_real_f("ns_river", ns_river, ier)
write(*,'("ns_river : ", f12.3)') ns_river
write(*,*)

!read(1,*) riv_thresh
call cg_iric_read_integer_f("riv_thresh", riv_thresh, ier)
write(*,'("riv_thresh : ", i7)') riv_thresh
write(*,*)

!read(1,*) width_param_c
!read(1,*) width_param_s
call cg_iric_read_real_f("width_param_c", width_param_c, ier)
call cg_iric_read_real_f("width_param_s", width_param_s, ier)

!read(1,*) depth_param_c
!read(1,*) depth_param_s
call cg_iric_read_real_f("depth_param_c", depth_param_c, ier)
call cg_iric_read_real_f("depth_param_s", depth_param_s, ier)

!read(1,*) height_param
!read(1,*) height_limit_param
call cg_iric_read_real_f("height_param", height_param, ier)
call cg_iric_read_real_f("height_limit_param", height_limit_param, ier)

!read(1,*) rivfile_switch
call cg_iric_read_integer_f("rivfile_switch", rivfile_switch, ier)

!read(1,'(a)') widthfile
!read(1,'(a)') depthfile
!read(1,'(a)') heightfile
widthfile=''; depthfile=''; heightfile=''
if(rivfile_switch == 1)then
	call cg_iric_read_string_f("widthfile", widthfile, ier)
	call cg_iric_read_string_f("depthfile", depthfile, ier)
	call cg_iric_read_string_f("heightfile", heightfile, ier)
end if

if(rivfile_switch.eq.0) then
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

!read(1,*)
write(*,*)

!--------------------------------------------------
!hotstart用　初期条件
!--------------------------------------------------
!read(1,*) init_slo_switch, init_riv_switch, init_gw_switch, init_gampt_ff_switch
call cg_iric_read_integer_f("init_slo_switch", init_slo_switch, ier)
call cg_iric_read_integer_f("init_riv_switch", init_riv_switch, ier)
call cg_iric_read_integer_f("init_gw_switch", init_gw_switch, ier)
call cg_iric_read_integer_f("init_gampt_ff_switch", init_gampt_ff_switch, ier)

!read(1,"(a)") initfile_slo
!read(1,'(a)') initfile_riv
!read(1,'(a)') initfile_gw
!read(1,'(a)') initfile_gampt_ff

initfile_slo=''; initfile_riv=''; initfile_gw=''; initfile_gampt_ff=''
if(init_slo_switch == 1) call cg_iric_read_string_f("initfile_slo", initfile_slo, ier)
if(init_riv_switch == 1) call cg_iric_read_string_f("initfile_riv", initfile_riv, ier)
if(init_gw_switch == 1) call cg_iric_read_string_f("initfile_gw", initfile_gw, ier)
if(init_gampt_ff_switch == 1) call cg_iric_read_string_f("initfile_gampt_ff", initfile_gampt_ff, ier)

if(init_slo_switch.ne.0) write(*,'("initfile_slo : ", a)') trim(adjustl(initfile_slo))
if(init_riv_switch.ne.0) write(*,'("initfile_riv : ", a)') trim(adjustl(initfile_riv))
if(init_gw_switch.ne.0) write(*,'("initfile_gw : ", a)') trim(adjustl(initfile_gw))
if(init_gampt_ff_switch.ne.0) write(*,'("initfile_gampt_ff : ", a)') trim(adjustl(initfile_gampt_ff))

!read(1,*)
write(*,*)


!--------------------------------------------------
!hs, hr境界条件　→　未実装　要確認
!--------------------------------------------------
!read(1,*) bound_slo_wlev_switch, bound_riv_wlev_switch
!read(1,'(a)') boundfile_slo_wlev
!read(1,'(a)') boundfile_riv_wlev
!if(bound_slo_wlev_switch.ne.0) write(*,'("boundfile_slo_wlev : ", a)') trim(adjustl(boundfile_slo_wlev))
!if(bound_riv_wlev_switch.ne.0) write(*,'("boundfile_riv_wlev : ", a)') trim(adjustl(boundfile_riv_wlev))
!read(1,*)
write(*,"(a)") "bound_slo_wlev_switch, bound_riv_wlev_switch have not implemented yet."
write(*,*)


!--------------------------------------------------
!qs, qr境界条件　→　未実装　要確認
!--------------------------------------------------
!read(1,*) bound_slo_disc_switch, bound_riv_disc_switch
!read(1,'(a)') boundfile_slo_disc
!read(1,'(a)') boundfile_riv_disc
!if(bound_slo_disc_switch.ne.0) write(*,'("boundfile_slo_disc : ", a)') trim(adjustl(boundfile_slo_disc))
!if(bound_riv_disc_switch.ne.0) write(*,'("boundfile_riv_disc : ", a)') trim(adjustl(boundfile_riv_disc))
!read(1,*)
write(*,"(a)") "bound_slo_disc_switch, bound_riv_disc_switch have not implemented yet."
write(*,*)

!--------------------------------------------------
!土地利用条件　→　格子属性として与える
!--------------------------------------------------
!read(1,*) land_switch
!read(1,'(a)') landfile
!if(land_switch.eq.1) write(*,'("landfile : ", a)') trim(adjustl(landfile))
!
!read(1,*)
write(*,"(a)") "landuse can be set as the grid attribute."
write(*,*)

!--------------------------------------------------
!Dam条件　→　境界条件属性として与える
!--------------------------------------------------
!read(1,*) dam_switch
!read(1,'(a)') damfile
!if(dam_switch.eq.1) write(*,'("damfile : ", a)') trim(adjustl(damfile))
!read(1,*)
write(*,"(a)") "Dam condition can be set as the boundary condition."
write(*,*)

!--------------------------------------------------
!div条件　→　未実装　要確認
!--------------------------------------------------
!read(1,*) div_switch
!call cg_iric_read_integer_f("div_switch", div_switch, ier)
!read(1,'(a)') divfile
!if(div_switch.eq.1) write(*,'("divfile : ", a)') trim(adjustl(divfile))
!read(1,*)
write(*,"(a)") "div condition has not implemented yet."
write(*,*)

!--------------------------------------------------
!蒸発条件　→　格子属性として与える
!--------------------------------------------------
!read(1,*) evp_switch
!read(1,'(a)') evpfile
!read(1,*) xllcorner_evp
!read(1,*) yllcorner_evp
!read(1,*) cellsize_evp_x, cellsize_evp_y
!if( evp_switch .ne. 0 ) then
! write(*,'("evpfile : ", a)') trim(adjustl(evpfile))
! write(*,'("xllcorner_evp : ", f15.5)') xllcorner_evp
! write(*,'("yllcorner_evp : ", f15.5)') yllcorner_evp
! write(*,'("cellsize_evp_x : ", f15.5, " cellsize_evp_y : ", f15.5)') cellsize_evp_x, cellsize_evp_y
!endif
!read(1,*)
write(*,"(a)") "evp condition can be set as the grid attribute."
write(*,*)

!--------------------------------------------------
!河道断面データ　→　未実装　要確認
!--------------------------------------------------
!read(1,*) sec_length_switch
!read(1,'(a)') sec_length_file
!if(sec_length_switch.eq.1) write(*,'("sec_length : ", a)') trim(adjustl(sec_length_file))
!read(1,*)
!write(*,*)
!read(1,*) sec_switch
!read(1,'(a)') sec_map_file
!read(1,'(a)') sec_file
!if(sec_switch.eq.1) write(*,'("sec_map_file : ", a)') trim(adjustl(sec_map_file))
!if(sec_switch.eq.1) write(*,'("sec_file : ", a)') trim(adjustl(sec_file))
!read(1,*)
write(*,"(a)") "Cross section condition has not implemented yet."
write(*,*)


!--------------------------------------------------
!??? 　→　要確認
!--------------------------------------------------
!read(1,*) emb_switch
!read(1,'(a)') embrfile
!read(1,'(a)') embbfile
!if(emb_switch.eq.1) write(*,'("embrfile : ", a)') trim(adjustl(embrfile))
!if(emb_switch.eq.1) write(*,'("embbfile : ", a)') trim(adjustl(embbfile))
write(*,*)

!--------------------------------------------------
!計算結果出力・ファイル
!--------------------------------------------------
!read(1,*) outswitch_hs, outswitch_hr, outswitch_hg, outswitch_qr, outswitch_qu, outswitch_qv, &
!          outswitch_gu, outswitch_gv, outswitch_gampt_ff, outswitch_storage

!read(1,'(a)') outfile_hs
outfile_hs = ''
call cg_iric_read_integer_f("outswitch_hs", outswitch_hs, ier)
if(outswitch_hs == 1) call cg_iric_read_string_f("outfile_hs", outfile_hs, ier)

!read(1,'(a)') outfile_hr
outfile_hr = ''
call cg_iric_read_integer_f("outswitch_hr", outswitch_hr, ier)
if(outswitch_hr == 1) call cg_iric_read_string_f("outfile_hr", outfile_hr, ier)

!read(1,'(a)') outfile_hg
outfile_hg = ''
call cg_iric_read_integer_f("outswitch_hg", outswitch_hg, ier)
if(outswitch_hg == 1) call cg_iric_read_string_f("outfile_hg", outfile_hg, ier)

!read(1,'(a)') outfile_qr
outfile_qr = ''
call cg_iric_read_integer_f("outswitch_qr", outswitch_qr, ier)
if(outswitch_qr == 1) call cg_iric_read_string_f("outfile_qr", outfile_qr, ier)

!read(1,'(a)') outfile_qu
outfile_qu = ''
call cg_iric_read_integer_f("outswitch_qu", outswitch_qu, ier)
if(outswitch_qu == 1) call cg_iric_read_string_f("outfile_qu", outfile_qu, ier)

!read(1,'(a)') outfile_qv
outfile_qv = ''
call cg_iric_read_integer_f("outswitch_qv", outswitch_qv, ier)
if(outswitch_qv == 1) call cg_iric_read_string_f("outfile_qv", outfile_qv, ier)

!read(1,'(a)') outfile_gu
outfile_gu = ''
call cg_iric_read_integer_f("outswitch_gu", outswitch_gu, ier)
if(outswitch_gu == 1) call cg_iric_read_string_f("outfile_gu", outfile_gu, ier)

!read(1,'(a)') outfile_gv
outfile_gv = ''
call cg_iric_read_integer_f("outswitch_gv", outswitch_gv, ier)
if(outswitch_gv == 1) call cg_iric_read_string_f("outfile_gv", outfile_gv, ier)

!read(1,'(a)') outfile_gampt_ff
outfile_gampt_ff = ''
call cg_iric_read_integer_f("outswitch_gampt_ff", outswitch_gampt_ff, ier)
if(outswitch_gampt_ff == 1) call cg_iric_read_string_f("outfile_gampt_ff", outfile_gampt_ff, ier)

!read(1,'(a)') outfile_storage
outfile_storage = ''
call cg_iric_read_integer_f("outswitch_storage", outswitch_storage, ier)
if(outswitch_storage == 1) call cg_iric_read_string_f("outfile_storage", outfile_storage, ier)


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

!read(1,*)
write(*,*)

!--------------------------------------------------
!iRICの基本機能で対応
!--------------------------------------------------
!read(1,*) hydro_switch
!read(1,'(a)') location_file
!!location_file = trim(rri_dir)//location_file(3:len(location_file))
!
!if(hydro_switch .eq. 1) write(*,'("location_file : ", a)') trim(adjustl(location_file))
!
write(*,*)

!close(1)
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

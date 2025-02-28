subroutine RRI_Read

    use globals
    use dam_mod
    use tecout_mod
    use sediment_mod    !added for RSR model 20240724
    use runge_mod       !added for RSR model 20240724
    use RRI_iric
    use iric

    implicit none

    integer i, num_of_bound_point
    character*256 format_version

    integer :: ier
!--------------------------------------------------
!CGNSファイルを開く
!--------------------------------------------------
    icount = iargc()
    if (icount == 1) then
        call getarg(1, cgns_name)
    else
        write (*, "(a)") "You should specify an argument."
        stop
    end if
    call cg_iric_open(cgns_name, IRIC_MODE_MODIFY, cgns_f, ier)
    if (ier /= 0) stop "cg_iric_open failed"

!--------------------------------------------------
!RRI バージョン情報
!--------------------------------------------------
!read(1,'(a)') format_version
!call cg_iric_read_string(cgns_f, "rri_ver", format_version, ier)
!
!write(*,'("format_version : ", a)') trim(adjustl(format_version))
!if( format_version .ne. "Ver1_4_2 for iRIC" ) stop "This RRI model requires RRI_Input_Format_Ver1_4_2"
!write(*,*)

!Run Type
    call cg_iric_read_integer(cgns_f, "run_type", run_type, ier)

!--------------------------------------------------
!外部ファイル　→　格子属性として設定
!--------------------------------------------------
!read(1,'(a)') rainfile
    call cg_iric_read_string(cgns_f, "rainfile", rainfile, ier)

!read(1,'(a)') demfile
    call cg_iric_read_string(cgns_f, "demfile", demfile, ier)

!read(1,'(a)') accfile
    call cg_iric_read_string(cgns_f, "accfile", accfile, ier)

!read(1,'(a)') dirfile
    call cg_iric_read_string(cgns_f, "dirfile", dirfile, ier)

    write (*, '("rainfile : ", a)') trim(adjustl(rainfile))
    write (*, '("demfile : ", a)') trim(adjustl(demfile))
    write (*, '("accfile : ", a)') trim(adjustl(accfile))
    write (*, '("dirfile : ", a)') trim(adjustl(dirfile))
!
!read(1,*)
    write (*, *)

!--------------------------------------------------
!基本オプション
!--------------------------------------------------
!read(1,*) utm
call cg_iric_read_integer(cgns_f, "utm", utm, ier)
!    utm = 0

!read(1,*) eight_dir
call cg_iric_read_integer(cgns_f, "eight_dir", eight_dir, ier)
   ! eight_dir = 1
!    eight_dir = 0 !using 4 direction 20240808
    write (*, '("utm : ", i5)') utm
    write (*, '("eight_dir : ", i5)') eight_dir
    write (*, *)

!--------------------------------------------------
!時間条件
!--------------------------------------------------
!read(1,*) lasth
    call cg_iric_read_integer(cgns_f, "lasth", lasth, ier)

!read(1,*) dt
    call cg_iric_read_integer(cgns_f, "dt", dt, ier)

!read(1,*) dt_riv
    call cg_iric_read_integer(cgns_f, "dt_riv", dt_riv, ier)

!read(1,*) outnum
    call cg_iric_read_integer(cgns_f, "outnum", outnum, ier)

    write (*, '("lasth : ", i8)') lasth
    write (*, '("dt : ", i12)') dt
    write (*, '("dt_riv : ", i8)') dt_riv
    write (*, '("outnum : ", i8)') outnum
    write (*, *)

!--------------------------------------------------
!降雨条件　格子属性として設定
!--------------------------------------------------
!read(1,*) xllcorner_rain
!read(1,*) yllcorner_rain
!read(1,*) cellsize_rain_x, cellsize_rain_y

    call cg_iric_read_real(cgns_f, "xllcorner_rain", xllcorner_rain, ier)
    call cg_iric_read_real(cgns_f, "yllcorner_rain", yllcorner_rain, ier)
    call cg_iric_read_real(cgns_f, "cellsize_rain_x", cellsize_rain_x, ier)
    call cg_iric_read_real(cgns_f, "cellsize_rain_y", cellsize_rain_y, ier)

    write (*, '("xllcorner_rain : ", f15.5)') xllcorner_rain
    write (*, '("yllcorner_rain : ", f15.5)') yllcorner_rain
    write (*, '("cellsize_rain_x : ", f15.5, "  cellsize_rain_y : ", f15.5)') cellsize_rain_x, cellsize_rain_y

    write (*, *)

!--------------------------------------------------
!slopeパラメータ　→　格子セル属性　複合条件として設定
!--------------------------------------------------
!read(1,*) num_of_landuse
!call cg_iric_read_complex_count(cgns_f, "landuse_c", num_of_landuse, ier)
!RRI for iRIC version
!Number of Land use type is 3 as fixed.
    num_of_landuse = 5
!
    allocate (dif(num_of_landuse))
    allocate (ns_slope(num_of_landuse), soildepth(num_of_landuse))
    allocate (gammaa(num_of_landuse))
!
!dif
    call cg_iric_read_integer(cgns_f, "dif_1", dif(1), ier)
    call cg_iric_read_integer(cgns_f, "dif_2", dif(2), ier)
    call cg_iric_read_integer(cgns_f, "dif_3", dif(3), ier)
    call cg_iric_read_integer(cgns_f, "dif_4", dif(4), ier)
    call cg_iric_read_integer(cgns_f, "dif_4", dif(5), ier)

!ns_slope
    call cg_iric_read_real(cgns_f, "ns_slope_1", ns_slope(1), ier)
    call cg_iric_read_real(cgns_f, "ns_slope_2", ns_slope(2), ier)
    call cg_iric_read_real(cgns_f, "ns_slope_3", ns_slope(3), ier)
    call cg_iric_read_real(cgns_f, "ns_slope_4", ns_slope(4), ier)
    call cg_iric_read_real(cgns_f, "ns_slope_5", ns_slope(5), ier)

!solid depth
    call cg_iric_read_real(cgns_f, "soildepth_1", soildepth(1), ier)
    call cg_iric_read_real(cgns_f, "soildepth_2", soildepth(2), ier)
    call cg_iric_read_real(cgns_f, "soildepth_3", soildepth(3), ier)
    call cg_iric_read_real(cgns_f, "soildepth_4", soildepth(4), ier)
    call cg_iric_read_real(cgns_f, "soildepth_5", soildepth(5), ier)

!Porosity
    call cg_iric_read_real(cgns_f, "gammaa_1", gammaa(1), ier)
    call cg_iric_read_real(cgns_f, "gammaa_2", gammaa(2), ier)
    call cg_iric_read_real(cgns_f, "gammaa_3", gammaa(3), ier)
    call cg_iric_read_real(cgns_f, "gammaa_4", gammaa(4), ier)
    call cg_iric_read_real(cgns_f, "gammaa_5", gammaa(5), ier)

!

    write (*, '("num_of_landuse : ", i5)') num_of_landuse
    write (*, '("dif : ", 100i5)') (dif(i), i=1, num_of_landuse)
    write (*, '("ns_slope : ", 100f12.3)') (ns_slope(i), i=1, num_of_landuse)
    write (*, '("soildepth : ", 100f12.3)') (soildepth(i), i=1, num_of_landuse)
    write (*, '("gammaa : ", 100f12.3)') (gammaa(i), i=1, num_of_landuse)
!
!
!read(1,*)
    write (*, *)
!
    allocate (ksv(num_of_landuse), faif(num_of_landuse))
!
!ksv
    call cg_iric_read_real(cgns_f, "ksv_1", ksv(1), ier)
    call cg_iric_read_real(cgns_f, "ksv_2", ksv(2), ier)
    call cg_iric_read_real(cgns_f, "ksv_3", ksv(3), ier)
    call cg_iric_read_real(cgns_f, "ksv_4", ksv(4), ier)
    call cg_iric_read_real(cgns_f, "ksv_5", ksv(5), ier)

!Sf
    call cg_iric_read_real(cgns_f, "faif_1", faif(1), ier)
    call cg_iric_read_real(cgns_f, "faif_2", faif(2), ier)
    call cg_iric_read_real(cgns_f, "faif_3", faif(3), ier)
    call cg_iric_read_real(cgns_f, "faif_4", faif(4), ier)
    call cg_iric_read_real(cgns_f, "faif_5", faif(5), ier)

!
    write (*, '("ksv : ", 100e12.3)') (ksv(i), i=1, num_of_landuse)
    write (*, '("faif : ", 100f12.3)') (faif(i), i=1, num_of_landuse)
!
!read(1,*)
    write (*, *)
!
    allocate (ka(num_of_landuse), gammam(num_of_landuse), beta(num_of_landuse))
!
!ka
    call cg_iric_read_real(cgns_f, "ka_1", ka(1), ier)
    call cg_iric_read_real(cgns_f, "ka_2", ka(2), ier)
    call cg_iric_read_real(cgns_f, "ka_3", ka(3), ier)
    call cg_iric_read_real(cgns_f, "ka_4", ka(4), ier)
    call cg_iric_read_real(cgns_f, "ka_5", ka(5), ier)

!Unsat. porosity
    call cg_iric_read_real(cgns_f, "gammam_1", gammam(1), ier)
    call cg_iric_read_real(cgns_f, "gammam_2", gammam(2), ier)
    call cg_iric_read_real(cgns_f, "gammam_3", gammam(3), ier)
    call cg_iric_read_real(cgns_f, "gammam_4", gammam(4), ier)
    call cg_iric_read_real(cgns_f, "gammam_5", gammam(5), ier)

!beta
    call cg_iric_read_real(cgns_f, "beta_1", beta(1), ier)
    call cg_iric_read_real(cgns_f, "beta_2", beta(2), ier)
    call cg_iric_read_real(cgns_f, "beta_3", beta(3), ier)
    call cg_iric_read_real(cgns_f, "beta_4", beta(4), ier)
    call cg_iric_read_real(cgns_f, "beta_5", beta(5), ier)

    write (*, '("ka : ", 100e12.3)') (ka(i), i=1, num_of_landuse)
    write (*, '("gammam : ", 100f12.3)') (gammam(i), i=1, num_of_landuse)
    write (*, '("beta : ", 100f12.3)') (beta(i), i=1, num_of_landuse)

    do i = 1, num_of_landuse
        if (gammam(i) .gt. gammaa(i)) stop "gammag must be smaller than gammaa"
    end do
!
!read(1,*)
    write (*, *)
!
    allocate (ksg(num_of_landuse), gammag(num_of_landuse), kg0(num_of_landuse), fpg(num_of_landuse), rgl(num_of_landuse))
!
! ksg = 0.0 means that this model has not been implemented into iRIC version. ****
!
    ksg = 0.0; gammag = 0.0; kg0 = 0.0; fpg = 0.0; rgl = 0.0

!
    write (*, '("ksg : ", 100e12.3)') (ksg(i), i=1, num_of_landuse)
    write (*, '("gammag : ", 100f12.3)') (gammag(i), i=1, num_of_landuse)
    write (*, '("kg0 : ", 100e12.3)') (kg0(i), i=1, num_of_landuse)
    write (*, '("fpg : ", 100f12.3)') (fpg(i), i=1, num_of_landuse)
    write (*, '("rgl : ", 100e12.3)') (rgl(i), i=1, num_of_landuse)
!read(1,*)
    write (*, *)

!--------------------------------------------------
!河道パラメータ　→　wc,ws,dc,dsで指定する
!--------------------------------------------------
!read(1,*) ns_river
    call cg_iric_read_real(cgns_f, "ns_river", ns_river, ier)
    write (*, '("ns_river : ", f12.3)') ns_river
    write (*, *)

!read(1,*) riv_thresh
    call cg_iric_read_integer(cgns_f, "riv_thresh", riv_thresh, ier)
    write (*, '("riv_thresh : ", i7)') riv_thresh
    write (*, *)

!read(1,*) width_param_c
!read(1,*) width_param_s
    call cg_iric_read_real(cgns_f, "width_param_c", width_param_c, ier)
    call cg_iric_read_real(cgns_f, "width_param_s", width_param_s, ier)

!read(1,*) depth_param_c
!read(1,*) depth_param_s
    call cg_iric_read_real(cgns_f, "depth_param_c", depth_param_c, ier)
    call cg_iric_read_real(cgns_f, "depth_param_s", depth_param_s, ier)

!read(1,*) height_param
!read(1,*) height_limit_param
    call cg_iric_read_real(cgns_f, "height_param", height_param, ier)
    call cg_iric_read_integer(cgns_f, "height_limit_param", height_limit_param, ier)

!read(1,*) rivfile_switch
!call cg_iric_read_integer(cgns_f, "rivfile_switch", rivfile_switch, ier)
!rivfile_switchは利用しない
    rivfile_switch = 0

!read(1,'(a)') widthfile
!read(1,'(a)') depthfile
!read(1,'(a)') heightfile
    widthfile = ''; depthfile = ''; heightfile = ''
    if (rivfile_switch == 1) then
        call cg_iric_read_string(cgns_f, "widthfile", widthfile, ier)
        call cg_iric_read_string(cgns_f, "depthfile", depthfile, ier)
        call cg_iric_read_string(cgns_f, "heightfile", heightfile, ier)
    end if

    if (rivfile_switch .eq. 0) then
        write (*, '("width_param_c : ", f12.2)') width_param_c
        write (*, '("width_param_s : ", f12.2)') width_param_s
        write (*, '("depth_param_c : ", f12.2)') depth_param_c
        write (*, '("depth_param_s : ", f12.2)') depth_param_s
        write (*, '("height_param : ", f12.2)') height_param
        write (*, '("height_limit_param : ", i10)') height_limit_param
    else
        write (*, '("widthfile : ", a)') trim(adjustl(widthfile))
        write (*, '("depthfile : ", a)') trim(adjustl(depthfile))
        write (*, '("heightfile : ", a)') trim(adjustl(heightfile))
    end if

!read(1,*)
    write (*, *)

!--------------------------------------------------
!hotstart用　初期条件
!--------------------------------------------------
!read(1,*) init_slo_switch, init_riv_switch, init_gw_switch, init_gampt_ff_switch
    call cg_iric_read_integer(cgns_f, "init_slo_switch", init_slo_switch, ier)
    call cg_iric_read_integer(cgns_f, "init_riv_switch", init_riv_switch, ier)
    call cg_iric_read_integer(cgns_f, "init_gw_switch", init_gw_switch, ier)
    call cg_iric_read_integer(cgns_f, "init_gampt_ff_switch", init_gampt_ff_switch, ier)

!read(1,"(a)") initfile_slo
!read(1,'(a)') initfile_riv
!read(1,'(a)') initfile_gw
!read(1,'(a)') initfile_gampt_ff

    initfile_slo = ''; initfile_riv = ''; initfile_gw = ''; initfile_gampt_ff = ''
    if (init_slo_switch == 1) call cg_iric_read_string(cgns_f, "initfile_slo", initfile_slo, ier)
    if (init_riv_switch == 1) call cg_iric_read_string(cgns_f, "initfile_riv", initfile_riv, ier)
    if (init_gw_switch == 1) call cg_iric_read_string(cgns_f, "initfile_gw", initfile_gw, ier)
    if (init_gampt_ff_switch == 1) call cg_iric_read_string(cgns_f, "initfile_gampt_ff", initfile_gampt_ff, ier)

    if (init_slo_switch .ne. 0) write (*, '("initfile_slo : ", a)') trim(adjustl(initfile_slo))
    if (init_riv_switch .ne. 0) write (*, '("initfile_riv : ", a)') trim(adjustl(initfile_riv))
    if (init_gw_switch .ne. 0) write (*, '("initfile_gw : ", a)') trim(adjustl(initfile_gw))
    if (init_gampt_ff_switch .ne. 0) write (*, '("initfile_gampt_ff : ", a)') trim(adjustl(initfile_gampt_ff))

!read(1,*)
    write (*, *)

!--------------------------------------------------
!hs, hr境界条件　→　境界条件設定に実装
!--------------------------------------------------
    bound_slo_wlev_switch = 0; bound_riv_wlev_switch = 0
    call cg_iric_read_bc_count(cgns_f, "bound_hs", num_of_bound_point)
    if (num_of_bound_point > 0) bound_slo_wlev_switch = 1

    call cg_iric_read_bc_count(cgns_f, "bound_hr", num_of_bound_point)
    if (num_of_bound_point > 0) bound_riv_wlev_switch = 1

!read(1,*) bound_slo_wlev_switch, bound_riv_wlev_switch
!read(1,'(a)') boundfile_slo_wlev
!read(1,'(a)') boundfile_riv_wlev
!if(bound_slo_wlev_switch.ne.0) write(*,'("boundfile_slo_wlev : ", a)') trim(adjustl(boundfile_slo_wlev))
!if(bound_riv_wlev_switch.ne.0) write(*,'("boundfile_riv_wlev : ", a)') trim(adjustl(boundfile_riv_wlev))
!read(1,*)
!write(*,*)

!--------------------------------------------------
!qs, qr境界条件　→　境界条件設定に実装
!--------------------------------------------------
    bound_slo_disc_switch = 0; bound_riv_disc_switch = 0
    call cg_iric_read_bc_count(cgns_f, "bound_qs", num_of_bound_point)
    if (num_of_bound_point > 0) bound_slo_disc_switch = 1

    call cg_iric_read_bc_count(cgns_f, "bound_qr", num_of_bound_point)
    if (num_of_bound_point > 0) bound_riv_disc_switch = 1

!read(1,*) bound_slo_disc_switch, bound_riv_disc_switch
!read(1,'(a)') boundfile_slo_disc
!read(1,'(a)') boundfile_riv_disc
!if(bound_slo_disc_switch.ne.0) write(*,'("boundfile_slo_disc : ", a)') trim(adjustl(boundfile_slo_disc))
!if(bound_riv_disc_switch.ne.0) write(*,'("boundfile_riv_disc : ", a)') trim(adjustl(boundfile_riv_disc))
!read(1,*)
!write(*,*)

!--------------------------------------------------
!土地利用条件　→　格子属性として与える
!--------------------------------------------------
!read(1,*) land_switch
!read(1,'(a)') landfile
!if(land_switch.eq.1) write(*,'("landfile : ", a)') trim(adjustl(landfile))
!
!read(1,*)
!write(*,*)

!--------------------------------------------------
!Dam条件　→　境界条件設定に実装
!--------------------------------------------------
    dam_switch = 0
    call cg_iric_read_bc_count(cgns_f, "dam", dam_num)
    if (dam_num > 0) dam_switch = 1
    !for dam
    if(dam_switch>0)then
    write(*,'("Number of dam : ", i7)') dam_num
    endif
!read(1,*) dam_switch
!read(1,'(a)') damfile
!if(dam_switch.eq.1) write(*,'("damfile : ", a)') trim(adjustl(damfile))
!read(1,*)
!write(*,*)

!--------------------------------------------------
!div条件　→　境界条件設定に実装
!--------------------------------------------------
    div_id_max = 0
    call cg_iric_read_bc_count(cgns_f, "div", div_id_max)
    if (div_id_max > 0) div_switch = 1
!read(1,*) div_switch
!call cg_iric_read_integer(cgns_f, "div_switch", div_switch, ier)
!read(1,'(a)') divfile
!if(div_switch.eq.1) write(*,'("divfile : ", a)') trim(adjustl(divfile))
!read(1,*)
!write(*,"(a)") "div condition has not implemented yet."
!write(*,*)

!--------------------------------------------------
!蒸発条件　→
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
!write(*,"(a)") "evp condition can be set as the grid attribute."
!write(*,*)

!--------------------------------------------------
!河道断面データ　→　実装しない
!--------------------------------------------------
!read(1,*) sec_length_switch
!read(1,'(a)') sec_length_file        !これはセルごとに河道長を設定するファイル

!if(sec_length_switch.eq.1) write(*,'("sec_length : ", a)') trim(adjustl(sec_length_file))
!read(1,*)
!write(*,*)
!read(1,*) sec_switch
!read(1,'(a)') sec_map_file        !これはセルごとにsecfileのidを指定するファイル
!read(1,'(a)') sec_file                !これはsection形状のファイル　なければwc,ws,dc,dsで指定された値となる
!if(sec_switch.eq.1) write(*,'("sec_map_file : ", a)') trim(adjustl(sec_map_file))
!if(sec_switch.eq.1) write(*,'("sec_file : ", a)') trim(adjustl(sec_file))
!read(1,*)
!write(*,"(a)") "Cross section condition has not implemented yet."
!write(*,*)

!--------------------------------------------------
!??? 　→　要確認
!--------------------------------------------------
!read(1,*) emb_switch
!read(1,'(a)') embrfile
!read(1,'(a)') embbfile
!if(emb_switch.eq.1) write(*,'("embrfile : ", a)') trim(adjustl(embrfile))
!if(emb_switch.eq.1) write(*,'("embbfile : ", a)') trim(adjustl(embbfile))
!write(*,*)

!--------------------------------------------------
!計算結果出力・ファイル
!--------------------------------------------------
!read(1,*) outswitch_hs, outswitch_hr, outswitch_hg, outswitch_qr, outswitch_qu, outswitch_qv, &
!          outswitch_gu, outswitch_gv, outswitch_gampt_ff, outswitch_storage

!read(1,'(a)') outfile_hs
    outfile_hs = ''
    call cg_iric_read_integer(cgns_f, "outswitch_hs", outswitch_hs, ier)
    if (outswitch_hs == 1) call cg_iric_read_string(cgns_f, "outfile_hs", outfile_hs, ier)

!read(1,'(a)') outfile_hr
    outfile_hr = ''
    call cg_iric_read_integer(cgns_f, "outswitch_hr", outswitch_hr, ier)
    if (outswitch_hr == 1) call cg_iric_read_string(cgns_f, "outfile_hr", outfile_hr, ier)

!read(1,'(a)') outfile_hg
    outfile_hg = ''
    call cg_iric_read_integer(cgns_f, "outswitch_hg", outswitch_hg, ier)
    if (outswitch_hg == 1) call cg_iric_read_string(cgns_f, "outfile_hg", outfile_hg, ier)

!read(1,'(a)') outfile_qr
    outfile_qr = ''
    call cg_iric_read_integer(cgns_f, "outswitch_qr", outswitch_qr, ier)
    if (outswitch_qr == 1) call cg_iric_read_string(cgns_f, "outfile_qr", outfile_qr, ier)

!read(1,'(a)') outfile_qu
    outfile_qu = ''
    call cg_iric_read_integer(cgns_f, "outswitch_qu", outswitch_qu, ier)
    if (outswitch_qu == 1) call cg_iric_read_string(cgns_f, "outfile_qu", outfile_qu, ier)

!read(1,'(a)') outfile_qv
    outfile_qv = ''
    call cg_iric_read_integer(cgns_f, "outswitch_qv", outswitch_qv, ier)
    if (outswitch_qv == 1) call cg_iric_read_string(cgns_f, "outfile_qv", outfile_qv, ier)

!read(1,'(a)') outfile_gu
    outfile_gu = ''
    call cg_iric_read_integer(cgns_f, "outswitch_gu", outswitch_gu, ier)
    if (outswitch_gu == 1) call cg_iric_read_string(cgns_f, "outfile_gu", outfile_gu, ier)

!read(1,'(a)') outfile_gv
    outfile_gv = ''
    call cg_iric_read_integer(cgns_f, "outswitch_gv", outswitch_gv, ier)
    if (outswitch_gv == 1) call cg_iric_read_string(cgns_f, "outfile_gv", outfile_gv, ier)

!read(1,'(a)') outfile_gampt_ff
    outfile_gampt_ff = ''
    call cg_iric_read_integer(cgns_f, "outswitch_gampt_ff", outswitch_gampt_ff, ier)
    if (outswitch_gampt_ff == 1) call cg_iric_read_string(cgns_f, "outfile_gampt_ff", outfile_gampt_ff, ier)

!read(1,'(a)') outfile_storage
    outfile_storage = ''
    call cg_iric_read_integer(cgns_f, "outswitch_storage", outswitch_storage, ier)
    if (outswitch_storage == 1) call cg_iric_read_string(cgns_f, "outfile_storage", outfile_storage, ier)

    if (outswitch_hs .ne. 0) write (*, '("outfile_hs : ", a)') trim(adjustl(outfile_hs))
    if (outswitch_hr .ne. 0) write (*, '("outfile_hr : ", a)') trim(adjustl(outfile_hr))
    if (outswitch_hg .ne. 0) write (*, '("outfile_hg : ", a)') trim(adjustl(outfile_hg))
    if (outswitch_qr .ne. 0) write (*, '("outfile_qr : ", a)') trim(adjustl(outfile_qr))
    if (outswitch_qu .ne. 0) write (*, '("outfile_qu : ", a)') trim(adjustl(outfile_qu))
    if (outswitch_qv .ne. 0) write (*, '("outfile_qv : ", a)') trim(adjustl(outfile_qv))
    if (outswitch_gu .ne. 0) write (*, '("outfile_gu : ", a)') trim(adjustl(outfile_gu))
    if (outswitch_gv .ne. 0) write (*, '("outfile_gv : ", a)') trim(adjustl(outfile_gv))
    if (outswitch_gampt_ff .ne. 0) write (*, '("outfile_gampt_ff : ", a)') trim(adjustl(outfile_gampt_ff))
    if (outswitch_storage .ne. 0) write (*, '("outfile_storage : ", a)') trim(adjustl(outfile_storage))

!read(1,*)
    write (*, *)

 !------added for RSR model 20240724
    allocate( kgv(num_of_landuse), tg(num_of_landuse))

    call cg_iric_read_integer(cgns_f, "sed_switch", sed_switch, ier)
    write (*, '("sed_switch : ", i7)') sed_switch
    write (*, *)
    !In case sediment computation, read the following
    if (sed_switch >= 1)then
        call cg_iric_read_real(cgns_f, "t_beddeform_start", t_beddeform_start, ier)
        write (*, '("t_beddeform_start(hour): ", f12.1)') t_beddeform_start
        t_beddeform_start = t_beddeform_start*3600.
        call cg_iric_read_real(cgns_f, "s", s, ier)
        write (*, '("s : ", f12.2)') s
        call cg_iric_read_real(cgns_f, "grav", grav, ier)
        write (*, '("grav : ", f12.2)') grav
        call cg_iric_read_real(cgns_f, "t_Crit", t_Crit, ier)
        write (*, '("t_Crit : ", f12.3)') t_Crit
        call cg_iric_read_real(cgns_f, "kin_visc", kin_visc, ier)
        write (*, '("kin_visc : ", f12.3)') kin_visc
        call cg_iric_read_real(cgns_f, "karmans", karmans, ier)
        write (*, '("karmans : ", f12.1)') karmans
        call cg_iric_read_real(cgns_f, "lambda", lambda, ier)
        write (*, '("lambda : ", f12.1)') lambda
        write (*, *)
        call cg_iric_read_integer(cgns_f, "sed_type_switch", sed_type_switch, ier)
        write (*, '("sed_type_switch : ", i7)') sed_type_switch
        call cg_iric_read_real(cgns_f, "ds_river", ds_river, ier)
        write (*, '("ds_river : ", f12.5)') ds_river
        call cg_iric_read_integer(cgns_f, "iidt", iidt, ier)
        write (*, '("iidt : ", i7)') iidt
        call cg_iric_read_integer(cgns_f, "isedeq", isedeq, ier)
        write (*, '("isedeq : ", i7)') isedeq
        call cg_iric_read_integer(cgns_f, "isuseq", isuseq, ier)
        write (*, '("isuseq : ", i7)') isuseq
        write (*, *)
        !Regulations
!        call cg_iric_read_integer(cgns_f, "max_acc_0th_riv", max_acc_0th_riv, ier)
!        write (*, '("max_acc_0th_riv : ", i7)') max_acc_0th_riv
        call cg_iric_read_integer(cgns_f, "min_num_cell_link", min_num_cell_link, ier)
        write (*, '("min_num_cell_link : ", i7)') min_num_cell_link
        call cg_iric_read_real(cgns_f, "perosion", perosion, ier)
        write (*, '("perosion : ", f12.5)') perosion
        call cg_iric_read_real(cgns_f, "min_slope", min_slope, ier)
        write (*, '("min_slope : ", f12.5)') min_slope
        call cg_iric_read_real(cgns_f, "max_slope", max_slope, ier)
        write (*, '("max_slope : ", f12.2)') max_slope
        call cg_iric_read_real(cgns_f, "min_hr", min_hr, ier)
        write (*, '("min_hr : ", f12.5)') min_hr
        call cg_iric_read_real(cgns_f, "alpha_ss1", alpha_ss1, ier)
        write (*, '("alpha_ss1 : ", f12.5)') alpha_ss1
        call cg_iric_read_real(cgns_f, "alpha_ss2", alpha_ss2, ier)
        write (*, '("alpha_ss2 : ", f12.5)') alpha_ss2
        call cg_iric_read_real(cgns_f, "thresh_ss", thresh_ss, ier)
        write (*, '("thresh_ss : ", f12.5)') thresh_ss
        call cg_iric_read_real(cgns_f, "hr0", hr0, ier)
        write (*, '("Initial flow depth : ", f12.5)') hr0
        call cg_iric_read_real(cgns_f, "hs0", hs0, ier)
        write (*, '("Initial water depth of slope cells: ", f12.5)') hs0
        !call cg_iric_read_integer(cgns_f, "cut_overdepo_switch", cut_overdepo_switch, ier) !tentatively turned of 20250122
         !write (*, '("Enforcing sediment overflow: ", i7)') cut_overdepo_switch
        !Non-uniform
        if(sed_type_switch ==2)then
            call cg_iric_read_real(cgns_f, "Em", Em, ier)
            write (*, '("Em : ", f12.5)') Em
            call cg_iric_read_integer(cgns_f, "no_of_layers", no_of_layers, ier)
            write (*, '("no_of_layers : ", i7)') no_of_layers
            call cg_iric_read_integer(cgns_f, "Nl", Nl, ier)
            write (*, '("Nl : ", i7)') Nl
            call cg_iric_read_real(cgns_f, "Emc", Emc, ier)
            write (*, '("Emc : ", f12.5)') Emc
        end if

        call cg_iric_read_integer(cgns_f, "detail_console", detail_console, ier)

    end if

 !------Landslide and debris flow 20240724
    call cg_iric_read_integer(cgns_f, "debris_switch", debris_switch, ier)
    write (*, '("debris_switch : ", i7)') debris_switch
    if (debris_switch >= 1)then
        call cg_iric_read_real(cgns_f, "cohe", cohe, ier)
        write (*, '("cohesion: ", f12.2)') cohe
        call cg_iric_read_real(cgns_f, "pwc", pwc, ier)
        write (*, '("pwc: ", f12.2)') pwc
        call cg_iric_read_real(cgns_f, "pf", pf, ier)
        write (*, '("pf: ", f12.2)') pf
        call cg_iric_read_real(cgns_f, "b_mp", b_mp, ier)
        write (*, '("b_mp: ", f12.2)') b_mp
        call cg_iric_read_real(cgns_f, "d_mp_ini", d_mp_ini, ier)
        write (*, '("d_mp_ini: ", f12.2)') d_mp_ini
        call cg_iric_read_real(cgns_f, "L_rain_ini", L_rain_ini, ier)
        write (*, '("L_rain_ini: ", f12.2)') L_rain_ini
        call cg_iric_read_integer(cgns_f, "debris_end_switch", debris_end_switch, ier)
        write (*, '("debris_end_switch: ", i7)') debris_end_switch
        if(debris_end_switch==1)then
            call cg_iric_read_real(cgns_f, "T_bebris_off", T_bebris_off, ier)
            write (*, '("T_bebris_off(hour): ", f12.2)') T_bebris_off
            T_bebris_off = T_bebris_off*3600.
        else
            T_bebris_off = 1000000000.
        end if            
    end if

 !------Slope erosion 20240724
 !---modified for slope erosion
    call cg_iric_read_integer(cgns_f, "slo_sedi_cal_switch", slo_sedi_cal_switch, ier)
    write (*, '("slo_sedi_cal_switch : ", i7)') slo_sedi_cal_switch
    if(slo_sedi_cal_switch>0) then
        call cg_iric_read_integer(cgns_f, "slope_ero_switch", slope_ero_switch, ier)
        write (*, '("slope_ero_switch : ", i7)') slope_ero_switch
        if(slope_ero_switch>0)then
        nm_cell =9 !tentative
        allocate (B_gully_r(nm_cell+1), D_gully(nm_cell+1))
        do i = 0, nm_cell
        write(cm,'(i1)') i
        gullyB_label = 'B_gully_r_'//trim(cm)  
        gullyD_label = 'D_gully_'//trim(cm)  
        call cg_iric_read_real(cgns_f, gullyB_label, B_gully_r(i+1), ier)      
        call cg_iric_read_real(cgns_f, gullyD_label, D_gully(i+1), ier)     
        enddo
        endif

        call cg_iric_read_real(cgns_f, "modirate_to_slo_dt", modirate_to_slo_dt, ier)
        write (*, '("modirate_to_slo_dt: ", f12.2)') modirate_to_slo_dt
        dt_slo_sed = dt/modirate_to_slo_dt 
        call cg_iric_read_real(cgns_f, "surflowdepth", surflowdepth, ier)
        write (*, '("surflowdepth: ", f12.2)') surflowdepth
    else
    slope_ero_switch= 0    
    end if

 !------Driftwood 20240724
    call cg_iric_read_integer(cgns_f, "j_drf", j_drf, ier)
    write (*, '("j_drf : ", i7)') j_drf
    call cg_iric_read_real(cgns_f, "Wood_density", Wood_density, ier)
    write (*, '("Wood_density : ", f12.5)') Wood_density

 !------Advanced setting eps 20240724
        call cg_iric_read_real(cgns_f, "eps", eps, ier)
        write (*, '("eps: ", f12.5)') eps
        call cg_iric_read_real(cgns_f, "ddt_min_riv", ddt_min_riv, ier)
        write (*, '("ddt_min_riv: ", f12.5)') ddt_min_riv
        call cg_iric_read_real(cgns_f, "ddt_min_slo", ddt_min_slo, ier)
        write (*, '("ddt_min_slo: ", f12.5)') ddt_min_slo

 !------Advanced output setting 20240724
    outfile_test = ''
    call cg_iric_read_integer(cgns_f, "outswitch_test", outswitch_test, ier)
    if (outswitch_test == 1)then
        call cg_iric_read_string(cgns_f, "outfile_test", outfile_test, ier)
        write(*,'("outfile_test : ", a)') trim(adjustl(outfile_test))
    end if
!for slope erosion
    outfile_slope = ''
    call cg_iric_read_integer(cgns_f, "outswitch_slope", outswitch_slope, ier)
    if (outswitch_slope == 1)then
        call cg_iric_read_string(cgns_f, "outfile_slope", outfile_slope, ier)
        write(*,'("outfile_slope : ", a)') trim(adjustl(outfile_slope))
    end if

!---tentatively turn off the channel width expansion and over deposition cutting 20250122
    riv_wid_expan_switch = 0 
    cut_overdepo_switch = 0
!--------------------------------------------------
!iRICの基本機能で対応
!--------------------------------------------------
!read(1,*) hydro_switch
!read(1,'(a)') location_file
!!location_file = trim(rri_dir)//location_file(3:len(location_file))
!
!if(hydro_switch .eq. 1) write(*,'("location_file : ", a)') trim(adjustl(location_file))
!
    write (*, *)

!close(1)
    call cg_iric_close(cgns_f, ier)

! Parameter Check
    do i = 1, num_of_landuse
        if (ksv(i) .gt. 0.d0 .and. ka(i) .gt. 0.d0) &
            stop "Error: both ksv and ka are non-zero."
        if (gammam(i) .gt. gammaa(i)) &
            stop "Error: gammam must be smaller than gammaa."
    end do

! Set da, dm and infilt_limit
    allocate (da(num_of_landuse), dm(num_of_landuse), infilt_limit(num_of_landuse))
    da(:) = 0.d0
    dm(:) = 0.d0
    infilt_limit(:) = 0.d0
    do i = 1, num_of_landuse
        if (soildepth(i) .gt. 0.d0 .and. ksv(i) .gt. 0.d0) infilt_limit(i) = soildepth(i)*gammaa(i)
        if (soildepth(i) .gt. 0.d0 .and. ka(i) .gt. 0.d0) da(i) = soildepth(i)*gammaa(i)
        if (soildepth(i) .gt. 0.d0 .and. ka(i) .gt. 0.d0 .and. gammam(i) .gt. 0.d0) &
            dm(i) = soildepth(i)*gammam(i)
    end do

! if ksg(i) = 0.d0 -> no gw calculation
    gw_switch = 0
    do i = 1, num_of_landuse
        if (ksg(i) .gt. 0.d0) then
            gw_switch = 1
        else
            gammag(i) = 0.d0
            kg0(i) = 0.d0
            fpg(i) = 0.d0
            rgl(i) = 0.d0
        end if
    end do

end subroutine RRI_Read

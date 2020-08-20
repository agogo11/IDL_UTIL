PRO yzracmes, OutPath,$
  dateArray, overpassArray,$
  RadarPath = radarpath,$
  CEILPath = ceilpath,$
  METPath = metpath,$
  sondepath = sondepath,$
  Plotfig = plotfig,$
  plotsonde = plotsonde, $
  timeinterval = timeinterval,$
  readradar = readradar,$
  readmet = readmet,$
  readceil = readceil,$
  readsonde  =readsonde, $
  metcontent_list = metcontent_list, $
  radarcontent_list = radarcontent_list, $
  outputradar = outputradar,$
  outputmet = outputmet,$
  outputceil = outputceil,$
  outputsonde = outputsonde,$
  ind_st = ind_st,$
  ind_nd = ind_nd

  ;set default
  setdefaultvalue, outputmet, 1
  setdefaultvalue, outputceil, 1
  setdefaultvalue, outputradar, 1
  setdefaultvalue, outputsonde, 1

  setdefaultvalue, readmet, 1
  setdefaultvalue, readceil, 1
  setdefaultvalue, readradar, 1
  setdefaultvalue, readsonde, 1

  setdefaultvalue, metcontent_list, ['temperature_ambient','rh_ambient','pressure_ambient']
  setdefaultvalue, radarcontent_list, ['height','reflectivity']
  setdefaultvalue, timeinterval, 3.
  setdefaultvalue, ind_st, 0
  setdefaultvalue, ind_nd, N_ELEMENTS(dateArray) - 1

  ;---consts
  Lv = 2501. ; kj/kg
  cpd = 1005. ; J/(kg.K)

  ;set output name via specifying the suffix
  IF 60.*timeinterval GE 100. THEN BEGIN
    suffix = STRING(60.*timeinterval, format = '(I3)')
  ENDIF ELSE BEGIN
    suffix = '0' + STRING(60.*timeinterval, format = '(I2)')
  ENDELSE

  ;for conventional nomination
  date_cases = dateArray
  overpass_cases = overpassArray

  FOR icases = ind_st, ind_nd DO BEGIN

    PRINT, STRTRIM(icases + 1,1) + '/' + STRTRIM(ind_nd + 1,1) + ':' + $
      date_cases[icases], overpass_cases[icases]


    ;#1 Determine date and time of the cases----------------------------------------------
    Year = LONG(STRMID(date_cases[icases],0,4))
    Month = LONG(STRMID(date_cases[icases],4,2))
    Day = LONG(STRMID(date_cases[icases],6,2))
    Hour = LONG(STRMID(Overpass_cases[icases],0,2))
    Minute = LONG(STRMID(Overpass_cases[icases],2,2))

    ;Hour from the start of the day for each case
    Hour_SAT_Day = Hour + Minute/60.0

    ;Hour from the start of the month for each case
    Hour_SAT_Month = 24*(Day - 1) + Hour + Minute/60.0

    ;cal the +-1 day of the case
    JULDAY_cases=julday(month,day,year,hour,minute,0)
    caldat, JULDAY_cases - 1, MONTH1, DAY1, YEAR1, HOUR1, MIN1, SEC1
    caldat, JULDAY_cases + 1, MONTH3, DAY3, YEAR3, HOUR3, MIN3, SEC3
    date = yzcreatedate(month1,day1, year1, month3, day3, year3)


    ;#2 1st-round determination of the availability of data------------------------------
    ;we assume datasets are available
    statusmet = 1
    statusceil = 1
    statusradar = 1
    statussonde = 1

    ;radar
    RADARFile0=FILE_SEARCH(RADARPath+'*'+ date[0] +'*', count = count0)
    RADARFile1=FILE_SEARCH(RADARPath+'*'+ date[1] +'*', count = count1)
    RADARFile2=FILE_SEARCH(RADARPath+'*'+ date[2] +'*', count = count2)
    RADARFile = [RADARFile0,RADARFile1,RADARFile2]

    ind = WHERE(STRLEN(RADARFile) EQ 0, COMPLEMENT = ind_c, count)
    RADARFile = RADARFile[ind_c]
    nRadarFile = count1

    ;Met
    METFile0=FILE_SEARCH(metPath+'*'+date[0]+'*',count = count0)
    METFile1=FILE_SEARCH(metPath+'*'+date[1]+'*',count = count1)
    METFile2=FILE_SEARCH(metPath+'*'+date[2]+'*',count = count2)
    METFile = [METFile0,METFile1,METFile2]
    nMETFile = count1

    ;ceilometer
    CEILFile0=FILE_SEARCH(ceilPath+'*'+date[0]+'*',count = count0)
    CEILFile1=FILE_SEARCH(ceilPath+'*'+date[1]+'*',count = count1)
    CEILFile2=FILE_SEARCH(ceilPath+'*'+date[2]+'*',count = count2)
    CEILFile = [CEILFile0,CEILFile1,CEILFile2]
    nCEILFile = count1

    ;radiosonde
    SondeFile1=FILE_SEARCH(SondePath+'*' + date[0] + '*',count = count0)
    SondeFile2=FILE_SEARCH(SondePath+'*' + date[1] + '*',count = count1)
    SondeFile3=FILE_SEARCH(SondePath+'*' + date[2] + '*',count = count2)

    nSONDEFile = count0 + count1 + count2
    SONDEFile = [Sondefile1,Sondefile2,Sondefile3]

    ;remove unfounded data
    ind = WHERE(STRLEN(SONDEFile) EQ 0, COMPLEMENT = ind_c, count)
    SONDEFile = SONDEFile[ind_c]
    nSONDEFile = N_ELEMENTS(SONDEFile)

    julday_array = MAKE_ARRAY(nSONDEFile, /double)
    FOR isondefile = 0, nSONDEFile - 1 DO BEGIN
      a = strsplit(SONDEFile[isondefile],'.', /extract)
      b = strsplit(a[0],'\', /extract)
      sondeproductname = b[N_ELEMENTS(b) - 1]

      Yeartemp = LONG(STRMID(a[2],0,4))
      Monthtemp = LONG(STRMID(a[2],4,2))
      Daytemp = LONG(STRMID(a[2],6,2))
      Hourtemp = LONG(STRMID(a[3],0,2))
      Minutetemp = LONG(STRMID(a[3],2,2))
      julday_array[isondefile] = julday(monthtemp,daytemp,yeartemp,hourtemp,minutetemp,0) - julday(1,1,2000,0,0,0)
    ENDFOR
    JULDAY_cases = JULDAY_cases - julday(1,1,2000,0,0,0)

    ;find the nearest one
    so = SORT(ABS(julday_array- JULDAY_cases))
    sondefile_use = SONDEFile[so[0]]
    Dif_sonde_cases = (julday_array[so[0]]- JULDAY_cases)*24.d

    ;updata the status of data
    statusmet = nMETFile*readmet < 1
    statusceil = nCEILFile*readceil < 1
    statusradar = nRadarFile*readradar < 1
    statussonde = (Dif_sonde_cases LE 1.)*readsonde < 1
  
    ;#3 Start reading data-----------------------------------------------------------
    ;MET
    IF statusmet EQ 1 THEN BEGIN
      ind = WHERE(STRLEN(METFile) EQ 0, count,COMPLEMENT = indc)
      METFile = METFile[indc]

      ;get the data product name
      a = strsplit(METFile[0],'.', /extract)
      b = strsplit(a[0],'\', /extract)
      metproductname = b[N_ELEMENTS(b) - 1]

      ;read met
      nMETFile = N_ELEMENTS(METFile)
      FOR iFile = 0, nMETFile - 1 DO BEGIN

        status = yzreadarmcdf(METFile[iFile],bufferMET, $
          content_list = metcontent_list,/check_missingvalue)

        tag_names0 = TAG_NAMES(bufferMET)
        ind_time = WHERE(tag_names0 EQ 'TIME')
        ind_day = WHERE(tag_names0 EQ 'DAY')
        ind_temp = WHERE(tag_names0 EQ STRUPCASE(metcontent_list[0]))
        ind_rh = WHERE(tag_names0 EQ STRUPCASE(metcontent_list[1]))
        ind_pres = WHERE(tag_names0 EQ STRUPCASE(metcontent_list[2]))

        time_tmp = bufferMET.(ind_time[0]) + 24.*(bufferMET.day[0]- 1.)
        Ta_tmp = buffermet.(ind_temp[0])
        rh_tmp = buffermet.(ind_rh[0])
        pres_tmp = buffermet.(ind_pres[0])

        ni = N_ELEMENTS(time_tmp)
        LCL_tmp = FLTARR(ni)
        FOR i = 0, ni - 1 DO BEGIN
          LCL_tmp[i] = romplcl(100.*pres_tmp[i],Ta_tmp[i] + 273.15,rh_tmp[i]/100.)/1000.
        ENDFOR

        IF iFile EQ 0 THEN BEGIN
          time_day_met = time_tmp
          LCL_day = LCL_tmp
          Ta_day = Ta_tmp
          rh_day = rh_tmp
          pres_day = pres_tmp
        ENDIF ELSE BEGIN
          time_day_met = [time_day_met,time_tmp]
          LCL_day = [LCL_day, LCL_tmp]
          Ta_day = [Ta_day, Ta_tmp]
          rh_day = [rh_day, rh_tmp]
          pres_day = [pres_day, pres_tmp]
        ENDELSE

      ENDFOR

    ENDIF

    ;Ceilometer data
    IF statusceil EQ 1 THEN BEGIN
      ind = WHERE(STRLEN(CeilFile) EQ 0, count,COMPLEMENT = indc)
      CeilFile = CeilFile[indc]

      ;get the data product name
      a = strsplit(CeilFile[0],'.', /extract)
      b = strsplit(a[0],'\', /extract)
      ceilproductname = b[N_ELEMENTS(b) - 1]

      nCeilFile = N_ELEMENTS(CeilFile)
      FOR iFile = 0, nCeilFile - 1 DO BEGIN
        status = yzreadarmcdf(CEILFile[iFile],bufferCEIL, $
          content_list = ['first_cbh', 'second_cbh'],/check_missingvalue)

        time_tmp = bufferCEIL.time + 24.*(bufferCEIL.day[0] - 1.)
        first_cbh_tmp = bufferCEIL.first_cbh
        second_cbh_tmp = bufferCEIL.second_cbh

        IF iFile EQ 0 THEN BEGIN
          first_cbh_day = first_cbh_tmp
          second_cbh_day = second_cbh_tmp
          time_day_ceil = time_tmp
        ENDIF ELSE BEGIN
          first_cbh_day = [first_cbh_day, first_cbh_tmp]
          second_cbh_day = [second_cbh_day, second_cbh_tmp]
          time_day_ceil = [time_day_ceil,time_tmp]
        ENDELSE
      ENDFOR
    ENDIF

    ;radiosonde data
    IF statussonde EQ 1 THEN BEGIN
      yzreadarmcdf_old, sondefile_use, tdry, dp, alt, pres, wspd, deg, rh,$
        VNames = ['tdry','dp','alt','pres','wspd','deg','rh'],$
        MissingValues = [-9999.,-9999.,1,-9999.,-9999.,-9999.,-9999.]

      ;PRINT, N_ELEMENTS(alt[WHERE(alt LT 3000.)])
      IF N_ELEMENTS(alt[WHERE(alt LT 3000.)]) LT 70 THEN CONTINUE

      IF alt[0] GT 2000. THEN alt = alt - alt[0]

      alt = alt/1000.0
      tdry = tdry + 273.15

      E = 6.112*EXP(17.67*(tdry-273.15)/(tdry-273.15 + 243.5))
      E = E*rh/100.

      q = 1000.*0.622*E/(pres - E)
      ptdry = tdry*(1000./pres)^(0.286)
      eptdry = (tdry + (Lv/cpd)*(q/1000.))*(1000./pres)^(0.286)

      z_inv_temp = yztempinversion_sonde(alt, tdry, q, pres)

      zi_inv_max = z_inv_temp.zi_inv_max
      zi_inv_max0 = z_inv_temp.zi_inv_max0
      t_inv_max = z_inv_temp.t_inv_max

      IF zi_inv_max0 NE 999. THEN zi_inv_max = zi_inv_max0
      PRINT, ' '
    ENDIF

    ;radar data
    IF statusradar EQ 1 THEN BEGIN

      RADARFile_use = RADARFile

      ;get the data product name
      a = strsplit(RADARFile_use[0],'.', /extract)
      b = strsplit(a[0],'\', /extract)
      radarproductname = b[N_ELEMENTS(b) - 1]

      FOR ifile = 0, N_ELEMENTS(RADARFile_use) - 1 DO BEGIN
        status = yzreadarmcdf(RADARFile_use[ifile],bufferRadar, $
          content_list = radarcontent_list)

        tag_names0 = TAG_NAMES(bufferRadar)
        ind_time = WHERE(tag_names0 EQ 'TIME')
        ind_day = WHERE(tag_names0 EQ 'DAY')
        ind_hgt = WHERE(tag_names0 EQ STRUPCASE(radarcontent_list[0]))
        ind_ref = WHERE(tag_names0 EQ STRUPCASE(radarcontent_list[1]))

        range = (bufferRadar.(ind_hgt[0]))/1000.

        ref_tmp = bufferRadar.(ind_ref[0])
        time_tmp = bufferRadar.(ind_time[0]) + 24.*(bufferRadar.day[0] - 1.)

        IF ifile EQ 0 THEN BEGIN
          ref_day = ref_tmp
          time_day_radar = time_tmp
        ENDIF ELSE BEGIN
          time_day_radar = [time_day_radar, time_tmp]
          ref_day = [[ref_day],[ref_tmp]]
        ENDELSE
      ENDFOR
      
      ;double check the availability of radar data
      ind = where(time_day_radar gt Hour_sat_month - timeinterval/2. and time_day_radar le Hour_sat_month + timeinterval/2., count)
      If count eq 0 then statusradar = 0
    ENDIF

    ;#4 determine key quantities and output-----------------------------------------------------
    ;radar data
    IF statusradar EQ 1 THEN BEGIN
      OutputRef = selectdata_radar(Ref_day, range, time_day_radar, range[0], 3.5,$
        Hour_sat_month - timeinterval/2., Hour_sat_month + timeinterval/2.,$
        H_Output = H_plot, T_Output = T_output)

      ;#4.1 column maximum reflectivity
      Maxref = MAX(OutputRef, dimension = 1)
      botref = OutputRef[0,*]

      ;output
      OPENW,LUN,OutPath + radarproductname + '_Maxref_'$
        + date_cases[icases] + '_' + overpass_cases[icases] + '.txt',/GET_LUN
      PRINTF,LUN, '*time','max_ref [dBZ]','bottom_ref [dBZ]',$
        FORMAT= '(3A15)'
      FOR i = 0, N_ELEMENTS(Maxref) - 1 DO BEGIN
        PRINTF,LUN, T_output[i] - 24.*(day - 1), Maxref[i], botref[i],$
          FORMAT = '(3F15.4)'
      ENDFOR
      FREE_LUN, LUN

      ;#4.2 reflectivity vertical profile
      range_vp = 0.2 + 0.1*INDGEN(34)
      nrange_vp = N_ELEMENTS(range_vp)

      ref_groups = [-999., -35, -25, -15, 0, 999.]
      ref_groups_name = ['<-35 dBZ','-35~-25 dBZ','-25~-15 dBZ','-15~0 dBZ','>0 dBZ']
      ngroups = N_ELEMENTS(ref_groups) - 1

      ref_vp = FLTARR(ngroups, nrange_vp)
      temp_vp = FLTARR(nrange_vp)

      FOR irange = 0, nrange_vp - 2 DO BEGIN
        ;reflectivity
        ind = WHERE(H_plot GT range_vp[irange] AND H_plot LE range_vp[irange + 1], count)
        binnedtemp = VALUE_LOCATE(ref_groups, OutputRef[ind, *])

        temp = cghistogram(binnedtemp,location = ind, /NAN)
        ref_vp[ind,irange] = cghistogram(binnedtemp,/NAN)

        ;temperature from radiosonde
        ind = WHERE(alt GT range_vp[irange] AND alt LE range_vp[irange + 1], count)
        temp_vp[irange] = mean(tdry[ind],/NAN)
      ENDFOR

      ;output
      OPENW,LUN,OutPath + radarproductname + '_ref_vp_'$
        + date_cases[icases] + '_' + overpass_cases[icases] + '.txt',/GET_LUN
      PRINTF,LUN, '*altitude','<-35 dBZ','-35~-25 dBZ','-25~-15 dBZ','-15~0 dBZ','>0 dBZ', 'temp [K]',$
        FORMAT= '(7A15)'

      FOR i = 0, nrange_vp - 2 DO BEGIN
        PRINTF,LUN, range_vp[i],ref_vp[0, i],ref_vp[1, i],ref_vp[2, i],ref_vp[3, i],ref_vp[4, i],temp_vp[i],$
          FORMAT = '(7F15.4)'
      ENDFOR
      FREE_LUN, LUN

      PRINT, ' '
    ENDIF

    ;radiosonde data
    IF statussonde*readsonde NE 0 THEN BEGIN

      alt_smooth = SMOOTH(alt, 3)
      pres_smooth = SMOOTH(pres, 3)
      tdry_smooth = SMOOTH(tdry, 3)
      rh_smooth = SMOOTH(rh, 3)
      q_smooth = SMOOTH(q, 3)
      eptdry_smooth = SMOOTH(eptdry, 3)

      ;determine Z_0.7zinv_LCL and Z_150m_LCL in Eq.(5) of McGibbon and Bretherton, 2017, JAMES
      near = MIN(ABS(alt_smooth - 0.7*zi_inv_max), ind)
      LCL_07zinv = romplcl(100.*pres_smooth[ind],tdry_smooth[ind],rh_smooth[ind]/100.)/1000. + 0.7*zi_inv_max

      near = MIN(ABS(alt_smooth - 0.15), ind)
      LCL_150m = romplcl(100.*pres_smooth[ind],tdry_smooth[ind],rh_smooth[ind]/100.)/1000.

      ;determine q_bot and q_top in section 2.2 of Jones et al., 2011, ACP
      alt_temp = alt_smooth[WHERE(alt_smooth LT zi_inv_max)]
      q_temp = q_smooth[WHERE(alt_smooth LT zi_inv_max)]
      eptdry_temp = eptdry_smooth[WHERE(alt_smooth LT zi_inv_max)]

      temp = cgpercentiles(alt_temp, Percentiles=[0., 0.25, 0.75, 1.])

      z_100zinv = temp[3]
      z_75zinv = temp[2]
      z_25zinv = temp[1]
      z_0zinv = temp[0]

      qbot = mean(q_temp[WHERE(alt_temp GT z_0zinv AND alt_temp LT z_25zinv)])
      qtop = mean(q_temp[WHERE(alt_temp GT z_75zinv AND alt_temp LT z_100zinv)])

      eptdrybot = mean(eptdry_temp[WHERE(alt_temp GT z_0zinv AND alt_temp LT z_25zinv)])
      eptdrytop = mean(eptdry_temp[WHERE(alt_temp GT z_75zinv AND alt_temp LT z_100zinv)])

      PRINT, ' '

      nvar = 20
      str = {varname:' ', varvalue:-999., varunit:' '}
      sonde_output = REPLICATE(str, nvar)

      sonde_output[0].varname = 'dif_sonde_cases' & sonde_output[0].varvalue = dif_sonde_cases & sonde_output[0].varunit = 'hr'
      sonde_output[1].varname = 'zi_inv_max' & sonde_output[1].varvalue = zi_inv_max & sonde_output[1].varunit = 'km'
      sonde_output[2].varname = 'LCL_07zinv' & sonde_output[2].varvalue = LCL_07zinv & sonde_output[2].varunit = 'km'
      sonde_output[3].varname = 'LCL_150m' & sonde_output[3].varvalue = LCL_150m & sonde_output[3].varunit = 'km'
      sonde_output[4].varname = 'qbot' & sonde_output[4].varvalue = qbot & sonde_output[4].varunit = 'g/kg'
      sonde_output[5].varname = 'qtop' & sonde_output[5].varvalue = qtop & sonde_output[5].varunit = 'g/kg'
      sonde_output[6].varname = 'eptdrybot' & sonde_output[6].varvalue = eptdrybot & sonde_output[6].varunit = 'K'
      sonde_output[7].varname = 'eptdrytop' & sonde_output[7].varvalue = eptdrytop & sonde_output[7].varunit = 'K'
      sonde_output[8].varname = 'zi_inv_max0' & sonde_output[8].varvalue = zi_inv_max0 & sonde_output[8].varunit = 'km'
      sonde_output[9].varname = 't_inv_max' & sonde_output[9].varvalue = t_inv_max & sonde_output[9].varunit = 'deg'
      

      sonde_output = sonde_output[WHERE(sonde_output.varvalue NE -999.)]
      nvar = N_ELEMENTS(sonde_output)

      ;output the results
      OPENW,LUN,OutPath + sondeproductname + '_'  + date_cases[icases] + '_' + $
        overpass_cases[icases] + '.txt',/GET_LUN
      PRINTF,LUN,'*Date:', date_cases[icases], FORMAT= '(A-10, A-10)'
      PRINTF,LUN,'*Satellite Overpass:', overpass_cases[icases], FORMAT= '(A-25, A-10)'
      FOR ivar = 0, nvar - 1 DO BEGIN
        PRINTF,LUN,''
        PRINTF,LUN,sonde_output[ivar].varname + ' = ', sonde_output[ivar].varvalue, sonde_output[ivar].varunit, FORMAT= '(A-30, F-9.3, A6)'
      ENDFOR
      FREE_LUN, LUN


      IF KEYWORD_SET(plotsonde) THEN BEGIN
        cgps_open, OutPath + sondeproductname + '_'$
          + date_cases[icases] + '_' + overpass_cases[icases] + '.eps',/CMYK,$
          /nomatch,Font = 1,/Times, xsize = 12/2.54, ysize = 10/2.54, fontsize =14

        pos = cglayout([1,1], OXMargin=[5, 5], OYMargin=[5, 5])
        unit = !D.Y_CH_SIZE / FLOAT(!D.Y_SIZE)
        mycharsize = 1.3
        
        ind = where(alt lt 3.5)
        ;        y = alt[0:350]
        ;        x1 = ptdry[0: 350]
        ;        x2 = q[0:350]

        y = alt[ind]
        x1 = ptdry[ind]
        x2 = q[ind]

        myxrange1 = [MIN(x1), MAX(x1)]
        myxrange2 = [MIN(x2), MAX(x2)]
        myyrange = [0., 3.5]

        cgplot, 0, 0, color = 'black', position = pos[*,0], XStyle=9, /nodata, $
          /noerase, ytitle = 'Height [km]', xtitle = 'Potential Temp [K]',$
          charsize = mycharsize, xrange = myxrange1, yrange = myyrange, xtickinterval = 3

        ;bottom 25% od z_inv
        nx = 100.
        x = myxrange1[0] + ((myxrange1[1] - myxrange1[0])/nx)*INDGEN(nx + 1)

        low_error = z_0zinv + FLTARR(nx + 1)
        high_error = z_25zinv + FLTARR(nx + 1)

        cgcolorfill, [x, reverse(x), x[0]], $
          [high_error, reverse(low_error), high_error[0]], $
          Color='light grey'

        ;upper 25% of z_inv
        nx = 100.
        x = myxrange1[0] + ((myxrange1[1] - myxrange1[0])/nx)*INDGEN(nx + 1)

        low_error = z_75zinv + FLTARR(nx + 1)
        high_error = z_100zinv + FLTARR(nx + 1)

        cgcolorfill, [x, reverse(x), x[0]], $
          [high_error, reverse(low_error), high_error[0]], $
          Color='light grey'

        ;plot pot temp
        cgoplot, x1, y

        cgplots, myxrange1, [zi_inv_max,zi_inv_max], color = 'orange', thick = 3., linestyle = 2
        cgplots, myxrange1, [zi_inv_max0,zi_inv_max0], color = 'yellow', thick = 3., linestyle = 2
        cgplots, myxrange1, [LCL_07zinv,LCL_07zinv], color = 'black', thick = 3., linestyle = 2
        cgplots, myxrange1, [LCL_150m,LCL_150m], color = 'red', thick = 3., linestyle = 2

        text_lg = ['zi_inv','LCL_07zinv','LCL_150m']
        colors_lg = ['orange','black','red']
        linestyle_lg = [2,2,2]
        location_lg = [myxrange1[1] - 7., 0.9*myyrange[1]]

        cglegend, title = text_lg, linestyle = linestyle_lg, color = colors_lg, Location = location_lg, /data,$
          vspace = 0.8, charsize = 0.7*mycharsize, tcolors = colors_lg

        ;plot q
        cgaxis, XAxis=1, XRange=myxrange2, xtitle = 'Water vapor mixing ratio [g/kg]' ,/Save, color = 'red', charsize = mycharsize
        cgoplot, x2, y, color = 'red'

        cgps_close, /png
        PRINT, ' '
      ENDIF

    ENDIF

    IF statusceil EQ 1 THEN BEGIN
      Ind = WHERE(time_day_ceil GE Hour_SAT_Month - timeinterval/2. AND time_day_ceil LE Hour_SAT_Month + timeinterval/2., count)
      count_tot = count

      time_cbh_use = time_day_ceil[ind]
      first_cbh_use = first_cbh_day[ind]
      second_cbh_use = second_cbh_day[ind]

      ;output
      OPENW,LUN,OutPath + ceilproductname + '_histo_'$
        + date_cases[icases] + '_' + overpass_cases[icases] + '.txt',/GET_LUN
      PRINTF,LUN, '*time','first_cbh [km]','second_cbh [km]',$
        FORMAT= '(3A15)'
      FOR i = 0, N_ELEMENTS(time_cbh_use) - 1 DO BEGIN
        PRINTF,LUN, time_cbh_use[i] - 24.*(day - 1), first_cbh_use[i], second_cbh_use[i],$
          FORMAT = '(3F15.4)'
      ENDFOR
      FREE_LUN, LUN

      ind = WHERE(first_cbh_use GE 100 AND first_cbh_use LE 3500 AND FINITE(first_cbh_use) EQ 1, count)
      count_1stcld = count

      temp = first_cbh_use[ind]/1000.
      first_cbh_CF = 100.*count_1stcld/count_tot
      first_cbh_median = count NE 0 ? MEDIAN(temp[ind], /even): -999.
      first_cbh_mean = count NE 0 ? mean(temp[ind],/NAN): -999.
      first_cbh_Std = count NE 0 ? stddev(temp[ind],/NAN): -999.
      first_cbh_Skew = count NE 0 ? skewness(temp[ind],/NAN): -999.

      a = cgpercentiles(temp, Percentiles=[0., 0.2])
      first_cbh_bot20mean = mean(temp[WHERE(temp GE a[0] AND temp LE a[1])])


      ind = WHERE(second_cbh_use GE 100 AND second_cbh_use LE 3500 AND FINITE(second_cbh_use) EQ 1, count)
      count_2ndcld = count

      temp = second_cbh_use[ind]/1000.
      second_cbh_CF = 100.*count_2ndcld/count_tot
      second_cbh_median = count NE 0 ? MEDIAN(temp[ind], /even): -999.
      second_cbh_mean = count NE 0 ? mean(temp[ind],/NAN): -999.
      second_cbh_Std = count NE 0 ? stddev(temp[ind],/NAN): -999.
      second_cbh_Skew = count NE 0 ? skewness(temp[ind],/NAN): -999.

      nvar = 20
      str = {varname:' ', varvalue:-999., varunit:' '}
      temp_output = REPLICATE(str, nvar)

      temp_output[0].varname = 'first_cbh_CF' & temp_output[0].varvalue = first_cbh_CF & temp_output[0].varunit = '%'
      temp_output[1].varname = 'first_cbh_median' & temp_output[1].varvalue = first_cbh_median & temp_output[1].varunit = 'km'
      temp_output[2].varname = 'first_cbh_mean' & temp_output[2].varvalue = first_cbh_mean & temp_output[2].varunit = 'km'
      temp_output[3].varname = 'first_cbh_Std' & temp_output[3].varvalue = first_cbh_Std & temp_output[3].varunit = 'km'
      temp_output[4].varname = 'first_cbh_Skew' & temp_output[4].varvalue = first_cbh_Skew & temp_output[4].varunit = ' '
      temp_output[5].varname = 'first_cbh_bot20mean' & temp_output[5].varvalue = first_cbh_bot20mean & temp_output[5].varunit = 'km'
      temp_output[6].varname = 'second_cbh_CF' & temp_output[6].varvalue = second_cbh_CF & temp_output[6].varunit = '%'
      temp_output[7].varname = 'second_cbh_median' & temp_output[7].varvalue = second_cbh_median & temp_output[7].varunit = 'km'
      temp_output[8].varname = 'second_cbh_mean' & temp_output[8].varvalue = second_cbh_mean & temp_output[8].varunit = 'km'
      temp_output[9].varname = 'second_cbh_Std' & temp_output[9].varvalue = second_cbh_Std & temp_output[9].varunit = 'km'
      temp_output[10].varname = 'second_cbh_Skew' & temp_output[10].varvalue = second_cbh_Skew & temp_output[10].varunit = ' '

      ceil_output = TEMPORARY(temp_output[WHERE(temp_output.varvalue NE -999.)])
      nvar = N_ELEMENTS(ceil_output)

      ;output the results
      OPENW,LUN,OutPath + ceilproductname + '_'  + date_cases[icases] + '_' + $
        overpass_cases[icases] + '.txt',/GET_LUN
      PRINTF,LUN,'*Date:', date_cases[icases], FORMAT= '(A-10, A-10)'
      PRINTF,LUN,'*Satellite Overpass:', overpass_cases[icases], FORMAT= '(A-25, A-10)'
      FOR ivar = 0, nvar - 1 DO BEGIN
        PRINTF,LUN,''
        PRINTF,LUN,ceil_output[ivar].varname + ' = ', ceil_output[ivar].varvalue, ceil_output[ivar].varunit, FORMAT= '(A-30, F-9.3, A6)'
      ENDFOR
      FREE_LUN, LUN

    ENDIF

    IF statusmet EQ 1 THEN BEGIN
      Ind = WHERE(time_day_met GT Hour_SAT_Month - timeinterval/2. AND time_day_met LT Hour_SAT_Month + timeinterval/2., count)
      LCL_mean_met = mean(LCL_day[ind],/NAN)
      Ta_mean_met = mean(Ta_day[ind],/NAN)
      rh_mean_met = mean(rh_day[ind],/NAN)
      pres_mean_met = mean(pres_day[ind],/NAN)

      nvar = 20

      str = {varname:' ', varvalue:-999., varunit:' '}
      met_output = REPLICATE(str, nvar)

      met_output[0].varname = 'LCL_mean_met' & met_output[0].varvalue = LCL_mean_met & met_output[0].varunit = 'km'
      met_output[1].varname = 'Ta_mean_met' & met_output[1].varvalue = Ta_mean_met & met_output[1].varunit = 'C'
      met_output[2].varname = 'rh_mean_met' & met_output[2].varvalue = rh_mean_met & met_output[2].varunit = '%'
      met_output[3].varname = 'pres_mean_met' & met_output[3].varvalue = pres_mean_met & met_output[3].varunit = 'hPa'

      met_output = met_output[WHERE(met_output.varvalue NE -999.)]

      nvar = N_ELEMENTS(met_output)

      ;output the results
      OPENW,LUN,OutPath + metproductname + '_'  + date_cases[icases] + '_' + $
        overpass_cases[icases] + '.txt',/GET_LUN
      PRINTF,LUN,'*Date:', date_cases[icases], FORMAT= '(A-10, A-10)'
      PRINTF,LUN,'*Satellite Overpass:', overpass_cases[icases], FORMAT= '(A-25, A-10)'
      FOR ivar = 0, nvar - 1 DO BEGIN
        PRINTF,LUN,''
        PRINTF,LUN,met_output[ivar].varname + ' = ', met_output[ivar].varvalue, met_output[ivar].varunit, FORMAT= '(A-30, F-9.3, A6)'
      ENDFOR
      FREE_LUN, LUN
    ENDIF

    ;#5 Plot the data
    IF KEYWORD_SET(plotfig) AND statusradar NE 0 THEN BEGIN

      ;pre-process the data for plot
      time_day_radar_plot = time_day_radar - 24.*(day - 1)
      time_day_met_plot = time_day_met - 24.*(day - 1)
      time_day_ceil_plot = time_day_ceil - 24.*(day - 1)

      time_st = Hour_SAT_Day - timeinterval/2.
      time_nd = Hour_SAT_Day + timeinterval/2.

      cgps_open, OutPath + radarproductname + '_'+ date_cases[icases] + '_' + $
        overpass_cases[icases] + '.eps',/CMYK,$
        /nomatch,Font = 1,/Times, xsize = 18/2.54, ysize = 14/2.54, fontsize =14

      xunit = !D.X_CH_SIZE / FLOAT(!D.X_SIZE)
      yunit = 0.5*(!D.Y_CH_SIZE / FLOAT(!D.Y_SIZE))
      Mycharsize = 0.7

      pos = cglayout([1, 2], OXMargin=[3, 5], OYMargin=[3, 2], YGap=3, XGap = 0.5)

      p0 = pos[*,0]
      p1 = pos[*,1]

      unit = !D.Y_CH_SIZE / FLOAT(!D.Y_SIZE)

      ref_plot = selectdata_radar(ref_day, range, time_day_radar_plot, $
        range[0], 4.,$
        time_st, time_nd,$
        H_Output = H_plot, T_Output = T_plot)

      ;fig1
      plotradar, ref_plot, T_plot, H_plot, PlotVariable = 'ref',$
        Position = pos[*,0], MultiMargin = [0,2,0,4],$
        Mycharsize = Mycharsize, MyTcharsize = 0.5, title = 'W-band ARM Cloud Radar image', xtitle = 'Time [h]', TTitle = 'Reflectivity [dBZ]',$
        TPosition = [pos[2,0] + xunit,pos[1,0],pos[2,0] + 2.*xunit ,pos[3,0]]

      cgoplot, time_day_ceil_plot, first_cbh_day/1000., psym = 16, symsize = 0.12, color = 'black'
      cgoplot, time_day_ceil_plot, second_cbh_day/1000., psym = 16, symsize = 0.12, color = 'green'
      cgoplot, time_day_met_plot, lcl_day, psym = 16, symsize = 0.1, color = 'red'

      ;cgoplot, [Hour_SAT_Day,Hour_SAT_Day],[0., 3.], linestyle = 2, color = 'red'
      cgoplot, [Hour_SAT_Day - 0.25,Hour_SAT_Day - 0.25],[0., 4.], linestyle = 2, color = 'grey'
      cgoplot, [Hour_SAT_Day + 0.25,Hour_SAT_Day + 0.25],[0., 4.], linestyle = 2, color = 'grey'

      cgoplot, [Hour_SAT_Day - 0.25, Hour_SAT_Day + 0.25], [zi_inv_max0, zi_inv_max0], color = 'yellow', linestyle = 2
      cgoplot, [Hour_SAT_Day - 0.25, Hour_SAT_Day + 0.25], [zi_inv_max, zi_inv_max], color = 'orange', linestyle = 2
      
      ind = where(alt lt 3., count)
      ptdry_plot = ptdry[ind]
      q_plot = q[ind]
      scale_factor1 = (0.5/(MAX(ptdry_plot) - MIN(ptdry_plot)))
      scale_factor2 = (0.5/(MAX(q_plot) - MIN(q_plot)))
      add_offset = Hour_SAT_Day - 0.5/2.

      x1 = scale_factor1*(ptdry_plot - MIN(ptdry_plot)) + add_offset
      x2 = scale_factor2*(q_plot - MIN(q_plot)) + add_offset
      y = alt[ind]

      cgoplot, x1, y, color = 'yellow'
      cgoplot, x2, y, color = 'green'

      ;plot some calculated quantities
      xloc1 = p0[0]
      xloc2 = p0[0] + (p0[2] - p0[0])/3.
      xloc3 = p0[0] + 2.*(p0[2] - p0[0])/3.
      yinc = 1.5*yunit

      yloc = p0[1] - 7.*yunit
      FOR ivar = 0, N_ELEMENTS(ceil_output) - 1 DO BEGIN
        cgtext, xloc1, yloc, ceil_output[ivar].varname + ' = ' + STRTRIM(STRING(ceil_output[ivar].varvalue, format = '(f0.2)'),1) + ' ' + ceil_output[ivar].varunit, /normal, charsize = mycharsize
        yloc = yloc - yinc
      ENDFOR

      yloc = p0[1] - 7.*yunit
      FOR ivar = 0, N_ELEMENTS(sonde_output) - 1 DO BEGIN
        cgtext, xloc2, yloc, sonde_output[ivar].varname + ' = ' + STRTRIM(STRING(sonde_output[ivar].varvalue, format = '(f0.2)'),1) + ' ' + sonde_output[ivar].varunit, /normal, charsize = mycharsize
        yloc = yloc - yinc
      ENDFOR

      yloc = p0[1] - 7.*yunit
      FOR ivar = 0, N_ELEMENTS(met_output) - 1 DO BEGIN
        cgtext, xloc3, yloc, met_output[ivar].varname + ' = ' + STRTRIM(STRING(met_output[ivar].varvalue, format = '(f0.2)'),1) + ' ' + met_output[ivar].varunit, /normal, charsize = mycharsize
        yloc = yloc - yinc
      ENDFOR

      ;  DEVICE,/times, /BOLD
      ;  cgtext,pos[0,0] - 0.075, pos[3,0]+0.02, 'a',/normal, charsize = 0.7, FONT = 0
      ;  cgtext,pos[0,1] - 0.075, pos[3,1]+0.02, 'b',/normal, charsize = 0.7, FONT = 0
      ;  cgtext,pos[0,2] - 0.075, pos[3,2]+0.02, 'c',/normal, charsize = 0.7, FONT = 0
      ;  cgtext,pos[0,3] - 0.075, pos[3,3]+0.02, 'd',/normal, charsize = 0.7, FONT = 0

      cgps_close, /PNG, density = 1000
    ENDIF

    PRINT, ' '

  ENDFOR

END
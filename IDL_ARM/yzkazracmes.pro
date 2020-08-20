PRO YZKAZRACMES, OutPath,$
  dateArray, overpassArray,$
  kazrPath = kazrpath,$
  CEILPath = ceilpath,$
  METPath = metpath,$
  sondepath = sondepath,$
  Plotfig = plotfig,$
  plotsonde = plotsonde, $
  timeinterval = timeinterval,$
  readkazr = readkazr,$
  readmet = readmet,$
  readceil = readceil,$
  readsonde  =readsonde, $
  metcontent_list = metcontent_list, $
  outputkazr = outputkzar,$
  outputmet = outputmet,$
  outputceil = outputceil,$
  outputsonde = outputsonde,$
  ind_st = ind_st,$
  ind_nd = ind_nd

  ;set default
  SETDEFAULTVALUE, outputmet, 1
  SETDEFAULTVALUE, outputceil, 1
  SETDEFAULTVALUE, outputdl, 1
  SETDEFAULTVALUE, outputsonde, 1

  SETDEFAULTVALUE, readmet, 1
  SETDEFAULTVALUE, readceil, 1
  SETDEFAULTVALUE, readdl, 1
  SETDEFAULTVALUE, readsonde, 1

  SETDEFAULTVALUE, metcontent_list, ['temp_mean','rh_mean','atmos_pressure']
  ;SETDEFAULTVALUE, kazrcontent_list, ['height','reflectivity_best_estimate']
  SETDEFAULTVALUE, kazrcontent_list, ['range','reflectivity_copol']
  SETDEFAULTVALUE, timeinterval, 2.
  SETDEFAULTVALUE, ind_st, 0
  SETDEFAULTVALUE, ind_nd, N_ELEMENTS(dateArray) - 1

  ;---consts
  Lv = 2501. ; kj/kg
  cpd = 1005. ; J/(kg.K)
  g = 9.80665d ; gravity const.
  rho_ref       = 1.d ;  density of air (kg/m3)
  p_ref = 1000.;  reference  pressure, hPa

  ref_thre = -42 ;threshold for determining cloudy pixels in DL data
  SNR_thre = 0.005 ; threshold for removing noisy data in DL
  gap_thre = 2 ; gap (s) for identifying single clouds

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

    ;    ;Hour from the start of the month for each case
    ;    Hour_SAT_Month = 24*(Day - 1) + Hour + Minute/60.0

    ;cal the +-1 day of the case
    JULDAY_cases=JULDAY(month,day,year,hour,minute,0)
    CALDAT, JULDAY_cases - 1, MONTH1, DAY1, YEAR1, HOUR1, MIN1, SEC1
    CALDAT, JULDAY_cases + 1, MONTH3, DAY3, YEAR3, HOUR3, MIN3, SEC3
    date = YZCREATEDATE(month1,day1, year1, month3, day3, year3)

    ;#2 1st-round determination of the availability of data------------------------------
    ;we assume datasets are available
    statusmet = 1
    statusceil = 1
    statuskazr = 1
    statussonde = 1

    ;kazr
    kazrFile=FILE_SEARCH(kazrPath+'*kazr*'+ date[1] +'*', count = count1)
    nkazrFile = count1

    ;Met
    METFile=FILE_SEARCH(metPath+'*'+date[1]+'*',count = count1)
    nMETFile = count1

    ;ceilometer
    CEILFile=FILE_SEARCH(ceilPath+'*'+date[1]+'*',count = count1)
    nCEILFile = count1

    ;radiosonde
    SondeFile1=FILE_SEARCH(SondePath+'*' + date[0] + '*',count = count0)
    SondeFile2=FILE_SEARCH(SondePath+'*' + date[1] + '*',count = count1)
    SondeFile3=FILE_SEARCH(SondePath+'*' + date[2] + '*',count = count2)

    nSONDEFile = count0 + count1 + count2
    SONDEFile = [Sondefile1,Sondefile2,Sondefile3]

    IF nSONDEFile EQ 0 THEN BEGIN
      statussonde = 0
      Dif_sonde_cases = 999.
    ENDIF ELSE BEGIN
      ;remove unfounded data
      ind = WHERE(STRLEN(SONDEFile) EQ 0, COMPLEMENT = ind_c, count)
      SONDEFile = SONDEFile[ind_c]
      nSONDEFile = N_ELEMENTS(SONDEFile)

      julday_array = MAKE_ARRAY(nSONDEFile, /double)
      FOR isondefile = 0, nSONDEFile - 1 DO BEGIN
        a = STRSPLIT(SONDEFile[isondefile],'.', /extract)
        b = STRSPLIT(a[0],'\', /extract)
        sondeproductname = b[N_ELEMENTS(b) - 1]

        Yeartemp = LONG(STRMID(a[2],0,4))
        Monthtemp = LONG(STRMID(a[2],4,2))
        Daytemp = LONG(STRMID(a[2],6,2))
        Hourtemp = LONG(STRMID(a[3],0,2))
        Minutetemp = LONG(STRMID(a[3],2,2))
        julday_array[isondefile] = JULDAY(monthtemp,daytemp,yeartemp,hourtemp,minutetemp,0) - JULDAY(1,1,2000,0,0,0)
      ENDFOR
      JULDAY_cases = JULDAY_cases - JULDAY(1,1,2000,0,0,0)

      ;find the nearest one
      so = SORT(ABS(julday_array- JULDAY_cases))
      sondefile_use = SONDEFile[so[0]]
      Dif_sonde_cases = (julday_array[so[0]]- JULDAY_cases)*24.d
    ENDELSE

    ;updata the status of data
    statusmet = nMETFile*readmet < 1
    statusceil = nCEILFile*readceil < 1
    statuskazr = nkazrFile*readkazr < 1
    statussonde = (Dif_sonde_cases LE 1.)*readsonde < 1

    ;#3 Start reading data-----------------------------------------------------------
    ;MET
    IF statusmet EQ 1 THEN BEGIN
      ind = WHERE(STRLEN(METFile) EQ 0, count,COMPLEMENT = indc)
      METFile = METFile[indc]

      ;get the data product name
      a = STRSPLIT(METFile[0],'.', /extract)
      b = STRSPLIT(a[0],'\', /extract)
      metproductname = b[N_ELEMENTS(b) - 1]

      ;read met
      nMETFile = N_ELEMENTS(METFile)
      FOR iFile = 0, nMETFile - 1 DO BEGIN

        status = YZREADARMCDF(METFile[iFile],bufferMET, $
          content_list = metcontent_list,/check_missingvalue)

        tag_names0 = TAG_NAMES(bufferMET)
        ind_time = WHERE(tag_names0 EQ 'TIME')
        ind_day = WHERE(tag_names0 EQ 'DAY')
        ind_temp = WHERE(tag_names0 EQ STRUPCASE(metcontent_list[0]))
        ind_rh = WHERE(tag_names0 EQ STRUPCASE(metcontent_list[1]))
        ind_pres = WHERE(tag_names0 EQ STRUPCASE(metcontent_list[2]))

        time_tmp = bufferMET.(ind_time[0])
        Ta_tmp = buffermet.(ind_temp[0])
        rh_tmp = buffermet.(ind_rh[0])
        pres_tmp = 10.*buffermet.(ind_pres[0]) ;from kPa to hPa

        ni = N_ELEMENTS(time_tmp)
        LCL_tmp = FLTARR(ni)
        FOR i = 0, ni - 1 DO BEGIN
          LCL_tmp[i] = ROMPLCL(100.*pres_tmp[i],Ta_tmp[i] + 273.15,rh_tmp[i]/100.)/1000.
        ENDFOR

        ;Td_tmp = Ta_tmp - (100.0 - rh_tmp)/5.0
        ;LCL_old = 125.*(Ta_tmp - Td_tmp)/1000.0

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
      a = STRSPLIT(CeilFile[0],'.', /extract)
      b = STRSPLIT(a[0],'\', /extract)
      ceilproductname = b[N_ELEMENTS(b) - 1]

      nCeilFile = N_ELEMENTS(CeilFile)
      FOR iFile = 0, nCeilFile - 1 DO BEGIN
        status = YZREADARMCDF(CEILFile[iFile],bufferCEIL, $
          content_list = ['first_cbh', 'second_cbh'],/check_missingvalue)

        time_tmp = bufferCEIL.time
        first_cbh_tmp = bufferCEIL.first_cbh
        second_cbh_tmp = bufferCEIL.second_cbh

        IF iFile EQ 0 THEN BEGIN
          first_cbh_day = first_cbh_tmp/1000.
          second_cbh_day = second_cbh_tmp/1000.
          time_day_ceil = time_tmp
        ENDIF ELSE BEGIN
          first_cbh_day = [first_cbh_day, first_cbh_tmp/1000.]
          second_cbh_day = [second_cbh_day, second_cbh_tmp/1000.]
          time_day_ceil = [time_day_ceil,time_tmp]
        ENDELSE
      ENDFOR
    ENDIF

    ;radiosonde data
    IF statussonde EQ 1 THEN BEGIN
      YZREADARMCDF_OLD, sondefile_use, tdry, dp, alt, pres, wspd, deg, rh,$
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
    ENDIF

    ;kazr data
    IF statuskazr EQ 1 THEN BEGIN

      kazrFile_use = kazrFile

      ;get the data product name
      a = STRSPLIT(kazrFile_use[0],'.', /extract)
      b = STRSPLIT(a[0],'\', /extract)
      kazrproductname = b[N_ELEMENTS(b) - 1]

      FOR ifile = 0, N_ELEMENTS(kazrFile_use) - 1 DO BEGIN
        status = YZREADARMCDF(kazrFile_use[ifile],bufferkazr, $
          content_list = kazrcontent_list)

        tag_names0 = TAG_NAMES(bufferkazr)
        ind_time = WHERE(tag_names0 EQ 'TIME')
        ind_day = WHERE(tag_names0 EQ 'DAY')
        ind_hgt = WHERE(tag_names0 EQ STRUPCASE(kazrcontent_list[0]))
        ind_ref = WHERE(tag_names0 EQ STRUPCASE(kazrcontent_list[1]))

        range = (bufferkazr.(ind_hgt[0]))/1000.

        ref_tmp = bufferkazr.(ind_ref[0])
        time_tmp = bufferkazr.(ind_time[0])

        IF ifile EQ 0 THEN BEGIN
          ref_day = ref_tmp
          time_day_kazr = time_tmp
        ENDIF ELSE BEGIN
          ref_day = [[ref_day],[ref_tmp]]
          time_day_kazr = [time_day_kazr, time_tmp]
        ENDELSE
      ENDFOR

      ;double check the availability of radar data
      ind = WHERE(time_day_kazr GT Hour_sat_day - timeinterval/2. AND time_day_kazr LE Hour_sat_day + timeinterval/2., count)
      IF count EQ 0 THEN statuskazr = 0

    ENDIF

    ;#4 determine key quantities and output-----------------------------------------------------

    IF statusmet EQ 1 THEN BEGIN
      Ind = WHERE(time_day_met GT Hour_SAT_day - timeinterval/2. AND time_day_met LT Hour_SAT_day + timeinterval/2., count)
      LCL_mean_met = MEAN(LCL_day[ind],/NAN)
      Ta_mean_met = MEAN(Ta_day[ind],/NAN)
      rh_mean_met = MEAN(rh_day[ind],/NAN)
      pres_mean_met = MEAN(pres_day[ind],/NAN)

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

    IF statusceil EQ 1 THEN BEGIN
      Ind = WHERE(time_day_ceil GE Hour_SAT_Day - timeinterval/2. AND time_day_ceil LE Hour_SAT_Day + timeinterval/2., count)
      count_tot = count

      time_cbh_use = time_day_ceil[ind]
      first_cbh_use = first_cbh_day[ind]
      second_cbh_use = second_cbh_day[ind]

      ind = WHERE(first_cbh_use GE 0.1 AND first_cbh_use LE LCL_mean_met*1.2 AND FINITE(first_cbh_use) EQ 1, count)
      count_1stcld = count

      temp = first_cbh_use[ind]
      first_cbh_CF = 100.*count_1stcld/count_tot
      first_cbh_median = count NE 0 ? MEDIAN(temp[ind], /even): -999.
      first_cbh_mean = count NE 0 ? MEAN(temp[ind],/NAN): -999.
      first_cbh_Std = count NE 0 ? STDDEV(temp[ind],/NAN): -999.
      first_cbh_Skew = count NE 0 ? SKEWNESS(temp[ind],/NAN): -999.

      a = CGPERCENTILES(temp, Percentiles=[0., 0.25])
      first_cbh_bot25mean = MEAN(temp[WHERE(temp GE a[0] AND temp LE a[1])])

      nvar = 10
      str = {varname:' ', varvalue:-999., varunit:' '}
      temp_output = REPLICATE(str, nvar)

      temp_output[0].varname = 'first_cbh_CF' & temp_output[0].varvalue = first_cbh_CF & temp_output[0].varunit = '%'
      temp_output[1].varname = 'first_cbh_median' & temp_output[1].varvalue = first_cbh_median & temp_output[1].varunit = 'km'
      temp_output[2].varname = 'first_cbh_mean' & temp_output[2].varvalue = first_cbh_mean & temp_output[2].varunit = 'km'
      temp_output[3].varname = 'first_cbh_Std' & temp_output[3].varvalue = first_cbh_Std & temp_output[3].varunit = 'km'
      temp_output[4].varname = 'first_cbh_Skew' & temp_output[4].varvalue = first_cbh_Skew & temp_output[4].varunit = ' '
      temp_output[5].varname = 'first_cbh_bot25mean' & temp_output[5].varvalue = first_cbh_bot25mean & temp_output[5].varunit = 'km'

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

    ;kazr data
    IF statuskazr EQ 1 THEN BEGIN

      ;===STEP 1: computing cloud-base related quantities
      Myxrange = [Hour_SAT_day - timeinterval/2., Hour_SAT_day + timeinterval/2.]
      Myyrange = [first_cbh_bot25mean, 10.]

      ref_plot = SELECTDATA_RADAR(ref_day, range, time_day_kazr, $
        Myyrange[0], Myyrange[1], $
        Myxrange[0], Myxrange[1],$
        H_Output = H_plot, T_Output = T_plot)

      ;create arrays for cloud-base
      si = SIZE(ref_plot)
      nz = si[1]
      nt = si[2]

      t0 = T_plot
      cbh0 = MAKE_ARRAY(nt, /float, VALUE = -999.)
      cth0 = MAKE_ARRAY(nt, /float, VALUE = -999.)

      ;interpolate ceilo cbh to the kazr time series
      loc = VALUE_LOCATE(T_plot, time_cbh_use)
      cbh0[loc] = first_cbh_use
      ind = WHERE(FINITE(cbh0) EQ 0 OR cbh0 EQ -999., count, NCOMPLEMENT = ind_c)
      cbh0[ind] = -999.

      ;determine cloud top height
      FOR it = 0, nt - 1 DO BEGIN
        IF cbh0[it] EQ -999. THEN CONTINUE
        ind = WHERE(H_plot GE cbh0[it])
        z_tmp = H_plot[ind]
        ref_tmp = ref_plot[ind, it]

        ind_use = WHERE(ref_tmp GT ref_thre, count)
        IF count GT 0 THEN BEGIN
          FIND_GAP, ind_use, gap_thre, ind_use_st, ind_use_nd
          IF z_tmp[ind_use[ind_use_st[0]]] - cbh0[it] GT 0.3 OR $
            z_tmp[ind_use[ind_use_st[0]]] - first_cbh_bot25mean GT 0.3 THEN CONTINUE

          cth0[it] = z_tmp[ind_use[ind_use_nd[0]]]
        ENDIF ELSE BEGIN
          CONTINUE
        ENDELSE
        ;PRINT,''
      ENDFOR

      ;for plotting purpose
      ind = WHERE(cth0 NE -999)
      cth0_plot = cth0[ind]
      cbh0_plot = cbh0[ind]
      t0_plot = t0[ind]

      cth_mean = MEAN(cth0_plot)
      cth_max = MAX(cth0_plot)

      ;===STEP 2: output

      ;output1: column-integrated variables
      nvar = 20

      str = {varname:' ', varvalue:-999., varunit:' '}
      kazr_output = REPLICATE(str, nvar)

      kazr_output[0].varname = 'cth_mean' & kazr_output[0].varvalue = cth_mean & kazr_output[0].varunit = 'km'
      kazr_output[1].varname = 'cth_max' & kazr_output[1].varvalue = cth_max & kazr_output[1].varunit = 'km'

      kazr_output = kazr_output[WHERE(kazr_output.varvalue NE -999.)]

      nvar = N_ELEMENTS(kazr_output)

      ;output the results
      OPENW,LUN,OutPath + kazrproductname + '_'  + date_cases[icases] + '_' + $
        overpass_cases[icases] + '.txt',/GET_LUN
      PRINTF,LUN,'*Date:', date_cases[icases], FORMAT= '(A-10, A-10)'
      PRINTF,LUN,'*Satellite Overpass:', overpass_cases[icases], FORMAT= '(A-25, A-10)'
      FOR ivar = 0, nvar - 1 DO BEGIN
        PRINTF,LUN,''
        PRINTF,LUN,kazr_output[ivar].varname + ' = ', kazr_output[ivar].varvalue, kazr_output[ivar].varunit, FORMAT= '(A-30, F-9.3, A6)'
      ENDFOR
      FREE_LUN, LUN

      ;output2: cth detail
      OPENW,LUN,OutPath + kazrproductname + '_cth_'$
        + date_cases[icases] + '_' + overpass_cases[icases] + '.txt',/GET_LUN
      PRINTF,LUN, '*t0_plot','cth0_plot',$
        FORMAT= '(2A15)'

      FOR i = 0, N_ELEMENTS(cth0_plot) - 1 DO BEGIN
        PRINTF,LUN, t0_plot[i],cth0_plot[i],$
          FORMAT = '(2F15.4)'
      ENDFOR
      FREE_LUN, LUN

      ;===STEP 4: plot or not?
      IF KEYWORD_SET(plotfig) AND statuskazr NE 0 THEN BEGIN
        CGPS_OPEN,  OutPath + kazrproductname + '_'+ date_cases[icases] + '_' + $
          overpass_cases[icases] + '.eps',/CMYK,$
          /nomatch,Font = 1,/Times, xsize = 12/2.54, ysize = 8/2.54, fontsize =14

        pos = CGLAYOUT([1,1],XGap=7, YGap = 6.0, OXMargin=[4,5], OYMargin=[4,3])

        p0 = pos[*,0]
        ;p1 = pos[*,1]

        unit = !D.Y_CH_SIZE / FLOAT(!D.Y_SIZE)

        ;===========(a) Plot attenuated backscatter=============
        ref_plot = SELECTDATA_RADAR(ref_day, range, time_day_kazr, $
          range[0], MAX(cth0_plot) + 0.2 > first_cbh_bot25mean,$
          Myxrange[0], Myxrange[1],$
          H_Output = H_plot, T_Output = T_plot)

        PLOTRADAR, ref_plot, T_plot, H_plot, PlotVariable = 'ref',$
          Position = pos[*,0], MultiMargin = [0,2,0,4],$
          TTitle = 'Reflectivity [dbz]',$
          mycharsize = 1., myTcharsize = 0.9, $
          TPosition = [p0[2] + 0.2*unit,P0[1],P0[2]+0.4*unit, P0[3]],$
          title = 'Hcb = ' + STRING(first_cbh_bot25mean,format = '(F0.2)') + 'km'

        IF statusceil THEN CGOPLOT, time_day_ceil, first_cbh_day, psym = 16, symsize = 0.3, color = 'green'
        IF statusmet THEN CGOPLOT, t0_plot, cth0_plot, psym = 2, symsize = 0.3, color = 'purple'
        CGPS_CLOSE, /png

      ENDIF
    ENDIF
    PRINT, ' '
  ENDFOR

END
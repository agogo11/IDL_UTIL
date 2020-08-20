PRO YZDLICMES, OutPath,$
  dateArray, overpassArray,$
  DLPath = dlpath,$
  CEILPath = ceilpath,$
  METPath = metpath,$
  sondepath = sondepath,$
  ecorpath = ecorpath,$
  Plotfig = plotfig,$
  plotsonde = plotsonde, $
  timeinterval = timeinterval,$
  readdl = readdl,$
  readmet = readmet,$
  readceil = readceil,$
  readsonde  =readsonde, $
  readecor = readecor, $
  metcontent_list = metcontent_list, $
  outputdl = outputdl,$
  outputmet = outputmet,$
  outputceil = outputceil,$
  outputsonde = outputsonde,$
  outputecor = outputecor,$
  ind_st = ind_st,$
  ind_nd = ind_nd

  ;set default
  SETDEFAULTVALUE, outputmet, 1
  SETDEFAULTVALUE, outputceil, 1
  SETDEFAULTVALUE, outputdl, 1
  SETDEFAULTVALUE, outputsonde, 1
  SETDEFAULTVALUE, outputecor, 1

  SETDEFAULTVALUE, readmet, 1
  SETDEFAULTVALUE, readceil, 1
  SETDEFAULTVALUE, readdl, 1
  SETDEFAULTVALUE, readsonde, 1
  SETDEFAULTVALUE, readecor, 1

  SETDEFAULTVALUE, metcontent_list, ['temp_mean','rh_mean','atmos_pressure']
  SETDEFAULTVALUE, dlcontent_list, ['range','radial_velocity','intensity','attenuated_backscatter']
  SETDEFAULTVALUE, timeinterval, 2.
  SETDEFAULTVALUE, ind_st, 0
  SETDEFAULTVALUE, ind_nd, N_ELEMENTS(dateArray) - 1

  ;---consts
  Lv = 2501. ; kj/kg
  cpd = 1005. ; J/(kg.K)
  g = 9.80665d ; gravity const.
  rho_ref       = 1.d ;  density of air (kg/m3)
  p_ref = 1000.;  reference  pressure, hPa

  ref_thre = 10.^(-4.6) ;threshold for determining cloudy pixels in DL data
  SNR_thre = 0.005 ; threshold for removing noisy data in DL
  gap_thre = 20. ; gap (s) for identifying single clouds

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
    statusdl = 1
    statusecor = 1
    statussonde = 1

    ;dl
    dlFile=FILE_SEARCH(dlPath+ date[1] + '\*'+ date[1] +'*', count = count1)
    ndlFile = count1

    IF count1 EQ 0 THEN BEGIN
      dlFile=FILE_SEARCH(dlPath+'*'+ date[1] +'*', count = count1)
      ndlFile = count1
    ENDIF

    ;Met
    METFile=FILE_SEARCH(metPath+'*'+date[1]+'*',count = count1)
    nMETFile = count1

    ;ceilometer
    CEILFile=FILE_SEARCH(ceilPath+'*'+date[1]+'*',count = count1)
    nCEILFile = count1

    ;ecor
    ecorFile=FILE_SEARCH(ecorPath+'*'+date[1]+'*',count = count1)
    necorFile = count1

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
    statusecor = necorFile*readecor < 1
    statusceil = nCEILFile*readceil < 1
    statusdl = ndlFile*readdl < 1
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

    ;ecor data
    IF statusecor EQ 1 THEN BEGIN
      ind = WHERE(STRLEN(EcorFile) EQ 0, count,COMPLEMENT = indc)
      EcorFile = EcorFile[indc]

      ;get the data product name
      a = STRSPLIT(EcorFile[0],'.', /extract)
      b = STRSPLIT(a[0],'\', /extract)
      ecorproductname = b[N_ELEMENTS(b) - 1]

      nEcorFile = N_ELEMENTS(EcorFile)
      FOR iFile = 0, nEcorFile - 1 DO BEGIN
        status = YZREADARMCDF(EcorFile[iFile],bufferECOR, $
          content_list = ['sensible_heat_flux', 'latent_heat_flux'],/check_missingvalue)

        time_tmp = bufferECOR.time
        SHF_tmp = bufferECOR.sensible_heat_flux
        LHF_tmp = bufferECOR.latent_heat_flux

        IF iFile EQ 0 THEN BEGIN
          SHF_day = SHF_tmp
          LHF_day = LHF_tmp
          time_day_ecor = time_tmp
        ENDIF ELSE BEGIN
          SHF_day = [SHF_day, SHF_tmp]
          LHF_day = [LHF_day, LHF_tmp]
          time_day_ecor = [time_day_ecor,time_tmp]
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

    ;DL data
    IF statusdl EQ 1 THEN BEGIN

      DLFile_use = DLFile

      ;get the data product name
      a = STRSPLIT(DLFile_use[0],'.', /extract)
      b = STRSPLIT(a[0],'\', /extract)
      dlproductname = b[N_ELEMENTS(b) - 1]

      FOR ifile = 5, N_ELEMENTS(DLFile_use) - 1 DO BEGIN
        status = YZREADARMCDF(DLFile_use[ifile],bufferDL, $
          content_list = dlcontent_list)

        tag_names0 = TAG_NAMES(bufferDL)
        ind_time = WHERE(tag_names0 EQ 'TIME')
        ind_day = WHERE(tag_names0 EQ 'DAY')
        ind_hgt = WHERE(tag_names0 EQ STRUPCASE(dlcontent_list[0]))
        ind_w = WHERE(tag_names0 EQ STRUPCASE(dlcontent_list[1]))
        ind_SNR = WHERE(tag_names0 EQ STRUPCASE(dlcontent_list[2]))
        ind_ref = WHERE(tag_names0 EQ STRUPCASE(dlcontent_list[3]))

        range = (bufferDL.(ind_hgt[0]))/1000.

        w_tmp = bufferDL.(ind_w[0])
        SNR_tmp = bufferDL.(ind_SNR[0])
        ref_tmp = bufferDL.(ind_ref[0])
        time_tmp = bufferDL.(ind_time[0])

        ind = WHERE(SNR_tmp - 1 LT SNR_thre)
        w_tmp[ind] = !values.f_nan
        SNR_tmp = !Null

        IF ifile EQ 5 THEN BEGIN
          w_day = w_tmp
          ;SNR_day = SNR_tmp
          ref_day = ref_tmp
          time_day_DL = time_tmp
        ENDIF ELSE BEGIN
          time_day_dl = [time_day_dl, time_tmp]
          w_day = [[w_day],[w_tmp]]
          ;SNR_day = [[SNR_day],[SNR_tmp]]
          ref_day = [[ref_day],[ref_tmp]]
        ENDELSE
      ENDFOR

      ;double check the availability of radar data
      ind = WHERE(time_day_dl GT Hour_sat_day - timeinterval/2. AND time_day_dl LE Hour_sat_day + timeinterval/2., count)
      IF count EQ 0 THEN statusdl = 0

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

    IF statusecor EQ 1 THEN BEGIN
      Ind = WHERE(time_day_ecor GT Hour_SAT_day - timeinterval/2. AND time_day_ecor LT Hour_SAT_day + timeinterval/2., count)
      SHF_mean_ecor = MEAN(SHF_day[ind],/NAN)
      LHF_mean_ecor = MEAN(LHF_day[ind],/NAN)

      nvar = 20

      str = {varname:' ', varvalue:-999., varunit:' '}
      ecor_output = REPLICATE(str, nvar)

      ecor_output[0].varname = 'SHF_mean_ecor' & ecor_output[0].varvalue = SHF_mean_ecor & ecor_output[0].varunit = 'W/m2'
      ecor_output[1].varname = 'LHF_mean_ecor' & ecor_output[1].varvalue = LHF_mean_ecor & ecor_output[1].varunit = 'W/m2'

      ecor_output = ecor_output[WHERE(ecor_output.varvalue NE -999.)]

      nvar = N_ELEMENTS(ecor_output)

      ;output the results
      OPENW,LUN,OutPath + ecorproductname + '_'  + date_cases[icases] + '_' + $
        overpass_cases[icases] + '.txt',/GET_LUN
      PRINTF,LUN,'*Date:', date_cases[icases], FORMAT= '(A-10, A-10)'
      PRINTF,LUN,'*Satellite Overpass:', overpass_cases[icases], FORMAT= '(A-25, A-10)'
      FOR ivar = 0, nvar - 1 DO BEGIN
        PRINTF,LUN,''
        PRINTF,LUN,ecor_output[ivar].varname + ' = ', ecor_output[ivar].varvalue, ecor_output[ivar].varunit, FORMAT= '(A-30, F-9.3, A6)'
      ENDFOR
      FREE_LUN, LUN
    ENDIF

    ;    ;dl horizontal wind data (temporary on July 15, 2020)
    ;    ;======================================================================
    ;    File=FILE_SEARCH('D:\Data\SGP_dlprofwind\*' + date[1]+'*',count = count1)
    ;
    ;    status = YZREADARMCDF(File[0],buffer, $
    ;      content_list = ['height', 'wind_speed'],/check_missingvalue)
    ;
    ;    time_tmp = buffer.time
    ;    range_tmp = buffer.height/1000.
    ;    wspd_tmp = buffer.wind_speed
    ;
    ;    Mytrange = [Hour_SAT_day - timeinterval/2., Hour_SAT_day + timeinterval/2.]
    ;    Myzrange = [first_cbh_bot25mean - 0.1, first_cbh_bot25mean + 0.1]
    ;
    ;    indt = where(time_tmp ge Mytrange[0] and time_tmp le Mytrange[1], count)
    ;    countt = count
    ;    indz = where(range_tmp ge Myzrange[0] and range_tmp le Myzrange[1], count)
    ;    countz = count
    ;
    ;    a = wspd_tmp[indz[0]:indz[countz-1], indt[0]:indt[countt-1]]
    ;    wspd_cb = mean(a,/nan)
    ;
    ;    nvar = 1
    ;
    ;    str = {varname:' ', varvalue:-999., varunit:' '}
    ;    dlwind_output = REPLICATE(str, nvar)
    ;
    ;    dlwind_output[0].varname = 'wspd_cb' & dlwind_output[0].varvalue = wspd_cb & dlwind_output[0].varunit = 'm/s'
    ;
    ;    ;output the results
    ;    OPENW,LUN,OutPath + 'sgpdlprofwind4newsC1_'  + date_cases[icases] + '_' + $
    ;      overpass_cases[icases] + '.txt',/GET_LUN
    ;    PRINTF,LUN,'*Date:', date_cases[icases], FORMAT= '(A-10, A-10)'
    ;    PRINTF,LUN,'*Satellite Overpass:', overpass_cases[icases], FORMAT= '(A-25, A-10)'
    ;    FOR ivar = 0, nvar - 1 DO BEGIN
    ;      PRINTF,LUN,''
    ;      PRINTF,LUN,dlwind_output[ivar].varname + ' = ', dlwind_output[ivar].varvalue, dlwind_output[ivar].varunit, FORMAT= '(A-30, F-9.3, A6)'
    ;    ENDFOR
    ;    FREE_LUN, LUN
    ;
    ;    print,''
    ;    ;======================================================================

    ;dl data
    IF statusdl EQ 1 THEN BEGIN

      ;===STEP 1: computing cloud-base related quantities
      Myxrange = [Hour_SAT_day - timeinterval/2., Hour_SAT_day + timeinterval/2.]
      Myyrange = [first_cbh_bot25mean - 0.5 > range[3], first_cbh_bot25mean + 2.]

      w_plot = SELECTDATA_RADAR(W_day, range, time_day_dl, $
        Myyrange[0], Myyrange[1], $
        Myxrange[0], Myxrange[1],$
        H_Output = H_plot, T_Output = T_plot)

      ref_plot = SELECTDATA_RADAR(ref_day, range, time_day_dl, $
        Myyrange[0], Myyrange[1], $
        Myxrange[0], Myxrange[1],$
        H_Output = H_plot, T_Output = T_plot)

      ;create arrays for cloud-base
      si = SIZE(w_plot)
      nz = si[1]
      nt = si[2]

      t0 = T_plot
      cbh0 = MAKE_ARRAY(nt, /float, VALUE = -999.)
      w0 = MAKE_ARRAY(8, nt, /float, VALUE = -999.)
      ref0 = MAKE_ARRAY(8, nt, /float, VALUE = -999.)

      ;determine cloud-base varaibles
      FOR it = 0, nt - 1 DO BEGIN
        ref_tmp = ref_plot[*, it]
        w_tmp = w_plot[*, it]
        z_tmp = H_plot

        ind = WHERE(ref_tmp GT ref_thre, count)
        IF count GT 0 THEN BEGIN
          cbh0[it] = z_tmp[ind[0]]

          IF cbh0[it] LE first_cbh_bot25mean*1.3 THEN BEGIN
            IF ind[0] - 4 LT 0 OR ind[0] + 3 GT nz - 1 THEN CONTINUE
            w0[*,it] = w_tmp[ind[0] - 4: ind[0] + 3]
            ref0[*,it] = ref_tmp[ind[0] - 4: ind[0] + 3]
          ENDIF
        ENDIF
      ENDFOR

      ;compute single cloud and related properties
      ind_use = WHERE(cbh0 NE -999., count)
      IF count GT 0 THEN BEGIN
        FIND_GAP, ind_use, gap_thre, ind_use_st, ind_use_nd
      ENDIF

      ;create arrays related to single-cloud properties
      ncld1 = N_ELEMENTS(ind_use_st)
      duration_cld1 = MAKE_ARRAY(ncld1, /float, VALUE = -999.)
      cbh_cld1 = MAKE_ARRAY(ncld1, /float, VALUE = -999.)
      t_cld1 = MAKE_ARRAY(ncld1, /float, VALUE = -999.)
      w_cld1 = MAKE_ARRAY(ncld1, /float, VALUE = -999.)

      FOR icld = 0, ncld1 - 1 DO BEGIN
        cbh_cld1[icld] = MEAN(cbh0[ind_use[ind_use_st[icld]:ind_use_nd[icld]]])
        t_cld1[icld] = MEAN(T_plot[ind_use[ind_use_st[icld]:ind_use_nd[icld]]])
        duration_cld1[icld] = ind_use_nd[icld] - ind_use_st[icld] + 1

        tmp = w0[*,[ind_use[ind_use_st[icld]:ind_use_nd[icld]]]]
        ind = WHERE(tmp GT 0., count)
        IF count EQ 0 THEN BEGIN
          w_cld1[icld] = 0.
        ENDIF ELSE BEGIN
          w_cld1[icld] = YZAVERAGE(tmp[ind])
        ENDELSE
      ENDFOR

      ;for plotting purpose
      ind = WHERE(duration_cld1 GE 10)
      cbh_cld1_plot = cbh_cld1[ind]
      t_cld1_plot = t_cld1[ind]
      w_cld1_plot = w_cld1[ind]

      ind = WHERE(ref0 GT ref_thre)
      w_up_hist = w0[ind]
      w_down_hist = w0[0:3, *]
      w_down_hist = w_down_hist[WHERE(w_down_hist NE -999.)]

      w_final = [w_down_hist, w_up_hist]

      wb_1st_00 = MEAN(w_final[WHERE(w_final GT 0.)])
      wb_1st_01 = MEAN(w_final[WHERE(w_final GT 0.1)])
      wb_1st_05 = MEAN(w_final[WHERE(w_final GT 0.5)])
      wb_2nd_00 = YZCALW(w_final, 0, 0)
      wb_2nd_01 = YZCALW(w_final, 0, 0.1)
      wb_2nd_05 = YZCALW(w_final, 0, 0.5)

      percent_cb = 100.*N_ELEMENTS(w_final)/(8*nt)

      ;===STEP 2: computing sub-cloud quantities

      ;determin how many vertical layers and create arrays
      dlev = 0.2 ;km
      nlev = FLOOR((first_cbh_bot25mean - range[3])/dlev)

      level_vp = range[3] + dlev*INDGEN(nlev)
      ww_vp = MAKE_ARRAY(nlev, /float, VALUE = -999.)
      wskew_vp = MAKE_ARRAY(nlev, /float, VALUE = -999.)
      percent_vp = MAKE_ARRAY(nlev, /float, VALUE = -999.)

      FOR ilev = 0, nlev - 1 DO BEGIN
        w_plot = SELECTDATA_RADAR(w_day, range, Time_day_dl, $
          level_vp[ilev], level_vp[ilev] + dlev,$
          Myxrange[0], Myxrange[1])

        count_tot = N_ELEMENTS(w_plot)
        Index = WHERE(FINITE(w_plot) EQ 1, count)
        w_final = w_plot[index]
        count_use = count

        ww_vp[ilev] = MEAN((W_final - MEAN(W_final))^2.)
        wskew_vp[ilev] = SKEWNESS(W_final)
        percent_vp[ilev] = 100.*count_use/count_tot
      ENDFOR

      wwmax = MAX(ww_vp)

      ;compute vertically averaged w variance
      w_plot = SELECTDATA_RADAR(w_day, range, Time_day_dl, $
        range[3], first_cbh_bot25mean - 0.1,$
        Myxrange[0], Myxrange[1])

      Index = WHERE(FINITE(w_plot) EQ 1, count)
      w_final = w_plot[index]
      wwml = MEAN((W_final - MEAN(W_final))^2.)

      ;===STEP 3: output

      ;output1: column-integrated variables
      nvar = 20

      str = {varname:' ', varvalue:-999., varunit:' '}
      dl_output = REPLICATE(str, nvar)

      dl_output[0].varname = 'wb_1st_00' & dl_output[0].varvalue = wb_1st_00 & dl_output[0].varunit = 'm/s'
      dl_output[1].varname = 'wb_1st_01' & dl_output[1].varvalue = wb_1st_01 & dl_output[1].varunit = 'm/s'
      dl_output[2].varname = 'wb_2nd_00' & dl_output[2].varvalue = wb_2nd_00 & dl_output[2].varunit = 'm/s'
      dl_output[3].varname = 'wb_2nd_01' & dl_output[3].varvalue = wb_2nd_01 & dl_output[3].varunit = 'm/s'
      dl_output[4].varname = 'wwml' & dl_output[4].varvalue = wwml & dl_output[4].varunit = 'm2/s2'
      dl_output[5].varname = 'wwmax' & dl_output[5].varvalue = wwmax & dl_output[5].varunit = 'm2/s2'
      dl_output[6].varname = 'percent_cb' & dl_output[6].varvalue = percent_cb & dl_output[6].varunit = '%'
      dl_output[7].varname = 'wb_1st_05' & dl_output[7].varvalue = wb_1st_05 & dl_output[7].varunit = 'm/s'
      dl_output[8].varname = 'wb_2nd_05' & dl_output[8].varvalue = wb_2nd_05 & dl_output[8].varunit = 'm/s'

      dl_output = dl_output[WHERE(dl_output.varvalue NE -999.)]

      nvar = N_ELEMENTS(dl_output)

      ;output the results
      OPENW,LUN,OutPath + dlproductname + '_'  + date_cases[icases] + '_' + $
        overpass_cases[icases] + '.txt',/GET_LUN
      PRINTF,LUN,'*Date:', date_cases[icases], FORMAT= '(A-10, A-10)'
      PRINTF,LUN,'*Satellite Overpass:', overpass_cases[icases], FORMAT= '(A-25, A-10)'
      FOR ivar = 0, nvar - 1 DO BEGIN
        PRINTF,LUN,''
        PRINTF,LUN,dl_output[ivar].varname + ' = ', dl_output[ivar].varvalue, dl_output[ivar].varunit, FORMAT= '(A-30, F-9.3, A6)'
      ENDFOR
      FREE_LUN, LUN

      ;output2: profiles
      OPENW,LUN,OutPath + dlproductname + '_ww_vp_'$
        + date_cases[icases] + '_' + overpass_cases[icases] + '.txt',/GET_LUN
      PRINTF,LUN, '*altitude','ww', 'wskew', 'percent',$
        FORMAT= '(4A15)'

      FOR ilev = 0, nlev - 1 DO BEGIN
        PRINTF,LUN, level_vp[ilev] + 0.5*dlev,ww_vp[ilev],wskew_vp[ilev],percent_vp[ilev],$
          FORMAT = '(4F15.4)'
      ENDFOR
      FREE_LUN, LUN

      ;output3: indivitual cumuli
      so = SORT(duration_cld1)
      ncu = N_ELEMENTS(so)

      OPENW,LUN,OutPath + dlproductname + '_pixelcbw_'$
        + date_cases[icases] + '_' + overpass_cases[icases] + '.txt',/GET_LUN
      PRINTF,LUN, '*w','above_cld',$
        FORMAT= '(2A15)'

      FOR i = 0, n_elements(w_up_hist) - 1 DO BEGIN
        PRINTF,LUN, w_up_hist[i],1,$
          FORMAT = '(2F15.4)'
      ENDFOR
      
      FOR i = 0, N_ELEMENTS(w_down_hist) - 1 DO BEGIN
        PRINTF,LUN, w_down_hist[i],0,$
          FORMAT = '(2F15.4)'
      ENDFOR
      
      FREE_LUN, LUN
      
      ;output4: pixel-level vertical velocity at cloud bases
      so = SORT(duration_cld1)
      ncu = N_ELEMENTS(so)

      OPENW,LUN,OutPath + dlproductname + '_cumuli_'$
        + date_cases[icases] + '_' + overpass_cases[icases] + '.txt',/GET_LUN
      PRINTF,LUN, '*duration_cld1','cbh_cld1', 'w_cld1', 't_cld1',$
        FORMAT= '(4A15)'

      FOR i = 0, ncu - 1 DO BEGIN
        PRINTF,LUN, duration_cld1[so[i]],cbh_cld1[so[i]],w_cld1[so[i]],t_cld1[so[i]],$
          FORMAT = '(4F15.4)'
      ENDFOR
      FREE_LUN, LUN

      ;===STEP 4: plot or not?
      IF KEYWORD_SET(plotfig) AND statusdl NE 0 THEN BEGIN
        CGPS_OPEN,  OutPath + dlproductname + '_'+ date_cases[icases] + '_' + $
          overpass_cases[icases] + '.eps',/CMYK,$
          /nomatch,Font = 1,/Times, xsize = 20/2.54, ysize = 20/2.54, fontsize =14

        pos = CGLAYOUT([2,4],XGap=7, YGap = 6.0, OXMargin=[4,5], OYMargin=[4,3])

        p0 = pos[*,0]
        p1 = pos[*,1]
        p2 = pos[*,2]
        p3 = pos[*,3]
        p4 = pos[*,4]
        p5 = pos[*,5]
        p6 = pos[*,6]
        p7 = pos[*,7]

        unit = !D.Y_CH_SIZE / FLOAT(!D.Y_SIZE)

        ;===========(a) Plot attenuated backscatter=============
        ref_plot = SELECTDATA_RADAR(ref_day, range, time_day_dl, $
          range[0], 3.6,$
          Myxrange[0], Myxrange[1],$
          H_Output = H_plot, T_Output = T_plot)

        ref_plot = ALOG10(ref_plot)

        PLOTRADAR, ref_plot, T_plot, H_plot, PlotVariable = 'ATTENUATED BACKSCATTER',$
          Position = pos[*,0], MultiMargin = [0,2,0,4],$
          TTitle = 'Atten Backscatter',$
          mycharsize = 1., myTcharsize = 0.9, $
          TPosition = [p0[2] + 0.2*unit,P0[1],P0[2]+0.4*unit, P0[3]],$
          title = 'Hcb = ' + STRING(first_cbh_bot25mean,format = '(F0.2)') + 'km'

        IF statusceil THEN CGOPLOT, time_day_ceil, first_cbh_day, psym = 16, symsize = 0.1, color = 'green'
        IF statusmet THEN CGOPLOT, time_day_met, LCL_day, psym = 2, symsize = 0.1, color = 'red'

        ;===========(b) Plot w=============
        w_plot = SELECTDATA_RADAR(w_day, range, time_day_dl, $
          range[0], 3.6,$
          Myxrange[0], Myxrange[1],$
          H_Output = H_plot, T_Output = T_plot)

        PLOTRADAR, w_plot, T_plot, H_plot, PlotVariable = 'w',$
          Position = pos[*,1], MultiMargin = [0,2,0,4],$
          TTitle = 'w',$
          mycharsize = 1., myTcharsize = 0.9, $
          TPosition = [p1[2] + 0.2*unit,P1[1],P1[2]+0.4*unit, P1[3]],$
          title = 'Hcb = ' + STRING(first_cbh_bot25mean,format = '(F0.2)') + 'km'

        IF statusceil THEN CGOPLOT, time_day_ceil, first_cbh_day, psym = 16, symsize = 0.1, color = 'green'
        IF statusmet THEN CGOPLOT, time_day_met, LCL_day, psym = 2, symsize = 0.1, color = 'red'


        ;===========(c) Plot ww vertical profile=============
        CGPLOT, ww_vp, level_vp, color = 'red', $
          xtitle = 'w variance [m2/s2]', ytitle = 'Altitude [km]',$
          xrange = [0,1.1*MAX(ww_vp)], yrange = [0,1.1*MAX(level_vp)],$
          Position = pos[*,2], charsize = 1., /noerase, $
          title = 'wwmax = ' + STRING(MAX(ww_vp),format = '(F0.2)') + 'm2/s2, ' +$
          'wwml = ' + STRING(MEAN(ww_vp),format = '(F0.2)') + 'm2/s2'

        ;===========(d) Plot cloud base DL w=============
        Myyrange = [first_cbh_bot25mean - 0.5, first_cbh_bot25mean + 1.]

        w_plot = SELECTDATA_RADAR(W_day, range, time_day_dl, $
          Myyrange[0], Myyrange[1], $
          Myxrange[0], Myxrange[1],$
          H_Output = H_plot, T_Output = T_plot)

        ref_plot = SELECTDATA_RADAR(ref_day, range, time_day_dl, $
          Myyrange[0], Myyrange[1], $
          Myxrange[0], Myxrange[1],$
          H_Output = H_plot, T_Output = T_plot)

        PLOTRADAR, w_plot, T_plot, H_plot, PlotVariable = 'w',$
          Position = pos[*,3], MultiMargin = [0,2,0,4],$
          TTitle = 'w',$
          mycharsize = 1., myTcharsize = 0.9, $
          TPosition = [p3[2] + 0.2*unit,P3[1],P3[2]+0.4*unit, P3[3]],$
          title = 'wb_1st_00 = ' + STRING(wb_1st_00,format = '(F0.2)') + 'm/s'

        IF statusmet THEN CGOPLOT, time_day_met, LCL_day, psym = 2, symsize = 0.1, color = 'red'

        CGCONTOUR, TRANSPOSE(ref_plot), T_plot, H_plot, LEVELS = [ref_thre], thick = 2, LABEL = 0, /overplot

        ;===========(f) Plot w for single cloud=============
        CGPLOT, /noerase, t_cld1, w_cld1, psym = 16, symsize = 0.8, $
          xrange = myxrange, yrange = [0, 1.2*MAX(w_cld1)], $
          position = pos[*,5], charsize = 1.,$
          title = 'wb_1st_01 = ' + STRING(wb_1st_01,format = '(F0.2)') + 'm/s'

        CGPS_CLOSE, /png

      ENDIF
    ENDIF
    PRINT, ' '
  ENDFOR

END
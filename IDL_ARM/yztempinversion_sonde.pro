FUNCTION yztempinversion_sonde, alt, t, q, pres, maxalt = maxalt

  SetDefaultValue, maxalt, 3.

  zi_inv_max0 = 999.
  dpt_dz_max0 = 999.

  pt = t*(1000./pres)^(0.286)

  pt = SMOOTH(pt, 3)
  q = SMOOTH(q, 3)

  ind = WHERE(alt LT maxalt AND alt GT 0.1)
  alt_use = alt[ind]
  pt_use = pt[ind]
  t_use = t[ind]
  q_use = q[ind]
  pres_use = pres[ind]
  rh_use = yzq2rh(q_use/1000., pt_use, pres_use)
  dpt_dz_use = deriv(alt_use, pt_use)

  ind_dptdz = WHERE(dpt_dz_use GT 5., count)
  find_gap, ind_dptdz , 1, st, nd

  nlayer = N_ELEMENTS(st)
  deltapt = FLTARR(nlayer)

  FOR ilayer = 0, nlayer - 1 DO BEGIN
    temp = pt_use[ind_dptdz[st[ilayer]]:ind_dptdz[nd[ilayer]]]
    deltapt[ilayer] = temp[N_ELEMENTS(temp) - 1] - temp[0]
  ENDFOR

  ;The main inversion is determined as the altitude where the maximum dtheta/dz occurs
  dpt_dz_max = MAX(dpt_dz_use, ind_max)
  z_inv_max = alt_use[ind_max]

  ;The inversion-base height is the altitude with lowest temperature below the main inversion
  temp = MIN(t_use[0:ind_max], ind_mint)
  zi_inv_max = alt_use[ind_mint]
  t_inv_max = t_use[ind_mint]

  ;examine if the inversion caps the clouds via examining the RH
  rh_inv_max = rh_use[ind_mint]

  IF rh_inv_max LT 10. THEN BEGIN
    dpt_dz_use0 = dpt_dz_use[WHERE(alt_use LT zi_inv_max - 0.1)]
    alt_use0 = alt_use[WHERE(alt_use LT zi_inv_max - 0.1)]
    t_use0 = t_use[WHERE(alt_use LT zi_inv_max - 0.1)]

    ;determine the inversion-layer below
    dpt_dz_max0 = MAX(dpt_dz_use0, ind_max)
    temp = MIN(t_use0[0:ind_max], ind_mint)
    zi_inv_max0 = alt_use0[ind_mint]
  ENDIF

  ;Method from McGibbon and Bretherton, 2017, JAMES. Noted as MB17
  ind = WHERE(deltapt GT 2., count)
  IF count EQ 0 THEN BEGIN
    z_inv_MB17 = -999.
  ENDIF ELSE BEGIN
    pt_temp = pt_use[ind_dptdz[st[ind[0]]: nd[ind[0]]]]
    alt_temp = alt_use[ind_dptdz[st[ind[0]]: nd[ind[0]]]]
    pt_diff_temp = pt_temp - pt_temp[0]

    ind = WHERE(pt_diff_temp GE 2.)
    z_inv_MB17 = alt_temp[ind[0]]
  ENDELSE

  ;prepare for output
  inv_sonde = create_struct('zi_inv_max',zi_inv_max)
  inv_sonde = create_struct('z_inv_max',z_inv_max)
  inv_sonde = create_struct(inv_sonde,'dpt_dz_max',dpt_dz_max)
  inv_sonde = create_struct(inv_sonde,'zi_inv_max0',zi_inv_max0)
  inv_sonde = create_struct(inv_sonde,'dpt_dz_max0',dpt_dz_max0)
  inv_sonde = create_struct(inv_sonde,'z_inv_MB17',z_inv_MB17)
  inv_sonde = create_struct(inv_sonde,'t_inv_max',t_inv_max)

  RETURN, inv_sonde
END
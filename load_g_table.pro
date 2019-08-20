;---------------------------------------------------------------------
;
; Brief description:
;
; Load look up table of:
;
; Emission Lines: contribution function G [erg    cm^+3 sec^-1]      for
; user-specified ion and line.
;
; EUV Broadbands: TRF          function G [photon cm^+3 sec^-1 sr-1] for
; user-specified insrument and band.
;
; Note the different units of both tables: PH   in EUVBANDS, while ERG in EMISSIONLINES,
;                              as well as: sr-1 in EUVBANDS.
;
; The table is an IDL SAVE file for lines and TXT file for EUV bands, containing:
;
; For lines:
; 3D array:  G(Te, Ne, r) 
; 1D arrays: log(Te [K]), log10(Ne [cm-3]), r [Rsun],
; scalar: photosphere Teff [K].
;
; For EUV bands:
; 1D arrays: TRF(Te), log(Te [K])
;
; Tables must be stored in:
; tomroot/MultiTom/Emissivity/LookUp_Tables/
;
; INPUTS:
;
; To select a line set /line, and provide:
; ion_label: a string specifying the ion, possible values are:
; 'fexiii', ....
; line_wavelength: a string specifying the wavelength in A, possible
; values are: '10747', '10801'...
;
; To select a EUV band set /euvband, and provide:
; instrument_label = 'aia', 'euvi', 'eit'
; band_label = '171', '193', '195', '211', '284', '335'.
;
; OUTPUT:
;
; EMISSIONLINE: 1D arrays: N_e [cm^-3], T_e [K], r [Rsun]
;               2D array:  G(T_e,N_e,r) [ERG    cm^+3 sec^-1]
;
; EUVBANDS:     1D arrays: N_e [cm^-3], T_e [K].
;               1D array:  G(T_e)       [PHOTON cm^+3 sec^-1]
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

pro load_g_table,ion_label=ion_label,line_wavelength=line_wavelength,emissionline=emissionline,$
                 euvband=euvband,instrument_label=instrument_label,band_label=band_label

  common G_table, G, T_e, N_e, r, photT
  common directories, tomroot

  data_dir  = tomroot+'MultiTom/Emissivity_LookUp_Tables/'

  if keyword_set(emissionline) then begin
     file_name = 'G_function_'+ion_label+'_'+line_wavelength+'.save'
     restore,data_dir+file_name
     G     = emissivity            ; [ERG cm^+3 sec^-1]
     T_e   = 10.^temp              ; [K]
     N_e   = 10.^dens              ; [cm^-3]
     r     = rphot                 ; [Rsun]
     photT = radtemp               ; [K]
  endif

  if keyword_set(euvband) then begin
     xstring = ''
     file_name = 'TRF_function_'+instrument_label+'_'+band_label+'.txt'
     openr,1,data_dir+file_name
     for i=1,8 do readf,1,xstring
     x=0.
     Ntemp=0
     readf,1,x,Ntemp,x,x,x,x
     for i=1,3 do readf,1,xstring
     logTe =fltarr(Ntemp)
     TRF   =fltarr(Ntemp)
     logTe0=0.
     TRF0  =0.
     for itemp=0,Ntemp-1 do begin
        readf,1,logTe0,TRF0,x
        logTe[itemp]=logTe0
        TRF[itemp]=TRF0
     endfor
     close,1
; Note that TRF look-up table includes sr^-1 in its units, multiplying
; by 4pi below to put G in same units as for emission lines. Emissivity
; is divided back by 4*!pi in s_function.pro.
     G     = TRF*4.*!pi            ; [PHOTON cm^+3 sec^-1]
     T_e   = 10.^logTe             ; [K]

; Make a Ne 1D-array for posterior use in integrals. In a future
; version G will trully be a function of both (Te,Ne), even if weakly
; dependent on Ne.
     NNe   = 50
     logNe = 5. + (10.-5.) * findgen(NNe)/float(NNe-1)
     N_e   = 10.^logNe
  endif
  
  return
end

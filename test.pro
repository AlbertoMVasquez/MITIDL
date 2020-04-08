;---------------------------------------------------------------------
;
; Brief description:
;
; Wrapper routine to test suite of MIT codes.
;
; History:  V1.0, A.M. Vasquez, CLaSP, Spring-2018.
;           V1.5, F.A. Nuevo, IAFE, March-2020.
;           V1.6, A.M. Vasquez, IAFE, March-2020. 
;---------------------------------------------------------------------

pro test,Ne0=Ne0,Te0=Te0,euvband=euvband,emissionline=emissionline,$
         instrument_label=instrument_label,band_label=band_label,$
         ion_label=ion_label,line_wavelength=line_wavelength
  
  common constants, Rsun, kB, h, c
  common G_table, G, T_e, N_e, r, photT
  common directories, tomroot
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common dimensions, NTe, NNe
  common NT_limits, Ne0_Limits, Te0_Limits
  common tomographic_measurements, y0, y
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common measurement_errors,sig_WL,sig_y
  common index_measurement, i_measurement
;--------------------------------------------------------------------------------------------------------------------------
;                                     VALUES TO PLAY WITH
  
  ; Measurement values:
  y0 = 1.5e8                                ; KCOR Ne [cm^-3]
  y  = [15.e-10, 7.e-10, 250. , 500. ,300.] ; [e_1074,e_1079,FBE_171,FBE_193,FBE_211]
                                            ; in their respective usual units.
  i_mea_vec           =[0       ,0       ,1    ,1    ,1    ]
  ion_label_vec       =['fexiii','fexiii',''   ,''   ,''   ]
  line_wavelength_vec =['10747' ,'10801' ,''   ,''   ,''   ]
  instrument_label_vec=[''      ,''      ,'aia','aia','aia']
  band_label_vec      =[''      ,''      ,'171','193','211']
  
  ; Fractional error of each measurement:
  f_wl = 0.1
  f_y  = 0.1 + findgen(n_elements(y))

  ; Absolute error of each measurement:
  sig_WL = f_wl* y0
  sig_y  = f_y * y

  ; Test values for the parameters of the joint bivariate Te-Ne normal distribution:
  Tem        = 1.30e6 ; K
  SigTe      = 0.50e6 ; K
  Nem        = 1.75e8 ; cm^-3
  SigNe      = 0.50e8 ; cm^-3
  q          = 0.5

  ; Test values for coronal heliocentric height and iron abundance:
  r0         = 1.1    ; Rsun
  fip_factor = 1.0    ; Note that [Fe] = [Fe]_Feldman * fip_factor

  ; Parameter vector for both e_function and cost_function:
  parameters = [Nem, fip_factor, Tem, SigTe, SigNe, q]
  
  ; Test values for Ne and Te at which evaluate several functions:
  Ne0        = 2.50e8 ; cm^-3
  Te0        = 1.50e6 ; K

;--------------------------------------------------------------------------------------------------------------------------

  set_tomroot

  ; Default values for a few things
  if not keyword_set(i_measurement)    then i_measurement    = 0
  if not keyword_set(ion_label)        then ion_label        = 'fexiii' ; always use lowercase
  if not keyword_set(line_wavelength)  then line_wavelength  =  '10747' ; 5-character string with wavelength in A
  if not keyword_set(fip_factor)       then fip_factor       =     1.0  ; Note that [Fe] = [Fe]_Feldman * fip_factor
  if not keyword_set(instrument_label) then instrument_label =    'aia' ; always use lowercase
  if not keyword_set(band_label)       then band_label       =    '171'
  if not keyword_set(Ne0)              then Ne0              = 2.50e8   ; cm^-3
  if not keyword_set(Te0)              then Te0              = 1.50e6   ; K
  
  load_g_table,ion_label=ion_label,line_wavelength=line_wavelength,instrument_label=instrument_label,band_label=band_label

  help,T_e,N_e,G
  
  NTe = 1
  NNe = 1
  if (size(Te0))(0) eq 1 then NTe = (size(Te0))(1)
  if (size(Ne0))(0) eq 1 then NNe = (size(Ne0))(1)

; goto,ef
  
  print
  print,'Input values of Ne [cm^-3], Te [K], rad [Rsun], fip_factor:'
  print,Ne0,Te0,r0,fip_factor
  print
  print,'G-function value [erg(/PH) cm+3 sec-1]:'
  print,(g_function(Te0, Ne0))(0)
  print
  print,'s [erg(/PH) cm-3 sec-1 sr-1]:'
  print,(s_function(Ne0, Te0))(0)
  print
  print,'p [cm+3 K-1]:'
  print,p_function(Ne0, Te0)
  print
  print,'s*p [erg(/PH) sec-1 sr-1 K-1]:'
  print,sp_function(Ne0, Te0)
  print

  ef:
  
  Ne0_Limits = [min(N_e),max(N_e)]
  Te0_Limits = [min(T_e),max(T_e)]

  print,'e [erg sec-1 sr-1 K-1]:'
  print, e_function(parameters)
  print

  print,'e2[erg sec-1 sr-1 K-1]:'
  print, e2_function(parameters)
  print
  
  print, 'cost_function:'
  print, cost_function(parameters)
  print

  print,'grad_p'
  print,transpose(grad_p(Ne0,Te0))
  print

  print,'s*grad_P_i:'
  print,sgradp1_function(Ne0,Te0)
  print,sgradp2_function(Ne0,Te0)
  print,sgradp3_function(Ne0,Te0)
  print,sgradp4_function(Ne0,Te0)
  print,sgradp5_function(Ne0,Te0)
  print

  print,'grad_e'
  print,grad_e_function(parameters)
  print

  print,'grad_cost_function'
  print,grad_cost_function(parameters)
  print

  p   = parameters 
  dp  = 1.e-6* p

  dphi        = cost_function(p+dp) - cost_function(p)
  dphi_grad   = total(grad_cost_function(p)*dp)

  print,'Phi(p+dp) - Phi(p):',dphi
  print,'gradPhi *  dp:'     ,dphi_grad
  print,'relative difference:', abs( dphi - dphi_grad)/abs(dphi)
  stop
  return
end

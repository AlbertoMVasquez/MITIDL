;---------------------------------------------------------------------
;
; Brief description:
;
; Wrapper routine to test suite of codes.
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
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
  common tomographic_measurements, y0, y, measurement_type, i_measurement
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common weights,sig_WL,sig_v

  r0         = 1.2    ; Rsun
  fip_factor = 1.0    ; Note that [Fe] = [Fe]_Feldman * fip_factor
  Tem        = 1.75e6 ; K
  SigTe      = 0.50e6 ; K
  Nem        = 1.75e8 ; cm^-3
  SigNe      = 0.50e8 ; cm^-3
  q          = .5

  measurement_type = [1,2]
  if not keyword_set(euvband) AND not keyword_set(emissionline) then i_measurement = 0
  
  set_tomroot
  if     keyword_set(emissionline)     then i_measurement    = 0
  if     keyword_set(euvband)          then i_measurement    = 1
  if not keyword_set(ion_label)        then ion_label        = 'fexiii' ; always use lowercase
  if not keyword_set(line_wavelength)  then line_wavelength  =  '10747' ; 5-character string with wavelength in A
  if not keyword_set(fip_factor)       then fip_factor       =     1.0  ; Feldmand's Adundance Set value
  if not keyword_set(instrument_label) then instrument_label =    'aia' ; always use lowercase
  if not keyword_set(band_label)       then band_label       =    '171'
  
  load_g_table,ion_label=ion_label,line_wavelength=line_wavelength,instrument_label=instrument_label,band_label=band_label

  help,T_e,N_e,G
  
  if not keyword_set(Ne0) then Ne0=2.5e8
  if not keyword_set(Te0) then Te0=1.5e6

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

  parameters = [Nem, fip_factor, Tem, SigTe, SigNe, q]

  print,'e [erg sec-1 sr-1 K-1]:'
  print, e_function(parameters)
  print

  print,'e2[erg sec-1 sr-1 K-1]:'
  print, e2_function(parameters)
  print
  

  y0 = Nem * 1.1
  ;y  = [9.8e-10 , 280.5]
  ;measurement_type = [1,2]
  y = [15.e-10, 7.e-10, 250. , 500. ,300.]

  ; DUMMY TEST 
  ;y0=Nem
  ;y = [9.8808411e-10,   4.9120305e-10,       197.88246,       585.35546,       373.01257]
  
  f_wl = 0.1
  f    = 0.1 + y*0.

  sig_WL = f_wl* y0
  sig_v  = f   * y
  


  i_mea_vec=[0,0,1,1,1]
  ion_label_vec=       ['fexiii','fexiii','','','']
  line_wavelength_vec= ['10747','10801'  ,'','','']
  instrument_label_vec=['','','aia','aia','aia']
  band_label_vec=      ['','','171','193','211']

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
  dp  = 0.01* p

  dphi        = cost_function(p+dp) - cost_function(p)
  dphi_grad   = total(grad_cost_function(p)*dp)

  print,'Phi(p+dp) - Phi(p):',dphi
  print,'gradPhi *  dp:'     ,dphi_grad
  print,'relative difference:', abs( dphi - dphi_grad) /dphi
  
  return
end

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
  common tables,Te1,Te2,Te3,Te4,Te5,Ne1,Ne2,Ne3,Ne4,Ne5,G1,G2,G3,G4,G5,r1,r2
  common directories, tomroot
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common dimensions, NTe, NNe
  common NT_limits, Ne0_Limits, Te0_Limits
  common tomographic_measurements, y0, y
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common measurement_errors,sig_WL,sig_y
  common index_measurement, i_measurement
;---------------------------------------------------------------------------------------------------------------
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
  f_y  = 0.1 + fltarr(n_elements(y))

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
  ; The order chosen for its elements follows Rich's notes, right after Eq. (2)
  parameters = [Nem, fip_factor, Tem, SigTe, SigNe, q]
  
  ; Test values for Ne and Te at which evaluate several functions:
  Ne0        = 2.50e8 ; cm^-3
  Te0        = 1.50e6 ; K

;----------------------------------------------------------------------------------------------------------------
  set_tomroot

  load_tables


  ; Default values for a few things
                                            i_measurement    =       0
  if not keyword_set(ion_label)        then ion_label        = 'fexiii' ; always use lowercase
  if not keyword_set(line_wavelength)  then line_wavelength  =  '10747' ; 5-character string with wavelength in A
  if not keyword_set(fip_factor)       then fip_factor       =     1.0  ; [Fe] = [Fe]_Feldman * fip_factor
  if not keyword_set(instrument_label) then instrument_label =    'aia' ; always use lowercase
  if not keyword_set(band_label)       then band_label       =    '171'
  if not keyword_set(Ne0)              then Ne0              =  2.50e8  ; cm^-3
  if not keyword_set(Te0)              then Te0              =  1.50e6  ; K
  
  load_g_table,ion_label=ion_label,line_wavelength=line_wavelength,instrument_label=instrument_label,band_label=band_label

  help,T_e,N_e,G
  
  NTe = 1
  NNe = 1
  if (size(Te0))(0) eq 1 then NTe = (size(Te0))(1)
  if (size(Ne0))(0) eq 1 then NNe = (size(Ne0))(1)

  goto,ef
  
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

  goto,eg
  print,'e [erg sec-1 sr-1 K-1]:', e_function(parameters)
  print,'e2[erg sec-1 sr-1 K-1]:', e2_function(parameters)

  tstart     = systime(/seconds)
  print, 'cost_function: ', cost_function(parameters),'  Elapsed time:',systime(/seconds)-tstart
  
  print,'grad_p'
  print,transpose(grad_p_function(Ne0,Te0))
  print
    
  print,'s*grad_P_i:'
  print,[sgradp1_function(Ne0,Te0),sgradp1_function_new(Ne0,Te0)]
  print,[sgradp2_function(Ne0,Te0),sgradp2_function_new(Ne0,Te0)]
  print,[sgradp3_function(Ne0,Te0),sgradp3_function_new(Ne0,Te0)]
  print,[sgradp4_function(Ne0,Te0),sgradp4_function_new(Ne0,Te0)]
  print,[sgradp5_function(Ne0,Te0),sgradp5_function_new(Ne0,Te0)]
  print

  tstart     = systime(/seconds)  
  print,'grad_e_function: ',grad_e_function(parameters),'  Elapsed time:',systime(/seconds)-tstart
  print

  eg:
  tstart     = systime(/seconds)  
  print,'grad_cost_function:',grad_cost_function(parameters),'  Elapsed time:',systime(/seconds)-tstart
  print

  stop
  
  return
end

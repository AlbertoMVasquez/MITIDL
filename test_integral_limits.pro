;-----------------------------------------------------
; This routine does a comparison of the calculated
; emissivities (using INT2D) for a given set of parameters:
; [Nem, fip_factor, Tem, SigTe, SigNe, q]
; for two different integral limits.
;----------------------------------------------------


pro test_integral_limits

  common G_table, G, T_e, N_e, r, photT
  common tables,Te1,Te2,Te3,Te4,Te5,Ne1,Ne2,Ne3,Ne4,Ne5,G1,G2,G3,G4,G5,r1,r2
  common directories, tomroot
  common NT_limits, Ne0_Limits, Te0_Limits
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common index_measurement, i_measurement
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  
;---------------------------------------------------------------------------------------------------------------
;                                     VALUES TO PLAY WITH
  
  i_mea_vec           =[0       ,0       ,1    ,1    ,1    ]
  ion_label_vec       =['fexiii','fexiii',''   ,''   ,''   ]
  line_wavelength_vec =['10747' ,'10801' ,''   ,''   ,''   ]
  instrument_label_vec=[''      ,''      ,'aia','aia','aia']
  band_label_vec      =[''      ,''      ,'171','193','211']

  ; Test values for coronal heliocentric height and iron abundance:
  r0         = 1.1    ; Rsun
  fip_factor = 1.0    ; Note that [Fe] = [Fe]_Feldman * fip_factor

  ; Test values for the parameters of the joint bivariate Te-Ne normal distribution:
  Tem        = 1.30e6 ; K
  SigTe      = 0.50e6 ; K
  Nem        = 1.75e8 ; cm^-3
  SigNe      = 0.50e8 ; cm^-3
  q          = 0.5

  ; Parameter vector for both e_function and cost_function:
  ; The order chosen for its elements follows Rich's notes, right after Eq. (2)
  parameters = [Nem, fip_factor, Tem, SigTe, SigNe, q]
  

 ;----------------------------------------------------------------------------------------------------------------  
  
  set_tomroot
  load_tables
  
; Valores originales de Fede
   Ne0_Limits = [1.e6, 5.e9]
   Te0_Limits = [0.5e6,5.0e6]

   y1 = synth_y_values(parameters)

   print,'Ne limits:',Ne0_limits
   print,'Te limits:',Te0_limits
   print,'e_k values:',y1
   print
; Valores dinámicos razonables 
  Ne0_Limits = [max([min(Ne1),min(Ne2),min(Ne3),min(Ne4),min(Ne5)]),min([max(Ne1),max(Ne2),max(Ne3),max(Ne4),max(Ne5)])]
  Te0_Limits = [max([min(Te1),min(Te2),min(Te3),min(Te4),min(Te5)]),min([max(Te1),max(Te2),max(Te3),max(Te4),max(Te5)])]
  
  y2 = synth_y_values(parameters)

  print,'Ne limits:',Ne0_limits
  print,'Te limits:',Te0_limits
  print,'e_k values:',y2
  print

  print,'relative difference [%]',abs(y2-y1)/y2

  return
end




pro wrapper1
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q

 ; Test values for coronal heliocentric height and iron abundance:
  r0         = 1.11             ; Rsun

; Test values for the parameters of the joint bivariate Te-Ne normal distribution:
  Nem        = 0.61315012e8   ; cm^-3
  fip_factor = 1.0            ; Note that [Fe] = [Fe]_Feldman * fip_factor 
  Tem        = 1.4769326e6    ; K
  SigTe      = 0.3*Tem        ; K
  SigNe      = 0.3*Nem        ; cm^-3
  q          = 0.5  

  test_data,/riemann,/loguniform

  return
end


pro wrapper2
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q

 ; Test values for coronal heliocentric height and iron abundance:
  r0         = 1.21             ; Rsun

; Test values for the parameters of the joint bivariate Te-Ne normal distribution:
  Nem        = 0.27979984e8   ; cm^-3
  fip_factor = 1.0            ; Note that [Fe] = [Fe]_Feldman * fip_factor 
  Tem        = 1.5809454e6    ; K
  SigTe      = 0.3*Tem        ; K
  SigNe      = 0.3*Nem        ; cm^-3
  q          = 0.5  

  test_data,/riemann,/loguniform

  return
end



pro test_data,riemann=riemann,$
              uniform=uniform,$
              loguniform=loguniform,$
              lnuniform=lnuniform,$
              NNe_provided= NNe_provided,$
              NTe_provided= NTe_provided
             

  common G_table, G, T_e, N_e, r, photT
  common tables,Te1,Te2,Te3,Te4,Te5,Ne1,Ne2,Ne3,Ne4,Ne5,G1,G2,G3,G4,G5,r1,r2,r3,r4,r5
  common directories, tomroot
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common dimensions, NTe, NNe
  common NT_limits, Ne0_Limits, Te0_Limits
  common tomographic_measurements, y0, y
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common measurement_errors,sig_WL,sig_y
  common index_measurement, i_measurement
  common sk_over_fip_factor_array,sk_over_fip_factor
  common NT_arrays,Ne_array,Te_array,dNe_array,dTe_array,dTN  

  
;---------------------------------------------------------------------------------------------------------------
;                                     VALUES TO PLAY WITH
  
  i_mea_vec           =[0       ,0       ,1    ,1    ,1    ]
  ion_label_vec       =['fexiii','fexiii',''   ,''   ,''   ]
  line_wavelength_vec =['10747' ,'10801' ,''   ,''   ,''   ]
  instrument_label_vec=[''      ,''      ,'aia','aia','aia']
  band_label_vec      =[''      ,''      ,'171','193','211']

  
;----------------------------------------------------------------------------------------------------------------   
  set_tomroot
  load_tables
  
; Limits for the grid. These are DYNAMICAL as to adjust to future changes in G-tables:
  Ne0_Limits = [max([min(Ne1),min(Ne2),min(Ne3),min(Ne4),min(Ne5)]),min([max(Ne1),max(Ne2),max(Ne3),max(Ne4),max(Ne5)])]
  Te0_Limits = [max([min(Te1),min(Te2),min(Te3),min(Te4),min(Te5)]),min([max(Te1),max(Te2),max(Te3),max(Te4),max(Te5)])]

 ; restricted Ne and Te ranges 
 ;  Ne0_Limits = [1.0e6,5.0e9]
 ;  Te0_Limits = [0.5e6,5.0e6]
  
  if keyword_set(Riemann) then begin
      if not keyword_set(NNe_provided) then NNe_provided = 100
      if not keyword_set(NTe_provided) then NTe_provided = 100

      if keyword_set(   uniform) then make_grid,   /uniform,NNe_provided=NNe_provided,NTe_provided=NTe_provided
      if keyword_set(loguniform) then make_grid,/loguniform,NNe_provided=NNe_provided,NTe_provided=NTe_provided
      if keyword_set( lnuniform) then make_grid, /lnuniform,NNe_provided=NNe_provided,NTe_provided=NTe_provided
      if not keyword_set(uniform) and not keyword_set(loguniform) and not keyword_set(lnuniform) then begin
         print,'choose a grid (uniform, loguniform, or lnuniform)'
         return
      endif
      make_sk_over_fip_factor
   endif


   c= [Nem, fip_factor, Tem, SigTe, SigNe, q]*1.d
   y0 = double(Nem)     
   if     keyword_set(riemann) then y  = synth_y_values_cs(c)
   if NOT keyword_set(riemann) then y  = synth_y_values   (c)
   print,'synthetic values of y0 and y:',y0,y


  return
end

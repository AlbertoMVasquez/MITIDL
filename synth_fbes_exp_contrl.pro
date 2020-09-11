
pro wrapper

  synth_fbes_exp_contrl,/Riemann,/lnuniform,$
                        NNe_provided=50,NTe_provided=50
  return
end

;---------------------------------------------------------------------
; This routine calculate synthetic values of 171, 193 and 211 FBEs
; using a set of  parameters... 
;
;INPUT:

;
; KEYWORDS:
; Riemann: if keyword set the Riemann aproach is used.
; uniform : if keyword set the grid is uniform in Ne and Te
; loguniform: if keyword set the grid is uniform in log10Ne and
; log10Te
; lnuniform: if keyword set the grid is uniform in lnNe and lnTe.
; Also, use the Jacobian of the transformation to calculate dNe_array
; and dTe_array.
; NNe_provided: number of points in the Ne grid
; NTe_provided: number of points in the Te grid

pro synth_fbes_exp_contrl, Riemann=Riemann,$
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
 

; Test value for coronal heliocentric height:
 

  i_mea_vec           =[0       ,0       ,1    ,1    ,1    ]
  ion_label_vec       =['fexiii','fexiii',''   ,''   ,''   ]
  line_wavelength_vec =['10747' ,'10801' ,''   ,''   ,''   ]
  instrument_label_vec=[''      ,''      ,'aia','aia','aia']
  band_label_vec      =[''      ,''      ,'171','193','211']

  r0         = 1.11             ; Rsun
  Nm_demt = 0.61315012e8
  Tm_demt = 1.4769326e6
  wt_demt=  0.24922289e6

  fip_factor_v = [0.25,0.50,1.00]
  factor_sigma = [0.10,0.25,0.50]
  q_v          = [0.25,0.50,0.75]

  n1    = n_elements(fip_factor_v)
  n2    = n_elements(factor_sigma)
  n3    = n_elements(q_v)
  nband = 3
  x     = fltarr(n1,n2,n3,nband)
;----------------------------------------------------------------------------------------------------------------   


  for i1=0,n1-1 do begin
     for i2=0,n2-1 do begin
        for i3=0,n3-1 do begin

           set_tomroot
           load_tables
  
         ; Test values for the parameters of the joint bivariate Te-Ne normal distribution:
           Nem        = Nm_demt                ; cm^-3
           fip_factor = fip_factor_v(i1)       ; Note that [Fe] = [Fe]_Feldman * fip_factor 
           Tem        = Tm_demt                ; K
           SigTe      = factor_sigma(i2) * Tem ; K
           SigNe      = factor_sigma(i2) * Nem ; cm^-3
           q          = q_v         (i3)
  
         ; Parameter vector for both e_function and cost_function:
         ; The order chosen for its elements follows Rich's notes, right after Eq. (2)
           par_orig= [Nem, fip_factor, Tem, SigTe, SigNe, q]*1.d
           

        ;  Restricted Ne and Te ranges 
           Ne0_Limits = [1.0e6,5.0e9]
           Te0_Limits = [0.5e6,5.0e6]
  
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

          
           y  = synth_y_values_CS(par_orig) 
           x(i1,i2,i3,0)= y[2]
           x(i1,i2,i3,1)= y[3]
           x(i1,i2,i3,2)= y[4]
        endfor
     endfor
  endfor


  dir='~/Downloads/'
  datafiles=['x171','x193','x211']
  for i=0,nband-1 do begin
     openw,1,dir+datafiles(i)
     tmp = reform(x(*,*,*,i))
     writeu,1,tmp
     close,1
  endfor
 
  return
end

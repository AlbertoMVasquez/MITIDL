;----------------------------------------------------------------------
; Esta rutina sirve para testear la suite the codigos (_cs) que usan
; cuadratura simple para calcular las integrales dobles.
;
; CS: \int f(x,y) dx dy >  \Sum_{i,j} f(x_i,y_j) Dx Dy
;
; HISTORY
; V1.0 F.A. Nuevo, IAFE, April-2020
; V1.1 A.M. Vasquez, IAFE, April-2020
;
;----------------------------------------------------------------------

pro test_cuadratura

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
  common sk_array,sk       
  common NT_arrays,Ne_array,Te_array
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
  f_y  = 0.1 + findgen(n_elements(y))

  ; Absolute error of each measurement:
  sig_WL = f_wl* y0
  sig_y  = f_y * y
  
;----------------------------------------------------------------------------------------------------------------  
  ; Test values for coronal heliocentric height and iron abundance:
  r0         = 1.1    ; Rsun
  ; Initial value of fip factor to make the calculation of s_k matrix
  fip_factor = 1.0    ; Note that [Fe] = [Fe]_Feldman * fip_factor

  set_tomroot
  load_tables
  
  make_grid

 ; Calculate de S_k in the grid 
 ;make_sk_load_routine,sk
  make_sk             ,sk
;----------------------------------------------------------------------------------------------------------------
  ; Test values for the parameters of the joint bivariate Te-Ne normal distribution:
  Tem        = 1.30e6 ; K
  SigTe      = 0.50e6 ; K
  Nem        = 1.75e8 ; cm^-3
  SigNe      = 0.50e8 ; cm^-3
  q          = 0.5
  fip_factor = 1.1    ; Note that [Fe] = [Fe]_Feldman * fip_factor


  ; Parameter vector for both e_function and cost_function:
  ; The order chosen for its elements follows Rich's notes, right after Eq. (2)
  parameters = [Nem, fip_factor, Tem, SigTe, SigNe, q]


  
 ; make a comparison between the double integral calculated with INT2D and CS
  compare_integrals,parameters

 ; compare cost function and time calculated with the two SCHEMES :) 
 ; (no "SCENARIOS" Fede)
  tstart  = systime(/seconds)
  PHI1    = cost_function(parameters)
  tend    = systime(/seconds)
  print, 'cost_function [INT2D]:', PHI1, '   Elapsed time [sec]:', tend-tstart

  tstart  = systime(/seconds)
  PHI2    = cost_function_cs(parameters)
  tend    = systime(/seconds)
  print, 'cost_function (middle Riemann sum):', PHI2, '   Elapsed time [sec]:', tend-tstart

  print,'Relative diference [%]:',100.*abs(phi2-phi1)/phi1

  tstart = systime(/seconds)
  dphi1  = grad_cost_function(parameters)
  tend   = systime(/seconds)
  print,'grad_cost_function [INT2D]:',dphi1, '   Elapsed time [sec]:', tend-tstart

  tstart = systime(/seconds)
  dphi2  = grad_cost_function_cs(parameters)
  tend   = systime(/seconds)
  print,'grad_cost_function [middle Riemann sum]:',dphi2, '   Elapsed time [sec]:', tend-tstart

  print,'Relative diference [%]:', 100.* abs(dphi2-dphi1)/abs(dphi1)

  return
end



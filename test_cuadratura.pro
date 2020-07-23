;----------------------------------------------------------------------
; This routine tests the suite of codes (_cs) that uses CS to calculate
; the double integral. 
;
; CS: \int f(x,y) dx dy >  \Sum_{i,j} f(x_i,y_j) Dx Dy
;

; KEYWORDS;
; uniform : if keyword set the grid is uniform in Ne and Te
; loguniform: if keyword set the grid is uniform in log10Ne and
; log10Te
; lnuniform: if keyword set the grid is uniform in lnNe and lnTe.
; Also, use the Jacobian of the transformation to calculate dNe_array
; and dTe_array.
; NNe_provided: number of points in the Ne grid
; NTe_provided: number of points in the Te grid

; HISTORY
; V1.0 F.A. Nuevo, IAFE, April-2020
; V1.1 A.M. Vasquez, IAFE, May-2020
;----------------------------------------------------------------------

; test_cuadratura,/uniform,NNe_provided=50,NTe_provided=50
; test_cuadratura,/loguniform,NNe_provided=50,NTe_provided=50
; test_cuadratura,/lnuniform,NNe_provided=50,NTe_provided=50
; test_cuadratura,/lnuniform,NNe_provided=40,NTe_provided=20

pro test_cuadratura,uniform=uniform,lnuniform=lnuniform,loguniform=loguniform,$
                    NNe_provided=NNe_provided,NTe_provided=NTe_provided
  
  
  common tables,Te1,Te2,Te3,Te4,Te5,Ne1,Ne2,Ne3,Ne4,Ne5,G1,G2,G3,G4,G5,r1,r2,r3,r4,r5
  common directories, tomroot
  common dimensions, NTe, NNe
  common NT_limits, Ne0_Limits, Te0_Limits
  common tomographic_measurements, y0, y
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common measurement_errors,sig_WL,sig_y
  common index_measurement, i_measurement
  common sk_over_fip_factor_array,sk_over_fip_factor
  common NT_arrays,Ne_array,Te_array,dNe_array,dTe_array,dTN
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common constants, Rsun, kB, h, c
  common G_table, G, T_e, N_e, r, photT
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
 
 ; Test values for coronal heliocentric height and iron abundance:
  r0         = 1.1    ; Rsun

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
  
;----------------------------------------------------------------------------------------------------------------  
 
  set_tomroot
  load_tables

  
; Limits for the grid. These are DYNAMICAL as to adjust to future changes in G-tables:
  ;Ne0_Limits = [max([min(Ne1),min(Ne2),min(Ne3),min(Ne4),min(Ne5)]),min([max(Ne1),max(Ne2),max(Ne3),max(Ne4),max(Ne5)])]
  ;Te0_Limits = [max([min(Te1),min(Te2),min(Te3),min(Te4),min(Te5)]),min([max(Te1),max(Te2),max(Te3),max(Te4),max(Te5)])]

  ; restricted Ne and Te ranges
  Ne0_Limits = [1.0e6,5.0e9]
  Te0_Limits = [0.5e6,5.0e6]

  print,'INTEGRAL LIMITS:'
  print,'temp. range [K]   :',Te0_limits
  print,'dens. range [cm-3]:',Ne0_limits
  
  ; Make the Ne and Te grid 
  if not keyword_set(NNe_provided) then NNe_provided = 100
  if not keyword_set(NTe_provided) then NTe_provided = 100
  if keyword_set(   uniform) then make_grid,   /uniform,NNe_provided=NNe_provided,NTe_provided=NTe_provided
  if keyword_set(loguniform) then make_grid,/loguniform,NNe_provided=NNe_provided,NTe_provided=NTe_provided
  if keyword_set( lnuniform) then make_grid, /lnuniform,NNe_provided=NNe_provided,NTe_provided=NTe_provided
  if not keyword_set (uniform) and not keyword_set(loguniform) and not keyword_set (lnuniform) then begin
     print,'please choose a grid'
     return
  endif
; Compute S_k/fip_factor in the grid 
  make_sk_over_fip_factor
;----------------------------------------------------------------------------------------------------------------
  
 ; make a comparison between the double integral calculated with INT2D and CS
  compare_integrals,parameters

 ; compare cost function and time calculated with the two schemes
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



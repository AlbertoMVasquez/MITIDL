
;---------------------------------------------------------------------
pro wrapper,method
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common tomographic_measurements, y0, y
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  
  i_mea_vec           =[0       ,0       ,1    ,1    ,1    ]
  ion_label_vec       =['fexiii','fexiii',''   ,''   ,''   ]
  line_wavelength_vec =['10747' ,'10801' ,''   ,''   ,''   ]
  instrument_label_vec=[''      ,''      ,'aia','aia','aia']
  band_label_vec      =[''      ,''      ,'171','193','211']

  ;synthetic values of y0 and y:   1.5000000e+08   1.1396292e-09   5.1159609e-10       275.48046       781.21159       319.38916
  ; Todas estas variables se pasan por commons!
  y0 = double(1.5000000e+08)
  y  = double([1.1396292e-09,   5.1159609e-10,       275.48046,       781.21159,       319.38916])
  r0  = 1.1              ; Rsun
  fip_factor = 1.        ; necesario para pòder calcular s_k con s_function
  
  minimizador,min_method=method,/Riemann,/lnuniform,NNe_provided=50,NTe_provided=50

 
  return
end

;---------------------------------------------------------------------
; This routine calculate synthetic values of y0 and y using a set of
; parameters and minimize the cost_function. 
;
;INPUT:
; min_method: this variable selects the minimization method
     
; min_method=1: Downhill Simplex
; min_method=2: Powell          
; min_method=3: BFGS            
; min_method=4: Polak-Ribiere   
; min_method=5: constrained minimization

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

pro minimizador,min_method=min_method,$
                Riemann=Riemann,$
                uniform=uniform,$
                loguniform=loguniform,$
                lnuniform=lnuniform,$
                NNe_provided= NNe_provided,$
                NTe_provided= NTe_provided
  
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
  common sk_over_fip_factor_array,sk_over_fip_factor
  common NT_arrays,Ne_array,Te_array,dNe_array,dTe_array,dTN
  
  
  if keyword_set(Riemann) then begin 
     Phi_name     ='cost_function_cs'
     grad_Phi_name='grad_cost_function_cs'
  endif else begin
     Phi_name     ='cost_function'
     grad_Phi_name='grad_cost_function'
  endelse

  print,'Minimization of ',Phi_name
  print
  
  if not keyword_set(min_method) then begin
     print,'minimization method (min_method keyword) not selected:'
     print,'1: Downhill Simplex'
     print,'2: Powell          '
     print,'3: BFGS            '
     print,'4: Polak-Ribiere   '
     print,'5: minimization with constrains'
     return
  endif
 
;---------------------------------------------------------------------------------------------------------------
;                                     VALUES TO PLAY WITH
  


 ; heliocentric height and iron abundance:
   print,'heliocentric height:',r0,' Rsun'

  ;----------------------------------------------------------------------------------------------------------------   
  set_tomroot
  load_tables
  
; Limits for the grid. These are DYNAMICAL as to adjust to future changes in G-tables:
 ; Ne0_Limits = [max([min(Ne1),min(Ne2),min(Ne3),min(Ne4),min(Ne5)]),min([max(Ne1),max(Ne2),max(Ne3),max(Ne4),max(Ne5)])]
 ; Te0_Limits = [max([min(Te1),min(Te2),min(Te3),min(Te4),min(Te5)]),min([max(Te1),max(Te2),max(Te3),max(Te4),max(Te5)])]

 ; restricted Ne and Te ranges 
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



  ; Fractional error of each measurement:
   f_wl = 0.1
   f_y  = 0.1 + fltarr (n_elements(i_mea_vec))
  
  

   print,'values of y0 and y:',y0,y
   print

  ; Absolute error of each measurement:
   sig_WL = f_wl* y0 
   sig_y  = f_y * y  




; INITIAL GUESS:
   make_guess_ini,guess_ini,PHIguess

  ; pass Ne and Te to units of 10^8 cm-3 and MK
   change_units,guess_ini

  
   print,'------------------------------------------------------------------------------------------------'
   print,'          Nem          fip_factor         Tem             SigTe          SigNe           q      '
   print,'initial guess:'
   print,guess_ini
  
  
 
 
   ftol = 1.0d-4
   P = Guess_ini
   tstart     = systime(/seconds)

   if min_method eq 1 then begin
      print,'Downhill simplex Method'
     ;scale = [1.e8, 1., 1.e6, 1.e6, 1.e8, 1.]*0.5d
      scale = [1., 1., 1., 1., 1., 1.]*0.5d
      P = AMOEBA(ftol,scale=scale, P0 = guess_ini ,FUNCTION_VALUE=fval,function_name=Phi_name)
   endif
    
   if min_method eq 2 then begin
      print,'Powell  Method'
 
      xi = TRANSPOSE([[1., 0. ,0. ,0. ,0. ,0.],$
                      [0., 1. ,0. ,0. ,0. ,0.],$
                      [0., 0. ,1. ,0. ,0. ,0.],$
                      [0., 0. ,0. ,1. ,0. ,0.],$
                      [0., 0. ,0. ,0. ,1. ,0.],$
                      [0., 0. ,0. ,0. ,0. ,1.]])*1.d
      POWELL, P, xi, ftol, fmin, Phi_name
   endif
  
   if min_method eq 3  then begin
      print,'BFGS Method '
      DFPMIN, P, ftol, Fmin, Phi_name, Grad_Phi_name, /double, itmax=1000
   endif

   IF min_method eq 4 then begin
      print,'Polak-Ribiere Method '
      if not keyword_set(riemann) then begin
         command1='cp phi1.pro phi.pro'
         command2='cp gradphi1.pro gradphi.pro'
      endif else begin
         command1='cp phi2.pro phi.pro'
         command2='cp gradphi2.pro gradphi.pro'
      endelse
      print,command1 & spawn,command1
      print,command2 & spawn,command2
      print
      pr_min,guess_ini,out,phiv,ftol   
      P = OUT     
   endif

   IF min_method eq 5 then begin
      print,"Constrained Minimization"
      gcomp='cost_function_constr_min'
      xbnd=[[0.1, 0.5, 0.6, 0.01, 0.01, 0.01],[10., 2., 5., 2., 1., 0.95]]
      gbnd=[[0.1,0.],[10.,0.]]
     ;xbnd=[[1.e8, 0.5, 0.6e6, 0.01e6, 0.01e8, 0.01],[1.e9, 2., 5.e6, 2.e6, 1.e8, 0.95]]
     ;gbnd=[[1.e7,0.],[1.0e9,0.]]
      nobj=1
      CONSTRAINED_MIN, P, xbnd, gbnd, nobj, gcomp, inform
   ENDIF


  
  
   t_elapsed  = systime(/seconds)-tstart
  
   print,'------------------------------------------------------------------------------------------------'
   print,'          Nem          fip_factor         Tem             SigTe          SigNe           q      '
   print,'initial guess:'
   print,guess_ini
   print,'obtainded minimum:'
   print,p
   print
   
   print,'Phi(guess):', cost_function_cs(guess_ini)
   print,'Phi(min):',   cost_function_cs(P)
   print

  print,'Elapsed time:',t_elapsed
  print

  ysynth  = synth_y_values_cs(P) & ysynth  = [P(0),ysynth]
  yv      = [y0, y]
  score = mean(abs((yv - ysynth)/yv))
  print,'Score:',score

  return
end



pro change_units,guess,direction=direction
; IF direction= 1:  cm-3, K      > 10^8cm-3, MK
; IF direction=-1:  10^8cm-3, MK > cm-3, K
  common tomographic_measurements, y0, y
  common measurement_errors,sig_WL,sig_y
  common NT_arrays,Ne_array,Te_array,dNe_array,dTe_array,dTN
  

  if not keyword_set(direction) then direction=1

  if direction eq 1 then begin
 ; -----------------------------------------------------
     y0 = y0/1.d8               ; WL measure [10^8 cm-3]              
     sig_WL = sig_WL/1.d8       ; error of WL measure [10^8 cm-3]
  ;Initial Guess in new units 
     Guess[0]= Guess[0]/1.d8    ; Nm [10^8 cm-3]
     Guess[2]= Guess[2]/1.d6    ; Tm [MK]
     Guess[3]= Guess[3]/1.d6    ; sigT [MK]
     Guess[4]= Guess[4]/1.d8    ; sigN [10^8 cm-3]
  ; Ne and Te grid 
     Ne_array = Ne_array /1.d8
     Te_array = Te_array /1.d6 
     dNe_array= dNe_array/1.d8
     dTe_array= dTe_array/1.d6
     dTN      = dTN/1.d6 /1.d8
  endif

  if direction eq -1 then begin
     y0 = y0*1.d8               ; WL measure [10^8 cm-3]              
     sig_WL = sig_WL*1.d8       ; error of WL measure [10^8 cm-3]
  ;Initial Guess in new units 
     Guess[0]= Guess[0]*1.d8    ; Nm [10^8 cm-3]
     Guess[2]= Guess[2]*1.d6    ; Tm [MK]
     Guess[3]= Guess[3]*1.d6    ; sigT [MK]
     Guess[4]= Guess[4]*1.d8    ; sigN [10^8 cm-3]
  ; Ne and Te grid 
     Ne_array = Ne_array*1.d8
     Te_array = Te_array*1.d6 
     dNe_array= dNe_array*1.d8
     dTe_array= dTe_array*1.d6
     dTN      = dTN*1.d6*1.d8
  endif

; -----------------------------------------------------

  

  return
end

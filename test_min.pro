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

;---------------------------------------------------------------------
;test_min,min_method=1,/Riemann,/uniform
;test_min,min_method=1,/Riemann,/loguniform 

pro test_min,min_method=min_method,$
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
     return
  endif
 
;---------------------------------------------------------------------------------------------------------------
;                                     VALUES TO PLAY WITH
  
  i_mea_vec           =[0       ,0       ,1    ,1    ,1    ]
  ion_label_vec       =['fexiii','fexiii',''   ,''   ,''   ]
  line_wavelength_vec =['10747' ,'10801' ,''   ,''   ,''   ]
  instrument_label_vec=[''      ,''      ,'aia','aia','aia']
  band_label_vec      =[''      ,''      ,'171','193','211']


  

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
  par_orig= [Nem, fip_factor, Tem, SigTe, SigNe, q]*1.d
  
;----------------------------------------------------------------------------------------------------------------   
  set_tomroot
  load_tables
  
; Limits for the grid. These are DYNAMICAL as to adjust to future changes in G-tables:
  Ne0_Limits = [max([min(Ne1),min(Ne2),min(Ne3),min(Ne4),min(Ne5)]),min([max(Ne1),max(Ne2),max(Ne3),max(Ne4),max(Ne5)])]
  Te0_Limits = [max([min(Te1),min(Te2),min(Te3),min(Te4),min(Te5)]),min([max(Te1),max(Te2),max(Te3),max(Te4),max(Te5)])]

 ; restricted Ne and Te ranges 
 ; Ne0_Limits = [1.0e6,5.0e9]
 ; Te0_Limits = [0.5e6,5.0e6]

  ; synthetic values of y0 and y
  ;y0 = 1.3*Nem
  y0 =     Nem
  y = synth_y_values(par_orig)

  ; Fractional error of each measurement:
  f_wl = 0.1
  f_y  = 0.1 + findgen(n_elements(i_mea_vec))

  
  ; Absolute error of each measurement:
  sig_WL = f_wl* y0
  sig_y  = f_y * y
  

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

  goto,skiptest
  Nc=10
  cmax = 1.5
  cmin = 0.5
  c = cmin+(cmax-cmin)*findgen(Nc+1)/float(Nc)
  phi1v = dblarr(Nc+1)
  phi2v = dblarr(Nc+1)
  for ic=0,Nc do begin
     phi1v[ic] = cost_function   (c[ic]*par_orig)
     phi2v[ic] = cost_function_cs(c[ic]*par_orig)
  endfor
  window,0
  plot,c,phi1v
  oplot,c,phi2v,psym=4
  stop
  skiptest:


  ftol = 1.0e-4
  Guess_ini = 0.8d * par_orig
  ;Guess_ini = [0.8d,1.6d,0.9d,0.5,1.2d,0.5d] * par_orig
  P = Guess_ini
  tstart     = systime(/seconds)

  if min_method eq 1 then begin
     print,'Downhill simplex Method'
     scale = [1.e8, 1., 1.e6, 1.e6, 1.e8, 1.]*0.5d
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
     DFPMIN, P, ftol, Fmin, Phi_name, Grad_Phi_name
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

  
  
  t_elapsed  = systime(/seconds)-tstart
  
  print,'initial guess:'
  print,guess_ini
  print,'obtainded minimum:'
  print,p
  print,'original values:'
  print,par_orig
  print
  print,'relative difference:',abs((P-par_orig)/par_orig)
  print

  print,'Phi(guess):', cost_function(guess_ini)
  print,'Phi(min):',   cost_function(P)
  print,'Phi(orig):',  cost_function(par_orig)
  print

  print,'Elapsed time:',t_elapsed
  
  
  return
end

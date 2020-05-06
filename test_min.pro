;---------------------------------------------------------------------


;---------------------------------------------------------------------
pro test_min,min_method=min_method
  
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
  common sk_over_fip_factor_array,sk_over_fip_factor
  common NT_arrays,Ne_array,Te_array,dNe_array,dTe_array,dTN
  
  if not keyword_set(min_method) then begin
     print,'minimization method (min_method keyword) not selected:'
     print,'1: Downhill Simplex'
     print,'2: Powell          '
     print,'3: BFGS            '
     print,'4: Polak-Ribiere   '
     return
  endif
  

  set_tomroot
  i_mea_vec           = [0,0,1,1,1]
  ion_label_vec       = ['fexiii','fexiii','','','']
  line_wavelength_vec = ['10747','10801'  ,'','','']
  instrument_label_vec= ['','','aia','aia','aia']
  band_label_vec      = ['','','171','193','211']

  
  ; Original values of parameters used to calculate synth emmisivities
  r0         = 1.2              ; Rsun                                                                                 
  fip_factor = 1.0              ; Note that [Fe] = [Fe]_Feldman * fip_factor        
  Tem        = 1.75e6           ; K                                               
  SigTe      = 0.50e6           ; K                                                                      
  Nem        = 1.75e8           ; cm^-3                                                                              
  SigNe      = 0.50e8           ; cm^-3                                                                         
  q          =.5
  par_orig   = [Nem, fip_factor, Tem, SigTe, SigNe, q]
  
  ; integral limits
  Ne0_Limits = [1.e6 ,5.0e9]
  Te0_Limits = [0.5e6,5.0e6]
  ; synthetic values of y0 and y
  y0= Nem
  y = synth_y_values(par_orig)
  
  ; Fractional error of each measurement:
  f_wl = 0.1
  f_y  = 0.1 + findgen(n_elements(y))

  ; Absolute error of each measurement:
  sig_WL = f_wl* y0
  sig_y  = f_y * y

   
  ; Ne and Te grid
  make_grid,/uniform,NNe_provided=80,NTe_provided=100 
  make_sk_over_fip_factor
  
  ; initial guess
  guess_ini= 0.8*par_orig

  ;print,cost_function(par_orig)
  ;print,cost_function(guess_ini)
  ;return

  ftol = 1.0e-4
  P = Guess_ini

  tstart     = systime(/seconds)

  if min_method eq 1 then begin
     print,'Downhill simplex Method'
     scale = [1.e8, 1., 1.e6, 1.e6, 1.e8, 1.]
     P = AMOEBA(ftol,scale=scale, P0 = guess_ini ,FUNCTION_VALUE=fval,function_name='cost_function_cs')
  endif
  

  
  if min_method eq 2 then begin
     print,'Powell  Method'
     xi = TRANSPOSE([[1., 0. ,0. ,0. ,0. ,0.],$
                     [0., 1. ,0. ,0. ,0. ,0.],$
                     [0., 0. ,1. ,0. ,0. ,0.],$
                     [0., 0. ,0. ,1. ,0. ,0.],$
                     [0., 0. ,0. ,0. ,1. ,0.],$
                     [0., 0. ,0. ,0. ,0. ,1.]])
     POWELL, P, xi, ftol, fmin, 'cost_function_cs'
  endif
  
  if min_method eq 3  then begin
     print,'BFGS Method '
     DFPMIN, P, ftol, Fmin, 'cost_function_cs', 'grad_cost_function_cs'
  endif

  if min_method eq 4 then begin
     PR_min,guess_ini,out,phiv
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
  
  print,'Phi(guess):', cost_function(guess_ini)
  print,'Phi(min):',   cost_function(P)
  print,'Phi(orig):',  cost_function(par_orig)
  print

  print,'Elapsed time:',t_elapsed
  
  
  return
end

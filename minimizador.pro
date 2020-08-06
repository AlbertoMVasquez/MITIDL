
pro minimizador,Phi_name,grad_phi_name,guess_ini,P,min_method=min_method
  
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
  

   ftol = 1.0d-4
   P = Guess_ini


   if min_method eq 1 then begin
  ;    print,'Downhill simplex Method'
     ;scale = [1.e8, 1., 1.e6, 1.e6, 1.e8, 1.]*0.5d
      scale = [1., 1., 1., 1., 1., 1.]*0.5d
      P = AMOEBA(ftol,scale=scale, P0 = guess_ini ,FUNCTION_VALUE=fval,function_name=Phi_name)
      if P eq -1 then P = fltarr(6) -666.
   endif
    
   if min_method eq 2 then begin
   ;   print,'Powell  Method'
 
      xi = TRANSPOSE([[1., 0. ,0. ,0. ,0. ,0.],$
                      [0., 1. ,0. ,0. ,0. ,0.],$
                      [0., 0. ,1. ,0. ,0. ,0.],$
                      [0., 0. ,0. ,1. ,0. ,0.],$
                      [0., 0. ,0. ,0. ,1. ,0.],$
                      [0., 0. ,0. ,0. ,0. ,1.]])*1.d
      POWELL, P, xi, ftol, fmin, Phi_name
   endif
  
   if min_method eq 3  then begin
   ;   print,'BFGS Method '
      DFPMIN, P, ftol, Fmin, Phi_name, Grad_Phi_name, /double, itmax=1000
   endif

   IF min_method eq 4 then begin
   ;   print,'Polak-Ribiere Method '
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
   return
end



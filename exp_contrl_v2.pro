
pro wrapper
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common NT_limits, Ne0_Limits, Te0_Limits

  ;This vectors define the order of the measurement
  i_mea_vec           =[0       ,0       ,1    ,1    ,1    ]
  ion_label_vec       =['fexiii','fexiii',''   ,''   ,''   ]
  line_wavelength_vec =['10747' ,'10801' ,''   ,''   ,''   ]
  instrument_label_vec=[''      ,''      ,'aia','aia','aia']
  band_label_vec      =[''      ,''      ,'171','193','211']
  ;File name of x-tomographic products of all instruments 
  ;(need to be consistent with above)
  xfiles   = ['',$
              '',$
              '',$
              '',$
              '',$
              '']

 file_par_in = ''

 
; Restricted Ne and Te ranges 
  Ne0_Limits = [1.0e6,5.0e9]
  Te0_Limits = [0.5e6,5.0e6]

  method=4 ; Polak-Ribiere method

  file_demt='ldem2198.out'
  file_out ='mit.out'
  dir      ='~/Downloads/'
  exp_contrl_v2,xfiles,/Riemann,min_method=method,$
                file_demt=file_demt,file_out=file_out,dir=dir,$
                /lnuniform,NNe_provided=50,NTe_provided=50

  return
end


;================================================================
; Master code to Multi-Instrument Tomography (MIT)

;===============================================================

pro exp_contrl_v2,rmin,rmax,xfiles,min_method=min_method,riemann=riemann,$
        uniform=uniform,loguniform=loguniform,lnuniform=lnuniform,$
        NNe_provided=NNe_provided,NTe_provided=NTe_provided,$
        file_demt=file_demt,file_out=file_out,dir=dir

  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common NT_limits, Ne0_Limits, Te0_Limits
  common NT_arrays,Ne_array,Te_array,dNe_array,dTe_array,dTN
  common directories, tomroot
  common tomographic_measurements, y0, y  
  common measurement_errors,sig_WL,sig_y
  common sk_over_fip_factor_array,sk_over_fip_factor


  if not keyword_set(min_method) then begin
     print,'minimization method (min_method keyword) not selected:'
     print,'1: Downhill Simplex'
     print,'2: Powell          '
     print,'3: BFGS            '
     print,'4: Polak-Ribiere   '
     print,'5: minimization with constrains'
     return
  endif
  
  if keyword_set(Riemann) then begin 
     Phi_name     ='cost_function_cs'
     grad_Phi_name='grad_cost_function_cs'
  endif else begin
     Phi_name     ='cost_function'
     grad_Phi_name='grad_cost_function'
  endelse
  
  print,'Minimization of ',Phi_name
  print

  if min_method eq 1 then print,'Downhill simplex Method'
  if min_method eq 2 then print,'Powell method'
  if min_method eq 3 then print,'BFGS method'
  IF min_method eq 4 then begin
     PRINT,'Polak-Ribiere method'
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
  endif
  IF min_method eq 5 then print,"Constrained Minimization"

  
  
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
  endif


  restore,dir_par_in+file_par_in
  sz = size(par_in)
  nr = sz(1)
  nt = sz(2)
  np = sz(3)


; Load tomographies into memory
  xread,file=xfiles(0),nr=nr,nt=nt,np=np,map=KCOR
  xread,file=xfiles(1),nr=nr,nt=nt,np=np,map=comp1074
  xread,file=xfiles(2),nr=nr,nt=nt,np=np,map=comp1079
  xread,file=xfiles(3),nr=nr,nt=nt,np=np,map=FBE171
  xread,file=xfiles(4),nr=nr,nt=nt,np=np,map=FBE193
  xread,file=xfiles(5),nr=nr,nt=nt,np=np,map=FBE211
 


  restore,'~/Downloads/ldem_mit.out'
  ndat=(size(demt_A))(4)
  nm_demt_array=reform(demt_A(*,*,*,ndat-4))
  tm_demt_array=reform(demt_A(*,*,*,ndat-3))
  wt_demt_array=reform(demt_A(*,*,*,ndat-2))

  Npar      = 6
  Ndat      = 6
  parA      = dblarr(Nr,Nt,Np,npar) -666.
  yA        = dblarr(Nr,Nt,Np,ndat) 
  ysynthA   = dblarr(Nr,Nt,Np,ndat) 
;--------------------------------


   set_tomroot
   load_tables
   fip_factor = 1.0
   make_sk_over_fip_factor

   change_units_grid
   

   irad1=0                      
   irad2=Nr-1                   
   Dirad=1
   ilat1=0
   ilat2=nt-1
   Dilat=1
   ilon1=0
   ilon2=np-1
   Dilon=1
  ;--------------------------
  
   for ir =irad1,irad2,Dirad do begin
      for ith=ilat1,ilat2,Dilat do begin
         for ip =ilon1,ilon2,Dilon do begin
          ; Load y0 and y in the voxel
            y0 = kcor(ir,ith,ip) /1.e8
            y  =[comp1074(ir,ith,ip),comp1079(ir,ith,ip),$
                 FBE171(ir,ith,ip),FBE193(ir,ith,ip),FBE211(ir,ith,ip)]
          ; Fractional error of each measurement:
            f_wl = 0.1 
            f_y  = 0.1 + fltarr (n_elements(i_mea_vec))
          ; Absolute error of each measurement:
            sig_WL = f_wl* y0 
            sig_y  = f_y * y  

            if min([y0,y]) le 0. then goto,skipmin ; SKIP ZDAs

            ;INITIAL GUESS:
            ;make_guess_ini_new_units,guess_ini,PHIguess
            nm_demt = nm_demt_array(ir,ith,ip)
            tm_demt = tm_demt_array(ir,ith,ip)
            wt_demt = wt_demt_array(ir,ith,ip)
            make_guess_ini_with_demt,nm_demt,tm_demt,wt_demt,guess_ini,PHIguess


           ;MINIMIZATION BLOCK
            tstart     = systime(/seconds)
            minimizador,phi_name,grad_phi_name,guess_ini,P,min_method=min_method
            t_elapsed  = systime(/seconds)-tstart
            print,'Elapsed time:',t_elapsed
            parA(ir,ith,ip,*) = P ; save the parameters vector in a 3D array
            ysynth  = synth_y_values_cs(P) & ysynth  = [P(0),ysynth]
            yv      = [y0, y]
            score = mean(abs((yv - ysynth)/yv))
            print,'Score R:',score
            print,'R_k    :',abs((yv - ysynth)/yv)
            
            skipmin:             
         endfor                 ; IP  loop
      endfor                    ; ITH loop
   endfor                       ; IR loop

;--------------------------------

   save,filename=dir+fileout,rat,lat,lon,parA,yA,ysynthA

  
   RETURN
END

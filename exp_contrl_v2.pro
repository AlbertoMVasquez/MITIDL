; COMENTARIOS:
; si flag_noise=0 EL OUTPUT se guarda en /data1/DATA/MIT/dir_out/sin_ruido/ 
; si flag_noise=1 EL OUTPUT se guarda en /data1/DATA/MIT/dir_out/con_ruido/ 
; para modificar ftol "global" > L18 minimizador.pro
; para modificar TOL de dbrent > L8 linmin.pro

pro wrapper
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common NT_limits, Ne0_Limits, Te0_Limits

;--------------------------------------------------------

 ;Suffixs to make the mame of x-tomographic products 
 ;used as input file
  ;noise_suffix = 'sin_ruido'      & flag_noise=0
  noise_suffix = 'con_ruido_0.01'  & flag_noise=1
  exp_suffix   = ''               
 ;exp_suffix   = 'exp_B'

 ;Name of output file  
  file_out ='mit_exp_contrl_'+noise_suffix+exp_suffix+'.out';'mit_exp_contrl.out'
 ;Directory where is writed the output file
  dir_out  ='exp_contrl_v2_ftol1e-8TOL1e-6'
 
 
 ;Use this flag to select the minimization method
 ;method=1  ; Downhill-Simplex method (AMOEBA)
 ;method=2  ; Powell method
 ;method=3  ; BFGS Method
  method=4  ; Polak-Ribiere method 
;------------------------------------------------------  


 ;This vectors define the order of the measurement
  i_mea_vec           =[0       ,0       ,1    ,1    ,1    ]
  ion_label_vec       =['fexiii','fexiii',''   ,''   ,''   ]
  line_wavelength_vec =['10747' ,'10801' ,''   ,''   ,''   ]
  instrument_label_vec=[''      ,''      ,'aia','aia','aia']
  band_label_vec      =[''      ,''      ,'171','193','211']


 ;File names of x-tomographic products of all instruments 
 ;In /data1/tomography/bindata/
 ;(need to be consistent with above)
  xfiles      =['xkcor',$
                'x.comp1074',$
                'x.comp1079',$
                'x.aia171',$
                'x.aia193',$
                'x.aia211']+'_exp_contrl_'+noise_suffix+exp_suffix+'.out'
  
 ;Input parameters used to calculate the synthetic x-tomographic products
  file_par_in ='param_input_exp_contrl_'  +noise_suffix+exp_suffix+'.out' 
 ;LDEM file to use as initial guess in the minimization
  file_demt   ='LDEM_AIA3_MIT_exp_contrl_'+noise_suffix+exp_suffix+'.sav'

 ;Restricted Ne and Te ranges 
  Ne0_Limits = [1.0e6,5.0e9]
  Te0_Limits = [0.5e6,5.0e6]
  
  exp_contrl_v2,xfiles,/Riemann,min_method=method,$
                file_demt=file_demt,file_out=file_out,dir_out=dir_out,$
                file_par_in=file_par_in,flag_noise=flag_noise,$
                /lnuniform,NNe_provided=50,NTe_provided=50
   return
end


;================================================================
; This routine calculate the MIT inversion using the 
; synthetic x-tomographic products (see synth_x_exp_contrl.pro)
;===============================================================

pro exp_contrl_v2,xfiles,min_method=min_method,riemann=riemann,$
                  uniform=uniform,loguniform=loguniform,lnuniform=lnuniform,$
                  NNe_provided=NNe_provided,NTe_provided=NTe_provided,$
                  file_demt=file_demt,file_out=file_out,dir_out=dir_out,$
                  file_par_in=file_par_in,flag_noise=flag_noise

  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common NT_limits, Ne0_Limits, Te0_Limits
  common NT_arrays,Ne_array,Te_array,dNe_array,dTe_array,dTN
  common directories, tomroot
  common tomographic_measurements, y0, y  
  common measurement_errors,sig_WL,sig_y
  common sk_over_fip_factor_array,sk_over_fip_factor
  common tables,Te1,Te2,Te3,Te4,Te5,Ne1,Ne2,Ne3,Ne4,Ne5,G1,G2,G3,G4,G5,r1,r2,r3,r4,r5
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common units,ne_unit,te_unit
  
  tstart_whole_exp     = systime(/seconds) ; to compute the time of all experiments

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

  dir_par_in = '/data1/tomography/bindata/' ; directory where is the input parameters file
  restore,dir_par_in+file_par_in ; restore the input parameters file
  sz = size(par_in) 
  nr = sz(1) 
  nt = sz(2)
  np = sz(3)


; Load synthetic tomographies into memory
  xread,file=xfiles(0),nr=nr,nt=nt,np=np,map=KCOR
  xread,file=xfiles(1),nr=nr,nt=nt,np=np,map=comp1074
  xread,file=xfiles(2),nr=nr,nt=nt,np=np,map=comp1079
  xread,file=xfiles(3),nr=nr,nt=nt,np=np,map=FBE171
  xread,file=xfiles(4),nr=nr,nt=nt,np=np,map=FBE193
  xread,file=xfiles(5),nr=nr,nt=nt,np=np,map=FBE211
 

; Restore DEMT file to use as initial guess
  restore,'/data1/DATA/ldem_files/'+file_demt 
  ndat=(size(demt_A))(4)
  nm_demt_array=reform(demt_A(*,*,*,ndat-4))
  tm_demt_array=reform(demt_A(*,*,*,ndat-3))
  wt_demt_array=reform(demt_A(*,*,*,ndat-2))

  Npar      = 6 ; number of parameters
  Ndat      = 6 ; number of measures
  parA      = dblarr(Nr,Nt,Np,npar) -666. ;output parameters
  yA        = dblarr(Nr,Nt,Np,ndat)       
  ysynthA   = dblarr(Nr,Nt,Np,ndat) 
  scoreR    = dblarr(Nr,Nt,Np)
  scoreRk   = dblarr(Nr,Nt,Np,ndat) 
;--------------------------------


   set_tomroot
   ; load G tables in memory
   load_tables
   ; make the sk array
   fip_factor = 1.0
   make_sk_over_fip_factor
   
   ; change units of Ne and Te grid
   ; to e8 cm-3 and MK.
   load_units
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
  ; Triple lazo de experimentos
  ;--------------------------
   for ir =irad1,irad2,Dirad do begin
      for ith=ilat1,ilat2,Dilat do begin
         for ip =ilon1,ilon2,Dilon do begin
          ; Load y0 and y in the voxel
            y0 = kcor(ir,ith,ip) /ne_unit
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
            nm_demt = nm_demt_array(ir,ith,ip)/ne_unit
            tm_demt = tm_demt_array(ir,ith,ip)/te_unit
            wt_demt = wt_demt_array(ir,ith,ip)/te_unit
            make_guess_ini_with_demt,nm_demt,tm_demt,wt_demt,guess_ini,PHIguess
            
            ;Input parameters in the voxel
            P_in=par_in(ir,ith,ip,*)/[ne_unit, 1., te_unit, te_unit, ne_unit, 1.]
            

            ;==========================================
            print,'Minimization...'
           ;MINIMIZATION BLOCK
            tstart     = systime(/seconds)
            minimizador,phi_name,grad_phi_name,guess_ini,P,min_method=min_method
            t_elapsed  = systime(/seconds)-tstart
            print,'Elapsed time:',t_elapsed
            parA(ir,ith,ip,*) = P ; save the parameters vector in a 3D array
            ysynth  = synth_y_values_cs(P) & ysynth  = [P(0),ysynth]
            yv      = [y0, y]
            Rk      = abs((yv - ysynth)/yv) ; score R_k
            score   = mean(Rk)              ; score R
            ; save the scores R and R_k in a 3D array
            scoreR (ir,ith,ip  ) = score
            scoreRk(ir,ith,ip,*) = Rk
            ;==========================================
            print,'Score R:',score
            print,'R_k    :',Rk
            print
            print,'in:',TRANSPOSE(P_in)
            print,'(out-in)/in:',(P-P_in)/P_in
            ;stop

            skipmin:             
         endfor                 ; IP  loop
      endfor                    ; ITH loop
   endfor                       ; IR loop
;--------------------------------

   par_out = parA

  
   dir_root ='/data1/DATA/MIT/'
   if flag_noise eq 0 then dir_out=dir_out+'/sin_ruido/'
   if flag_noise eq 1 then dir_out=dir_out+'/con_ruido/'
   spawn,'mkdir -p '+dir_root+dir_out

   save,filename=dir_root+dir_out+file_out,par_in,par_out,scoreR,scoreRk

   print
   print
   t_elapsed_whole_exp  = systime(/seconds)-tstart_whole_exp
   print,'Total elapsed time (whole experiment):',t_elapsed_whole_exp

   print
   print,'Results are in: ',dir_root+dir_out+file_out


   RETURN
END

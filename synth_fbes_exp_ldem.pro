
pro synth_fbes_exp_ldem,noise=noise,exp_suffix=exp_suffix

  f_suffix='_0.1'
  f       =  0.1

  if not keyword_set(noise) then noise_suffix='sin_ruido'
  if     keyword_set(noise) then noise_suffix='con_ruido'+f_suffix
  if not keyword_set(exp_suffix) then exp_suffix=''


; AIA channels
  instrument_label_vec=['aia','aia','aia']
  band_label_vec      =['171','193','211']
  nband               = n_elements(band_label_vec)  
; Temperature grid
  Tmin = 0.5 ; MK
  Tmax = 5.0 ; MK
  NT   = 500
  dT   = (Tmax - Tmin)/NT ; MK
  T    = Tmin + findgen(NT)*dT + dT/2 ; MK

; ARRAY of Nband X NT  with the Nband TRFs in the Temp. Grid 
  Q_k =  fltarr(nband,NT)
  for iband = 0, nband-1 do begin
     instrument_label = instrument_label_vec(iband)
     band_label       = band_label_vec      (iband)
     read_and_interpol_TRF,instrument_label,band_label,T,TRF
     Q_k (iband,*) = TRF
  endfor
  
 ; Parameter space to explore in the experiment
  npar       = 3
  Ne0_v      = [1.0e8,3.0e8,5.0e8]         & n1 = n_elements(Ne0_v) 
  Tc_v       = [1.1,1.5,2.0]               & n2 = n_elements(Tc_v) 
  factor_sigT= [0.1,0.25,0.5]              & n3 = n_elements(factor_sigT) 
  
  ; ARRAY WITH THE SYNTH FBEs
  x_fbe     = fltarr(n1,n2,n3,nband)
  ; ARRAY WITH THE INPUT PARAMETERS
  par_in    = fltarr(n1,n2,n3,npar)


  for i = 0, n1-1 do begin
     for j= 0, n2-1 do begin
        for k = 0, n3-1 do begin
           
           Ne0 = Ne0_v (i)
           Tc  = Tc_v  (j)
           sigT= factor_sigT(k)*Tc
           
           LDEM = (Ne0^2/sqrt(2.d0*!pi)/sigT)*exp(-0.5d0*((T-Tc)/sigT)^2)
           ;plot,T,LDEM,xtitle='T [MK]',ytitle='LDEM [cm-6 MK]'
           ;STOP
           FBE  = (dT*LDEM)##Q_k
           IF Keyword_set(noise) then begin
              f_v  = f + fltarr (nband)
              FBE  = FBE  * (1.0 + f_v *randomn(seed,nband))
           endif
           ;FBE2 = FBE*0. & FBE3 = FBE*0.
           ;for iband = 0,nband-1 do begin
           ;   FBE2[iband] = total(reform(Q_k(iband,*))*LDEM *dT)
           ;   FBE3[iband] = INT_TABULATED(T,reform(Q_k(iband,*))*LDEM)
           ;endfor
           par_IN  (i,j,k,*) = [Ne0,Tc,sigT]
           x_fbe   (i,j,k,*) = FBE
           ;stop
        endfor
     endfor
  endfor
  


  dir = '/data1/tomography/bindata/'
  datafiles=['x.aia171','x.aia193','x.aia211']$
            +'_exp_contrl_ldem_'+noise_suffix+exp_suffix+'.out'
  for i=0,nband-1 do begin
     openw,1,dir+datafiles(i)
     tmp = reform(x_fbe(*,*,*,i))
     writeu,1,tmp
     close,1
  endfor

  file_par_input = 'param_input'+'_exp_contrl_ldem_'+noise_suffix+exp_suffix+'.out'
  save,filename=dir+file_par_input,par_IN


  return
end


pro read_and_interpol_TRF,instrument_label,band_label,T,TRF
  
  data_dir = '/data1/tomography/MITIDL/Emissivity_LookUp_Tables/'
  xstring = ''
  file_name = 'TRF_function_'+instrument_label+'_'+band_label+'.txt'
  openr,1,data_dir+file_name
  for i=1,8 do readf,1,xstring
  x=0.
  Ntemp=0
  readf,1,x,Ntemp,x,x,x,x
  for i=1,3 do readf,1,xstring
  logTe =fltarr(Ntemp)
  TRF   =fltarr(Ntemp)
  logTe0=0.
  TRF0  =0.
  for itemp=0,Ntemp-1 do begin
     readf,1,logTe0,TRF0,x
     logTe[itemp]=logTe0
     TRF[itemp]=TRF0
  endfor
  close,1
  T_e   = 10.^logTe             ; [K]

  TRF_interp  = interpol( TRF , T_e, T*1.e6)
  ;plot,T_e/1.e6,TRF,xrange=[0.5,5.0]
  ;oplot,T,TRF_interp,psym=4
  ;stop
  TRF = TRF_interp

  return
end

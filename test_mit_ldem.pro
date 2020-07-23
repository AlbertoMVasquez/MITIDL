 

pro test_mit_ldem
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
  
  i_mea_vec           =[0       ,0       ,1    ,1    ,1    ]
  ion_label_vec       =['fexiii','fexiii',''   ,''   ,''   ]
  line_wavelength_vec =['10747' ,'10801' ,''   ,''   ,''   ]
  instrument_label_vec=[''      ,''      ,'aia','aia','aia']
  band_label_vec      =[''      ,''      ,'171','193','211']


 ; Test values for coronal heliocentric height and iron abundance:
  r0         = 1.1              ; Rsun

 ; Test values for the parameters of the joint bivariate Te-Ne normal distribution:
  Nem        = 1.75e8           ; cm^-3
  fip_factor = 1.1              ; Note that [Fe] = [Fe]_Feldman * fip_factor 
  Tem        = 1.30e6           ; K
  SigTe      = 0.50e6           ; K
  SigNe      = 0.50e8           ; cm^-3
  q          = 0.5
  

;----------------------------------------------------------------------------------------------------------------   
  set_tomroot
  load_tables
  
; Limits for the grid. These are DYNAMICAL as to adjust to future changes in G-tables:
  Ne0_Limits = [max([min(Ne1),min(Ne2),min(Ne3),min(Ne4),min(Ne5)]),min([max(Ne1),max(Ne2),max(Ne3),max(Ne4),max(Ne5)])]
  Te0_Limits = [max([min(Te1),min(Te2),min(Te3),min(Te4),min(Te5)]),min([max(Te1),max(Te2),max(Te3),max(Te4),max(Te5)])]

 ; restricted Ne and Te ranges 
 ; Ne0_Limits = [1.0e6,5.0e9]
 ; Te0_Limits = [0.5e6,5.0e6]


 ; synthetic values of y0 and y,     
 ; using the fractional errors above to simulate (normal) uncertainty of measurement.
   
   


  nr =3
  nt =3
  np =3
  x = fltarr(nr,nt,np,3)
  
  Nmmin= 1.0e8
  Nmmax= 3.0e8
  Tmmin= 1.0e6
  Tmmax= 2.0e6
  sTmin= 0.5e6
  sTmax= 1.0e6

  Nm_vector = Nmmin + (Nmmax-Nmmin)*findgen(nr)/float(nr-1) 
  Tm_vector = Tmmin + (Tmmax-Tmmin)*findgen(nt)/float(nt-1) 
  sT_vector = sTmin + (sTmax-sTmin)*findgen(np)/float(np-1) 
  
  NmA = fltarr(nr,nt,np)
  TmA = fltarr(nr,nt,np)
  STA = fltarr(nr,nt,np)

  for ir=0,nr-1 do begin
     for it=0,nt-1 do begin
        for ip=0,np-1 do begin
           param= [Nm_vector(ir), fip_factor, Tm_vector(it), ST_vector(ip), SigNe, q]
           y0 = Nm_vector(ir)         ;* (1.0+f_wl*randomn(seed,1))
           y  = synth_y_values(param) ;* (1.0+f_y *randomn(seed,5))
           x(ir,it,ip,0)= y[2]
           x(ir,it,ip,1)= y[3]
           x(ir,it,ip,2)= y[4]
           NmA(ir,it,ip)= Nm_vector(ir)
           TmA(ir,it,ip)= Tm_vector(it)
           STA(ir,it,ip)= ST_vector(iP)
        endfor
     endfor
  endfor

  dir='~/Downloads/'
  datafiles=['x171','x193','x211']
  nband = 3
  for i=0,nband-1 do begin
     openw,1,dir+datafiles(i)
     tmp = reform(x(*,*,*,i))
     writeu,1,tmp
     close,1
  endfor
  par_file='param.out'
  save,filename=dir+par_file,NmA,TmA,STA

  return
end

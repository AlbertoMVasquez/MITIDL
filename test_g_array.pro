;---------------------------------------------------------
; This routine is a wrapper to execute the routines:
; - load_G_array
; - load_sk_array
;--------------------------------------------------------

pro test_g_array

  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common directories, tomroot

  
  set_tomroot
  i_mea_vec           =[0       ,0       ,1    ,1    ,1    ]
  ion_label_vec       =['fexiii','fexiii',''   ,''   ,''   ]
  line_wavelength_vec =['10747' ,'10801' ,''   ,''   ,''   ]
  instrument_label_vec=[''      ,''      ,'aia','aia','aia']
  band_label_vec      =[''      ,''      ,'171','193','211']

  
; Craete a r-Ne-Te 3D grid to interpolate the G-functions 
  NNE=60
  NTe=60
  nmin=1.e6
  nmax=1.e9
  tmin=0.5e6
  tmax=5.0e6
  dNe=(nmax -nmin)/NNe
  dTe=(Tmax -Tmin)/NTe
  fip_factor = 1.

  Nr=26
  dr = 0.01
  r_array  = 1. + dr * findgen(Nr)    + dr/2
  r_array  = [1.1]              ; array de un solo elemento
  Ne_array = Nmin + dNe* findgen(NNe) + dNe/2
  Te_array = tmin + dTe* findgen(NTe) + dTe/2
  
  
  load_G_array,Ne_array,Te_array,r_array,G_A
  load_sk_array,fip_factor,Ne_array,Te_array,r_array,sk_A
  stop

  return
end

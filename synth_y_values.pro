


function synth_y_values,parameters
  
  common NT_limits, Ne0_Limits, Te0_Limits
  common tomographic_measurements, y0, y
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common measurement_errors,sig_WL,sig_y
  common index_measurement, i_measurement
  

  set_tomroot
  Ne0_limits=[1.e5,1.e10]
  Te0_limits=[0.5e6,1.e7]  

  Nem        = parameters[0]
  fip_factor = parameters[1]
  Tem        = parameters[2]
  SigTe      = parameters[3]
  SigNe      = parameters[4]
  q          = parameters[5]

  
  measurement_type = [1,2]

  
  i_mea_vec=[0,0,1,1,1]
  ion_label_vec=       ['fexiii','fexiii','','','']
  line_wavelength_vec=  ['10747','10801'  ,'','','']
  instrument_label_vec=['','','aia','aia','aia']
  band_label_vec=      ['','','171','193','211']
  
  y_synth = i_mea_vec*0d
  M          = n_elements(y_synth)


 
 
  for k = 0, M-1 do begin   
     i_measurement=i_mea_vec(k)  
     ion_label = ion_label_vec(k)
     line_wavelength=line_wavelength_vec(k)
     instrument_label=instrument_label_vec(k)
     band_label = band_label_vec (k)
     if measurement_type[i_measurement] eq 1 then begin
        load_g_table,ion_label=ion_label,line_wavelength=line_wavelength
     endif
     if measurement_type[i_measurement] eq 2 then begin
        load_g_table,instrument_label=instrument_label,band_label=band_label
     endif
     y_synth(k) = e_function(parameters) 
  endfor
  
  return, y_synth
end

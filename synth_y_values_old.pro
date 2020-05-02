
; This routine calculates the synthetic values of emissivities for a
; given set of paramaters.


function synth_y_values_old,parameters
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common index_measurement, i_measurement
  

  Nem        = parameters[0]
  fip_factor = parameters[1]
  Tem        = parameters[2]
  SigTe      = parameters[3]
  SigNe      = parameters[4]
  q          = parameters[5]

  
  M          = n_elements(i_mea_vec)
  y_synth    = dblarr(M) 

 
  for k = 0, M-1 do begin   
     i_measurement    =            i_mea_vec(k)  
     ion_label        =        ion_label_vec(k)
     line_wavelength  =  line_wavelength_vec(k)
     instrument_label = instrument_label_vec(k)
     band_label       =       band_label_vec(k)
     load_g_table,ion_label=ion_label,line_wavelength=line_wavelength,instrument_label=instrument_label,band_label=band_label
     y_synth[k] = e_function(parameters) 
  endfor
  
  return, y_synth
end

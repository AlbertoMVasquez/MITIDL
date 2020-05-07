;-----------------------------------------------------
; This routine does a comparison of the calculated
; emissivities using INT2D and simple Riemann cuadrature,
; respectively. 
;
;INPUT:
; parameters: a 1D array of 6 elements: [Nem, fip_factor, Tem, SigTe, SigNe, q]
;
;
;-----------------------------------------------------

pro compare_integrals,parameters
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common index_measurement, i_measurement
  
  M = n_elements(i_mea_vec)
  
    for k = 0, M-1 do begin   
       
     i_measurement    =            i_mea_vec(k)  
     ion_label        =        ion_label_vec(k)
     line_wavelength  =  line_wavelength_vec(k)
     instrument_label = instrument_label_vec(k)
     band_label       =       band_label_vec(k)
     load_g_table,ion_label=ion_label,line_wavelength=line_wavelength,instrument_label=instrument_label,band_label=band_label

       print,'Relative difference for the emissivity of index k='+string(k)+':'+string(100.*(e_function(parameters) - e_function_cs(k,parameters))/e_function(parameters))+' %'

    endfor
    i_measurement=0
    
  return
end

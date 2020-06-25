; Brief description:
;
; This routine calculates the synthetic values of emissivities for a
; given set of paramaters using the Riemann aproach.
;
; ARGUMENT:
; parameters: a 1D array of 6 elements: [Nem, fip_factor, Tem, SigTe, SigNe, q]
;
;
; OUTPUTS:
; y_synth:  an M-element 1D vector containing M synthetic tomograhic measurements in
; the voxel, possibly including: EUV FBEs, CoMP/UCoMP line emissivities. 
;
; History:  V1.0, V1.0, F.A. Nuevo, IAFE, May 2020.
;---------------------------------------------------------------------


function synth_y_values_cs,parameters
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  ;common index_measurement, i_measurement
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q  

  Nem        = parameters[0]
  fip_factor = parameters[1]
  Tem        = parameters[2]
  SigTe      = parameters[3]
  SigNe      = parameters[4]
  q          = parameters[5]

  M          = n_elements(i_mea_vec)
  y_synth    = dblarr(M) 
 
        
  for k=0,M-1 DO  y_synth[k] = e_function_cs(k,parameters) 
  
  return, y_synth
end



  
  
  

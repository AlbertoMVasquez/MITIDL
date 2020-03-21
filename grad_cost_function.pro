;---------------------------------------------------------------------
;
; Brief description:
;
; Gradient of the Cost funtion to use in the minimization.
;
; Argument:
; parameters: a 1D array of 6 elements: [Nem, fip_factor, Tem, SigTe, SigNe, q]
;
; Parameters in COMMON BLOCKS:
;
; common tomographic_measurements: contains:
;     - y0: white-light tomography electron density of the voxel.
;     - y:  an M-element 1D vector containing M tomograhic measurements in
;           the voxel, possibly including: EUV FBEs, CoMP/UCoMP line emissivities.
;     - measurement_type: am M-element 1D vector containing the type
;       of measurement of the corresponding element in array y,
;       possible values are: 1, for CoMP/UCoMP line emissivity,
;                            2, for EUV FBE.
;
; To-be-done: There are yet no weighting factors 1/SIGMA² in each term
;             of the cost function.
;
; OUTPUTS:
; Value of the gradient of the cost function for the given values 
; of the inputs and the parameters.
;
; History:  V1.0, Federico Nuevo, IAFE, March-2020.
;      
;---------------------------------------------------------------------
function grad_cost_function, parameters
  common NT_limits, Ne0_Limits, Te0_Limits
  common tomographic_measurements, y0, y, measurement_type, i_measurement
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common weights,sig_WL,sig_v

  RESULT     = parameters * 0d
  Nem        = parameters[0]
  fip_factor = parameters[1]
  Tem        = parameters[2]
  SigTe      = parameters[3]
  SigNe      = parameters[4]
  q          = parameters[5]
  
  M          = n_elements(y)
  for k = 0, M-1 do begin
     
     i_measurement   = i_mea_vec           (k)  
     ion_label       = ion_label_vec       (k)
     line_wavelength = line_wavelength_vec (k)
     instrument_label= instrument_label_vec(k)
     band_label      = band_label_vec      (k)
     
     if measurement_type[i_measurement] eq 1 then begin
        load_g_table,ion_label=ion_label,line_wavelength=line_wavelength
     endif
     if measurement_type[i_measurement] eq 2 then begin
        load_g_table,instrument_label=instrument_label,band_label=band_label
     endif
     RESULT = RESULT + 2*(e_function(parameters) - y[k]) /sig_v[k]^2   * grad_e_function(parameters)
  endfor
  result(0) = result(0) + 2*(Nem-y0)/sig_WL^2
  

  return, RESULT
end

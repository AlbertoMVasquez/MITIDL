;---------------------------------------------------------------------
;
; Brief description:
;
; Gradient of the Cost funtion to use in the minimization. Use CS to
; calculate the double integrals
;
; ARGUMENT:
; parameters: a 1D array of 6 elements: [Nem, fip_factor, Tem, SigTe, SigNe, q]
;
; Parameters in COMMON BLOCKS:
;
; common tomographic_measurements: contains:
;     - y0: white-light tomography electron density of the voxel.
;     - y:  an M-element 1D vector containing M tomograhic measurements in
;           the voxel, possibly including: EUV FBEs, CoMP/UCoMP line emissivities.
;     - sig_WL, sig_y: error of the y0 and y, respectively. Used as
;       weights in the cost function. 
;
;
; OUTPUT:
; Value of the gradient of the cost function (vector of 6 components) for the given values 
; of the inputs and the parameters.
;
; History:  V1.0, Federico Nuevo, IAFE, April-2020.
;      
;---------------------------------------------------------------------
function gradphi, parameters
  common tomographic_measurements, y0, y
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common measurement_errors,sig_WL,sig_y
  common index_measurement, i_measurement
  
  RESULT     = parameters * 0d
  Nem        = parameters[0]  
  q          = parameters[5]  
  M          = n_elements(y)
  for k = 0, M-1 do begin
     RESULT = RESULT + 2*(e_function_cs(k,parameters) - y[k]) /sig_y[k]^2   * grad_e_function_cs(k,parameters)
  endfor
  result(0) = result(0) + 2*(Nem-y0)/sig_WL^2
  result(5) = result(5) + df_penalty_q(q)
  return, RESULT
end

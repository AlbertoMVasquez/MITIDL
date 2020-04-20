;---------------------------------------------------------------------
;
; Brief description:
;
; Cost funtion to be minimazied in a each voxel of the tomographic grid.
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
; Value of the function for the given values of the inputs and the parameters.
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;           V1.1 la versión 1.0 no andaba, varias cosas ... 
;---------------------------------------------------------------------
function cost_function_cs, parameters
  
  common measurement_errors,sig_WL,sig_y
  common tomographic_measurements, y0, y
  

  Nem        = parameters[0]
  M          = n_elements(y)
  
  RESULT = (Nem-y0)^2/sig_WL^2
  for k = 0, M-1 do begin   
     RESULT = RESULT + (e_function_cs(k,parameters) - y[k])^2/sig_y[k]^2
     
  endfor

  return, RESULT
  end
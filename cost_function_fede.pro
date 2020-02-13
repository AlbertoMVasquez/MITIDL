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
;
;---------------------------------------------------------------------
function cost_function_fede, parameters
  common NT_limits, Ne0_Limits, Te0_Limits
  common tomographic_measurements, y0, y, measurement_type, i_measurement
  Nem        = parameters[0]
  fip_factor = parameters[1]
  Tem        = parameters[2]
  SigTe      = parameters[3]
  SigNe      = parameters[4]
  q          = parameters[5]
  M          = n_elements(y)
  ;--------------------------------SOME TESTS--------------------------------------------------------
  for i_measurement = 0, M-1 do begin
     if measurement_type[i_measurement] ne 1 AND measurement_type[i_measurement] ne 2 then begin
        print,'Wrong measurement type for y-element #',i_measurement
        STOP
     endif
     print, i_measurement, measurement_type[i_measurement], e_function_fede(parameters), y[i_measurement]
  endfor
  ;--------------------------------------------------------------------------------------------------
  RESULT = (Nem-y0)^2
  for i_measurement = 0, M-1 do $
  RESULT = RESULT + (e_function_fede(parameters) - y[i_measurement])^2
  return, RESULT
end

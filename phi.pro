;---------------------------------------------------------------------
;
; Brief description:
;
; Cost funtion to be minimizied in a each voxel of the tomographic grid.
; This version (_cs) uses CS to calculate the emissivities
;
; CS: \int f(x,y) dx dy >  \Sum_{i,j} f(x_i,y_j) Dx Dy
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
; common measurement_errors,sig_WL,sig_y: contains:     
;     - sig_wl, sig_y: error of y0 and y, respectively. used as
;       weights in the cost function.
;    
;                          
; OUTPUT:
; Value of the cost function for the given values of the inputs and the parameters.
;
; History:  V1.0, F.A. Nuevo, IAFE, April-2020.
;
;
;---------------------------------------------------------------------
function phi, parameters  
  common measurement_errors,sig_WL,sig_y
  common tomographic_measurements, y0, y
  RESULT     = 0d
  Nem        = parameters[0]
  q          = parameters[5]
  M          = n_elements(y) 
  RESULT = (Nem-y0)^2/sig_WL^2
  for k = 0, M-1 do RESULT = RESULT + (e_function_cs(k,parameters) - y[k])^2/sig_y[k]^2

  ;RESULT = RESULT + f_penalty_q(q)

  return, RESULT
end

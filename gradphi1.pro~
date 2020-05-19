;---------------------------------------------------------------------
;
; Brief description:
;
; Gradient of the Cost funtion to use in the minimization.
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
;     - imea_vec: a vector with the i_measurement of the M tomographic measurements
;       (i_measurement: scalar the type of measurement (0=line, 1=FBE). 
;     - sig_WL,sig_y: error of the y0 and y, respectively. Used as
;       weights in the cost function. 

; OUTPUT:
; Value of the gradient (vector of 6 component) of the cost function 
; for the given values of the inputs and the parameters.
;
; History:  V1.0, F.A. Nuevo, IAFE, March-2020.
;           V1.1, A.M. Vasquez, HOME, April-2020
;      
;---------------------------------------------------------------------
function grad_cost_function, parameters
  common tomographic_measurements, y0, y
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common measurement_errors,sig_WL,sig_y
  common index_measurement, i_measurement
  common G_table, G, T_e, N_e, r, photT
  common tables,Te1,Te2,Te3,Te4,Te5,Ne1,Ne2,Ne3,Ne4,Ne5,G1,G2,G3,G4,G5,r1,r2
  
  RESULT     = parameters * 0d
  Nem        = parameters[0]
  fip_factor = parameters[1]
  Tem        = parameters[2]
  SigTe      = parameters[3]
  SigNe      = parameters[4]
  q          = parameters[5]
  
  result(0)  = result(0) + 2*(Nem-y0)/sig_WL^2

  M          = n_elements(y)  
  for k = 0, M-1 do begin
     i_measurement = i_mea_vec(k)  
     CASE k of
     0: BEGIN
        G   = G1
        T_e = Te1
        N_e = Ne1
        r   = r1
     END
     1: BEGIN
        G   = G2
        T_e = Te2
        N_e = Ne2
        r   = r2
     END
     2: BEGIN
        G   = G3
        T_e = Te3
        N_e = Ne3
     END
     3: BEGIN
        G   = G4
        T_e = Te4
        N_e = Ne4
     END
     4: BEGIN
        G   = G5
        T_e = Te5
        N_e = Ne5
     END
     ENDCASE     
     RESULT = RESULT + (2*(e_function(parameters) - y[k])/sig_y[k]^2)*grad_e_function(parameters)
  endfor
  
  return, RESULT
end

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
;           V2.0, Federico A. Nuevo, IAFE, March 2020
;           El indice del loop pasa a ser k (i_meausurement es para
;           evaluar measurement_type)
;           Carga g_table para para medición 
;---------------------------------------------------------------------
function cost_function, parameters
  common NT_limits, Ne0_Limits, Te0_Limits
  common tomographic_measurements, y0, y, measurement_type, i_measurement
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common weights,sig_WL,sig_v
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common G_table, G, T_e, N_e, r, photT

  Nem        = parameters[0]
  M          = n_elements(y)
  
  RESULT = (Nem-y0)^2/sig_WL^2
  
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
     
     RESULT = RESULT + (e_function(parameters) - y[k])^2/sig_v[k]^2
  endfor

  return, RESULT
  end

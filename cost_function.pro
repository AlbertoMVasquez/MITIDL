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
;     - i_measurement: scalar the type of measurement (0=line, 1=FBE)
;       of measurement of the corresponding element in array y,
;       possible values are: 1, for CoMP/UCoMP line emissivity,
;                            2, for EUV FBE.
;
; OUTPUTS:
; Value of the function for the given values of the inputs and the parameters.
;
; History:  V1.0, A.M. Vasquez, CLaSP, Spring-2018.
;           V2.0, F.A. Nuevo, IAFE, March 2020.
;                 Load correct g_table to compute appropriate e_function for each term.
;                 Also introduced the weight factors 1/sigma².
;           V2.1, A.M. Vasquez, IAFE, March-2020.
;                 Only one call to load_g_table was needed (there were
;                 two calls in V2.0).
;                 Also, "measurement_type" is not needed anymore.
;                 Now i_measurement indicates the type of measurement,
;                 provided to load_g_table through a dedicated common block,
;                 ultimately needed by function_g.pro
;---------------------------------------------------------------------
function cost_function, parameters
  common tomographic_measurements, y0, y
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common measurement_errors,sig_WL,sig_y
  common index_measurement, i_measurement
  Nem    = parameters[0]
  M      = n_elements(y)  
  RESULT = (Nem - y0)^2/sig_WL^2  
  for k = 0, M-1 do begin
     i_measurement    =            i_mea_vec(k)  
     ion_label        =        ion_label_vec(k)
     line_wavelength  =  line_wavelength_vec(k)
     instrument_label = instrument_label_vec(k)
     band_label       =       band_label_vec(k)
     load_g_table,ion_label=ion_label,line_wavelength=line_wavelength,instrument_label=instrument_label,band_label=band_label
     ; The previous call to the load_g_table must be replaced here
     ; by selection of the proper G table and associated arrays
     ; from the arrays stored in a new common block named "G_tables"
     ; which should contain all tables. The selected arrays are 
     ; to be named with as the variables of the current "G_table"
     ; common block, hence passing those to "g_function.pro".
     RESULT = RESULT + (e_function(parameters) - y[k])^2/sig_y[k]^2
  endfor
  return, RESULT
  end

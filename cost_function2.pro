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
function cost_function2, parameters
  common tomographic_measurements, y0, y
  common measurement_errors,sig_WL,sig_y
  common G_table, G, T_e, N_e, r, photT
  common tables,Te1,Te2,Te3,Te4,Te5,Ne1,Ne2,Ne3,Ne4,Ne5,G1,G2,G3,G4,G5
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec 
  common index_measurement, i_measurement

  Nem    = parameters[0]
  M      = n_elements(y)  
  RESULT = (Nem - y0)^2/sig_WL^2  
  
  for k = 0, M-1 do begin
     i_measurement    =            i_mea_vec(k)  


     if k eq 0 then begin
        G   = G1
        T_e = Te1
        N_e = Ne1
     endif
     if k eq 1 then begin
        G   = G2
        T_e = Te2
        N_e = Ne2
     endif
     if k eq 2 then begin
        G   = G3
        T_e = Te3
        N_e = Ne3
     endif
     if k eq 3 then begin
        G   = G4
        T_e = Te4
        N_e = Ne4
     endif
     if k eq 4 then begin
        G   = G5
        T_e = Te5
        N_e = Ne5
     endif
  
   
     RESULT = RESULT + (e_function(parameters) - y[k])^2/sig_y[k]^2
    
  endfor
  return, RESULT
  end

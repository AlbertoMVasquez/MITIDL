; Brief description:
;
; This routine calculates the synthetic values of emissivities for a
; given set of paramaters.
;
; ARGUMENT:
; parameters: a 1D array of 6 elements: [Nem, fip_factor, Tem, SigTe, SigNe, q]
;
;
; OUTPUTS:
; y:  an M-element 1D vector containing M synthetic tomograhic measurements in
; the voxel, possibly including: EUV FBEs, CoMP/UCoMP line emissivities. 
;
; History:  V1.0, V1.0, F.A. Nuevo, IAFE, May 2020.
;---------------------------------------------------------------------


function synth_y_values,parameters
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common index_measurement, i_measurement
  common tables,Te1,Te2,Te3,Te4,Te5,Ne1,Ne2,Ne3,Ne4,Ne5,G1,G2,G3,G4,G5,r1,r2
  common G_table, G, T_e, N_e, r, photT
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q  

  Nem        = parameters[0]
  fip_factor = parameters[1]
  Tem        = parameters[2]
  SigTe      = parameters[3]
  SigNe      = parameters[4]
  q          = parameters[5]

  M          = n_elements(i_mea_vec)
  y_synth    = dblarr(M) 
 
  FOR k = 0, M-1 do begin   
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
     y_synth[k] = e_function(parameters) 
  ENDFOR
  return, y_synth
end



  
  
  

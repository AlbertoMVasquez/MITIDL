;---------------------------------------------------------------------
;
; Brief description: 
; This routine computes the S_k(Ne,Te,r) functions (divided by fip_factor)
; for a fixed Ne and Te grid (evaluated in Ne_array and Te_array) and
; given r0 value. 
;
; INPUTS 
;  (in common NT_ARRAYS):
;  Te_array: 1D array of NTe elements with the electron temperature in units of [K]
;  Ne_array: 1D array of NNe elements with the electron density     in units of [cm-3]
;  (in common parameters):
;  r0: heliocentric height where the s_k is evaluated. 

; OUTPUT: sk_over_fip_factor, An 3D array  of M x NTe x NNe elements
; with S_k(Ne,Te, r) evaluated in Te_array and Ne_array.
 
; 
; HISTORY
; V1.0 F.A. Nuevo, IAFE, April-2020
; V1.1 A.M. Vasquez, IAFE, April-2020
;      Load G tables from memory.
;---------------------------------------------------------------------
pro make_sk_over_fip_factor
  common sk_over_fip_factor_array,sk_over_fip_factor  
  common NT_arrays,Ne_array,Te_array,dNe_array,dTe_array,dTN
  common dimensions, NTe, NNe
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common index_measurement, i_measurement
  common tables,Te1,Te2,Te3,Te4,Te5,Ne1,Ne2,Ne3,Ne4,Ne5,G1,G2,G3,G4,G5,r1,r2,r3,r4,r5
  common G_table, G, T_e, N_e, r, photT
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  M  = n_elements(i_mea_vec)
  sk_over_fip_factor = dblarr (M, NTe , NNe)
 
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
        r   = r3 
     END
     3: BEGIN
        G   = G4
        T_e = Te4
        N_e = Ne4
        r   = r4
     END
     4: BEGIN
        G   = G5
        T_e = Te5
        N_e = Ne5
        r   = r5
     END
     ENDCASE
     sk_over_fip_factor(k,*,*) = s_function(Ne_array,Te_array) / fip_factor
  endfor
  return
end


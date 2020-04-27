;---------------------------------------------------------------------
;
; Brief description:
;
; Load  all the tables of the G(Ne,Te) functions:
; Gi, Tei, Nei (i=1,...,5)
; the arrays are in common tables
;
; IMPORTANT NOTE:
;It stores G tables in a specific/rigid order:
;1: CoMP-1074
;2: CoMP-1079
;3: AIA-171
;4: AIA-193
;5: AIA-211
;
; For lines:
; 3D array:  G(Te, Ne, r) 
; 1D arrays: log10(Te [K]), log10(Ne [cm-3]), r [Rsun],
; scalar: photosphere Teff [K].
;
; For EUV bands:
; 1D arrays: TRF(Te), log10(Te [K])
;
; Tables must be stored in:
; tomroot/MITIDL/Emissivity/LookUp_Tables/
;
; 
; History:  V1.0, Federico A. Nuevo, IAFE, April-2020.
;

;---------------------------------------------------------------------
pro  load_tables
  
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common index_measurement, i_measurement
  common G_table, G, T_e, N_e, r, photT
  common directories, tomroot
  common tables,Te1,Te2,Te3,Te4,Te5,Ne1,Ne2,Ne3,Ne4,Ne5,G1,G2,G3,G4,G5,r1,r2
   
  M      = n_elements(i_mea_vec)  
 
  for k = 0, M-1 do begin
     i_measurement    =            i_mea_vec(k)  
     ion_label        =        ion_label_vec(k)
     line_wavelength  =  line_wavelength_vec(k)
     instrument_label = instrument_label_vec(k)
     band_label       =       band_label_vec(k)
     load_g_table,ion_label=ion_label,line_wavelength=line_wavelength,instrument_label=instrument_label,band_label=band_label
     CASE k of
     0: BEGIN
        G1  = G
        Te1 = T_e
        Ne1 = N_e
        r1  = r
     END
     1: BEGIN
        G2  = G
        Te2 = T_e
        Ne2 = N_e
        r2  = r
     END
     2: BEGIN
        G3  = G
        Te3 = T_e
        Ne3 = N_e 
     END
     3: BEGIN
        G4  = G
        Te4 = T_e
        Ne4 = N_e 
     END
     4: BEGIN
        G5  = G
        Te5 = T_e
        Ne5 = N_e 
     END
     ENDCASE
  endfor
  return
  end

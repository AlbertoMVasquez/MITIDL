;---------------------------------------------------------------------
;
; Brief description: 
; This routine computes the S_k(Ne,Te) functions for a fixed Ne and Te
; grid. 
; Es necesario que se hallan creado los arrays 1D Ne_array y Te_array
; Output: sk array de M x NTe x NNe
; 
; History: 
;           V1.0, F.A. Nuevo, IAFE, April-2020.
;                      
;---------------------------------------------------------------------


pro make_sk_load_routine,sk
  
  common NT_arrays,Ne_array,Te_array
  common dimensions, NTe, NNe
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common index_measurement, i_measurement
 
  M  = n_elements(i_mea_vec)
  sk= dblarr (M, NTe , NNe)
 
  for k = 0, M-1 do begin   

     i_measurement    =            i_mea_vec(k)  
     ion_label        =        ion_label_vec(k)
     line_wavelength  =  line_wavelength_vec(k)
     instrument_label = instrument_label_vec(k)
     band_label       =       band_label_vec(k)
     load_g_table,ion_label=ion_label,line_wavelength=line_wavelength,instrument_label=instrument_label,band_label=band_label
     
     sk(k,*,*) = s_function (Ne_array,Te_array)

  endfor
  
  i_measurement=0

  return
end


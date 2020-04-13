;---------------------------------------------------------------------

pro  load_tables
  
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common index_measurement, i_measurement
  common G_table, G, T_e, N_e, r, photT
  common directories, tomroot
  common tables,TeCoMP,NeCoMP,TeEUV,NeEUV,G1,G2,G3,G4,G5
  
 
  M      = n_elements(i_mea_vec)  
 
  for k = 0, M-1 do begin
     i_measurement    =            i_mea_vec(k)  
     ion_label        =        ion_label_vec(k)
     line_wavelength  =  line_wavelength_vec(k)
     instrument_label = instrument_label_vec(k)
     band_label       =       band_label_vec(k)
     load_g_table,ion_label=ion_label,line_wavelength=line_wavelength,instrument_label=instrument_label,band_label=band_label
          
     if k eq 0  then begin
                     G1=G
        TeCoMP=T_e
        NeCoMP=N_e 
     endif
     if k eq 1 then  G2=G
     if k eq 2 then begin
                     G3=G
        TeEUV=T_e
        NeEUV=N_e
     endif
     if k eq 3 then  G4=G
     if k eq 4 then  G5=G
     
     i_measurement=0
    
  endfor
  return
  end

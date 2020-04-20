;---------------------------------------------------------------------
pro  load_tables2
  
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common index_measurement, i_measurement
  common G_table, G, T_e, N_e, r, photT
  common directories, tomroot
  common tables,T1,T2,T3,T4,T5,N1,N2,N3,N4,N5,G1,G2,G3,G4,G5
   
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
        T1=T_e
        N1=N_e 
     endif
     if k eq 1  then begin
        G2=G
        T2=T_e
        N2=N_e 
     endif
     if k eq 2  then begin
        G3=G
        T3=T_e
        N3=N_e 
     endif
     if k eq 3  then begin
        G4=G
        T4=T_e
        N4=N_e 
     endif
     if k eq 4  then begin
        G5=G
        T5=T_e
        N5=N_e 
     endif
     i_measurement=0
    
  endfor
  return
  end

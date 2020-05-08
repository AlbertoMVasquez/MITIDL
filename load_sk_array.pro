;---------------------------------------------------------------------
;
; Brief description:
;
; Load  all the tables of the G(Ne,Te,r) functions:
; Gi, Tei, Nei (i=1,...,5) in one array and calculate directly
; an array with s_k(Ne,Te,r)
;
;
;
;INPUT: 1D arrays: Te_A, Ne_A, r_A 
; where the original G tables are interpolated
;OUTPUT: 4D array:  sk_A(k,Te, Ne, r) 
;
;
; Tables must be stored in:
; tomroot/MITIDL/Emissivity/LookUp_Tables/
;
; 
; History:  V1.0, Federico A. Nuevo, IAFE, April-2020.
;

;---------------------------------------------------------------------
pro  load_sk_array,fip_factor,Ne_A,Te_A,r_A,sk_A
  
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common index_measurement, i_measurement
  common G_table, G, T_e, N_e, r, photT
  common directories, tomroot
  
   
  M      = n_elements(i_mea_vec)  
  NTe    = (size(Te_A))(1)
  NNe    = (size(Ne_A))(1)
  Nr     = (size(r_A ))(1)
  sk_A    = dblarr(M,NTe,NNe,Nr)
  
  for k = 0, M-1 do begin
     i_measurement    =            i_mea_vec(k)  
     ion_label        =        ion_label_vec(k)
     line_wavelength  =  line_wavelength_vec(k)
     instrument_label = instrument_label_vec(k)
     band_label       =       band_label_vec(k)
     load_g_table,ion_label=ion_label,line_wavelength=line_wavelength,instrument_label=instrument_label,band_label=band_label

     CASE i_measurement OF
        0: BEGIN
        
           for iTe=0,NTe-1 do begin
              for iNe=0,NNe-1 do begin
                 for ir=0,Nr-1 do begin
                    sk_A[k,iTe,iNe,ir] = Ne_A[iNe]^2 * findval3d_function(G,T_e,N_e,r,Te_A[iTe],Ne_A[iNe],r_A[ir])
                 endfor
              endfor
           endfor
        END
        1: BEGIN
           
           G_atTe0       = interpol( G    , T_e, Te_A)
           for iNe=0,NNe-1 do begin
              for ir=0,Nr-1 do begin
                 sk_A[k,*,iNe,ir] = Ne_A[iNe]^2 * G_atTe0
              endfor
           endfor
        END
     ENDCASE
  endfor

   sk_A= sk_A / 4./!pi

  return
end

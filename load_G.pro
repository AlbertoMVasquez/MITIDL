;---------------------------------------------------------------------
;
; Brief description:
;
; Load look up table of contribution function G for user-specified Ion
; and line. The table is an IDL SAVE file containing:
; G [erg cm^+3 sec^-1], log(Te [K]), log10(Ne [cm-3]), rad [Rsun],
; and must be stored in:
; tomroot/tomography/MultiTom/Emissivity/LookUp_Tables/
;
; INPUTS:
; ion_label: a string specifying the ion, possible values are:
; 'fexiii', ....
;
; line_wavelength: a string specifying the wavelength in A, possible
; values are: '10747'.
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

pro load_G,ion_label=ion_label,line_wavelength=line_wavelength
  common G_table,G,T_e,N_e,r,photT
  common directories,tomroot
  
  data_dir  = tomroot+'MultiTom/Emissivity_LookUp_Tables/'
  file_name = 'G_function_'+ion_label+'_'+line_wavelength+'.save'

  restore,data_dir+file_name

  G     = emissivity ; erg cm^+3 sec^-1 ; Note: G(Te, Ne, r)
  N_e   = 10.^dens   ; cm^-3
  T_e   = 10.^temp   ; K
  r     = rphot      ; Rsun
  photT = radtemp    ; K

  return
end

;---------------------------------------------------------------------
;
; Brief description:
;
; Wrapper routine to test suite of codes.
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

pro test,Ne0=Ne0,Te0=Te0,euvband=euvband,emissionline=emissionline
  common constants,Rsun,kB,h,c
  common G_table,G,T_e,N_e,r,photT
  common directories,tomroot
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common dimensions,NTe,NNe

  r0=1.2
  fip_factor=1.
  Tem=1.75e6
  Nem=1.75e8
  SigTe=0.5e6
  SigNe=0.5e8
  q=0.
  
  set_tomroot
  if not keyword_set(ion_label)        then ion_label        = 'fexiii' ; always use lowercase
  if not keyword_set(line_wavelength)  then line_wavelength  =  '10747' ; string with wavelength in A
  if not keyword_set(fip_factor)       then fip_factor       =     1.0  ; Feldmand's Adundance Set value
  if not keyword_set(instrument_label) then instrument_label =    'aia' ; always use lowercase
  if not keyword_set(band_label)       then band_label       =    '171'

  load_g_table,ion_label=ion_label,line_wavelength=line_wavelength,instrument_label=instrument_label,band_label=band_label,emissionline=emissionline,euvband=euvband

  if not keyword_set(Ne0) then Ne0=2.5e8
  if not keyword_set(Te0) then Te0=1.5e6

  NTe = 1
  NNe = 1
  if (size(Te0))(0) eq 1 then NTe = (size(Te0))(1)
  if (size(Ne0))(0) eq 1 then NNe = (size(Ne0))(1)

  print
  print,'Input values of Ne [cm^-3], Te [K], rad [Rsun], fip_factor:'
  print,Ne0,Te0,r0,fip_factor
  print
  print,'G-function value [erg cm+3 sec-1]:'
  print,(g_function(Te0, Ne0,euvband=euvband,emissionline=emissionline))(0)
  print
  print,'s [erg cm-3 sec-1 sr-1]:'
  print,(s_function(Ne0, Te0,euvband=euvband,emissionline=emissionline))(0)
  print
  print,'p [cm+3 K-1]:'
  print,p_function(Ne0, Te0)
  print
  print,'s*p [erg sec-1 sr-1 K-1]:'
  print,sxp_function(Ne0, Te0,emissionline=emissionline,euvband=euvband)
  print

  Ne0_Limits = [min(N_e),max(N_e)]
  Te0_Limits = [min(T_e),max(T_e)]
  print,'e [erg sec-1 sr-1 K-1]:'
  print, e_function( Ne0_Limits , Te0_Limits, emissionline=emissionline ,euvband=euvband)
  print
  print,'e2[erg sec-1 sr-1 K-1]:'
  print, e2_function( Ne0_Limits , Te0_Limits )
  print
  
  return
end

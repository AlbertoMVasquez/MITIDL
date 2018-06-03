;---------------------------------------------------------------------
;
; Brief description:
;
; Wrapper routine to test suite of codes.
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

pro test,Ne0=Ne0,Te0=Te0;,r0=r0,fip_factor=fip_factor,line_wavelength=line_wavelength
  common constants,Rsun,kB,h,c
  common G_table,G,T_e,N_e,r,photT
  common directories,tomroot
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q

  r0=1.2
  fip_factor=1.
  Tem=1.75e6
  Nem=1.75e8
  SigTe=0.5e6
  SigNe=0.5e8
  q=0.5

  set_tomroot
  if not keyword_set(ion_label)       then ion_label       = 'fexiii' ; always use lowercase
  if not keyword_set(line_wavelength) then line_wavelength =  '10747' ; string with wavelength in A
  if not keyword_set(fip_factor)      then fip_factor      =     1.0  ; Feldmand's Adundance Set value

  load_g_table,ion_label=ion_label,line_wavelength=line_wavelength

  print
  print,'Input values of Ne [cm^-3], Te [K], rad [Rsun], fip_factor:'
  print,Ne0,Te0,r0,fip_factor
  print
  print,'G-function value [erg cm+3 sec-1]:'
  print,(g_function(Te0, Ne0))(0)
  print
  print,'s [erg cm-3 sec-1 sr-1]:'
  print,(s_function(Ne0, Te0))(0)
  print
  print,'p [cm+3 K-1]:'
  print,p_function(Ne0, Te0)
  print
  print,'s*p [erg sec-1 sr-1 K-1]:'
  print,sxp_function(Ne0, Te0)
  print

  Ne0_Limits = [1.e7, 1.e9]
  Te0_Limits = [1.e6, 3.e6]
  print,'s*p [erg sec-1 sr-1 K-1]:'
  print, e_function( Ne0_Limits , Te0_Limits )
  print
  
  return
end

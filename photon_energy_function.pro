;---------------------------------------------------------------------
;
; Brief description:
;
; This functions returns the energy [erg] of a photon of wavelength
; line_wavelength.
;
; INPUTS:
; line_wavelength: ; string with wavelength in A (see routine
; load_G_table.pro for valid values).
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

function photon_energy_function,line_wavelength
  common constants,Rsun,kB,h,c
  if wavelength_line eq '10747' then lambda0 = 10747.e-8 ; cm
  E0 = h*c/lambda_0                                      ; erg
  return,E0
end


; Note: 1 W m-2 nm-1 = 1.e+2 erg s-1 cm-2 A-1

pro convert_units_comp

  Rsun = 6.96e10 ; cm
  r     = [1.11, 1.21] ; Rsun

; Median tomographic values in units of [1.e-6 Bsun*A]
  x_comp_1074 = [3.01564,1.09300]
  x_comp_1079 = [1.05946,1.01825] ; These values will probably be corrected.
;  x_comp_1074 = [3.01406,1.13503]
;  x_comp_1079 = [1.52335,0.865927] ; These values will probably be corrected.


; ETR Values for Bsun
  Bsun_1075 = 6.1768E-01 * 1.e+2 ; erg s-1 cm-2 A-1
  Bsun_1080 = 6.2250E-01 * 1.e+2 ; erg s-1 cm-2 A-1

; Median tomographic values in units of [erg s-1 cm-3 sr-1]
  x_comp_1074 = x_comp_1074 * Bsun_1075 / (4.*!pi) / rsun 
  x_comp_1079 = x_comp_1079 * Bsun_1080 / (4.*!pi) / rsun 

  print
  print,'CoMP tomographic emissivities [erg s-1 cm-3 sr-1]'
  print,'      r[Rs]    E_1074       E_1079'
  for ir=0,n_elements(r)-1 do print,r[ir],x_comp_1074[ir],x_comp_1079[ir]
  print
  
  return
end

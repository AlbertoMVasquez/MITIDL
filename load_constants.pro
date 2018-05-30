;---------------------------------------------------------------------
;
; Brief description:
;
; This routine loads useful solar and universal constants
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

pro load_constants
common constants,Rsun,kB,h,c
  Rsun = 6.957e+10        ; cm
  kB   = 1.38064852e-16   ; erg K^-1
  h    = 6.62607004e10-27 ; cm^2 gr sec^-1
  c    = 299792458.e+2    ; cm sec^-1 
end

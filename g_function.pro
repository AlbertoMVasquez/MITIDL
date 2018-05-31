;---------------------------------------------------------------------
;
; Brief description:
;
; This function returns:
;
;  - For lines: the contribution function G provided values for Ne, Te and r.
;
;  - For EUV bands: the temperature response function TRF provided values for Ne, Te.
;
; In both cases the routine tri-linearly interpolates G/TRF into (Ne,Te,r).
;
; INPUTS: values for Ne, Te, rad
;
; OUTPUTS:
;
; G [erg cm+3 sec-1]
; or
; TRF [units?]
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

function g_function, Te0, Ne0, r0
  common G_table,G,T_e,N_e,r,photT

  goto,tri_linear_interpolator
 
  ; Rough 1st-approach: Assign closest value from look-up table
  fTe = abs(Te0-T_e) & iTe = median(where(fTe eq min(fTe)))
  fNe = abs(Ne0-N_e) & iNe = median(where(fNe eq min(fNe)))
  fr  = abs(r0 -r  ) & ir  = median(where(fr  eq min(fr )))
  if iTe eq -1 or iNe eq -1 or ir eq -1 then begin
     print,'Can not assign value to G_function'
     stop
  endif
  G_value = G(iTe,iNe,ir)
  goto,exit
  
  tri_linear_interpolator:
  ; Tri-linear interpolator
  ; Rename dimensions: Te->X, Ne->Y, rad->Z
  DATA_ARRAY = G
  xa = T_e & x0 = Te0
  ya = N_e & y0 = Ne0
  za = r   & z0 = r0

  G_value = findval3D_function(DATA_ARRAY,xa,ya,za,x0,y0,z0)

  exit:
  return,G_value
end

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
 
  ; Assign closest value from look-up table
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
  ; Rename dimensions: Ne->X, Te->Y, rad->Z
  xa = N_e & x0 = Ne0
  ya = T_e & y0 = Te0
  za = r   & z0 = r0
  ; Find the x-planes that surround x0
  ixA = max(where(xa le x0))
  if ixA eq -1 then begin
     print,'Ne value is out of range' ; we may need something better to treat this case.
     stop
  endif
  ixB=ixA+1
  ; Define two 2D y-z arrays at x-planes that surround x0. 
  DATA_ARRAY_xA = reform(G(*,ixA,*))
  DATA_ARRAY_xB = reform(G(*,ixB,*))
  ; Bi-linearly interpolate in xA and xB planes.   
  G_xA = findval2D_function( DATA_ARRAY_xA ,ya ,za , y0, z0 )
  G_xB = findval2D_function( DATA_ARRAY_xB ,ya ,za , y0, z0 )
  ; Linearly interpolate the value og G along x, betwee xA and xB planes.
  G_value = G_xA + (G_xB-G_xA)*(x0-xa(ixA))/(xa(ixB)-xa(ixA))

  exit:
  return,G_value
end

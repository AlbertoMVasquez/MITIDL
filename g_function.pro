;---------------------------------------------------------------------
;
; Brief description:
;
; This function returns:
;
;  - For lines: the contribution function G provided values for Ne, Te and r.
;
;  - For EUV bands: the temperature response function TRF provided
;    values for Ne, Te.
;
;  - For f=G/TRF, it also returns the derivatives:
;    df/dTe, df/dNe, df/dr, evaluated at (Te0,Ne0,r0)
;
; In both cases the routine tri-linearly interpolates G/TRF into (Te0,Ne0,r0).
;
; INPUTS: values for Te0, Ne0, r0
;
; OUTPUTS:
;
; RESULT = [f, df/dTe, df/dNe, df/dr]
;
; where f is:
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
  RESULT = findval3D_function(G,T_e,N_e,r,Te0,Ne0,r0)
  return, RESULT
end

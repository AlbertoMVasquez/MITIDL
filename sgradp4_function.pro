;---------------------------------------------------------------------
;
; Brief description:
;
; This function computes the product of the s and dp/dsigN functions.
; dp/dsigN is calculated with grad_p_function.
;
; INPUTS: Ne0, Te0.
; OUTPUT: The value of s*dp/sigN
;
; History:  V1.0, F.A. Nuevo, IAFE, March-2020.
;           V1.1, A.M. Vasquez, HOME, April-2020.
;
;---------------------------------------------------------------------
function sgradp4_function, Ne0, Te0
  s =  s_function(Ne0,Te0) ;*0. + 1. ;ACTIVATE to make s=1.
  gradP = grad_p_function(Ne0, Te0)
  gradp4=reform(gradP(*,*,3)) ; dP/dsigN
  return, s*gradp4
end

;---------------------------------------------------------------------
;
; Brief description:
;
; This function computes the product of the s and dp/dNm functions.
; dp/dNm is calculated with grad_p_function.
;
; INPUTS: Ne0, Te0.
;
;
; History:  V1.0, Federico A. Nuevo, IAFE, March-2020.
;
;---------------------------------------------------------------------
function sgradp1_function, Ne0, Te0
; 's' defined in the next line is the emissivity/or/FBE as a function of Te 
  s =  s_function(Ne0,Te0) ;*0. + 1. ;ACTIVATE to make s=1.
  
  gradP = grad_p_function(Ne0, Te0)
  gradp1= reform(gradP(*,*,0))
  RESULT = s*gradp1
  return, RESULT
end

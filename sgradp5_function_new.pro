;---------------------------------------------------------------------
;
; Brief description:
;
; This function computes the product of the s and dp/dq functions.
; dp/dq is calculated directly (without using grad_p_function)
;
; INPUTS: Ne0, Te0.
;
;
; History:  V1.0, Federico A. Nuevo, IAFE, April-2020.
;
;---------------------------------------------------------------------


function sgradp5_function_new,Ne0,Te0
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q


; s' defined in the next line is the emissivity/or/FBE as a function of Te 
  s =  s_function(Ne0,Te0)      ;*0. + 1. ;ACTIVATE to make s=1.
  
  expT    = ((Te0-Tem)/SigTe)
  expN    = ((Ne0-Nem)/SigNe)
  expT2   = expT^2
  expN2   = expN^2
  expTN   = expT * expN
  p_value = (1./(2.*!pi*sigTe*sigNe*sqrt(1.-q^2)))*$
            exp( - (1./2./(1.-q^2))*( expT2 + expN2 - 2.*q*expTN ) )
  
  gradp5  = p_value/(1.-q^2) * ( q  - q/(1.-q^2) * ( expT2 + expN2 - 2.*q*expTN )+ expTN ) ; dP/dq
  
  RESULT = s*gradp5

  

  return,RESULT
end

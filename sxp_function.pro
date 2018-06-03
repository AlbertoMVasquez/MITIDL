;---------------------------------------------------------------------
;
; Brief description:
;
; This function computes the product of the s and p functions.
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

function sxp_function, Ne0, Te0
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  s = (s_function(Te0,Ne0,r0,fip_factor))(0)
  p =  p_function(Te0,Ne0,Tem,Nem,SigTe,SigNe,q)
  return, s*p
end

;---------------------------------------------------------------------
;
; Brief description:
;
; This function computes the product of the s and p functions.
;
; INPUTS: Ne0, Te0.
;
; Note that here Te0 is a 1D-array, as this fucntion is called from
; int_2D to integrate it. See IMPORTANT NOTE in routine e_function.pro.
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

function sxp_function, Ne0, Te0
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  s = reform( (s_function(Ne0,Te0))(*,0,0) );*0. + 1. ;ACTIVATE to make s=1.
  p =          p_function(Ne0,Te0)
  return, s*p
end

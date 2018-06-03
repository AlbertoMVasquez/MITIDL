;---------------------------------------------------------------------
;
; Brief description:
;
; This function computes the product of the s and p2 functions.
;
; INPUTS: Z1, Z2
;
; Note that here Z2 is a 1D-array, as this fucntion is called from
; int_2D to integrate it. See IMPORTANT NOTE in routine e2_function.pro.
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

function sxp2_function, Z1, Z2
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  Ne0 = Nem + SigNe*Z1
  Te0 = Tem + SigTe*( q*Z1 + sqrt(1.-q^2)*Z2)
  s  = reform( ( s_function(Ne0,Te0))(*,0,0) );*0. + 1. ;ACTIVATE to make s=1.
  p  = (1./2./!pi)*exp(-0.5*(Z1^2+Z2^2)) 
  return, s*p
end

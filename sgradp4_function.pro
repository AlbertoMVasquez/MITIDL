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
function sgradp4_function, Ne0, Te0
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common tomographic_measurements, y0, y, measurement_type, i_measurement
; 's' defined in the next line is the emissivity/or/FBE as a function of Te 
  s =  s_function(Ne0,Te0) ;*0. + 1. ;ACTIVATE to make s=1.
  gradP = grad_p(Ne0, Te0)
  RESULT = s*(gradP(3))
  return, RESULT
end

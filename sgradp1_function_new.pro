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
function sgradp1_function_new, Ne0, Te0
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common tomographic_measurements, y0, y, measurement_type, i_measurement
  
; 's' defined in the next line is the emissivity/or/FBE as a function of Te 
  s =  s_function(Ne0,Te0) ;*0. + 1. ;ACTIVATE to make s=1.
  
  ;gradP = grad_p(Ne0, Te0)
  ;gradp1= reform(gradP(*,*,0))
  expT    = ((Te0-Tem)/SigTe)^2
  expN    = ((Ne0-Nem)/SigNe)^2
  expTN   = (Te0-Tem)*(Ne0-Nem)/(SigTe*SigNe)
  p_value = (1./(2.*!pi*sigTe*sigNe*sqrt(1.-q^2)))*$
            exp( - (1./2./(1.-q^2))*( expT + expN - 2.*q*expTN ) )
  
  gradp1  = p_value * ( 1./(1.-q^2)/sigNe * ( (Ne0-Nem)/SigNe + q * (Te0-Tem)/SigTe )) ; dP/dNm
  
  RESULT = s*gradp1
  return, RESULT
end

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
function sgradp5_function_new, Ne0, Te0
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common tomographic_measurements, y0, y, measurement_type, i_measurement
  
; 's' defined in the next line is the emissivity/or/FBE as a function of Te 
  s =  s_function(Ne0,Te0) ;*0. + 1. ;ACTIVATE to make s=1.
  
  
  expT    = ((Te0-Tem)/SigTe)^2
  expN    = ((Ne0-Nem)/SigNe)^2
  expTN   = (Te0-Tem)*(Ne0-Nem)/(SigTe*SigNe)
  p_value = (1./(2.*!pi*sigTe*sigNe*sqrt(1.-q)))*$
            exp( - (1./2./(1.-q^2))*( expT + expN - 2.*q*expTN ) )
  
  gradp5  = p_value * ( q/(1.-q^2) -q/(1.-q^2)^2/2.*( expT + expN - 2.*q*expTN )-1./(1.-q^2)*expTN) ; dP/dq
  
  RESULT = s*gradp5
  return, RESULT
end

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
;           V1.1. F.A. Nuevo.
;           V1.2. A.M.Vasquez.
;                 Eliminated unnecesary common blocks.
;---------------------------------------------------------------------
function sp_function, Ne0, Te0
; 's' defined in the next line is the emissivity/or/FBE as a function of Te 
  s = s_function(Ne0,Te0)
 ;s = 1. ; activate this line to make s=1.
  p = p_function(Ne0,Te0)
  RESULT = s*p
  return, RESULT
end

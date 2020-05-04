;---------------------------------------------------------------------
;
; Brief description:
;
; This function computes the emissivity e of a line/band in a voxel.
;
; INPUTS: values for:
; Z1_Limits..
;
; OUTPUTS:
;
; emissivity e in units of:
;
; For lines: [erg cm-3 sec-1 sr-1]
; or
; For EUV bands:
;
; IMPORTANT NOTE: In its current implementation, this function
; computes the 2D integral using the INT_2D.PRO routine, assuming
; order=0, which means that the integration order is dy-dx, i.e the
; independent variable is x (external integral), and the dependent
; variable is y, where here x=Ne and y=Te.
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------
function e2_function, parameters
  common Ylimits, Y_Limits
  common NT_limits, Ne0_Limits, Te0_Limits
  Nem        = parameters[0]
  fip_factor = parameters[1]
  Tem        = parameters[2]
  SigTe      = parameters[3]
  SigNe      = parameters[4]
  q          = parameters[5]
  Z1_Limits  = (Ne0_Limits-Nem)/SigNe
   Y_Limits  = (Te0_Limits-Tem)/SigTe
  RESULT = INT_2D('sp2_function',Z1_Limits,'z2_limits',96,/double,order=0)
  return, RESULT
end

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

function e2_function, Ne0_Limits, Te0_Limits
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common Ylimits,Y_Limits
  Z1_Limits = (Ne0_Limits-Nem)/SigNe
   Y_Limits = (Te0_Limits-Tem)/SigTe
  RESULT = INT_2D('sxp2_function',Z1_Limits,'z2_limits',96,/double,order=0)
  return, RESULT
end

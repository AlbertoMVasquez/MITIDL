;---------------------------------------------------------------------
;
; Brief description:
;
; This function computes the emissivity e of a line/band in a voxel.
;
; INPUTS: values for:
; Ne0_Limits, Te0_Limits.
;
; OUTPUTS:
;
; emissivity e in units of:
;
; For lines: [erg cm-3 sec-1 sr-1]
; or
; For EUV bands:
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

function e_function, Ne0_Limits, Te0_Limits
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  RESULT = INT_2D('sxp_function',Ne0_Limits,'Te_Limits',96,/double)
  return, RESULT
end

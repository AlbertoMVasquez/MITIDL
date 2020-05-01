;---------------------------------------------------------------------;
; Brief description:
;
; This function computes the emissivity e of a line/band in a voxel.
;
; ARGUMENT: 
; parameters: a 1D array of 6 elements: [Nem, fip_factor, Tem, SigTe, SigNe, q]
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
;           v1.1, test that is the same dxdy-order than dydx-order
;---------------------------------------------------------------------
function e_function, parameters,order=order
  common NT_limits, Ne0_Limits, Te0_Limits
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  Nem        = parameters[0]
  fip_factor = parameters[1]
  Tem        = parameters[2]
  SigTe      = parameters[3]
  SigNe      = parameters[4]
  q          = parameters[5]
  if not keyword_set(order) then $
     RESULT = INT_2D('sp_function',Ne0_Limits,'te_limits',96,/double,order=0) ; dydx-order
  if     keyword_set(order) then $
     RESULT = INT_2D('sp_function',Te0_Limits,'ne_limits',96,/double,order=1) ; dxdy-order 
 return, RESULT
end

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
function e_function_cs, k, parameters
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common NT_arrays,Ne_array,Te_array
  common sk_array,sk
 

  Nem        = parameters[0]
  fip_factor = parameters[1]
  Tem        = parameters[2]
  SigTe      = parameters[3]
  SigNe      = parameters[4]
  q          = parameters[5]
 
  Ne0= Ne_array
  Te0= Te_array
  dNe = Ne0(1)-Ne0(0)
  dTe = Te0(1)-Te0(0)

  tmp = sk(k,*,*) * p_function_loop(Ne0,Te0)
  result = total(tmp)* dNe * dTe


  return, RESULT
end

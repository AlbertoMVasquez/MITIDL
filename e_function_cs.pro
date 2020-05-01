;---------------------------------------------------------------------;
; Brief description:
;
; This function computes the emissivity e_k for a given k-th line/band in a voxel.
; The routine uses CS to calculate the double integrals in the
; definition of e_k.
;
; CS: \int f(x,y) dx dy >  \Sum_{i,j} f(x_i,y_j) Dx Dy
;
; ARGUMENT: 
; k         : the index of the emissivity
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
; History:  V1.0, Federico A. Nuevo, IAFE, April-2020.
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
  ; Ne and Te grid
  Ne0= Ne_array
  Te0= Te_array
  dNe = Ne0(1)-Ne0(0)
  dTe = Te0(1)-Te0(0)
  ; CS : \Sum_{i,j} f(x_i,y_j) Dx Dy 
  tmp = fip_factor*sk(k,*,*)*p_function_loop(Ne0,Te0)
  result = total(tmp)*dNe*dTe
  return, RESULT
end

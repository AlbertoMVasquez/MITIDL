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
; History:  V1.0, F.A. Nuevo, IAFE, April-2020.
;           V1.1, A.M. Vasquez, IAFE, April-2020.
;                 Simplified.
;---------------------------------------------------------------------
function e_function_cs, k, parameters
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common NT_arrays,Ne_array,Te_array,dNe_array,dTe_array,dTN
  common sk_over_fip_factor_array,sk_over_fip_factor    

  Nem        = parameters[0]
  fip_factor = parameters[1]
  Tem        = parameters[2]
  SigTe      = parameters[3]
  SigNe      = parameters[4]
  q          = parameters[5]
 
  tmp = fip_factor*reform(sk_over_fip_factor(k,*,*))*p_function_cs(Ne_array,Te_array) ; array NTe x NNe
  ;result = Sum_{i,j} f(x_i,y_j) dx_i dy_j 
  return, total( tmp * dTN ) ; = total(dTe_array*(tmp#dNe_array))
end

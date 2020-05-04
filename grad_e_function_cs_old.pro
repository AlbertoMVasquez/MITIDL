;---------------------------------------------------------------------;
; Brief description:
;
; This function computes the gradient of the emissivity respect to the
; parameters using CS to calculate the double integral.

;
; ARGUMENT: 
; K         : index of the emissivity
; parameters: a 1D array of 6 elements: [Nem, fip_factor, Tem, SigTe, SigNe, q]
;
;
;
; History:  V1.0, Federico A. Nuevo, IAFE, April-2020.
;
;---------------------------------------------------------------------
function grad_e_function_cs,k, parameters
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common NT_limits, Ne0_Limits, Te0_Limits
  common sk_array,sk
  common NT_arrays,Ne_array,Te_array,dNe_array,dTe_array

  
  Nem        = parameters[0]
  fip_factor = parameters[1]
  Tem        = parameters[2]
  SigTe      = parameters[3]
  SigNe      = parameters[4]
  q          = parameters[5]
  
  RESULT=parameters*0d

  ; Ne and Te grid
  Ne0= Ne_array
  Te0= Te_array
  dNe = dNe_array
  dTe = dTe_array

  ; P gradient
  gradP=grad_p_function_cs(Ne0,Te0)
  
  ; Calculate the double integrals with CS 
  Result(0) = total ( fip_factor * sk(k,*,*) * gradP (*,*,0) ) * dNe * dTe
  Result(1) = e_function_cs(k,parameters)/fip_factor
  Result(2) = total ( fip_factor * sk(k,*,*) * gradP (*,*,1) ) * dNe * dTe
  Result(3) = total ( fip_factor * sk(k,*,*) * gradP (*,*,2) ) * dNe * dTe
  Result(4) = total ( fip_factor * sk(k,*,*) * gradP (*,*,3) ) * dNe * dTe
  Result(5) = total ( fip_factor * sk(k,*,*) * gradP (*,*,4) ) * dNe * dTe  

  return, RESULT
end

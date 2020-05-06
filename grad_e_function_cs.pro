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
;           V1.1, A.M. Vasquez, IAFE, April-2020.
;                 Simplified.
;
;---------------------------------------------------------------------
function grad_e_function_cs,k, parameters
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common NT_limits, Ne0_Limits, Te0_Limits
  common NT_arrays,Ne_array,Te_array,dNe_array,dTe_array,dTN
  common sk_over_fip_factor_array,sk_over_fip_factor    
  
  Nem        = parameters[0]
  fip_factor = parameters[1]
  Tem        = parameters[2]
  SigTe      = parameters[3]
  SigNe      = parameters[4]
  q          = parameters[5]
  
  RESULT=parameters*0d

  n_par = n_elements(parameters)

  ; P gradient
  gradP=grad_p_function_cs(Ne_array,Te_array)
  
  ; Calculate the double integrals with CS 
  ;Result(0) = total ( dTe_array * fip_factor * ((reform(sk_over_fip_factor(k,*,*)) * reform(gradP (*,*,0)) )  #  dNe_array))
  ;Result(1) = e_function_cs(k,parameters)/fip_factor
  ;Result(2) = total ( dTe_array * fip_factor * ((reform(sk_over_fip_factor(k,*,*)) * reform(gradP (*,*,1)) )  #  dNe_array))
  ;Result(3) = total ( dTe_array * fip_factor * ((reform(sk_over_fip_factor(k,*,*)) * reform(gradP (*,*,2)) )  #  dNe_array))
  ;Result(4) = total ( dTe_array * fip_factor * ((reform(sk_over_fip_factor(k,*,*)) * reform(gradP (*,*,3)) )  #  dNe_array))
  ;Result(5) = total ( dTe_array * fip_factor * ((reform(sk_over_fip_factor(k,*,*)) * reform(gradP (*,*,4)) )  #  dNe_array))

  Result(0) = total( fip_factor * reform(sk_over_fip_factor(k,*,*)) * reform(gradP (*,*,0)) * dTN )
  Result(1) = e_function_cs(k,parameters)/fip_factor
  Result(2) = total( fip_factor * reform(sk_over_fip_factor(k,*,*)) * reform(gradP (*,*,1)) * dTN )
  Result(3) = total( fip_factor * reform(sk_over_fip_factor(k,*,*)) * reform(gradP (*,*,2)) * dTN )
  Result(4) = total( fip_factor * reform(sk_over_fip_factor(k,*,*)) * reform(gradP (*,*,3)) * dTN )
  Result(5) = total( fip_factor * reform(sk_over_fip_factor(k,*,*)) * reform(gradP (*,*,4)) * dTN )
  
  return, RESULT
end

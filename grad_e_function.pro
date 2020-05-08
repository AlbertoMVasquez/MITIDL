;---------------------------------------------------------------------;
; Brief description:
;
; This function computes the gradient of the emissivity of a line/band in a voxel.
;
; ARGUMENT: 
; parameters: a 1D array of 6 elements: [Nem, fip_factor, Tem, SigTe, SigNe, q]
;
; OUTPUTS:
;
; gradient of the emissivity (vector of 6 elements) respect to parameters
;
;
; History:  V1.0, F.A. Nuevo, IAFE, March-2020.
;
;---------------------------------------------------------------------
function grad_e_function, parameters
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common NT_limits, Ne0_Limits, Te0_Limits

  Nem        = parameters[0]
  fip_factor = parameters[1]
  Tem        = parameters[2]
  SigTe      = parameters[3]
  SigNe      = parameters[4]
  q          = parameters[5]

  nodes=96
  result=parameters*0d
  
     result(0) = INT_2D('sgradp1_function',Ne0_Limits,'te_limits',nodes,/double,order=0)
     result(1) = e_function(parameters)/fip_factor
     result(2) = INT_2D('sgradp2_function',Ne0_Limits,'te_limits',nodes,/double,order=0)
     result(3) = INT_2D('sgradp3_function',Ne0_Limits,'te_limits',nodes,/double,order=0)
     result(4) = INT_2D('sgradp4_function',Ne0_Limits,'te_limits',nodes,/double,order=0)
     result(5) = INT_2D('sgradp5_function',Ne0_Limits,'te_limits',nodes,/double,order=0)
 
  return, RESULT
end

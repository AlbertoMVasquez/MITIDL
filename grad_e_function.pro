;---------------------------------------------------------------------;
; Brief description:
;
; This function computes the gradient of the emissivity respect to the parameters.
;
; ARGUMENT: 
; parameters: a 1D array of 6 elements: [Nem, fip_factor, Tem, SigTe, SigNe, q]
;
;
;
; History:  V1.0, Federico A. Nuevo, IAFE, March-2020.
;
;---------------------------------------------------------------------
function grad_e_function, parameters
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common NT_limits, Ne0_Limits, Te0_Limits
  common tomographic_measurements, y0, y, measurement_type, i_measurement
  
  Nem        = parameters[0]
  fip_factor = parameters[1]
  Tem        = parameters[2]
  SigTe      = parameters[3]
  SigNe      = parameters[4]
  q          = parameters[5]

  nodes=48
  result=parameters*0d
  new = 1

  if new eq 0 then begin
     result(0) = INT_2D('sgradp1_function',Ne0_Limits,'te_limits',nodes,/double,order=0)
     result(1) = e_function(parameters)/fip_factor
     result(2) = INT_2D('sgradp2_function',Ne0_Limits,'te_limits',nodes,/double,order=0)
     result(3) = INT_2D('sgradp3_function',Ne0_Limits,'te_limits',nodes,/double,order=0)
     result(4) = INT_2D('sgradp4_function',Ne0_Limits,'te_limits',nodes,/double,order=0)
     result(5) = INT_2D('sgradp5_function',Ne0_Limits,'te_limits',nodes,/double,order=0)
  endif
  if new eq 1 then begin
     result(0) = INT_2D('sgradp1_function_new',Ne0_Limits,'te_limits',nodes,/double,order=0)
     result(1) = e_function(parameters)/fip_factor
     result(2) = INT_2D('sgradp2_function_new',Ne0_Limits,'te_limits',nodes,/double,order=0)
     result(3) = INT_2D('sgradp3_function_new',Ne0_Limits,'te_limits',nodes,/double,order=0)
     result(4) = INT_2D('sgradp4_function_new',Ne0_Limits,'te_limits',nodes,/double,order=0)
     result(5) = INT_2D('sgradp5_function_new',Ne0_Limits,'te_limits',nodes,/double,order=0)
  endif
  if new eq 2 then begin
     result(0) = INT_2D('sgradp1',Ne0_Limits,'te_limits',nodes,/double,order=0)
     result(1) = e_function(parameters)/fip_factor
     result(2) = INT_2D('sgradp2',Ne0_Limits,'te_limits',nodes,/double,order=0)
     result(3) = INT_2D('sgradp3',Ne0_Limits,'te_limits',nodes,/double,order=0)
     result(4) = INT_2D('sgradp4',Ne0_Limits,'te_limits',nodes,/double,order=0)
     result(5) = INT_2D('sgradp5',Ne0_Limits,'te_limits',nodes,/double,order=0)
  endif

  
  return, RESULT
end

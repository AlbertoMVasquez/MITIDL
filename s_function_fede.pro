;---------------------------------------------------------------------
;
; Brief description:
;
; This function returns:
;
;  - For lines: s is the line emissivity function evaluated at the
;    provided values Ne0, Te0 and r0, using the provided function G.
;          
;  - For EUV bands: s is the FBE function evaluated at the
;    provided values for Ne0 and Te0, using the provided function TRF.
;
; G and TRF are tri-linearly interpolated onto the grid (Ne,Te,r).
;
;
; INPUTS: values for Ne, Te, rad and fip_factor (default value is 1).
;
; Note: fip_factor=1 means using the abundance set assumed in
; computing the loop-up table for G/TRF.
;
; OUTPUT:
;
; RESULT = [s0]
;
; With the s0 value returned in the following units:
;
; For lines: Emissivity [erg cm-3 sec-1 sr-1]
; or
; For EUV bands: FBE    [ph cm-3 sec-1 sr-1]
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;           V1.1, elimino derivadas, reescribo para usar un solo loop 
;---------------------------------------------------------------------

function s_function_fede, Ne0, Te0
  common G_table, G, T_e, N_e, r, photT
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common dimensions, NTe, NNe
  common tomographic_measurements, y0, y, measurement_type, i_measurement
 

  NTe = 1
  NNe = 1
  if (size(Te0))(0) eq 1 then NTe = (size(Te0))(1)
  if (size(Ne0))(0) eq 1 then NNe = (size(Ne0))(1)

; Linearly interpolate G from the look-up table, and compute its derivatives:
  RESULT_s = dblarr(NTe,NNe)
; calcula G(Ne,Te) en Te0 y Ne0 arrays  
  RESULT_g = g_function_fede(Te0,Ne0) ; [erg/PH cm+3 sec-1]; esto ya es un array de NTe X NNe
  for iTe=0,NTe-1 do begin
     ; Compute emissivity assuming isotropic emission:
     result_s[iTe,*]=(fip_factor/4./!pi)*Ne0^2*result_g[iTe,*] ; [erg/PH cm-3 sec-1 sr-1]
  endfor

  return, RESULT_s
end

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
; Derivatives ds/dTe, ds/dNe, ds/dr are also returned.
;
; INPUTS: values for Ne, Te, rad and fip_factor (default value is 1).
;
; Note: fip_factor=1 means using the abundance set assumed in
; computing the loop-up table for G/TRF.
;
; OUTPUT:
;
; RESULT = [s0, ds_dTe, ds_dNe, ds_dr]
;
; With the s0 value returned in the following units:
;
; For lines: Emissivity [erg cm-3 sec-1 sr-1]
; or
; For EUV bands: FBE    [erg cm-3 sec-1 sr-1]
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

function s_function, Ne0, Te0
  common G_table, G, T_e, N_e, r, photT
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common dimensions, NTe, NNe
  common tomographic_measurements, y0, y, measurement_type, i_measurement
  
; Set default fip_factor:
  if not keyword_set(fip_factor) then fip_factor = 1.0 
stop
; Linearly interpolate G from the look-up table, and compute its derivatives:
  RESULT_g = g_function(Te0,Ne0)
stop
  RESULT_s = dblarr(NTe,NNe,4)
  for iTe=0,NTe-1 do begin
  for iNe=0,NNe-1 do begin

  G0     = RESULT_g(iTe,iNe,0) ; [erg cm+3 sec-1     ]
  dG_dTe = RESULT_g(iTe,iNe,1) ; [erg cm+3 sec-1 K-1 ]
  dG_dNe = RESULT_g(iTe,iNe,2) ; [erg cm+6 sec-1     ]
  dG_dr  = RESULT_g(iTe,iNe,3) ; [erg cm+3 sec-1 Rs-1]
  
; Compute emissivity and its derivatives, assuming isotropic emission:
  s0     = (fip_factor/4./!pi)*Ne0^2*G0                   ; [erg cm-3 sec-1 sr-1     ]
  ds_dTe = (fip_factor/4./!pi)*Ne0^2*dG_dTe               ; [erg cm-3 sec-1 K-1  sr-1]
  ds_dNe = (fip_factor/4./!pi)*(2.*Ne0*G0 + Ne0^2*dG_dNe) ; [erg      sec-1      sr-1]
  ds_dr  = (fip_factor/4./!pi)*Ne0^2*dG_dr                ; [erg cm-3 sec-1 Rs-1 sr-1]

  RESULT_s[iTe,iNe,*] = [s0, ds_dTe, ds_dNe, ds_dr]
  endfor
  endfor
  
  return, RESULT_s
end

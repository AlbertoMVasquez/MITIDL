;---------------------------------------------------------------------
;
; Brief description:
;
; This function returns:
;
;  - For lines: the emissivity (s) provided values for Ne, Te and r,
;               and the contribution function G.
;
;  - For EUV bands: the FBE (s) provided values for Ne, Te, and the
;    temperature response TRF
;
; In both cases the routine tri-linearly interpolates G/TRF into
; (Ne,Te,r).
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
; With s0 having the following units:
;
; For lines: Emissivity [erg cm-3 sec-1 sr-1]
; or
; For EUV bands: FBE [units?]
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

function s_function,Te0,Ne0,r0,fip_factor,derivatives=derivatives
  common G_table,G,T_e,N_e,r,photT

; Set default fip_factor:
  if not keyword_set(fip_factor) then fip_factor = 1.0 

; Tri-linearly interpolate G from the look-up table, and compute its
; derivatives:
  RESULT_g = g_function(Te0,Ne0,r0) 
  G0     = RESULT_g(0) ; [erg cm+3 sec-1     ]
  dG_dTe = RESULT_g(1) ; [erg cm+3 sec-1 K-1 ]
  dG_dNe = RESULT_g(2) ; [erg cm+6 sec-1     ]
  dG_dr  = RESULT_g(3) ; [erg cm+3 sec-1 Rs-1]
  
; Compute emissivity and its derivatives, assuming isotropic emission:
  s0     = (fip_factor/4./!pi)*Ne0^2*G0                   ; [erg cm-3 sec-1 sr-1     ]
  ds_dTe = (fip_factor/4./!pi)*Ne0^2*dG_dTe               ; [erg cm-3 sec-1 K-1  sr-1]
  ds_dNe = (fip_factor/4./!pi)*(2.*Ne0*G0 + Ne0^2*dG_dNe) ; [erg      sec-1      sr-1]
  ds_dr  = (fip_factor/4./!pi)*Ne0^2*dG_dr                ; [erg cm-3 sec-1 Rs-1 sr-1]

  return, [s0, ds_dTe, ds_dNe, ds_dr]
end

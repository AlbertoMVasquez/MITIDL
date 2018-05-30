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
; Derivatives: ds/dTe and ds/dNe are also returned if the keyword
; /derivatives is set in the calling sequence.
;
; INPUTS: values for Ne, Te, rad and fip_factor (default value is 1).
;
; Note: fip_factor=1 means using the abundance set assumed in
; computing the loop-up table for G/TRF.
;
; OUTPUTS:
;
; Emissivity [erg cm-3 sec-1 sr-1]
; or
; FBE [units?]
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

function s_function,Te0=Te0,Ne0=Ne0,r0=r0,fip_factor=fip_factor,derivatives=derivatives
  common G_table,G,T_e,N_e,r,photT

  ; Set default fip_factor:
  if not keyword_set(fip_factor) then fip_factor = 1.0 

  ; Tri-linearly interpolate G from the look-up table:
  G0 = G_function(Te0,Ne0,r0) ; [erg cm+3 sec-1]

  ; Compute emissivity assmuing isotropic emission:
  s0 = fip_factor * Ne0^2 * G0 / (4.*!pi) ; [erg cm-3 sec-1 sr-1]

  return,s0
end

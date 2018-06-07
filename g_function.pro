;---------------------------------------------------------------------
;
; Brief description:
;
; This function returns:
;
;  - For lines: the contribution function G provided values for Ne, Te and r.
;
;  - For EUV bands: the temperature response function TRF provided
;    values for Ne, Te.
;
;  - For f=G/TRF, it also returns the derivatives:
;    df/dTe, df/dNe, df/dr, evaluated at (Te0,Ne0,r0)
;
; In both cases the routine tri-linearly interpolates G/TRF into (Te0,Ne0,r0).
;
; INPUTS: values for Te0, Ne0, r0
;
; OUTPUTS:
;
; RESULT = [f, df/dTe, df/dNe, df/dr]
;
; where f is:
;
; G [erg cm+3 sec-1]
; or
; TRF [units?]
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

function g_function, Te0, Ne0, emissionline=emissionline, euvband=euvband
  common G_table,G,T_e,N_e,r,photT
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common dimensions,NTe,NNe  
  
  NTe = 1
  NNe = 1
  if (size(Te0))(0) eq 1 then NTe = (size(Te0))(1)
  if (size(Ne0))(0) eq 1 then NNe = (size(Ne0))(1)
  RESULT = dblarr(NTe,NNe,4)
if keyword_set(emissionline) then begin
  for iTe=0,NTe-1 do begin
     for iNe=0,NNe-1 do begin
        RESULT[iTe,iNe,*] = findval3d_function(G,T_e,N_e,r,Te0[iTe],Ne0[iNe],r0)
     endfor  
  endfor
endif
if keyword_set(euvband) then begin
   for iNe=0,NNe-1 do begin      
      RESULT[*,iNe,0] = interpol(G, T_e, Te0)
      dG_dTe = 0.
      RESULT[*,iNe,1:3] = [dG_dTe,0.,0.]
      ;; NEED TO COMPUTE dG_dTe above!
      ;; Maybe EASIER to compute G(Te,Ne), and make a new routine to
      ;; compute both derivatives...
   endfor  
endif
  
  return, RESULT
end

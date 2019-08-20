;---------------------------------------------------------------------
;
; Brief description:
;
; This function returns:
;
;  - For lines: the contribution function G provided values for Ne0, Te0 and r0.
;
;  - For EUV bands: the temperature response function TRF provided values for Ne0 and Te0.
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
; G   [ERG    cm^+3 sec^-1]
; or
; TRF [PHOTON cm^+3 sec^-1]
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------
function g_function, Te0, Ne0
  common G_table, G, T_e, N_e, r, photT
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common dimensions, NTe, NNe  
  common tomographic_measurements, y0, y, measurement_type, i_measurement
  NTe = 1
  NNe = 1
  if (size(Te0))(0) eq 1 then NTe = (size(Te0))(1)
  if (size(Ne0))(0) eq 1 then NNe = (size(Ne0))(1)
  RESULT = dblarr(NTe,NNe,4)
  CASE measurement_type[i_measurement] OF
  1: BEGIN
  print,'Selected EMISSION LINE in g_function.pro'
  for iTe=0,NTe-1 do begin
     for iNe=0,NNe-1 do begin
        RESULT[iTe,iNe,*] = findval3d_function(G,T_e,N_e,r,Te0[iTe],Ne0[iNe],r0)
     endfor  
  endfor
  END
  2: BEGIN
  print,'Selected EUV BAND in g_function.pro'
   for iNe=0,NNe-1 do begin
     ;print,'Selected EUVBAND in g_function.pro'
      dG_dTe       = deriv(T_e,G)
      G_atTe0      = interpol( G    , T_e, Te0)
      dG_dTe_atTe0 = interpol(dG_dTe, T_e, Te0)
      RESULT[*,iNe,0] = G_atTe0
      RESULT[*,iNe,1] = dG_dTe_atTe0
      ; Note that dG_dNe and dG_dr are set to ZERO.
      ; In a next version G will be a function of (Te,Ne), even if
      ; very weakly dependent on Ne, will compute dG_dNe too.
   endfor
   END
  return, RESULT
end

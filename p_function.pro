;---------------------------------------------------------------------
;
; Brief description:
;
; For values of Te0 and Ne0, this function returns the value of the
; bivariate normal joint temperature-density probability distribution
; with user-provided parameters: Tem, Nem, SigTe, SigNe, q.
;
; INPUTS:
;  Te0: float with the electron temperature       in units of [K]
;  Ne0: float with the electron density           in units of [cm-3]
;  Tem: float with the mean electron temperature  in units of [K]
;  Nem: float with the mean electron density      in units of [cm-3]
;SigTe: float with the electron temperature StDev in units of [K]
;SigNe: float with the electron density     StDev in units of [cm-3]
;    q: temperature-density correlation dimensionless coefficient:
;       0 means no-correlation, 1 is full-correlation.
;
; OUTPUT:
;       value of the probablility function p, normalized to volume 1
;       over [-Infty,+Infty] space for both Te and Ne.
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;           V1.1, Federio A. Nuevo, IAFE, Mach-2020
;                 The normalization factor of P lacked the square
;                 power in q. 
;---------------------------------------------------------------------
function p_function, Ne0, Te0
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  eps=1.e-8
  if abs(q-1.) le eps then begin
     print,'Ne-Te correlation can not be 1.'
     stop
  endif
  expT2    = ((Te0-Tem)/SigTe)^2
  expN2    = ((Ne0-Nem)/SigNe)^2
  expTN    = (Te0-Tem)*(Ne0-Nem)/(SigTe*SigNe)
  p_value  = (1./(2.*!pi*sigTe*sigNe*sqrt(1.-q^2)))*$
             exp( - (1./2./(1.-q^2))*( expT2 + expN2 - 2.*q*expTN ) )  
  return,p_value
end

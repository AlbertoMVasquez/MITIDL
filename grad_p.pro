;---------------------------------------------------------------------
;
; Brief description:
;
; For values of Te0 and Ne0, this function returns the value of the
; gradient of the bivariate normal joint temperature-density probability distribution
; respect to: Tem, Nem, SigTe, SigNe, q.
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
; 
;
; History:  V1.0, Federico Nuevo, 15-02-20.
;
;---------------------------------------------------------------------

function grad_p, Ne0, Te0
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q


  grad=dblarr(5)

  eps=1.e-8
  if abs(q-1.) le eps then begin
     print,'Ne-Te correlation can not be 1.'
     stop
  endif
  expT    = ((Te0-Tem)/SigTe)^2
  expN    = ((Ne0-Nem)/SigNe)^2
  expTN   = (Te0-Tem)*(Ne0-Nem)/(SigTe*SigNe)
  p_value = (1./(2.*!pi*sigTe*sigNe*sqrt(1.-q)))*$
            exp( - (1./2./(1.-q^2))*( expT + expN - 2.*q*expTN ) )  

  grad(0) = p_value * ( 1./(1.-q^2)/sigNe * ( (Ne0-Nem)/SigNe + q * (Te0-Tem)/SigTe )) ; dP/dNm

  grad(1) = p_value * ( 1./(1.-q^2)/sigTe * ( (Te0-Tem)/SigTe + q * (Ne0-Nem)/SigNe )) ; dP/dTm

  grad(2) = p_value * ( 1./(1.-q^2)/sigNe * ( expN + q * expTN - ( 1.-q^2)) ; dP/sigN
 
  grad(3) = p_value * ( 1./(1.-q^2)/sigTe * ( expT + q * expTN - ( 1.-q^2)) ; dP/sigT

  grad(4) = p_value * ( q/(1.-q^2) -q/(1.-q^2)^2/2.*( expT + expN - 2.*q*expTN )-1./(1.-q^2)*expTN) ; dP/dq

  return,grad
end

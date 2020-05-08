;---------------------------------------------------------------------
;
; Brief description:
;
; For values of Te0 and Ne0, this function returns the value of the
; gradient of the bivariate normal joint temperature-density
; probability distribution respect to: Tem, Nem, SigTe, SigNe, q.
;
; INPUTS:
;  Te0: float with the electron temperature       in units of [K]
;  Ne0: float with the electron density           in units of [cm-3]
;  (in COMMON parameters) 
;  Tem: float with the mean electron temperature  in units of [K]
;  Nem: float with the mean electron density      in units of [cm-3]
;SigTe: float with the electron temperature StDev in units of [K]
;SigNe: float with the electron density     StDev in units of [cm-3]
;    q: temperature-density correlation dimensionless coefficient:
;       0 means no-correlation, 1 is full-correlation.
;
; OUTPUT: An array of NTe x NNe x5 components with the gradient of P
;         respect to the parameters vector.
;
; NOTE: Te0 can be an 1D array and Ne0 an scalar, or the opposite.
;       But, is not possible that both be arrays.
;
; History:  V1.0, F.A. Nuevo,   15-02-20.
;           V1.1, A.M. Vasquez, 22-04-20.
;---------------------------------------------------------------------

function grad_p_function, Ne0, Te0
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  
  eps=1.e-8
  if abs(q-1.) le eps then begin
     print,'Ne-Te correlation can not be 1.'
     stop
  endif 

  NTe = 1
  NNe = 1
  if (size(Te0))(0) eq 1 then NTe = (size(Te0))(1)
  if (size(Ne0))(0) eq 1 then NNe = (size(Ne0))(1)

  grad=dblarr(NTe,NNe,5)
  
  expT    = (Te0-Tem)/SigTe
  expN    = (Ne0-Nem)/SigNe
  expT2   = expT^2
  expN2   = expN^2
  expTN   = expT*expN
  
  p_value = (1./(2.*!pi*sigTe*sigNe*sqrt(1.-q^2)))*$
            exp( - (1./2./(1.-q^2))*(expT2 + expN2 - 2.*q*expTN) )  

; Correct NaNs due to INFINITY exp(expTN)
  index_nan = where(finite(exp((1./(1.-q^2))*q*expTN)) ne 1)
  if index_nan(0) ne -1 then p_value(index_nan) = 0.

  grad(*,*,0) = p_value * (1./sigNe)*(expN  - q*expT)                             ; \propto dP/dNm
  grad(*,*,1) = p_value * (1./sigTe)*(expT  - q*expN)                             ; \propto dP/dTm
  grad(*,*,2) = p_value * (1./sigTe)*(expT2 - q*expTN - (1.-q^2))                 ; \propto dP/dsigT
  grad(*,*,3) = p_value * (1./sigNe)*(expN2 - q*expTN - (1.-q^2))                 ; \propto dP/dsigN
  grad(*,*,4) = p_value * (q + expTN - (q/(1.-q^2))*(expN2 + expT2 - 2.*q*expTN)) ; \propto dP/dq
  return, (1/(1-q^2))*grad
end

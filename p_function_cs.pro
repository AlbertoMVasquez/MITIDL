;---------------------------------------------------------------------
;
; Brief description:
;
; For values of 1D arrays  of Te0 and Ne0, this function returns the value of the
; bivariate normal joint temperature-density probability distribution
; with user-provided parameters: Tem, Nem, SigTe, SigNe, q.
;
; INPUTS:
;  Te0: 1D array with the electron temperature       in units of [K]
;  Ne0: 1D array with the electron density           in units of [cm-3]

;  Tem: float with the mean electron temperature  in units of [K]
;  Nem: float with the mean electron density      in units of [cm-3]
;SigTe: float with the electron temperature StDev in units of [K]
;SigNe: float with the electron density     StDev in units of [cm-3]
;    q: temperature-density correlation dimensionless coefficient:
;       0 means no-correlation, 1 is full-correlation.
;
; OUTPUT:
;       value of the probablility function p in a 2D array
;       of NTe X NNe, evaluated in Ne0 and Te0. 
;
; History:  V1.0, Federio A. Nuevo, IAFE, May-2020
;           V1.1, Alberto M. Vásquez, IAFE, May-2020
;           Correct NaNs due to INFINITY exp(expTN)
;         
;---------------------------------------------------------------------
function p_function_cs, Ne0, Te0
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  eps=1.e-8
  if abs(q-1.) le eps then begin
     print,'Ne-Te correlation can not be 1.'
     stop
  endif
  expT2    = double(((Te0-Tem)/SigTe)^2)
  expN2    = double(((Ne0-Nem)/SigNe)^2)
  expTN    = double(((Te0-Tem)/sigTe)#((Ne0-Nem)/sigNe))                       ; array NTe x NNe
  p_value  = (1./(2.*!pi*sigTe*sigNe*sqrt(1.-q^2))) * $                        ; Scalar
             (exp(-(1./2./(1.-q^2))*expT2) # exp(-(1./2./(1.-q^2))*expN2)) * $ ; array NTe x NNe
             exp((1./(1.-q^2))*q*expTN)                                        ; array NTe x NNe

; Correct NaNs due to INFINITY exp(expTN)
  index_nan = where(finite(exp((1./(1.-q^2))*q*expTN)) ne 1)
  if index_nan(0) ne -1 then p_value(index_nan) = 0.

  return,p_value
end

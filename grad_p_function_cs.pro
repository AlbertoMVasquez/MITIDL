;---------------------------------------------------------------------
;
; Brief description:
;
; For values of Te0 and Ne0, this function returns the value of the
; gradient of the bivariate normal joint temperature-density
; probability distribution respect to: Tem, Nem, SigTe, SigNe, q.
; The routine use a double loop in Ne and Te.
; Esta rutina es necesaria si se calculan las integrales dobles
; usando cuadratura simple (e_function_cs)
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
; History:  V1.0, Federico Nuevo, 15-02-20.
;
;---------------------------------------------------------------------

function grad_p_function_cs, Ne0, Te0
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

  expT     = double(((Te0-Tem)/SigTe))
  expN     = double(((Ne0-Nem)/SigNe))
  expT2    = double(((Te0-Tem)/SigTe)^2)
  expN2    = double(((Ne0-Nem)/SigNe)^2)
  expTN    = double(((Te0-Tem)/sigTe)#((Ne0-Nem)/sigNe))
  p_value  = (1./(2.*!pi*sigTe*sigNe*sqrt(1.-q^2)))*$
             exp( - (1./2./(1.-q^2))* expT2 ) # exp( - (1./2./(1.-q^2))* expN2 ) * $
             exp((1./(1.-q^2))*q*expTN ) 

  PexpN =dblarr(NTe,NNe)
  PexpT =dblarr(NTe,NNe)
  PexpN2=dblarr(NTe,NNe)
  PexpT2=dblarr(NTe,NNe)

  for iTe=0,NTe-1 do begin
     pexpN (iTe,*) = P_value(iTe,*) * expN
     pexpN2(iTe,*) = P_value(iTe,*) * expN2
  endfor

  for iNe=0,NNe-1 do begin
     pexpT (*,iNe) = P_value(*,iNe) * expT
     pexpT2(*,iNe) = P_value(*,iNe) * expT2
  endfor

  pexpTN = p_value * expTN

  
  grad(*,*,0) = (1./sigNe)*(pexpN  - q*pexpT)                                       ; \propto dP/dNm
  grad(*,*,1) = (1./sigTe)*(pexpT  - q*pexpN)                                       ; \propto dP/dTm
  grad(*,*,2) = (1./sigTe)*(pexpT2 - q*pexpTN - (1.-q^2)*p_value)                   ; \propto dP/dsigT  
  grad(*,*,3) = (1./sigNe)*(pexpN2 - q*pexpTN - (1.-q^2)*p_value)                   ; \propto dP/dsigN
  grad(*,*,4) = (q*p_value + pexpTN - (q/(1.-q^2))*(pexpN2 + pexpT2 - 2.*q*pexpTN)) ; \propto dP/dq

  return, (1./(1-q^2))*grad
end

;---------------------------------------------------------------------
;
; Brief description:
;
; For values of Te0 and Ne0, this function returns the value of the
; gradient of the bivariate normal joint temperature-density
; probability distribution respect to: Tem, Nem, SigTe, SigNe, q.
;
; This routine is necesary if the double integrals are calculated 
; using CS (e_function_cs). 
;
; INPUTS:
;  Te0: 1D array of NTe elements with the electron temperature in units of [K]
;  Ne0: 1D array of NNe elements with the electron density     in units of [cm-3]

;  (in common parameters)
;  Tem: float with the mean electron temperature  in units of [K]
;  Nem: float with the mean electron density      in units of [cm-3]
;SigTe: float with the electron temperature StDev in units of [K]
;SigNe: float with the electron density     StDev in units of [cm-3]
;    q: temperature-density correlation dimensionless coefficient
;       
; OUTPUT:
;       value of the gradient of the probabilility function p in a 
;       3D array of NTe X NNe x 5, evaluated in the 1D arrays: Ne0 and Te0.
;
; History:  V1.0, Federico Nuevo, 01-05-20.
;           V1.1, A.M. Vasquez,   06-05-20.
;                 Removed use of loops.
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

  if abs(q) lt 1 then begin

  expT     = double(((Te0-Tem)/SigTe))
  expN     = double(((Ne0-Nem)/SigNe))
  expT2    = double(((Te0-Tem)/SigTe)^2)
  expN2    = double(((Ne0-Nem)/SigNe)^2)
  expTN    = double(((Te0-Tem)/sigTe)#((Ne0-Nem)/sigNe))
  p_value  = (1./(2.*!pi*sigTe*sigNe*sqrt(1.-q^2)))*$
             exp( - (1./2./(1.-q^2))* expT2 ) # exp( - (1./2./(1.-q^2))* expN2 ) * $
             exp((1./(1.-q^2))*q*expTN )

  ; Correct NaNs due to INFINITY exp(expTN)
  index_nan = where(finite(exp((1./(1.-q^2))*q*expTN)) ne 1)
  if index_nan(0) ne -1 then p_value(index_nan) = 0.


;-----------------------------------------------------------------------
; Fede, aqui te dejo el codigo original tuyo y mi versión SIN LOOPS:
; una vez entendido y verificado, borra la parte de loops.
; si hay otros códigos donde pudieras hacer lo mismo, do it.
; Notá que SIN LOOPS no precisás ni siquiera definir los 4 arrays 2D,
; que se crean en la operación matemática misma. 
; Para confirmar, vos  podés hacer una corrida de NNe=7, NTe=5
; (siempre usa dos N diferentes), activa el codigo viejo y pone un
; STOP justo antes del nuevo. Imprimi PexpN, PexpT, PexpN2 y PexpT2
; Y luego impirmi las cuatro expresiones que yo codeé. verás que son
; IDENTICAS.
goto,skip_code_with_loops
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
skip_code_with_loops:
;-----------------------------------------------------------------------
     pexpN  = P_value*((1.+dblarr(NTe))#expN )
     pexpN2 = P_value*((1.+dblarr(NTe))#expN2)
     pexpT  = P_value*(expT #(1.+dblarr(NNe)))
     pexpT2 = P_value*(expT2#(1.+dblarr(NNe)))

  pexpTN = p_value * expTN
  
  grad(*,*,0) = (1./sigNe)*(pexpN  - q*pexpT)                                       ; \propto dP/dNm
  grad(*,*,1) = (1./sigTe)*(pexpT  - q*pexpN)                                       ; \propto dP/dTm
  grad(*,*,2) = (1./sigTe)*(pexpT2 - q*pexpTN - (1.-q^2)*p_value)                   ; \propto dP/dsigT  
  grad(*,*,3) = (1./sigNe)*(pexpN2 - q*pexpTN - (1.-q^2)*p_value)                   ; \propto dP/dsigN
  grad(*,*,4) = (q*p_value + pexpTN - (q/(1.-q^2))*(pexpN2 + pexpT2 - 2.*q*pexpTN)) ; \propto dP/dq

  grad =  (1./(1-q^2))*grad
endif
  

  return, grad
  
end

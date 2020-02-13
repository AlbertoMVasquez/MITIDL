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
;           V1.1, adapto la rutina para que Ne0 y Te0 puedan ser arrays
;---------------------------------------------------------------------

function p_function_fede, Ne0, Te0
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

  p_value=dblarr(NTe,NNe)
  for iNe=0,NNe-1 do begin
     for iTe=0,NTe-1 do begin
        expT    = ((Te0(iTe)-Tem)/SigTe)^2
        expN    = ((Ne0(iNe)-Nem)/SigNe)^2
        expTN   = (Te0(iTe)-Tem)*(Ne0(iNe)-Nem)/(SigTe*SigNe)
        p_value(iTe,iNe) = (1./(2.*!pi*sigTe*sigNe*sqrt(1.-q)))*$
                           exp( - (1./2./(1.-q^2))*( expT + expN - 2.*q*expTN ) )  
     ENDFOR
  ENDFOR
  return,p_value
end

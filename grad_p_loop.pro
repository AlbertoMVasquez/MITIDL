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

function grad_p_loop, Ne0, Te0
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
  
  for iNe=0,NNe-1 do begin
     for iTe=0,NTe-1 do begin


        expT    = ((Te0(iTe)-Tem)/SigTe)
        expN    = ((Ne0(iNe)-Nem)/SigNe)
        expT2   = expT^2
        expN2   = expN^2
        expTN   = expT * expN
        p_value = (1./(2.*!pi*sigTe*sigNe*sqrt(1.-q^2)))*$
                  exp( - (1./2./(1.-q^2))*( expT2 + expN2 - 2.*q*expTN ) )  

        grad(iTe,iNe,0) = p_value * ( 1./(1.-q^2)/sigNe * ( expN  + q * expT ) ) ; dP/dNm
        
        grad(iTe,iNe,1) = p_value * ( 1./(1.-q^2)/sigTe * ( expT  + q * expN ) ) ; dP/dTm
        
        grad(iTe,iNe,2) = p_value * ( 1./(1.-q^2)/sigTe * ( expT2 + q * expTN - ( 1.-q^2))) ; dP/sigT
        
        grad(iTe,iNe,3) = p_value * ( 1./(1.-q^2)/sigNe * ( expN2 + q * expTN - ( 1.-q^2))) ; dP/sigN
        
        grad(iTe,iNe,4) = p_value * ( q/(1.-q^2) -q/(1.-q^2)^2 * ( expT2 + expN2 - 2.*q*expTN )-1./(1.-q^2)*expTN ) ; dP/dq

     endfor
   endfor

        return,grad
     end

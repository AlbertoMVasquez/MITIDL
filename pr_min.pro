pro pr_min,guess,out,phiv
  P=guess
  ftol=1.d-4                    ; convergence tolerance on the cost function value 
  retry:
  frprmn,P,ftol,iter,fret       ; entra en frprmn (Polak-Ribiere)
  if fret eq -666. then begin
     ftol=ftol*5.
     if ftol ge 1.d-1 then begin
        P   = 0.*P+1.
        fret= +999.
        goto,fin
     endif
     goto,retry
  endif
  fin:
  out=P                         ; out: parameter array that minimized PHI
  PHIv=fret                     ; the value PHI(out)
  return
end

;---------------------------------------------------------------
;  Polak-Ribiere minimization (Numerical Recipes, p. 416):
;---------------------------------------------------------------
pro frprmn,p,ftol,iter,fret
  ITMAX=1000
  EPS=1.d-10
; This thing in NumRecipes: "fa=func(a,df=xi)" changes to these next 2 lines:
  fp =     PHI(p)
  xi = gradPHI(p)
  g  =-xi
  h  =  g
  xi =  h
  for its=1,ITMAX do begin
     iter=its
     linmin,p,xi,fret
;NEW CODE NEXT LINE
     if fret eq -666. then return
     if 2.*abs(fret-fp) le ftol*(abs(fret)+abs(fp)+EPS) then return
     fp=fret                    ; = PHI(p)
     xi=gradPHI(p)
     gg=total(g^2)
     dgg=total((xi+g)*xi)
     if gg eq 0.0 then return   ; Unlikely, but if Gradient is exactly zero then we are done
     gam=dgg/gg
     g=-xi
     h=g+gam*h
     xi=h
  endfor
  print,'FRPRMN maximum iterations exceeded'
  return
end

;---------------------------------------------------------------
;  Minimization along a line (p. 412):
;---------------------------------------------------------------
pro linmin,p,xi,fret
  common f1com,pcom,xicom
  common iteracionesdbrent,iter
;
  TOL=1.d-4;1.d-7                     ; fractional precision of 1D minimization = Sqrt(machine-precision)
;
  pcom=p
  xicom=xi
  ax=0.d0
  xx=1.d0
  mnbrak,ax,xx,bx,fa,fx,fb
  fret=dbrent(ax,xx,bx,TOL,xmin)

;NEW CODE NEXT LINE
  if fret eq -666. then return

  xi=xmin*xi
  p=p+xi
  return
end

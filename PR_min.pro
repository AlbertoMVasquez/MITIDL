

pro PR_min,guess,out,phiv
  P=guess
  ftol=1.d-4              ; convergence tolerance on the cost function value 
  retry:
  frprmn,P,ftol,iter,fret ; entra en frprmn (Polak-Ribiere)
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
  out=P      ; out: parameter array that minimized PHI
  PHIv=fret  ; the value PHI(out)
  return
end

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
  fp=fret ; = PHI(p)
  xi=gradPHI(p)
  gg=total(g^2)
  dgg=total((xi+g)*xi)
  if gg eq 0.0 then return ; Unlikely, but if Gradient is exactly zero then we are done
  gam=dgg/gg
  g=-xi
  h=g+gam*h
  xi=h
endfor
print,'FRPRMN maximum iterations exceeded'
return
end

;---------------------------------------------------------------
;  Minimization along a line (p. 412):
;---------------------------------------------------------------
pro linmin,p,xi,fret
common f1com,pcom,xicom
common iteracionesdbrent,iter
;
TOL=1.d-7; fractional precision of 1D minimization = Sqrt(machine-precision)
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
;---------------------------------------------------------------
;  Bracketing a minimum along a line (Numerical Recipes, p. 393):
;---------------------------------------------------------------
pro mnbrak,ax,bx,cx,fa,fb,fc
GOLD=1.618034
GLIMIT=100.
TINY=1.d-20
;
fa=f1dim(ax)
fb=f1dim(bx)
if fb gt fa then begin     ; Switch roles so we can go downhill from a to b.
  dum=ax
  ax=bx
  bx=dum
  dum=fb
  fb=fa
  fa=dum
endif
cx=bx+GOLD*(bx-ax)         ; First guess of c.
fc=f1dim(cx)
while (fb ge fc) do begin
  r=(bx-ax)*(fb-fc)
  q=(bx-cx)*(fb-fa)
  dqr=max([abs(q-r),TINY])
  if q-r lt 0. then dqr=-dqr
  u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*dqr)
  ulim=bx+GLIMIT*(cx-bx)
  if (bx-u)*(u-cx) gt 0. then begin ; Parabolic u is between b and c. Try it.
    fu=f1dim(u)
    if fu lt fc then begin          ; Got a minimum between b and c.
      ax=bx
      fa=fb
      bx=u
      fb=fu
      return
    endif
    if fu gt fb then begin          ; Got a minimum between a and u.
      cx=u
      fc=fu
      return
    endif
    u=cx+GOLD*(cx-bx)   ; Parabolic fit was no use. Use default magnification.
    fu=f1dim(u)
  endif else begin
    if (cx-u)*(u-ulim) gt 0. then begin  ; Parabolic fit between c and limit.
      fu=f1dim(u)
      if fu lt fc then begin
        bx=cx
        cx=u
        u=cx+GOLD*(cx-bx)
        fb=fc
        fc=fu
        fu=f1dim(u)
      endif
    endif else begin
      if (u-ulim)*(ulim-cx) ge 0. then begin  ; Limit parabolic u to maximum.
        u=ulim
        fu=f1dim(u)
      endif else begin       ; Reject parabolic u, use default magnification.
        u=cx+GOLD*(cx-bx)
        fu=f1dim(u)
      endelse
    endelse
  endelse
  ax=bx                      ; Eliminate oldest point and continue.
  bx=cx
  cx=u
  fa=fb
  fb=fc
  fc=fu
endwhile
return
end
;---------------------------------------------------------------
;  Search minimum along a line using derivatives (p. 400):
;---------------------------------------------------------------
function dbrent,ax,bx,cx,tol,xmin
common iteracionesdbrent,iter
;
ITMAX=1000
ZEPS=1.d-10
;
a=min([ax,cx])
b=max([ax,cx])
v=bx
w=v
x=v
e=0.d0
dx=1
fx=f1dim(x)
dx=df1dim(x)
fv=fx  &  dv=dx
fw=fx  &  dw=dx
for iter=1,ITMAX do begin
    xm=0.5d0*(a+b)
    tol1=tol*abs(x)+ZEPS
    tol2=2.d0*tol1
;   print,'dbrent iter#:',iter,'..... convergence:',abs(x-xm)/(tol2-0.5d0*(b-a))                 
    if abs(x-xm) le (tol2-0.5d0*(b-a)) then begin
      xmin=x
      return,fx
    endif
    if abs(e) gt tol1 then begin
      d1=2.d0*(b-a)
      d2=d1
      if dw ne dx then d1=(w-x)*dx/(dx-dw)
      if dv ne dx then d2=(v-x)*dx/(dx-dv)
      u1=x+d1
      u2=x+d2
      ok1=((a-u1)*(u1-b) gt 0.) and (dx*d1 le 0.)
      ok2=((a-u2)*(u2-b) gt 0.) and (dx*d2 le 0.)
      olde=e
      e=d
      if not (ok1 or ok2) then goto,L1
      if (ok1 and ok2) then begin
        if abs(d1) lt abs(d2) then d=d1 else d=d2
      endif else begin
        if ok1 then d=d1 else d=d2
      endelse
      if abs(d) gt abs(0.5d0*olde) then goto,L1
      u=x+d
      if ((u-a) lt tol2) or ((b-u) lt tol2) then begin
        if (xm-x) ge 0. then d=tol1 else d=-tol1
      endif
      goto,L2
    endif
L1: if dx ge 0. then e=a-x else e=b-x
    d=0.5d0*e
L2: if abs(d) ge tol1 then begin
      u=x+d
      fu=f1dim(u)
    endif else begin
      if d ge 0. then u=x+tol1 else u=x-tol1
      fu=f1dim(u)
      if fu gt fx then begin
        xmin=x
        return,fx
      endif
    endelse
    du=df1dim(u) ;!!!
    if fu le fx then begin
      if u ge x then a=x else b=x
      v=w
      fv=fw
      dv=dw
      w=x
      fw=fx
      dw=dx
      x=u
      fx=fu
      dx=du
    endif else begin
      if u lt x then a=u else b=u
      if (fu le fw) or (w eq x) then begin
        v=w
        fv=fw
        dv=dw
        w=u
        fw=fu
        dw=du
      endif else begin
        if (fu le fv) or (v eq x) or (v eq w) then begin
          v=u
          fv=fu
          dv=du
        endif
      endelse
    endelse
endfor

print,'Error: dbrent exceeded maximum iterations. Retry.'

;ORIGINAL CODE NEXT TWO LINES
;xmin=x
;return,fx

;NEW CODE NEXT TWO LINES
 xmin= -666.
return,-666.

end
;---------------------------------------------------------------
;  Search minimum along a line using derivatives (p. 400):
;---------------------------------------------------------------
function dbrent_new,ax,bx,cx,tol,xmin
common iteracionesdbrent,iter
;
ITMAX=100
ZEPS=1.d-10
;
a=min([ax,cx])
b=max([ax,cx])
v=bx
w=v
x=v
e=0.d0
dx=1
fx=f1dim(x)
dx=df1dim(x)
fv=fx  &  dv=dx
fw=fx  &  dw=dx
for iter=1,ITMAX do begin
    xm=0.5d0*(a+b)
    tol1=tol*abs(x)+ZEPS
    tol2=2.d0*tol1
;   print,'dbrent iter#:',iter,'..... convergence:',abs(x-xm)/(tol2-0.5d0*(b-a))                 
    if abs(x-xm) le (tol2-0.5d0*(b-a)) then begin
      xmin=x
      return,fx
    endif
    if abs(e) gt tol1 then begin
      d1=2.d0*(b-a)
      d2=d1
      if dw ne dx then d1=(w-x)*dx/(dx-dw)
      if dv ne dx then d2=(v-x)*dx/(dx-dv)
      u1=x+d1
      u2=x+d2
      ok1=((a-u1)*(u1-b) gt 0.) and (dx*d1 le 0.)
      ok2=((a-u2)*(u2-b) gt 0.) and (dx*d2 le 0.)
      olde=e
      e=d
      if not (ok1 or ok2) then goto,L1
      if (ok1 and ok2) then begin
        if abs(d1) lt abs(d2) then d=d1 else d=d2
      endif else begin
        if ok1 then d=d1 else d=d2
      endelse
      if abs(d) gt abs(0.5d0*olde) then goto,L1
      u=x+d
      if ((u-a) lt tol2) or ((b-u) lt tol2) then begin
        if (xm-x) ge 0. then d=tol1 else d=-tol1
      endif
      goto,L2
    endif
L1: if dx ge 0. then e=a-x else e=b-x
    d=0.5d0*e
L2: if abs(d) ge tol1 then begin
      u=x+d
      fu= f1dim(u)
      du=df1dim(u)
    endif else begin
      if d ge 0. then u=x+tol1 else u=x-tol1
      fu= f1dim(u)
      du=df1dim(u)
      if fu gt fx then begin
        xmin=x
        return,fx
      endif
    endelse
    if fu le fx then begin
      if u ge x then a=x else b=x
      v=w
      fv=fw
      dv=dw
      w=x
      fw=fx
      dw=dx
      x=u
      fx=fu
      dx=du
    endif else begin
      if u lt x then a=u else b=u
      if (fu le fw) or (w eq x) then begin
        v=w
        fv=fw
        dv=dw
        w=u
        fw=fu
        dw=du
      endif else begin
        if (fu le fv) or (v eq x) or (v eq w) then begin
          v=u
          fv=fu
          dv=du
        endif
      endelse
    endelse
endfor
print,'Error: dbrent exceeded maximum iterations'
xmin=x
return,fx
end
;---------------------------------------------------------------
;  Evaluate function along a line (p. 413):
;---------------------------------------------------------------
function f1dim,x
common f1com,pcom,xicom
return,PHI(pcom+x*xicom)
end
;---------------------------------------------------------------
;  Evaluate derivative along a line (p. 417):
;---------------------------------------------------------------
function df1dim,x
common f1com,pcom,xicom
return,total( (gradPHI(pcom+x*xicom)) * xicom)
end

;-----------------End of the code-----------------------------------
;---------The rest is not used anymore------------------------------





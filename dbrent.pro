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
     du=df1dim(u)               ;!!!
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

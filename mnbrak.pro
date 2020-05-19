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

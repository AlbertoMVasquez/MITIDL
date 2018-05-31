;---------------------------------------------------------------------
;
; Brief description:
;
; This function returns the trilinearly interpolated value of
; F(xa,ya,za) into the point (x0,y0,z0)
;
; It also returns its 3 partial derivatives at that point.
;
; INPUTS:
; F:          3-dimensional array.
; xa, ya, za: three 1-dimensional vectors with the grid values.
; x0, y0, z0: three floats indicating the point where to interpolate.
;
; OUTPUT:
;
; RESULT = [F, dF/dx, dF/dy, dF/dz]
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

function findval3d_function, F, xa, ya, za, x0, y0, z0
  ; Find the x-planes that surround x0
  ixA = max(where(xa le x0))
  if ixA eq -1 then begin
     print,'x0 value is out of range' ; we may need something better to treat this case.
     stop
  endif
  ixB=ixA+1
  ; Define two 2D y-z arrays at x-planes that surround x0. 
  F_xA = reform(F(ixA,*,*))
  F_xB = reform(F(ixB,*,*))
  ; Bi-linearly interpolate in xA and xB planes.   
  D_xA = findval2D_function( F_xA ,ya ,za , y0, z0 )
  D_xB = findval2D_function( F_xB ,ya ,za , y0, z0 )
  ; Linearly interpolate the value og G along x, betwee xA and xB planes.
  F_value = D_xA + (D_xB-D_xA)*(x0-xa(ixA))/(xa(ixB)-xa(ixA))

  dF_dx = 0.
  dF_dy = 0.
  dF_dz = 0.

  return, [F_value, dF_dx, dF_dy, dF_dz]
end

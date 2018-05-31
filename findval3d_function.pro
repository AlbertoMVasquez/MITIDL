;---------------------------------------------------------------------
;
; Brief description:
;
; This function returns the trilinearly interpolated value of
; DATA_ARRAY(xa,ya,za) into the point (x0,y0,z0)
;
; INPUTS:
; DATA_ARRAY: 3-dimensional array.
; xa, ya, za: three 1-dimensional vectors with the grid values.
; x0, y0, z0: three floats indicating the point where to interpolate to.
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

function findval3d_function, DATA_ARRAY, xa, ya, za, x0, y0, z0
  ; Find the x-planes that surround x0
  ixA = max(where(xa le x0))
  if ixA eq -1 then begin
     print,'x0 value is out of range' ; we may need something better to treat this case.
     stop
  endif
  ixB=ixA+1
  ; Define two 2D y-z arrays at x-planes that surround x0. 
  DATA_ARRAY_xA = reform(DATA_ARRAY(ixA,*,*))
  DATA_ARRAY_xB = reform(DATA_ARRAY(ixB,*,*))
  ; Bi-linearly interpolate in xA and xB planes.   
  D_xA = findval2D_function( DATA_ARRAY_xA ,ya ,za , y0, z0 )
  D_xB = findval2D_function( DATA_ARRAY_xB ,ya ,za , y0, z0 )
  ; Linearly interpolate the value og G along x, betwee xA and xB planes.
  D_value = D_xA + (D_xB-D_xA)*(x0-xa(ixA))/(xa(ixB)-xa(ixA))
  return,D_value
end

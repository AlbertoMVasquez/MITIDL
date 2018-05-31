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

; Find the x-indices that surround x0
  ixA = max(where(xa le x0))
  if ixA eq -1 then begin
     print,'x0 value is out of range' ; we may need something better to treat this case.
     stop
  endif
  ixB=ixA+1
; Find the y-indices that surround y0
  iyA = max(where(ya le y0))
  if iyA eq -1 then begin
     print,'y0 value is out of range' ; we may need something better to treat this case.
     stop
  endif
  iyB=iyA+1
; Find the z-indices that surround z0
  izA = max(where(za le z0))
  if iyA eq -1 then begin
     print,'z0 value is out of range' ; we may need something better to treat this case.
     stop
  endif
  izB=izA+1

; Define two 2D y-z arrays at x-planes that surround x0:
  F_xA = reform(F(ixA,*,*))
  F_xB = reform(F(ixB,*,*))
; Bi-linearly interpolate in xA and xB planes:
  D_xA = findval2D_function( F_xA ,ya ,za , y0, z0 )
  D_xB = findval2D_function( F_xB ,ya ,za , y0, z0 )

; Define two 2D x-z arrays at y-planes that surround y0:
  F_yA = reform(F(*,iyA,*))
  F_yB = reform(F(*,iyB,*))
; Bi-linearly interpolate in yA and yB planes:
  D_yA = findval2D_function( F_yA ,xa ,za , x0, z0 )
  D_yB = findval2D_function( F_yB ,xa ,za , x0, z0 )

; Define two 2D x-y arrays at z-planes that surround z0:
  F_zA = reform(F(*,*,izA))
  F_zB = reform(F(*,*,izB))
; Bi-linearly interpolate in zA and zB planes:
  D_zA = findval2D_function( F_zA ,xa ,ya , x0, y0 )
  D_zB = findval2D_function( F_zB ,xa ,ya , x0, y0 )

; Linearly interpolate the value of G along x, between xA and xB planes:
  F_value_x = D_xA + (D_xB-D_xA)*(x0-xa(ixA))/(xa(ixB)-xa(ixA))

; Linearly interpolate the value of G along y, between yA and yB planes:
  F_value_y = D_yA + (D_yB-D_yA)*(y0-ya(iyA))/(ya(iyB)-ya(iyA))

; Linearly interpolate the value of G along z, between zA and zB planes:
  F_value_z = D_zA + (D_zB-D_zA)*(z0-za(izA))/(za(izB)-za(izA))

  eps=1.e-8
  if abs(1.-F_value_y/F_value_x) gt eps or abs(1.-F_value_z/F_value_x) gt eps then begin
     print,'3-linear interpolation failed.'
     stop
  endif

; Assign one value to F
  F_value = F_value_x

; Compute derivatives:
  dF_dx = (D_xB - D_xA)/(xa(ixB)-xa(ixA))
  dF_dy = (D_yB - D_yA)/(ya(iyB)-ya(iyA))
  dF_dz = (D_zB - D_zA)/(za(izB)-za(izA))

  return, [F_value, dF_dx, dF_dy, dF_dz]
end

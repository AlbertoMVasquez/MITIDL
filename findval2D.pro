;---------------------------------------------------------------------
;
; Brief description:
;
; This function returns the bilinearly interpolated value of
; DATA_ARRAY(ya,za) into the point (y0,z0)
;
; INPUTS:
; DATA_ARRAY: 2-dimensional array.
; ya, za: two 1-dimensional vectors with the grid values.
; y0, z0: two floats indicating the point where to interpolate to.
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

function findval2D, DATA_ARRAY ,ya ,za , y0, z0
iyA=max(where(ya le y0))
izA=max(where(za le z0))
Df=0.
if iyA eq -1 or izA eq -1 then begin
   print,'Te and/or r out of range.'
   stop
endif
iyB=iyA+1
izB=izA+1
if iyA eq n_elements(ya)-1 then iyB=iyA
if izA eq n_elements(za)-1 then izB=izA
D1=DATA_ARRAY(iyA,izA) 
D2=DATA_ARRAY(iyB,izA) 
D4=DATA_ARRAY(iyA,izB) 
D5=DATA_ARRAY(iyB,izB)
if iyA lt iyB AND izA lt izB then begin
 D3=D1+(D2-D1)*(y0-yA(iyA))/(ya(iyB)-ya(iyA))
 D6=D4+(D5-D4)*(y0-yA(iyA))/(ya(iyB)-ya(iyA))
 Df=D3+(D6-D3)*(z0-zA(izA))/(za(izB)-za(izA))
endif
if iyA lt iyB AND izA eq izB then Df=D1+(D2-D1)*(y0-ya(iyA))/(ya(iyB)-ya(iyA))
if iyA eq iyB AND izA lt izB then Df=D1+(D4-D1)*(z0-za(izA))/(za(izB)-za(izA))
if iyA eq iyB AND izA eq izB then Df=D1
return,Df
end

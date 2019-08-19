;---------------------------------------------------------------------
;
; Brief description:
;
; This code reads the tomographic solution "x" for all instruments.
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

pro xread,dir=dir,file=file,nr=nr,nt=nt,np=np,map=map
if not keyword_set(dir) then dir = '/data1/tomography/bindata/'
map=fltarr(nr,nt,np)
openr,1,dir+file
readu,1,map
close,1
return
end

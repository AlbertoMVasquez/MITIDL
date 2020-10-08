

  pro plot_scatter_and_historatio_eps,x,y,$
                                  dir=dir,filename=filename,xsuffix=xsuffix,$
                                  ysuffix=ysuffix,titulo=titulo,xr=xr,yr=yr,$
                                  historange=historange

    ;Coyote library path
    !PATH = Expand_Path('+~/idlfiles/coyote/') + ':' + !PATH
    ;Create custom made symbol (psym=8) for scatter plots
    N=25
    A = FINDGEN(N) * (!PI*2/float(N-1))
    f=2.
    USERSYM, COS(A)/f, SIN(A)/f,/FILL

    if not keyword_set(dir) then dir='~/Downloads/'
    if not keyword_set(filename) then filename='scatter'
    if not keyword_set(xsuffix) then xsuffix='x'
    if not keyword_set(ysuffix) then ysuffix='y'
    if not keyword_set(titulo) then titulo='Scatter Plot'
    if not keyword_set(xr) then xr=[min(x),max(x)]
    if not keyword_set(yr) then yr=[min(y),max(y)]
    ps1,dir+filename+'.eps',0
    DEVICE,/INCHES,XSIZE=5.7,SCALE_FACTOR=1
    DEVICE,/INCHES,YSIZE=8.5,SCALE_FACTOR=1
    !P.CHARTHICK=6
    !p.charsize=2
    !p.multi=[0,1,2]


    ;================Scatter Plot=====================================================
    plot,x,y,psym=3,xr=xr,yr=yr,xstyle=1,ystyle=1,xtitle=xsuffix,ytitle=ysuffix,$
         title=titulo,color=0,background=255,xthick=6,ythick=6
    loadct,12
    oplot,[0.,100],[0,100],linestyle=0,color=100,th=4
    xyouts,[.7],[.64],'!4q!3='+[strmid(string(correlate(x,y)),6,4)],/normal,col=[100],charthick=[6]
    loadct,0

    ;================Histogram y/x=====================================================   
    
    f  =y/x
    i=where(finite(f) ne 0)
    f = f(i)
    med=median(f)
    mea=mean(f)
    std=stddev(f)
    valmin= 0.5
    valmax= 1.5
    if keyword_set(historange) then begin
       valmin=historange(0)
       valmax=historange(1)
    endif
    i= where (f gt valmax)
    if i(0) ne -1 then f(i)=valmax
    i= where (f lt valmin)
    if i(0) ne -1 then f(i)=valmin
    binsize=(valmax-valmin)/1.e2*3.

    cgHistoplot, f, binsize=binsize, xtitle='Ratio '+ysuffix+'/'+xsuffix, title='',$
                 xrange=[valmin-binsize,valmax+binsize],/frequency,/fill,xstyle=1
    
    loadct,12
    oplot,[1,1],[0,10],th=6,color=100
   ;nstat=3
    nstat=2
    i0med=6
    i0mea=6
    i0std=6
    if med lt 1. then i0med=5
    if mea lt 1. then i0mea=5
    if std lt 1. then i0std=5
   ;xyouts,0.5+fltarr(nstat),.4-.05*findgen(nstat),['Median =','  Mean =','  Stdev =']+$
   ;       [strmid(string(med),i0med,5),strmid(string(mea),i0mea,5),strmid(string(std),i0std,5)],$
   ;       charthick=6,color=100+fltarr(nstat),/normal
   xyouts,0.5+fltarr(nstat),.4-.05*findgen(nstat),['Mean =','Stdev =']+$
          [strmid(string(mea),i0mea,5),strmid(string(std),i0std,5)],$
          charthick=6,color=100+fltarr(nstat),/normal
   loadct,0
   ps2

   return
end
  


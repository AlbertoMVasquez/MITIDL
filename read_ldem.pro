
pro read_ldem,file,ldem=ldem,dem=dem,gauss1=gauss1,$
              doublegauss=doublegauss,newpar4=newpar4,newpar5=newpar5,$
              gauss2c=gauss2c,glc=glc,gauss3c=gauss3c,gauss7par=gauss7par,$
              equalizer=equalizer,dgfw=dgfw,g2cfw=g2cfw,dir=dir
  common results_images,ima,sima,ra,pa,ya,za,em,npx
  common comunes,tm,wt,nband,demc,PHI,parametrizacion,Tmin,Tmax,nr,nth,np,rad,lat,lon,lambda,WTc
  common results_tomo,tfbe,sfbe,N_e
  common demc_extra,demc2,demc3
  common loss_rate,Er
  common fixed_width,sigma
  common fixed_parameter_equalizer,xmax,sigma_v,demv
  sca  =0
  Tmin =0.
  Tmax =0.
  nr   =0
  nth  =0
  np   =0
  npx  =0
  nband=0L
  ndat =0L
  nx   =0
  ny   =0

  if keyword_set(ldem) then resultado = 'LDEM'
  if keyword_set( dem) then resultado =  'DEM'

  if keyword_set(gauss1) then begin
     parametrizacion='gauss1'
     npar=3
  endif

if keyword_set(doublegauss) then begin
   parametrizacion='gauss2'
   npar=6
endif

if keyword_set(dgfw)         then begin
   parametrizacion='gauss2FW'
   npar=5
endif

if keyword_set(g2cfw)       then begin
   parametrizacion='gauss2cFW'
   npar=4
endif

if keyword_set(gauss2c) then begin
   parametrizacion = 'gauss2c'
   npar=5
endif

if keyword_set(equalizer) then begin
  parametrizacion = 'equalizer'
  npar = 0
endif

if keyword_set(glc) then begin
   parametrizacion = 'glc'
   npar=5
endif

if keyword_set(gauss3c) then begin
   parametrizacion = 'gauss3c'
   npar=6
endif

if keyword_set(gauss7par) then begin
   parametrizacion = 'gauss7par'
   npar=7
endif

if keyword_set(newpar4) then begin
   parametrizacion='newpar4'
   npar=4
endif

if keyword_set(newpar5) then begin
   parametrizacion='newpar5'
   npar=5
endif


 if not keyword_set(dir) then  dir ='/data1/DATA/ldem_files/'
  
 
  openr,1,dir+file
  readu,1,resultado
  readu,1,parametrizacion
  readu,1,sca,Tmin,Tmax,nr,nth,np,npx,nband,ndat
  if keyword_set(dem) then begin
     readu,1,nx,ny
     if nx eq npx and ny eq npx then print,'this is a full image DEM file'
  endif
  if parametrizacion eq 'equalizer' then begin
     xmax=dblarr(nband)
     sigma_v=dblarr(nband)
     readu,1,xmax,sigma_v
  endif

  
  if keyword_set(ldem) then DEMT_A = dblarr(nr ,nth,np,ndat)
  if keyword_set( dem) then DEMT_A = dblarr(nx ,ny    ,ndat)
    
  readu,1,DEMT_A
  close,1

 if keyword_set(ldem) then begin
    TFBE = DEMT_A(*,*,*,0    :  nband-1)
    SFBE = DEMT_A(*,*,*,nband:2*nband-1)
    i=0

    if parametrizacion eq 'gauss1' OR parametrizacion ne 'lorentz' then begin
         DEMc = DEMT_A(*,*,*,    2*nband  )
         Tc   = DEMT_A(*,*,*,    2*nband+1)
         i    = 2
    endif
    if parametrizacion eq 'gauss2' then begin
         DEMc  = DEMT_A(*,*,*,   2*nband  )
         DEMc2 = DEMT_A(*,*,*,   2*nband+1)
         i     = 2 
    endif
    if parametrizacion eq 'gauss2FW' or parametrizacion eq 'gauss2cFW' then begin
         DEMc  = DEMT_A(*,*,*,   2*nband  )
         DEMc2 = DEMT_A(*,*,*,   2*nband+1)
         sigma = DEMT_A(7,0,0,   2*nband+2)
         i     = 3 
    endif
    if parametrizacion eq 'gauss2c' then begin
         DEMc  = DEMT_A(*,*,*,   2*nband  )
         DEMc2 = DEMT_A(*,*,*,   2*nband+1)
         i     = 2
    endif
    if parametrizacion eq 'glc' then begin
         DEMc  = DEMT_A(*,*,*,   2*nband  )
         DEMc2 = DEMT_A(*,*,*,   2*nband+1)
         i     = 2
      endif
    if parametrizacion eq 'gauss3c' then begin
         DEMc  = DEMT_A(*,*,*,   2*nband  )
         DEMc2 = DEMT_A(*,*,*,   2*nband+1)
         DEMc3 = DEMT_A(*,*,*,   2*nband+2)
         i     = 3
    endif
    if parametrizacion eq 'gauss7par' then begin
         DEMc  = DEMT_A(*,*,*,   2*nband  )
         DEMc2 = DEMT_A(*,*,*,   2*nband+1)
         DEMc3 = DEMT_A(*,*,*,   2*nband+2)
         i     = 3
    endif
    if  parametrizacion eq 'equalizer' then begin
         npar = nband
         DEMv = DEMT_A(*,*,*,    2*nband:2*nband+npar-1  )
         DEMc = reform(DEMv(*,*,*,0))
        ;DEMc = total(DEMv,4)
         i    = npar
    endif



  lambda = DEMT_A(*,*,*,2*nband+i:2*nband+i+npar-1)
     PHI = DEMT_A(*,*,*,          2*nband+i+npar  )
    ftol = DEMT_A(*,*,*,          2*nband+i+npar+1)
      i0 = 2 ; set to ZERO if PHI and ftol above are commented out.
     N_e = DEMT_A(*,*,*,          2*nband+i+npar+i0  )
      Tm = DEMT_A(*,*,*,          2*nband+i+npar+i0+1)
      WT = DEMT_A(*,*,*,          2*nband+i+npar+i0+2)
      Er = DEMT_A(*,*,*,          2*nband+i+npar+i0+3)
   endif

 if keyword_set( dem) then begin
    if parametrizacion eq 'gauss1' then begin
       IMA = DEMT_A(*,*,0    :  nband-1)
       SIMA = DEMT_A(*,*,nband:2*nband-1)
       DEMc = DEMT_A(*,*,      2*nband  )
       Tc = DEMT_A(*,*,      2*nband+1)
       i=0
       lambda = DEMT_A(*,*,2*nband+2+i:2*nband+2+i+npar-1)
       PHI = DEMT_A(*,*,            2*nband+2+i+npar  )
       ftol = DEMT_A(*,*,            2*nband+2+i+npar+1)
       i0 = 2                   ; set to ZERO if PHI and ftol above are commented out.
       EM = DEMT_A(*,*,            2*nband+2+i+npar  +i0)
       Tm = DEMT_A(*,*,            2*nband+2+i+npar+1+i0)
       WT = DEMT_A(*,*,            2*nband+2+i+npar+2+i0)
       Er = DEMT_A(*,*,            2*nband+2+i+npar+3+i0)
       ra = DEMT_A(*,*,            2*nband+2+i+npar+4+i0)
       pa = DEMT_A(*,*,            2*nband+2+i+npar+5+i0)
       ya = DEMT_A(*,*,            2*nband+2+i+npar+6+i0)
        za = DEMT_A(*,*,            2*nband+2+i+npar+7+i0)
     endif
    if parametrizacion eq 'gauss2FW' then begin
       IMA = DEMT_A(*,*,0    :  nband-1)
       SIMA = DEMT_A(*,*,nband:2*nband-1)
       DEMc  = DEMT_A(*,*,      2*nband  )
       DEMc2 = DEMT_A(*,*,      2*nband+1)
       sigma = DEMT_A(0,0,      2*nband+2)
       i=1
       lambda = DEMT_A(*,*,2*nband+2+i:2*nband+2+i+npar-1)
       PHI = DEMT_A(*,*,            2*nband+2+i+npar  )
       i0=2
       EM = DEMT_A(*,*,            2*nband+2+i+npar  +i0)
       Tm = DEMT_A(*,*,            2*nband+2+i+npar+1+i0)
       WT = DEMT_A(*,*,            2*nband+2+i+npar+2+i0)
       Er = DEMT_A(*,*,            2*nband+2+i+npar+3+i0)
       ra = DEMT_A(*,*,            2*nband+2+i+npar+4+i0)
       pa = DEMT_A(*,*,            2*nband+2+i+npar+5+i0)
       ya = DEMT_A(*,*,            2*nband+2+i+npar+6+i0)
       za = DEMT_A(*,*,            2*nband+2+i+npar+7+i0)
     endif

    if  parametrizacion eq 'equalizer' then begin
        npar = nband
        IMA = DEMT_A(*,*,0    :  nband-1)
        SIMA = DEMT_A(*,*,nband:2*nband-1)
        DEMv = DEMT_A(*,*,    2*nband:2*nband+npar-1  )
        DEMc = reform(DEMv(*,*,0))
        i    =npar
        lambda = DEMT_A(*,*,2*nband+i:2*nband+i+npar-1)
        PHI = DEMT_A(*,*,            2*nband+i+npar  )
        ftol = DEMT_A(*,*,           2*nband+i+npar+1)
        i0=2
        EM = DEMT_A(*,*,            2*nband+i+npar  +i0)
        Tm = DEMT_A(*,*,            2*nband+i+npar+1+i0)
        WT = DEMT_A(*,*,            2*nband+i+npar+2+i0)
        Er = DEMT_A(*,*,            2*nband+i+npar+3+i0)
        ra = DEMT_A(*,*,            2*nband+i+npar+4+i0)
        pa = DEMT_A(*,*,            2*nband+i+npar+5+i0)
        ya = DEMT_A(*,*,            2*nband+i+npar+6+i0)
        za = DEMT_A(*,*,            2*nband+i+npar+7+i0)
     endif
  endif

  Nrad=nr
  Nlat=nth
  Nlon=np
  dt = 2.
  dr = 0.01
  rad =   1. + dr/2. + dr * findgen(Nrad)
  lon =   0. + dt/2. + dt * findgen(Nlon)
  lat = -90. + dt/2. + dt * findgen(Nlat)
 ; print,rad

 L=0
 if Tmax gt 3.4e6 and Tmax lt 3.6e6 then L=171
 if Tmax gt 4.9e6 and Tmax lt 5.1e6 then L=256
 if Tmax gt 14.e6 and Tmax lt 16.e6 then L=850
 if Tmax gt 9.e6  and Tmax lt 11.e6 then L=540
 if L eq 0 then STOP
 resol_temp= (tmax - tmin)/L/1.e6
 WTc       = resol_temp*1. ; MK

return
end

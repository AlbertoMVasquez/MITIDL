
pro wrapper


  ;============================
  ;noise_suffix='sin_ruido'
   noise_suffix='con_ruido_0.05'
  exp_suffix  =''
  ;exp_suffix  ='_exp_B'
  ;===========================

  dir_par_input = '/data1/tomography/bindata/'
  file_par_input= 'param_input_exp_contrl_'+noise_suffix+exp_suffix+'.out'

  restore,dir_par_input+file_par_input
  sz = size (par_in)


  Tmin = 0.5                       ; MK
  Tmax = 3.5                       ; MK
  L    = 171   &   Lstring='_171' ; Number of Temp bins
  nr   = sz(1)                   ; Number of tomographic radial bins
  nth  = sz(2)                   ; Number of latitudinal bins
  np   = sz(3)                  ; Number of longitudinal bins
  npx  = 1024                     ; Number of Image pixels, for DEM.
  box  = [400,620,240,380]        ; box of pixeles where do DEM  

   bandsindexes = [2,3,4]   &  Expstring='_AIA3_MIT_exp_contrl_'+noise_suffix+exp_suffix

   finalname = 'LDEM'+Expstring
 
   dir    = '/data1/tomography/bindata/' ;Tomography directory
   ;dir    = '~/Downloads/' ;Tomography directory

   parastring='_gauss1'

   datafiles=[ '',$
               '',$
               'x.aia171_exp_contrl_',$
               'x.aia193_exp_contrl_',$
               'x.aia211_exp_contrl_',$
               '' ]+noise_suffix+exp_suffix+'.out'   &   CRstring='_mit'
   
; The TRF
   file_ioneq='chianti.ioneq'
   file_abund='sun_coronal_1992_feldman_ext.abund'
   dirqkl='/data1/work/dem/TRF/'
   qklfiles=dirqkl+['Qkl_094_','Qkl_131_','Qkl_171_','Qkl_193_','Qkl_211_','Qkl_335_']+$
             file_ioneq+'_'+file_abund+$
            '.AIA-lastoptimizedWL-photons-Abund-1e-3-ALL-withCONTINUUM_Ne1E08_C71.extended.interpolated.out'

   

   suffix = 'v2'+CRstring+'_'+file_ioneq+'_'+file_abund+Expstring+Lstring
   
   binfactor=2
   fact     =1
   ldem_mit,Tmin,Tmax,datafiles,qklfiles,bandsindexes,suffix,dir,L,nr,nth,np,npx,binfactor,fact,$
           /aia,/norm_median,/ldem,/gauss1,/single,/linear,finalname=finalname;,oneheight=6 

   return
end





pro wrapper_old

  Tmin = 0.5                       ; MK
  Tmax = 3.5                       ; MK
  L    = 171   &   Lstring='_171' ; Number of Temp bins
  nr   = 3                       ; Number of tomographic radial bins
  nth  = 3                       ; Number of latitudinal bins
  np   = 3                      ; Number of longitudinal bins
  npx  = 1024                     ; Number of Image pixels, for DEM.
  box  = [400,620,240,380]        ; box of pixeles where do DEM  

   bandsindexes = [2,3,4]   &  Expstring='_AIA3_MIT_exp_contrl'

   finalname = 'LDEM'+Expstring
 
   dir    = '/data1/tomography/bindata/' ;Tomography directory
   dir    = '~/Downloads/' ;Tomography directory

   parastring='_gauss1'

   datafiles=[ '',$
               '',$
               'x171',$
               'x193',$
               'x211',$
               '' ]   &   CRstring='_mit'
   
; The TRF
   file_ioneq='chianti.ioneq'
   file_abund='sun_coronal_1992_feldman_ext.abund'
   dirqkl='/data1/work/dem/TRF/'
   qklfiles=dirqkl+['Qkl_094_','Qkl_131_','Qkl_171_','Qkl_193_','Qkl_211_','Qkl_335_']+$
             file_ioneq+'_'+file_abund+$
            '.AIA-lastoptimizedWL-photons-Abund-1e-3-ALL-withCONTINUUM_Ne1E08_C71.extended.interpolated.out'

   

   suffix = 'v2'+CRstring+'_'+file_ioneq+'_'+file_abund+Expstring+Lstring
   
   binfactor=2
   fact     =1
   ldem_mit,Tmin,Tmax,datafiles,qklfiles,bandsindexes,suffix,dir,L,nr,nth,np,npx,binfactor,fact,$
           /aia,/norm_median,/ldem,/gauss1,/single,/linear,finalname=finalname;,oneheight=6 

   return
end
;----------------------------------------------------------------------------------------------------------------------
; INPUTS      :
; Tminima     : The minimum temperature grid in K.
; Tmaxima     : The maximum temperature grid in K.
; datafiles   : The image or tomography files to do DEM or DEMT.
; qklfiles:   : The TRF files.
; bandsindexes: The indexes of band to use in the DEM invertion.
; suffix      : The suffix to use in the name of the DEM file.
; dir         : The directory where are the image or tomography files.
; L           : The temperature bins.
; nr          : The radial bins of the tomographic grid.
; nth         : The latitudinal bins of the tomographic grid.
; np          : The longitudinal bins of the tomographic grid.
; npx         : The number of pixel of the image.
; binfactor   : The factor to rebin the image if flag rebin is used.
; fact        : This variable is from a old version of this code and
;               is not use now.

; KEYWORDS    :
; box         : The box of pixeles where do the DEM invertion [ixa,ixb,iya,iyb]. If
;               you dont set this keyword the code use the full image.
; sustrato    : A AIA-3/1-normal DEM or (LDEM) results to use as initial guess
;               for AIA-4/2-normal 
; finalname   : If set the final name of the LDEM file can be choosen
; oneheight   : If set the DEMT is calculated in one height

; FLAGS       : 
; DEM         : If set a 2D DEM is calculated.
; LDEM        : If set a 3D DEMT is calculated.
; euvi        : Set if use EUVI data.
; aia         : Set if use AIA data.
; gauss1      : If set a 1-normal parametric model is used
; doublegauss : If set a 2-normal parametric model (with six free
;               parameters) is used.
; gauss2c     : If set a 2-normal parametric model (the centroid
;               of the second normal distribution fixed) with five
;               free parameters.
; dgfw       : If set a 2-normal parametric model (the width                                                                                                        
;              of the second normal distribution fixed) with five free
;              parameters.
   
;---------------------------------------------------------------------------------------------------------------
pro LDEM_MIT,Tminima,Tmaxima,datafiles,qklfiles,bandsindexes,suffix,dir,L,nr,nth,np,npx,binfactor,fact,$
            ldem=ldem,dem=dem,linear=linear,loga=loga,custom=custom,expe=expe,$
            norm_ones=norm_ones,norm_mean=norm_mean,norm_median=norm_median,norm_fbe=norm_fbe,$
            gauss1=gauss1,doublegauss=doublegauss,multiple=multiple,single=single,rebin=rebin,euvi=euvi,aia=aia,$
            newpar4=newpar4,newpar5=newpar5,gauss2c=gauss2c,glc=glc,lorentz=lorentz,gauss3c=gauss3c,gauss7par=gauss7par,$
            equalizer=equalizer,dgfw=dgfw,g2cfw=g2cfw,qgfw=qgfw,triplegauss=triplegauss,$
            box=box,sustrato=sustrato,finalname=finalname,oneheight=oneheight,loops=loops
  
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
common weights,Ck
common parametrizacion,parametrization
common paraloss,temp,loss
common intent,intento
common gauss_cond,T335,X335,Q335_max,lambda_inicial,DEMc1,DEMc2,DEMc3,T131,X131,Q131_max,SIGMA,X284,T284
common read_sustrato_gauss2c,lambda_inicial_A,demc1_A
common read_sustrato_gauss3c,lambda_inicial_Array,demc1_Array,demc2_Array
common centroid_equal,xmax,sigma_v,DEMcv,factor

 ;Set factor variable
 factor=fact

 ;Change to Kelvin T-range
 Tmin   = Tminima*1.e6
 Tmax   = Tmaxima*1.e6
 nTbins = L

 ;Set tiny
  tiny=1.d-100

; Keep only the needed datafiles and Qkl, define nband
 datafiles =  datafiles(bandsindexes)
  qklfiles =   qklfiles(bandsindexes)
     nband = n_elements(qklfiles)

 if keyword_set(aia ) then instrument='aia'
 if keyword_set(euvi) then instrument='euvi'

 print,'  T-range [K]: ',tmin,tmax
 print,' Band indexes: ',bandsindexes
 print,'        Tbins: ',L
 print,'   Data files: '
 print,transpose( datafiles)
 print,'    Qkl files: '
 print,transpose(  qklfiles)
 if keyword_set(box) then print,'range of pixels where do DEM:',box
 if keyword_set(loops) then print,'Bright AR loops DEM analysis:',loops

 if keyword_set(loga  ) then begin
    sca = 0 ; X-grid regular in log(T)
    gridstring='_log'
 endif
 if keyword_set(linear) then begin
    sca = 1 ; X-grid regular in     T
    gridstring='_lin'
 endif
 if keyword_set(custom) then begin
    sca = 2 ; X-grid custom made
    gridstring='_custom'
 endif

 read_Qkl,qklfiles,logT,kernel,nband
 Interpol_Qkl,logT,kernel,nband,nTbins,sca,Tmin,Tmax,Qkl,X,DX,L,logaT,instrument

;----------------------------------------------------------------
 filename='loss_func_chianti.ioneq_sun_coronal.abund.new.out' ; calculado con el CHIANTI reinstalado
 filename='loss_func_chianti.ioneq_sun_coronal.abund.out'
;filename='loss_func_arnaud_raymond.ioneq_sun_coronal.abund.out'
 N=10000.
 read_loss,'/data1/work/dem/TRF/'+filename,temp,loss
 Interpol_loss,N,temp,loss
;----------------------------------------------------------------

;Determine number of data values to store for each voxel according to chosen parametrization:

 if keyword_set(gauss1) then begin
  nsingleT = 2
  npar     = 3
  nPHI     = 2
  nmom     = 3
  nloss    = 1
  ndat     = 2 * nband + nsingleT + npar + nPHI + nmom + nloss
  parametrization = 'gauss1'
  endif

;+
if keyword_set(lorentz) then begin
  nsingleT = 2
  npar     = 3
  nPHI     = 2
  nmom     = 3
  nloss    = 1
  ndat     = 2 * nband + nsingleT + npar + nPHI + nmom + nloss
  parametrization = 'lorentz'
 endif
;-
 if keyword_set(doublegauss) then begin
  nsingleT = 2 
  npar     = 6
  nPHI     = 2
  nmom     = 3
  nloss    = 1
  ndat     = 2 * nband  + nsingleT + npar + nPHI + nmom + nloss
  parametrization = 'gauss2'
 endif

 if keyword_set(g2cfw) then begin
  nsingleT = 2 
  npar     = 4
  nPHI     = 2
  nmom     = 3
  nloss    = 1
  ndat     = 2 * nband  + nsingleT + npar + nPHI + nmom + nloss + 1
  parametrization = 'gauss2cFW'
 endif

;+
if keyword_set(gauss2c) or keyword_set(glc) or keyword_set(dgfw) then begin
  nsingleT = 2
  npar     = 5
  nPHI     = 2
  nmom     = 3
  nloss    = 1
  ndat     = 2 * nband  + nsingleT + npar + nPHI + nmom + nloss
endif 
if keyword_set(gauss2c) then parametrization = 'gauss2c'
if keyword_set(glc) then parametrization = 'glc'
if keyword_set(dgfw) then parametrization= 'gauss2FW' 
if keyword_set(dgfw) then ndat=ndat+1 
;-

;FEDE--------------------------------------------------------------
if keyword_set(qgfw) then begin
  nsingleT = 4
  npar     = 7
  nPHI     = 2
  nmom     = 3
  nloss    = 1
  ndat     = 2 * nband  + nsingleT + npar + nPHI + nmom + nloss 
  parametrization='gauss4FW'
endif

if keyword_set(triplegauss) then begin
  nsingleT = 3
  npar     = 8
  nPHI     = 2
  nmom     = 3
  nloss    = 1
  ndat     = 2 * nband  + nsingleT + npar + nPHI + nmom + nloss 
  parametrization='gauss3FW'
endif
;--------------------------------------------------------------------

;+
if keyword_set(gauss3c) then begin
  nsingleT = 3
  npar     = 6
  nPHI     = 2
  nmom     = 3
  nloss    = 1
  ndat     = 2 * nband  + nsingleT + npar + nPHI + nmom + nloss
  parametrization = 'gauss3c'
endif

if keyword_set(gauss7par) then begin
  nsingleT = 3
  npar     = 7
  nPHI     = 2
  nmom     = 3
  nloss    = 1
  ndat     = 2 * nband  + nsingleT + npar + nPHI + nmom + nloss
  parametrization = 'gauss7par'
endif
;-


if keyword_set(newpar4) then begin
  nsingleT = 2
  npar     = 4
  nPHI     = 2
  nmom     = 3
  nloss    = 1
  ndat     = 2 * nband + nsingleT + npar + nPHI + nmom + nloss
  parametrization = 'newpar4'
endif

if keyword_set(newpar5) then begin
  nsingleT = 2
  npar     = 5
  nPHI     = 2
  nmom     = 3
  nloss    = 1
  ndat     = 2 * nband + nsingleT + npar + nPHI + nmom + nloss
  parametrization = 'newpar5'
endif

if keyword_set(equalizer) then begin
  nsingleT = nband
  npar     = nband
  nPHI     = 2
  nmom     = 3
  nloss    = 1
  ndat     = 2 * nband + nsingleT + npar + nPHI + nmom + nloss
  parametrization = 'equalizer'
endif

;For Images also record the pixels' coordinates
 if keyword_set(dem) then ndat = ndat + 4
 

 suffix=suffix+'_'+parametrization+gridstring

;Define the arrays to store the data and the results
; At the end FBE_A is a 4-dimensional array containing NBAND
; tomographic reconstructed FBEs
       ;np = 2*nth ; NOW, np is an imput of the program..
if keyword_set(ldem) then begin
    FBE_A = fltarr(nr,nth,np,nband)
   DEMT_A = dblarr(nr,nth,np,ndat )
      tmp = fltarr(nr,nth,np)
   for i=0,nband-1 do begin
     openr,1,dir+datafiles(i)
     readu,1,tmp
     FBE_A(*,*,*,i)=tmp
     close,1
   endfor
endif
;+
if keyword_set(gauss2c) or keyword_set(glc) or keyword_set(doublegauss) $
   or keyword_set(dgfw) or keyword_set(g2cfw) then begin
;=======================================================>
   if keyword_set(ldem) then read_lambda_inicial_3band,sustrato=sustrato,/ldem
   if keyword_set( dem) then read_lambda_inicial_3band,sustrato=sustrato,/dem
;=======================================================<
  centroide_gaussiana_condicionada,nband,X,Qkl,Tmin
endif
;-

;+
if keyword_set(gauss3c) or $
   keyword_set(gauss7par) or $
   keyword_set(qgfw) or $
   keyword_set(triplegauss) then begin
 
   read_lambda_inicial_4band
   centroide_gaussiana_condicionada,nband,X,Qkl,Tmin
endif

if keyword_set(equalizer) then begin
centroide_equalizer,nband,X,Qkl,Tmin
endif
;-

; If doing DEM, FBE_A is a 3-dimensional array containing NBAND
; Images.
if keyword_set( dem) then begin
   if not keyword_set(loops) then begin
      if keyword_set(rebin) then npx=npx/binfactor
      FBE_A = fltarr(npx,npx,nband)
      for i=0,nband-1 do begin
         mreadfits,dir+datafiles(i),hdr,ima
         if NOT keyword_set(rebin) then FBE_A(*,*,i) = ima
         if keyword_set(rebin) then begin
            FBE_A(*,*,i)=rebin(ima,npx,npx)
            hdr.CRPIX1 = 1 + (hdr.CRPIX1 -1) * (hdr.NAXIS1/binfactor - 1) / (hdr.NAXIS1-1)
            hdr.CRPIX2 = 1 + (hdr.CRPIX2 -1) * (hdr.NAXIS2/binfactor - 1) / (hdr.NAXIS2-1)
            hdr.NAXIS1 = hdr.NAXIS1 / binfactor
            hdr.NAXIS2 = hdr.NAXIS2 / binfactor
            hdr.CDELT1 = hdr.CDELT1 * binfactor
            hdr.CDELT2 = hdr.CDELT2 * binfactor
         endif
      endfor
      if keyword_set(rebin) then begin
         p=where(FBE_A le 0.)
         if p(0) ne -1 then FBE_A(p)=-999.
      endif  
      if keyword_set(euvi) then ANGLE = -hdr.crota2
      if keyword_set(aia ) then ANGLE = -hdr.crota2 
      if ANGLE ne 0. then rotate,FBE_A,hdr,nband,ANGLE
      if keyword_set(euvi) then compute_grid,hdr,ra,pa,ya,za,/euvi
      if keyword_set(aia ) then compute_grid,hdr,ra,pa,ya,za,/aia
      nx=npx
      ny=npx
      if keyword_set(box)  then begin  
         ixa=box(0)
         ixb=box(1)
         iya=box(2)
         iyb=box(3)
         FBE_A=FBE_A(ixa:ixb,iya:iyb,*)
         ra=ra(ixa:ixb,iya:iyb)
         pa=pa(ixa:ixb,iya:iyb)
         ya=ya(ixa:ixb,iya:iyb)
         za=za(ixa:ixb,iya:iyb)
         nx=ixb-ixa+1
         ny=iyb-iya+1
      endif 
   endif else begin
      ;======LAPAPA===================
      dirloop='/data1/work/loops/'
      restore,dirloop+loops
      rloop  = sqrt(xloop^2+yloop^2)
      nw     =(size(cube))(2)
      nseg   =(size(cube))(3)
      nx     =1;fix(nw)
      ny     =fix(nseg)
      FBE_A  =f_loop(*,*,bandsindexes)
      ;===============================
   endelse
   DEMT_A = dblarr(nx,ny,ndat ) 
endif


;Read and store tomographic FBE or images Intensities, always call it FBE

 ;******************************************************************
 ;For LDEM::
 ;CORRECT FBE UNITS TO ACCOUNT FOR W* USED IN RECONSTRUCTION
 if keyword_set(ldem) then FBE_A = FBE_A ;/ 0.0696
 ;For  DEM::
 ;CORRECT  THE IMAGE in a factor 1x10^12 
  if keyword_set( dem) then FBE_A = FBE_A * 1.e12
 ;******************************************************************

;Store the global mean, median and min values of all FBE>0
    Nones       = fltarr(nband) + 1.
    GLOBAL_MEAN = FLTARR(nband)
  GLOBAL_MEDIAN = FLTARR(nband)
     GLOBAL_MIN = FLTARR(nband)
 for ib=0,nband-1 do begin
     if keyword_set(ldem) then tmp=reform(FBE_A(*,*,*,ib))
     if keyword_set( dem) then tmp=reform(FBE_A(*,*  ,ib))
     p=where(tmp gt 0.)
     global_mean(ib)=  mean(tmp(p))
   global_median(ib)=median(tmp(p))
      global_min(ib)=   min(tmp(p))
 endfor

;Choose weights for cost function:
 flag_norm_fbe=0
 if keyword_set(norm_mean)   then begin
    Ck =   GLOBAL_MEAN   &   print,'weights: MEANs'
    suffix=suffix+'_Norm-mean'
 endif
 if keyword_set(norm_median) then begin
    Ck = GLOBAL_MEDIAN   &   print,'weights: MEDIANs'
    suffix=suffix+'_Norm-median'
 endif
 if keyword_set(norm_ones)   then begin
    Ck = Nones           &   print,'no weights'
    suffix=suffix+'_Norm-one'
 endif
 if keyword_set(norm_fbe)    then begin
    flag_norm_fbe=1      &   print,'weights: local FBE or Intensity'
    suffix=suffix+'_Norm-fbe'
 endif

;Choose multiple or single starting point
 if keyword_set(single)   then begin
    multiple_flag=0
    suffix=suffix+'_singlStart'
 endif
 if keyword_set(multiple) then begin 
    multiple_flag=1
    suffix=suffix+'_multiStart' 
 endif

;-----------------------------------------------------------------------------------------
;Call test experiment module
if keyword_set(expe) then begin
                          box=0
   if keyword_set(b) then box=1
   Xc=0.
  ;experiment,sca,Xc,logaT,Q094,Q131,Q171,Q193,Q211,Q335,Tmin,noiselevel,box,exper,Ck
   experiment,sca,Xc,logaT,Qkl,noiselevel,box,exper,Ck
   goto,the_end
endif
;-----------------------------------------------------------------------------------------

print,'Full Output Suffix: '+suffix

;--------------------------------
; TOMOGRAPHY FBE
if keyword_set(ldem) then begin
   irad1=0                      ;1.025 Rsun
   irad2=Nr-1;10                ;1.255 Rsun
   if keyword_set(oneheight) then begin
      irad1=oneheight
      irad2=irad1
   endif
   Dirad=1
   ilat1=0
   ilat2=nth-1
   Dilat=1
   ilon1=0
   ilon2=np-1
   Dilon=1
  ;--------------------------
   goto,skiptesting
   ilat1 = 5                    ;0
   ilat2 = ilat1                ;nth-1
   Dilat = 1                    ;1
   ilon1 = 0
   ilon2 = np-1
   Dilon = 5                    ;1
   skiptesting:
   for ir =irad1,irad2,Dirad do begin
      for ith=ilat1,ilat2,Dilat do begin
         for ip =ilon1,ilon2,Dilon do begin
            goto,computation_block
            continue_tomography:
         endfor                 ; IP  loop
      endfor                    ; ITH loop
   endfor                       ; IR loop
endif                           ; END TOMOGRAPHY
;--------------------------------

;--------------------------------
; IMAGE DEM
if keyword_set(dem) then begin
  ix1=0
  ix2=nx-1
  iy1=0
  iy2=ny-1
    ;--------------------------
    ;for testing-purposes only:
     goto,skiptesting2
     ix1=0
     ix2=npx-1
     iy1=511
     iy2=512
    skiptesting2:
    ;--------------------------
     for iy=iy1,iy2    do begin  
        for ix=ix1,ix2    do begin
           goto,computation_block
           continue_image:
        endfor                  ; iy loop
     endfor                     ; ix loop
  endif                         ; END IMAGE
;--------------------------------

;Save results and finish
 goto,save

;--------------------------------------------------------------------------------------------------
computation_block:
 if keyword_set(ldem) then FBE=reform( FBE_A(ir,ith,ip,*) )
 if keyword_set( dem) then FBE=reform( FBE_A(ix,iy,    *) )

  if min(FBE) le 0. then begin
    DEMc=-999.
      Tc=-999.

     if keyword_set(gauss2c) or keyword_set(glc) or keyword_set(doublegauss) $
        or keyword_set(dgfw) or keyword_set(g2cfw) then begin
    DEMC1=-999.
    DEMc2=-999.
     endif


     if keyword_set(gauss3c) or keyword_set(gauss7par) then begin
    DEMC1=-999.
    DEMc2=-999.
    DEMc3=-999.
     endif
     
     if keyword_set(equalizer) then begin
    DEMCv= -999. + dblarr(nband)
    endif



     out=fltarr(npar)-1.
     goto,skipvoxelanalysis
  endif


if keyword_set(gauss2c) or keyword_set(glc) or keyword_set(doublegauss) $
   or keyword_set(dgfw) or keyword_set(g2cfw) then begin 
;=======================================================>
  if keyword_set(ldem) then begin
  lambda_inicial=reform(lambda_inicial_A(ir,ith,ip,*))
  DEMc1         =demc1_A                (ir,ith,ip  )
  endif
  if keyword_set(dem) then begin
  lambda_inicial=reform(lambda_inicial_A(ix,iy,*))
  DEMc1         =demc1_A                (ix,iy  )
  endif
;=======================================================<
  DEMc2         =FBE(nband-1)/Q335_max
  DEMc3         =FBE(0)/Q131_max ; mientras no este 94
endif
;-

;+
if keyword_set(gauss3c) or keyword_set(gauss7par) then begin 
  lambda_inicial=reform(lambda_inicial_Array(ir,ith,ip,*))
  DEMc1         =demc1_Array         (ir,ith,ip  )
  DEMc2         =demc2_Array         (ir,ith,ip  )
  DEMc3         =FBE(0)/Q131_max ; mientras no este 94
if nband eq 6 then DEMc3         =FBE(1)/Q131_max
endif
;-


if keyword_set(equalizer) then begin
DEMcv=dblarr(nband)
 for i=0,nband-1 do begin
 DEMcv(i) = FBE(i)/max(reform(Qkl(i,*)))
 endfor
endif


  DEMc=0.d
    Tc=0.
  new_guess=fltarr(npar)-666.
  intento=0

  if flag_norm_fbe eq 1 then Ck=FBE

  normalize,FBE,Qkl,nband

;+

  
  if keyword_set(dgfw) and nband eq 5 then begin
     ;r2=where(X*Tmin/1.e6 ge 2. and X*Tmin/1.e6 le 3.)
     ;X_2=X(r2)
     ;dx_2=dx(r2)
     FBE_2=FBE([nband-4,nband-3])
     Qkl_2=Qkl([nband-4,nband-3],*)
     ;ntbins_2=n_elements(x_2)
     ;singletempsolution,sca,Tmin,FBE_2,logaT,Qkl_2,DEMc_2,Tc,Xc,nband,nTbins_2,X_2,dX_2
     second_gaussian_guess,sca,Tmin,FBE_2,logaT,Qkl_2,DEMc2,Tc,Xc,nband,nTbins,X,dX,sigma
  ;demc2=demc_2
  endif
  

  if keyword_set(gauss1) then $
  singletempsolution,sca,Tmin,FBE,logaT,Qkl,DEMc,Tc,Xc,nband,nTbins,X,dX
  if keyword_set(lorentz)then $
  singletempsolution,sca,Tmin,FBE,logaT,Qkl,DEMc,Tc,Xc,nband,nTbins,X,dX
  if keyword_set(newpar4)then $
  singletempsolution,sca,Tmin,FBE,logaT,Qkl,DEMc,Tc,Xc,nband,nTbins,X,dX
  if keyword_set(newpar5)then $
  singletempsolution,sca,Tmin,FBE,logaT,Qkl,DEMc,Tc,Xc,nband,nTbins,X,dX

  if keyword_set(doublegauss)then DEMc=DEMc1
  if keyword_set(dgfw)       then DEMc=DEMc1
  if keyword_set(g2cfw)      then DEMc=DEMc1
  if keyword_set(gauss2c)    then DEMc=DEMc1
  if keyword_set(glc)        then DEMc=DEMc1
  if keyword_set(gauss3c)    then DEMc=DEMc1
  if keyword_set(gauss7par)  then DEMc=DEMc1
  if keyword_set(equalizer)  then DEMc=total(DEMcv)
;-


  
  retry:
  retry_multiple:
 ;print,'flag: ',multiple_flag

  if keyword_set(gauss1)  then begin
    if multiple_flag eq 0 then start_point,new_guess,Xc,nn,new_guess_a,npar
    if multiple_flag eq 1 then start_point,new_guess,Xc,nn,new_guess_a,npar,/multiple
  endif
;+
 if keyword_set(lorentz)  then begin
    if multiple_flag eq 0 then start_point_lorentz,new_guess,Xc,nn,new_guess_a,npar
    if multiple_flag eq 1 then start_point_lorentz,new_guess,Xc,nn,new_guess_a,npar,/multiple
  endif
;-
  if keyword_set(newpar4) then begin
    if multiple_flag eq 0 then start_point_newpar4,new_guess,Xc,nn,new_guess_a,npar
    if multiple_flag eq 1 then start_point_newpar4,new_guess,Xc,nn,new_guess_a,npar,/multiple
  endif

  if keyword_set(newpar5) then begin
    if multiple_flag eq 0 then start_point_newpar5,new_guess,Xc,nn,new_guess_a,npar
    if multiple_flag eq 1 then start_point_newpar5,new_guess,Xc,nn,new_guess_a,npar,/multiple
  endif

  if keyword_set(doublegauss) then begin
     if multiple_flag eq 0 then start_point_double,new_guess,Xc,nn,new_guess_a
     if multiple_flag eq 1 then start_point_double,new_guess,Xc,nn,new_guess_a,/multiple
  endif

  if keyword_set(dgfw) then begin
     if multiple_flag eq 0 then start_point_dgfw,new_guess,Xc,nn,new_guess_a
     if multiple_flag eq 1 then start_point_dgfw,new_guess,Xc,nn,new_guess_a,/multiple
  endif

  if keyword_set(g2cfw) then begin
     if multiple_flag eq 0 then start_point_g2cfw,new_guess,Xc,nn,new_guess_a
     if multiple_flag eq 1 then start_point_g2cfw,new_guess,Xc,nn,new_guess_a,/multiple
  endif
;+
  if keyword_set(gauss2c) or keyword_set(glc) then begin
     if multiple_flag eq 0 then start_point_2c,new_guess,Xc,nn,new_guess_a
     if multiple_flag eq 1 then start_point_2c,new_guess,Xc,nn,new_guess_a,/multiple
  endif

 if keyword_set(gauss3c) then begin
     if multiple_flag eq 0 then start_point_3c,new_guess,Xc,nn,new_guess_a
     if multiple_flag eq 1 then start_point_3c,new_guess,Xc,nn,new_guess_a,/multiple
  endif

if keyword_set(gauss7par) then begin
     if multiple_flag eq 0 then start_point_7par,new_guess,Xc,nn,new_guess_a
     if multiple_flag eq 1 then start_point_7par,new_guess,Xc,nn,new_guess_a,/multiple
  endif

if keyword_set(equalizer) then begin
     if multiple_flag eq 0 then start_point_eq,new_guess,Xc,nn,new_guess_a
     if multiple_flag eq 1 then start_point_eq,new_guess,Xc,nn,new_guess_a,/multiple
  endif


;-
  out_A=fltarr(nn,npar)
  PHI_A=fltarr(nn)
 ftol_A=fltarr(nn)

  for i=0,nn-1 do begin
   minPHI,smplx,sca,out,dem_syn,reform(new_guess_a(i,*)),PHIv,ftol
   out_A(i,*)=out
   PHI_A(i  )=PHIv
  ftol_A(i  )=ftol
  endfor

 ; Select lambda that corresponds to global minimum of PHI_A:
   p=(where(PHI_A eq min(PHI_A)))(0) 
   out=reform(out_A(p,*))
   ftol=ftol_A(p)

;  if multiple_flag eq 1 then STOP

   Residual1 = total((1.-SynthFBE(out)/FBE)^2)/float(nband)
   PHIfinal2 = min(PHI_A)/float(nband)
   PHIfinal  = PHIfinal2

; if nn eq 9 and min(PHI_A) gt 0.01 then out=-999.+fltarr(npar)
   
   goto,skip_multiple
   if keyword_set(gauss1) then begin
      if multiple_flag eq 0 and PHIfinal gt 1.e-3 then begin
         multiple_flag=1
         goto,retry
      endif
      if multiple_flag eq 1 then multiple_flag=0
   endif
  skip_multiple:

goto,skipsocotroco
  if keyword_set(newpar5) then begin
    ;if intento eq 0 AND min(PHI_A) gt 1.e-3 then begin
     if intento eq 0 AND Residual1  gt 1.e-3 then begin
        intento=1 
        print,'Residual1=',Residual1,'... reintentando'
        goto,retry
     endif
     if intento eq 1 then intento=0
       print,'Residual1=',Residual1,'.. es el valor mejorado'
  endif
skipsocotroco:

  unnormalize,fbe,qkl,nband

  skipvoxelanalysis:

;Compute LDEM-moments and Radiative Loss rate:
  if demc ne -999. then begin
        T_X = X * Tmin ; Temperatures in KELVIN:
        Ne2 = TOTAL(               demc*dX*FDEM(X,out) )
        Tm  = TOTAL(  T_X        * demc*dX*FDEM(X,out) ) / Ne2 ; K
        WT2 = TOTAL( (T_X-Tm)^2  * demc*dX*FDEM(X,out) ) / Ne2 ; K
        WT  = Sqrt(WT2)
        R=total(abs(1-SynthFBE(out)/FBE))/nband ; score R

        if keyword_set(ldem) then density_product = Sqrt(Ne2) ; Local density
        if keyword_set( dem) then density_product =      Ne2  ; Total LOS EM
        if keyword_set(loops)then w_cm=loop_par(iy,6)
        if keyword_set(loops)then density_product = Sqrt(Ne2/w_cm)
        Er  = 0. 
        if keyword_set(ldem) then Er = TOTAL( Loss_func(T_X) * demc*dX*FDEM(X,out) ) ; erg cm-3 sec-1
     endif else begin
        density_product = -999.
        Tm = -999.
        WT = -999.
        Er = -999.
        PHI_A =[-999.]
        PHIfinal = -999.
        ftol = -999.
        R    = -999.
     endelse

  
   
   ; For Lloveras et al (2020), We want minimization with multiple start-point if
   ; the LDEM have R > 0.01 and Tm ~ 0.5 MK. 
     ;goto,skip_multiple2
     if keyword_set(gauss1) and demc ne -999. then begin
        if multiple_flag eq 0 and R gt 1.e-2 and Tm lt 0.6e6 then begin
           multiple_flag=1
           normalize,FBE,Qkl,nband
           goto,retry_multiple
        endif
        if multiple_flag eq 1 then multiple_flag=0
     endif
     skip_multiple2:


 ; Store all results in single output array
 
 if keyword_set(gauss1 )    then results = [FBE,SynthFBE(out),DEMc,Tc,        out,PHIfinal,ftol,density_product,Tm,WT,Er]
 if keyword_set(lorentz)    then results = [FBE,SynthFBE(out),DEMc,Tc,        out,PHIfinal,ftol,density_product,Tm,WT,Er]
 if keyword_set(newpar4)    then results = [FBE,SynthFBE(out),DEMc,Tc,        out,PHIfinal,ftol,density_product,Tm,WT,Er]
 if keyword_set(newpar5)    then results = [FBE,SynthFBE(out),DEMc,Tc,        out,PHIfinal,ftol,density_product,Tm,WT,Er]
 if keyword_set(doublegauss)then results = [FBE,SynthFBE(out),DEMc1,DEMc2,    out,PHIfinal,ftol,density_product,Tm,WT,Er]
 if keyword_set(dgfw)       then results = [FBE,SynthFBE(out),DEMc1,DEMc2,sigma,    out,PHIfinal,ftol,density_product,Tm,WT,Er]
 if keyword_set(g2cfw)      then results = [FBE,SynthFBE(out),DEMc1,DEMc2,sigma,    out,PHIfinal,ftol,density_product,Tm,WT,Er]
 if keyword_set(gauss2c)    then results = [FBE,SynthFBE(out),DEMc1,DEMc2,          out,PHIfinal,ftol,density_product,Tm,WT,Er]
 if keyword_set(glc    )    then results = [FBE,SynthFBE(out),DEMc1,DEMc2,          out,PHIfinal,ftol,density_product,Tm,WT,Er]
 if keyword_set(gauss3c)    then results = [FBE,SynthFBE(out),DEMc1,DEMc2,DEMc3,    out,PHIfinal,ftol,density_product,Tm,WT,Er]
 if keyword_set(gauss7par)  then results = [FBE,SynthFBE(out),DEMc1,DEMc2,DEMc3,    out,PHIfinal,ftol,density_product,Tm,WT,Er]
 if keyword_set(equalizer)  then results = [FBE,SynthFBE(out),DEMcv            ,    out,PHIfinal,ftol,density_product,Tm,WT,Er]
 ;-
 if keyword_set(ldem) then DEMT_A(ir,ith,ip,*) =  results
 if keyword_set( dem) and not keyword_set(loops) then DEMT_A(ix,iy    ,*) = [results,ra(ix,iy),pa(ix,iy),ya(ix,iy),za(ix,iy)]
 if keyword_set( dem) and     keyword_set(loops) then DEMT_A(ix,iy    ,*) = [results,xloop(iy),yloop(iy),rloop(iy),sloop(iy)]
if keyword_set(ldem) then begin
 if demc ne -999. then $
 print,ir,ith,ip,reform([density_product/1.e8,Tm/1.e6,WT/1.e6, R])
 if demc eq -999. then $
 print,ir,ith,ip,reform([density_product     ,Tm     ,WT     ])
endif

if keyword_set(dem) then begin
 if demc ne -999. then $
 print,ix,iy,reform([density_product,Tm/1.e6,WT/1.e6, R])
 if demc eq -999. then $
 print,ix,iy,reform([density_product       ,Tm     ,WT    ])
endif

;CELDA
 goto,skipcelda
;goto,skipgraph
 window,1
 plot,x*Tmin/1.e6,qkl(1,*)/max(qkl),/nodata,ystyle=1,xstyle=1,xtitle='T [MK]',yr=[0,1]
for ib=0,nband-1 do oplot,x*Tmin/1.e6,qkl(ib,*)/max(qkl);,psym=1
                    oplot,x*Tmin/1.e6,fdem(x,out)/max(fdem(x,out)),psym=4
skipgraph:
inte=dblarr(nband)
for ib=0,nband-1 do inte(ib)=int_tabulated(X,DEMc*FDEM(X,out)*Qkl(ib,*))
print,'     Sums/Inte:',((DEMc*dX*FDEM(X,out))##Qkl)/inte
print,' synth/tom FBE:',SynthFBE(out)/FBE
print,'        FBE/Ck:',1.d*FBE/Ck
print,'      PHIfinal:',PHIfinal
;print,'      residual:',residual1
print,'          ftol:',ftol
print,'        lambda:',out
stop
;if PHIfinal gt 1.e-1 then stop
skipcelda:

if keyword_set(ldem) then goto,continue_tomography
if keyword_set( dem) then goto,continue_image
;--------------------------------------------------------------------------------------------------

save:

  if keyword_set(ldem) then output_file = 'LDEM.'+suffix
  if keyword_set( dem) then output_file =  'DEM.'+suffix
  if keyword_set( dem) and keyword_set(loops) then output_file =  'DEM.loops.'+suffix

  if keyword_set(ldem) then resultado = 'LDEM'
  if keyword_set( dem) then resultado =  'DEM'
  
  if keyword_set(finalname) then output_file=finalname

  openw,1,'/data1/DATA/ldem_files/'+output_file
  writeu,1,resultado
  writeu,1,parametrization
  writeu,1,sca,Tmin,Tmax,nr,nth,np,npx,nband,ndat
;;---------------------------FEDE----------------------------
 if keyword_set( dem) then writeu,1,nx,ny
;;---------------------------FEDE----------------------------
 if keyword_set(equalizer) then  writeu,1,xmax,sigma_v
 writeu,1,DEMT_A
  close,1
 print,'Processing finished. All results are in: ',output_file

 ; save the results also in .sav file
 save,filename='/data1/DATA/ldem_files/'+output_file+'.sav',DEMT_A
 

the_end:
 return
end            

;--------------------------------------------------------------------------------
pro read_Qkl,qklfiles,logT,kernel,nband

x=''
N_e=0.
logem=0.
inst_fwhm=0.
Ntemp=0.
w1=0.
w2=0.
wbin=0.
Nlambda=0.
psi=0.

psi_const=fltarr(nband)
tt=0.
ii=0.

for ip=0,nband-1 do begin
openr,1,qklfiles(ip)
for i=1,4 do readf,1,x
readf,1,n_e,logem
readf,1,x
readf,1,Inst_fwhm, NTemp, w1, w2, wbin, Nlambda
readf,1,x
readf,1,psi
psi_const(ip)=psi
readf,1,x
if ip eq 0 then begin
kernel=fltarr(nband,ntemp)
logT  =fltarr(ntemp)
endif
for i=0,ntemp-1 do begin
 readf,1,tt,ii
 if ip eq 0 then logT(i)=tt
 kernel(ip,i)=ii
endfor
close,1
endfor
kernel=kernel*1.d0

return
end

;+------------------------------------------------------------------

;=======================================================>
pro read_lambda_inicial_3band,sustrato=sustrato,dem=dem,ldem=ldem
;=======================================================<
common read_sustrato_gauss2c,lambda_inicial_A,demc1_A
sca_aux  =0
Tmin_aux =0.
Tmax_aux =0.
nr_aux   =0
nth_aux  =0
np_aux   =0
npx_aux  =0
nband_aux=0L
ndat_aux =0L
nx_aux=0
ny_aux=0
;=======================================================>
if keyword_set(ldem) then resultado_aux = 'LDEM'
if keyword_set( dem) then resultado_aux = 'DEM'
;=======================================================<
parametrizacion_aux='gauss1'
npar_aux=3
dir ='/data1/DATA/ldem_files/'

if not keyword_set(sustrato) then begin
   print,'el codigo necesita un archivo de sustrato'
   stop
endif
;=======================================================>
file =sustrato
;=======================================================<

  print,'archivo de sustrato: ',file

  openr,1,dir+file
  readu,1,resultado_aux
  readu,1,parametrizacion_aux
  readu,1,sca_aux,Tmin_aux,Tmax_aux,nr_aux,nth_aux,np_aux,npx_aux,nband_aux,ndat_aux
  if keyword_set(dem) then readu,1,nx_aux,ny_aux
;=======================================================>
  if keyword_set(ldem) then  DEMT_A = dblarr(nr_aux  ,nth_aux,np_aux,ndat_aux) 
  if keyword_set( dem) then  DEMT_A = dblarr(nx_aux  ,ny_aux        ,ndat_aux) 
  readu,1,DEMT_A
  close,1
  i=0
  if keyword_set(ldem) then begin
  lambda_inicial_A = DEMT_A(*,*,*,2*nband_aux+2+i:2*nband_aux+2+i+npar_aux-1)
  DEMc1_A          = DEMT_A(*,*,*,2*nband_aux  )
  endif
  if keyword_set( dem) then begin
  lambda_inicial_A = DEMT_A(*,*  ,2*nband_aux+2+i:2*nband_aux+2+i+npar_aux-1)
  DEMc1_A          = DEMT_A(*,*  ,2*nband_aux  )
  endif
;=======================================================<

;lambda_inicial_A(*,*,*,0) = sqrt(abs( lambda_inicial_A(*,*,*,0)) ) ;
;ya estan modificadas las parametrizaciones
return
end
;-

;+------------------------------------------------------------------

pro read_lambda_inicial_3band_exp
common read_sustrato_gauss2c,lambda_inicial_A,demc1_A

sca_aux  =0
Tmin_aux =0.
Tmax_aux =0.
nc_aux   =0
nw_aux   =0
nband_aux=0L
ndat_aux =0L

;if keyword_set(ldem) then 
resultado_aux = 'LDEM_simulated'
;if keyword_set( dem) then resultado_aux = 'DEM'


parametrizacion_aux='gauss1'
 
npar_aux=3

 ;dir ='/data1/work/dem/'                                                                                                                                   
  dir ='./ldem_files/'

  
 file ='LDEM.v2_simulated_chianti.ioneq_sun_coronal.abund_AIA-3band_128_gauss1_lin_Norm-median_singlStart'

  openr,1,dir+file
  readu,1,resultado_aux
  readu,1,parametrizacion_aux
  readu,1,sca_aux,Tmin_aux,Tmax_aux,nc_aux,nw_aux,nband_aux,ndat_aux
 DEMT_A = dblarr(nc_aux,nw_aux,ndat_aux) 
 readu,1,DEMT_A
  close,1
  
  
  
  


i=0
lambda_inicial_A = DEMT_A(*,*,2*nband_aux+2+i:2*nband_aux+2+i+npar_aux-1)
 DEMc1_A         = DEMT_A(*,*,      2*nband_aux  )

;lambda_inicial_A(*,*,*,0) = sqrt(abs( lambda_inicial_A(*,*,*,0)) ) ;
;ya estan modificadas las parametrizaciones
return
end
;-




;+------------------------------------------------------------------

pro read_lambda_inicial_4band
common read_sustrato_gauss3c,lambda_inicial_Array,demc1_Array,demc2_Array

sca_aux  =0
Tmin_aux =0.
Tmax_aux =0.
nr_aux   =0
nth_aux  =0
np_aux   =0
npx_aux  =0
nband_aux=0L
ndat_aux =0L

;if keyword_set(ldem) then 
resultado_aux = 'LDEM'
;if keyword_set( dem) then resultado_aux = 'DEM'


parametrizacion_aux='gauss3c'
   
npar_aux=5

; dir ='/data1/work/dem/'                                                                                                                     
  dir ='./ldem_files/'
  
  file='LDEM.v2_cr2099_l1.0_chianti.ioneq_sun_coronal_1992_feldman_ext.abund_AIA-4band_256_gauss2c_lin_Norm-median_singlStart'
 ;file='LDEM.v2_cr2099_l0.75_chianti.ioneq_sun_coronal_1992_feldman_ext.abund_AIA-4band-newQkl_1075Rsun_256_gauss2c_lin_Norm-median_singlStart'
 ;file='LDEM.v2_cr2099_l1.0_chianti.ioneq_sun_coronal.abund_AIA-4band-oldQkl_1075Rsun_256_gauss2c_lin_Norm-median_singlStart'
 ;file='LDEM.v2_cr2099_l1.0_chianti.ioneq_sun_coronal_1992_feldman_ext.abund_AIA-4band_256_gauss2c_lin_Norm-median_singlStart'
 ;file='LDEM.v2_cr2099_l0.75_chianti.ioneq_sun_coronal.abund_AIA-4band_256_gauss2c_lin_Norm-median_singlStart' ;CR-2099
 ;file='LDEM.v2_cr2106_l0.75_chianti.ioneq_sun_coronal.abund_AIA-4band_256_gauss2c_lin_Norm-median_singlStart' ;CR-2106
 ;file='LDEM.v2_cr2107_l0.75_chianti.ioneq_sun_coronal.abund_AIA-4band_256_gauss2c_lin_Norm-median_singlStart' ;CR-2107 
 
print,'archivo de sustrato:',file

  openr,1,dir+file
  readu,1,resultado_aux
  readu,1,parametrizacion_aux
  readu,1,sca_aux,Tmin_aux,Tmax_aux,nr_aux,nth_aux,np_aux,npx_aux,nband_aux,ndat_aux
 DEMT_A = dblarr(nr_aux ,nth_aux,np_aux,ndat_aux) 
 readu,1,DEMT_A
  close,1
i=0
lambda_inicial_Array = DEMT_A(*,*,*,2*nband_aux+2+i:2*nband_aux+2+i+npar_aux-1)
 DEMc1_Array         = DEMT_A(*,*,*,      2*nband_aux     )
 DEMc2_Array         = DEMT_A(*,*,*,      2*nband_aux + 1 )

return
end

;+------------------------------------------------------------------

pro centroide_gaussiana_condicionada,nband,X,Qkl,Tmin
common gauss_cond,T335,X335,Q335_max,lambda_inicial,DEMc1,DEMc2,DEMc3,T131,X131,Q131_max,SIGMA,X284,T284

if nband eq 4 or nband eq 6 then i=2
if nband eq 5 then begin
   i=3
   Q284=reform(Qkl(nband-2,*))
   Q284_max = max(Q284)
   p        = where(Q284 eq Q284_max)
   X284     = X(p(0))
   T284     = Tmin*x284
endif
Q335                    = reform(Qkl(nband-1,*))
Q211                    = reform(Qkl(nband-i,*))
Q131                    = reform(Qkl(0,*)) ;ojo! si no usamos todavia 94
Q171                    = reform(Qkl(1,*)) ;ojo! si no usamos todavia 94
;FEDE------------------------
if nband eq 6 then begin 
Q131 = reform(Qkl(1,*))
Q094 = reform(Qkl(0,*))
endif
;----------------------------
r = where ( X gt (0.5e6)/tmin )
Q335aux  = Q335 (r)
Q211aux  = Q211 (r)
X335aux  = x (r)
Q335_max = max(Q335aux)
Q211_max = max(Q211aux)
Q171_max = max(Q171)
Q131_max = max(Q131)
p1       = where(Q335aux eq Q335_max)
p2       = where(Q211aux eq Q211_max)
s        = where(Q171 eq Q171_max)
q        = where(Q131 eq Q131_max)
X335     = X335aux(p1(0))
X211     = X335aux(p2(0))
X171     = X(s(0))
X131     = X(q(0))
T335     = Tmin*X335
T211     = Tmin*X211
T171     = Tmin*X171
T131     = Tmin*X131
factor=2.*sqrt(2*alog(2.))
; ancho termico fijo para la segunda gaussiana
sigmaT= (T335-T211)/factor                               ;K
if nband eq 5 then sigmaT= ((T335+T284)*0.5-T211)/factor ;K 
;sigmaT= (T171-T131)/factor  ;K
 sigma = sigmaT/Tmin
; si quiero centroide de gaussiana para 131 en el maximo de
; temperaturas altas
 print,'Temperatura 335 A [MK]:',T335/1.e6
 print,'sigma T 335 A     [MK]:',sigmaT/1.e6
 stop
goto,skiptemp
;FEDE------------------------------------
r = where ( X gt (10^6.5)/tmin )
Xaux = x (r)
Q131aux = q131 (r)
Q094aux = q094 (r)
Q131aux_max = max(Q131aux)
Q094_max = max(Q094aux)
q1        = where(Q131aux eq Q131aux_max)
q2        = where(Q094aux eq Q094_max)
X131     = Xaux(q1(0))
X094     = Xaux(q2(0))
T131     = Tmin*X131
T094     = Tmin*X094
FWHM131  = 6.05e6
FWHM094  = 5.39e6
sigma131 = (X131-X094)/factor;FWHM131/factor/Tmin
sigma094 = (X131-X094)/factor;FWHM094/factor/Tmin

 print,'Temperatura 131 A [MK]:',T131/1.e6
 print,'sigma T 131 A     [MK]:',sigma131*Tmin/1.e6
 print,'Temperatura 094 A [MK]:',T094/1.e6
 print,'sigma T 094 A     [MK]:',sigma094*Tmin/1.e6
;FEDE------------------------------------
skiptemp:
; las lineas de abajo corresponden a gauss3c
;sigma_log=.05
;SIGMA = x131*abs(1.-10.^sigma_log)
;print,'Sigma LogT:',sigma_log
 goto,pedro
   area=1.
      y=area*exp(-0.5d0*((X-X131)/SIGMA131)^2) + 2*area*exp(-0.5d0*((X-X094)/SIGMA094)^2)
      y=area*exp(-0.5d0*((X-X131/2-X094/2)/SIGMA131)^2)
   Q171= reform(Qkl(1,*))
 plot,x*Tmin/1.e6,Q131/Q131aux_max
oplot,x*Tmin/1.e6,Q094/Q094_max
oplot,x*Tmin/1.e6,y,psym=4
stop
pedro:

return
end
;---------------------------------------------------------------------

pro centroide_equalizer,nband,X,Qkl,Tmin
common centroid_equal,xmax,sigma_v,demcv,factor

Xmax =dblarr(nband)
sigma_v=dblarr(nband)
sigma_logt= 0.05

for i=0,nband-1 do begin
p       = where ( reform(qkl(i,*)) eq max(reform(qkl(i,*))) )
xmax (i)= x(p)
;sigma_v(i)= x(i)*abs(1. -10.^sigma_logt)
endfor

factor=2.*sqrt(2*alog(2.))

;------albert-begin--------------------------------------------------
sigma_v(0      )=(xmax(      1)-xmax(      0))/factor
sigma_v(nband-1)=(xmax(nband-1)-xmax(nband-2))/factor
for i=1,nband-2 do sigma_v(i)=mean([xmax(i)-xmax(i-1),xmax(i+1)-xmax(i)])/factor


print,'centroides [MK]: ',Tmin*xmax/1.e6
print,'anchos [MK]: ',Tmin*sigma_v/1.e6

goto,skipgraph
!p.charsize=2
Tx=x*tmin/1.e6
Txmax=xmax*tmin/1.e6
sigtx=sigma_v*tmin/1.e6
                    plot,tx,qkl(0,*)
for i=0,nband-1 do oplot,tx,qkl(i,*)
for i=0,nband-1 do oplot,txmax(i)*[1,1],[0,1]
;for i=0,nband-1 do oplot,(txmax(i)-sigtx(i))*[1,1],[0,1],linestyle=2
;for i=0,nband-1 do oplot,(txmax(i)+sigtx(i))*[1,1],[0,1],linestyle=2
stop
skipgraph:
;------albert-end-------------------------------------------------------

return
end

;----------------------------------------------------------------------
pro read_loss,filename,temp,loss
; para generar el archivo.out con la loss function utilizar el archivo loss_rate.pro
x=''
N_temp=0
t=double(0)
lr=double(0)
openr,2,filename
readf,2,x
readf,2,N_temp
readf,2,x
temp=dblarr(N_temp)
loss=dblarr(N_temp)
for i=0,N_temp-1 do begin
readf,2,t,lr
temp(i)=t
loss(i)=lr
endfor
close,2
return
end
;------------------------------------------------------------------------------
pro interpol_loss,N,temp,loss
tmin=min(temp)
tmax=max(temp)
logTfina=alog10(tmin) + (alog10(tmax) - alog10(tmin))*(dindgen(N)/(N-1))
Tfina=10.^logTfina
lossfina=interpol(loss,temp,Tfina,/spline)
loss=lossfina; el CHIANTI da directamente la Radiative Loss Function
temp=tfina
return
end
;------------------------------------------------------------------------------
function loss_func,T0A  ; erg cm3 sec-1
common paraloss,temp,loss
nt=n_elements(T0A)
lossA=fltarr(nt)
for i=0,nt-1 do begin
 T0=T0A(i)
 dift=abs(Temp-T0)
 ind=where(dift eq min(dift))
 lossA(i)=loss(ind)
endfor
return,lossA
end

;--------------------------------------------------------------------------------
FUNCTION SynthFBE,lambda
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
return,(DEMc*dX*FDEM(X,lambda))##Qkl
end
;--------------------------------------------------------------------------------
FUNCTION PHI,lambda
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
common gauss_cond,T335,X335,Q335_max,lambda_inicial,DEMc1,DEMc2,DEMc3,T131,X131,Q131_max,SIGMA,x284,t284
return,total((FBE-(DEMc*dX*FDEM(X,lambda))##Qkl)^2)
end
;--------------------------------------------------------------------------------

;--------------------------------------------------------------------------------
FUNCTION FDEM,XX,lambda
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
common parametrizacion,parametrization
common gauss_cond,T335,X335,Q335_max,lambda_inicial,DEMc1,DEMc2,DEMc3,T131,X131,Q131_max,SIGMA,x284,t284
common centroid_equal,xmax,sigma_v,demcv,factor

if parametrization eq 'gauss1' then $
   return,(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((XX-lambda[1])/lambda[2])^2)

if parametrization eq 'lorentz' then $
   return,lambda[0]*lambda[2]/(!pi*( (XX-lambda[1])^2 + lambda[2]^2 ))

if parametrization eq 'newpar4' then $
   return,(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*abs((XX-lambda[1])/lambda[2])^(2*lambda[3]))

if parametrization eq 'newpar5' then $
   return,(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*abs((XX-lambda[1])/lambda[2])^(2*lambda[3]))*xx^lambda[4]

if parametrization eq 'gauss2' then begin
   FDEMa=(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((XX-lambda[1])/lambda[2])^2)
   FDEMb=(lambda[3]^2*(DEMc2/DEMc1)/sqrt(2.d0*!pi)/lambda[5])*exp(-0.5d0*((XX-lambda[4])/lambda[5])^2)
   return,FDEMa + FDEMb
endif

if parametrization eq 'gauss2FW' then begin
   FDEMa=(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((XX-lambda[1])/lambda[2])^2)
   FDEMb=(lambda[3]^2*(DEMc2/DEMc1)/sqrt(2.d0*!pi)/sigma)*exp(-0.5d0*((XX-lambda[4])/sigma)^2)
   return,FDEMa + FDEMb
endif
;FEDE----------------------------------------------------------------------------------------------
if parametrization eq 'gauss4FW' then begin
   FDEMa=(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((XX-lambda[1])/lambda[2])^2)
   FDEMb=(lambda[3]^2*(DEMc2/DEMc1)/sqrt(2.d0*!pi)/sigma)*exp(-0.5d0*((XX-lambda[4])/sigma)^2)
   FDEMc=(lambda[5]^2*(DEMc3/DEMc1)/sqrt(2.d0*!pi)/sigma131)*exp(-0.5d0*((XX-X131)/sigma131)^2)
   FDEMd=(lambda[6]^2*(DEMc4/DEMc1)/sqrt(2.d0*!pi)/sigma094)*exp(-0.5d0*((XX-X094)/sigma094)^2)
   return,FDEMa + FDEMb + FDEMc + FDEMd
endif

if parametrization eq 'gauss3FW' then begin
   FDEMa=(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((XX-lambda[1])/lambda[2])^2)
   FDEMb=(lambda[3]^2*(DEMc2/DEMc1)/sqrt(2.d0*!pi)/sigma)*exp(-0.5d0*((XX-lambda[4])/sigma)^2)
   FDEMc=(lambda[5]^2*(DEMc3/DEMc1)/sqrt(2.d0*!pi)/lambda[7])*exp(-0.5d0*((XX-lambda[6])/lambda[7])^2)
   return,FDEMa + FDEMb + FDEMc 
endif

;----------------------------------------------------------------------------------------------

if parametrization eq 'gauss2c' then begin
FDEMa=(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((XX-lambda[1])/lambda[2])^2)
FDEMb=(lambda[3]^2*(DEMc2/DEMc1)/sqrt(2.d0*!pi)/lambda[4])*exp(-0.5d0*((XX-X335     )/lambda[4])^2)
return,FDEMa + FDEMb
endif

if parametrization eq 'gauss2cFW' then begin
FDEMa=(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((XX-lambda[1])/lambda[2])^2)
FDEMb=(lambda[3]^2*(DEMc2/DEMc1)/sqrt(2.d0*!pi)/sigma)*exp(-0.5d0*((XX-X335)/sigma)^2)
return,FDEMa + FDEMb
endif

if parametrization eq 'glc' then begin
FDEMa=(lambda[0]^2/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((XX-lambda[1])/lambda[2])^2)
FDEMb=lambda[3]^2*lambda[4]/(!pi*( (XX-X335)^2 + lambda[4]^2 ))
return,FDEMa + FDEMb
endif

if parametrization eq 'gauss3c' then begin
   FDEMa=(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((XX-lambda[1])/lambda[2])^2)
   FDEMb=(lambda[3]^2*(DEMc2/DEMc1)/sqrt(2.d0*!pi)/lambda[4])*exp(-0.5d0*((XX-X335)/lambda[4])^2)
   FDEMc=(lambda[5]^2*(DEMc3/DEMc1)/sqrt(2.d0*!pi)/SIGMA)*exp(-0.5d0*((XX-X131)/SIGMA)^2)
   return,FDEMa + FDEMb + FDEMc
endif

if parametrization eq 'gauss7par' then begin
   FDEMa=(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((XX-lambda[1])/lambda[2])^2)
   FDEMb=(lambda[3]^2*(DEMc2/DEMc1)/sqrt(2.d0*!pi)/lambda[4])*exp(-0.5d0*((XX-X335)/lambda[4])^2)
   FDEMc=(lambda[5]^2*(DEMc3/DEMc1)/sqrt(2.d0*!pi)/lambda[6])*exp(-0.5d0*((XX-X131)/lambda[6])^2)
   return,FDEMa + FDEMb + FDEMc
endif

if parametrization eq 'equalizer' then begin
   FDEMv = dblarr(ntbins,nband)
   for i=0,nband-1 do FDEMv(*,i)= (lambda[i]^2/sqrt(2.d0*!pi)/sigma_v[i])*exp(-0.5d0*((XX-xmax[i])/sigma_v[i])^2)*(demcv[i]/demc)
   return,total(FDEMv,2)
endif
;-

end
;--------------------------------------------------------------------------------

;--------------------------------------------------------------------------------
FUNCTION gradPHI,lambda
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
common parametrizacion,parametrization
common gauss_cond,T335,X335,Q335_max,lambda_inicial,DEMc1,DEMc2,DEMc3,T131,X131,Q131_max,SIGMA,x284,T284
common centroid_equal,xmax,sigma_v,demcv,factor

if parametrization eq 'gauss1' then begin
 Sk   = (DEMc*dX*FDEM(X,lambda))                                          ##Qkl
 Sk_0 = Sk/lambda[0]
 Sk_1 = (DEMc*dX*FDEM(X,lambda)*  (X-lambda[1])/lambda[2]      /lambda[2])##Qkl
 Sk_2 = (DEMc*dX*FDEM(X,lambda)*(((X-lambda[1])/lambda[2])^2-1)/lambda[2])##Qkl
 return,-2.*[total((FBE-Sk)*Sk_0),total((FBE-Sk)*Sk_1),total((FBE-Sk)*Sk_2)]
endif
;+
if parametrization eq 'lorentz' then begin
 Sk   = (DEMc*dX*FDEM(X,lambda))                                                  ##Qkl
 Sk_0 =  Sk/lambda[0]
 Sk_1 = (DEMc*dX*2*!pi* FDEM(X,lambda)^2 *(X-lambda[1])/lambda[0]^2/lambda[2])    ##Qkl
 Sk_2 = (DEMc*dX*FDEM(X,lambda)*(1-2*!pi*lambda[2]*FDEM(X,lambda)/lambda[0]^2))   ##Qkl
 return,-2.*[total((FBE-Sk)*Sk_0),total((FBE-Sk)*Sk_1),total((FBE-Sk)*Sk_2)]
endif
;-
if parametrization eq 'gauss2' then begin
 FDEMa=(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((X-lambda[1])/lambda[2])^2) 
 FDEMb=(lambda[3]^2*(DEMc2/DEMc1)/sqrt(2.d0*!pi)/lambda[5])*exp(-0.5d0*((X-lambda[4])/lambda[5])^2)
 Sk   = (DEMc*dX* FDEM(X,lambda))                                          ##Qkl
 SkA  = (DEMc*dX*FDEMa)                                                    ##Qkl
 SkB  = (DEMc*dX*FDEMb)                                                    ##Qkl
 Sk_1 = SkA/lambda[0]
 Sk_2 = (DEMc*dX*FDEMa          *  (X-lambda[1])/lambda[2]      /lambda[2])##Qkl
 Sk_3 = (DEMc*dX*FDEMa          *(((X-lambda[1])/lambda[2])^2-1)/lambda[2])##Qkl 
 Sk_4 = 2.*SkB/lambda[3]
 Sk_5 = (DEMc*dX*FDEMb          *  (X-lambda[4])/lambda[5]      /lambda[5])##Qkl
 Sk_6 = (DEMc*dX*FDEMb          *(((X-lambda[4])/lambda[5])^2-1)/lambda[5])##Qkl 
return,-2.*[total((FBE-Sk)*Sk_1),total((FBE-Sk)*Sk_2),total((FBE-Sk)*Sk_3),total((FBE-Sk)*Sk_4),total((FBE-Sk)*Sk_5),total((FBE-Sk)*Sk_6)]
endif

if parametrization eq 'gauss2FW' then begin
 FDEMa=(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((X-lambda[1])/lambda[2])^2) 
 FDEMb=(lambda[3]^2*(DEMc2/DEMc1)/sqrt(2.d0*!pi)/sigma)*exp(-0.5d0*((X-lambda[4])/sigma)^2)
 Sk   = (DEMc*dX* FDEM(X,lambda))                                          ##Qkl
 SkA  = (DEMc*dX*FDEMa)                                                    ##Qkl
 SkB  = (DEMc*dX*FDEMb)                                                    ##Qkl
 Sk_1 = SkA/lambda[0]
 Sk_2 = (DEMc*dX*FDEMa          *  (X-lambda[1])/lambda[2]      /lambda[2])##Qkl
 Sk_3 = (DEMc*dX*FDEMa          *(((X-lambda[1])/lambda[2])^2-1)/lambda[2])##Qkl 
 Sk_4 = 2.*SkB/lambda[3]
 Sk_5 = (DEMc*dX*FDEMb          *  (X-lambda[4])/sigma      /sigma)##Qkl
return,-2.*[total((FBE-Sk)*Sk_1),total((FBE-Sk)*Sk_2),total((FBE-Sk)*Sk_3),total((FBE-Sk)*Sk_4),total((FBE-Sk)*Sk_5)]
endif

;FEDE---------------------------------------------------------------------------------------------
if parametrization eq 'gauss4FW' then begin
 FDEMa=(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((X-lambda[1])/lambda[2])^2) 
 FDEMb=(lambda[3]^2*(DEMc2/DEMc1)/sqrt(2.d0*!pi)/sigma   )*exp(-0.5d0*((X-lambda[4])/sigma)  ^2)
 FDEMc=(lambda[5]^2*(DEMc3/DEMc1)/sqrt(2.d0*!pi)/sigma131)*exp(-0.5d0*((X-X131    )/sigma131)^2)
 FDEMd=(lambda[6]^2*(DEMc4/DEMc1)/sqrt(2.d0*!pi)/sigma094)*exp(-0.5d0*((X-X094    )/sigma094)^2)
 Sk   = (DEMc*dX* FDEM(X,lambda))                                          ##Qkl
 SkA  = (DEMc*dX*FDEMa)                                                    ##Qkl
 SkB  = (DEMc*dX*FDEMb)                                                    ##Qkl
 SkC  = (DEMc*dX*FDEMc)                                                    ##Qkl
 SkD  = (DEMc*dX*FDEMd)                                                    ##Qkl
 Sk_1 = SkA/lambda[0]
 Sk_2 = (DEMc*dX*FDEMa          *  (X-lambda[1])/lambda[2]      /lambda[2])##Qkl
 Sk_3 = (DEMc*dX*FDEMa          *(((X-lambda[1])/lambda[2])^2-1)/lambda[2])##Qkl 
 Sk_4 = 2.*SkB/lambda[3]
 Sk_5 = (DEMc*dX*FDEMb          *  (X-lambda[4])/sigma      /sigma)##Qkl
 Sk_6 = 2.*SkC/lambda[5] 
 Sk_7 = 2.*SkD/lambda[6]
return,-2.*[total((FBE-Sk)*Sk_1),total((FBE-Sk)*Sk_2),total((FBE-Sk)*Sk_3),total((FBE-Sk)*Sk_4),total((FBE-Sk)*Sk_5),$
            total((FBE-Sk)*Sk_6),total((FBE-Sk)*Sk_7)]
endif

if parametrization eq 'gauss3FW' then begin
   FDEMa=(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((XX-lambda[1])/lambda[2])^2)
   FDEMb=(lambda[3]^2*(DEMc2/DEMc1)/sqrt(2.d0*!pi)/sigma)*exp(-0.5d0*((XX-lambda[4])/sigma)^2)
   FDEMc=(lambda[5]^2*(DEMc3/DEMc1)/sqrt(2.d0*!pi)/lambda[7])*exp(-0.5d0*((XX-lambda[6])/lambda[7])^2)
 Sk   = (DEMc*dX* FDEM(X,lambda))                                          ##Qkl
 SkA  = (DEMc*dX*FDEMa)                                                    ##Qkl
 SkB  = (DEMc*dX*FDEMb)                                                    ##Qkl
 SkC  = (DEMc*dX*FDEMc)                                                    ##Qkl
 Sk_1 = SkA/lambda[0]
 Sk_2 = (DEMc*dX*FDEMa          *  (X-lambda[1])/lambda[2]      /lambda[2])##Qkl
 Sk_3 = (DEMc*dX*FDEMa          *(((X-lambda[1])/lambda[2])^2-1)/lambda[2])##Qkl 
 Sk_4 = 2.*SkB/lambda[3]
 Sk_5 = (DEMc*dX*FDEMb          *  (X-lambda[4])/sigma      /sigma)##Qkl
 Sk_6 = SkC/lambda[5]
 Sk_7 = (DEMc*dX*FDEMC          *  (X-lambda[6])/lambda[7]      /lambda[7])##Qkl
 Sk_8 = (DEMc*dX*FDEMC          *(((X-lambda[6])/lambda[7])^2-1)/lambda[7])##Qkl
 
return,-2.*[total((FBE-Sk)*Sk_1),total((FBE-Sk)*Sk_2),total((FBE-Sk)*Sk_3),total((FBE-Sk)*Sk_4),total((FBE-Sk)*Sk_5),$
            total((FBE-Sk)*Sk_6),total((FBE-Sk)*Sk_7),total((FBE-Sk)*Sk_8)]
endif
;-------------------------------------------------------------------------------------------------------

if parametrization eq 'newpar4' then begin
 Sk   = (DEMc*dX*FDEM(X,lambda)                                                                                               )##Qkl
 Sk_0 =  Sk/lambda[0]
 Sk_1 = (DEMc*dX*FDEM(X,lambda)*  lambda[3]*abs((X-lambda[1])/lambda[2]+tiny)^(2*lambda[3]-2)*(X-lambda[1])/lambda[2]^2       )##Qkl
 Sk_2 = (DEMc*dX*FDEM(X,lambda)*( lambda[3]*abs((X-lambda[1])/lambda[2])^(2*lambda[3]) - 1 )/lambda[2]                        )##Qkl
 Sk_3 = (DEMc*dX*FDEM(X,lambda)*(-0.5d0)   *abs((X-lambda[1])/lambda[2])^(2*lambda[3])*alog(((X-lambda[1])/lambda[2])^2+tiny) )##Qkl
 return,-2.*[total((FBE-Sk)*Sk_0),total((FBE-Sk)*Sk_1),total((FBE-Sk)*Sk_2),total((FBE-Sk)*Sk_3)]
endif

if parametrization eq 'newpar5' then begin
 Sk   = (DEMc*dX*FDEM(X,lambda)                                                                                               )##Qkl
 Sk_0 =  Sk/lambda[0]
 Sk_1 = (DEMc*dX*FDEM(X,lambda)*  lambda[3]*abs((X-lambda[1])/lambda[2]+tiny)^(2*lambda[3]-2)*(X-lambda[1])/lambda[2]^2       )##Qkl
 Sk_2 = (DEMc*dX*FDEM(X,lambda)*( lambda[3]*abs((X-lambda[1])/lambda[2])^(2*lambda[3]) - 1 )/lambda[2]                        )##Qkl
 Sk_3 = (DEMc*dX*FDEM(X,lambda)*(-0.5d0)   *abs((X-lambda[1])/lambda[2])^(2*lambda[3])*alog(((X-lambda[1])/lambda[2])^2+tiny) )##Qkl
 Sk_4 = (DEMc*dX*FDEM(X,lambda)*alog(X)                                                                                       )##Qkl
 return,-2.*[total((FBE-Sk)*Sk_0),total((FBE-Sk)*Sk_1),total((FBE-Sk)*Sk_2),total((FBE-Sk)*Sk_3),total((FBE-Sk)*Sk_4)]
endif

;+
if parametrization eq 'gauss2c' then begin
 FDEMa=(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((X-lambda[1])/lambda[2])^2)
 FDEMb=(lambda[3]^2*(DEMc2/DEMc1)/sqrt(2.d0*!pi)/lambda[4])*exp(-0.5d0*((X-X335)/lambda[4])^2)
 Sk   = (DEMc*dX* FDEM(X,lambda))                                          ##Qkl
 SkA  = (DEMc*dX*FDEMa)                                                    ##Qkl
 SkB  = (DEMc*dX*FDEMb)                                                    ##Qkl
 Sk_1 = SkA/lambda[0]
 Sk_2 = (DEMc*dX*FDEMa *  (X-lambda[1])/lambda[2]      /lambda[2])         ##Qkl
 Sk_3 = (DEMc*dX*FDEMa *(((X-lambda[1])/lambda[2])^2-1)/lambda[2])         ##Qkl                                            
 Sk_4 = 2.*SkB/lambda[3]
 Sk_5 = (DEMc*dX*FDEMb*(((X-X335)/lambda[4])^2-1)/lambda[4])               ##Qkl
return,-2.*[total((FBE-Sk)*Sk_1),total((FBE-Sk)*Sk_2),total((FBE-Sk)*Sk_3),total((FBE-Sk)*Sk_4),total((FBE-Sk)*Sk_5)]
endif

if parametrization eq 'gauss2cFW' then begin
 FDEMa=(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((X-lambda[1])/lambda[2])^2) 
 FDEMb=(lambda[3]^2*(DEMc2/DEMc1)/sqrt(2.d0*!pi)/sigma)*exp(-0.5d0*((X-X335)/sigma)^2)
 Sk   = (DEMc*dX* FDEM(X,lambda))                                          ##Qkl
 SkA  = (DEMc*dX*FDEMa)                                                    ##Qkl
 SkB  = (DEMc*dX*FDEMb)                                                    ##Qkl
 Sk_1 = SkA/lambda[0]
 Sk_2 = (DEMc*dX*FDEMa          *  (X-lambda[1])/lambda[2]      /lambda[2])##Qkl
 Sk_3 = (DEMc*dX*FDEMa          *(((X-lambda[1])/lambda[2])^2-1)/lambda[2])##Qkl 
 Sk_4 = 2.*SkB/lambda[3]
return,-2.*[total((FBE-Sk)*Sk_1),total((FBE-Sk)*Sk_2),total((FBE-Sk)*Sk_3),total((FBE-Sk)*Sk_4)]
endif

if parametrization eq 'glc' then begin
 FDEMa=(lambda[0]^2/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((X-lambda[1])/lambda[2])^2)
 FDEMb=lambda[3]^2*lambda[4]/(!pi*( (X-X335)^2 + lambda[4]^2 ))
 Sk   = (DEMc*dX* FDEM(X,lambda))                                          ##Qkl
 SkA  = (DEMc*dX*FDEMa)                                                    ##Qkl
 SkB  = (DEMc*dX*FDEMb)                                                    ##Qkl
 Sk_1 = 2.*SkA/lambda[0]
 Sk_2 = (DEMc*dX*FDEMa *  (X-lambda[1])/lambda[2]      /lambda[2])         ##Qkl
 Sk_3 = (DEMc*dX*FDEMa *(((X-lambda[1])/lambda[2])^2-1)/lambda[2])         ##Qkl                                            
 Sk_4 = 2.*SkB/lambda[3]
 Sk_5 = (DEMc*dX*FDEMb*(1-2*!pi*lambda[4]*FDEMb/lambda[3]^2))              ##Qkl
return,-2.*[total((FBE-Sk)*Sk_1),total((FBE-Sk)*Sk_2),total((FBE-Sk)*Sk_3),total((FBE-Sk)*Sk_4),total((FBE-Sk)*Sk_5)]
endif

if parametrization eq 'gauss3c' then begin
 FDEMa=(lambda[0]/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((X-lambda[1])/lambda[2])^2)
 FDEMb=(lambda[3]^2*(DEMc2/DEMc1)/sqrt(2.d0*!pi)/lambda[4])*exp(-0.5d0*((X-X335)/lambda[4])^2)
 FDEMc=(lambda[5]^2*(DEMc3/DEMc1)/sqrt(2.d0*!pi)/SIGMA)*exp(-0.5d0*((X-X131)/SIGMA)^2)
 Sk   = (DEMc*dX* FDEM(X,lambda))                                          ##Qkl
 SkA  = (DEMc*dX*FDEMa)                                                    ##Qkl
 SkB  = (DEMc*dX*FDEMb)                                                    ##Qkl
 SkC  = (DEMc*dX*FDEMc)                                                    ##Qkl
 Sk_1 = SkA/lambda[0]
 Sk_2 = (DEMc*dX*FDEMa *  (X-lambda[1])/lambda[2]      /lambda[2])         ##Qkl
 Sk_3 = (DEMc*dX*FDEMa *(((X-lambda[1])/lambda[2])^2-1)/lambda[2])         ##Qkl                                            
 Sk_4 = 2.*SkB/lambda[3]
 Sk_5 = (DEMc*dX*FDEMb*(((X-X335)/lambda[4])^2-1)/lambda[4])               ##Qkl
 Sk_6 = 2.*SkC/lambda[5]                                                   
return,-2.*[total((FBE-Sk)*Sk_1),total((FBE-Sk)*Sk_2),total((FBE-Sk)*Sk_3),total((FBE-Sk)*Sk_4),total((FBE-Sk)*Sk_5),total((FBE-Sk)*Sk_6)]
endif

if parametrization eq 'gauss7par' then begin
 FDEMa=(lambda[0]^2/sqrt(2.d0*!pi)/lambda[2])*exp(-0.5d0*((X-lambda[1])/lambda[2])^2)
 FDEMb=(lambda[3]^2*(DEMc2/DEMc1)/sqrt(2.d0*!pi)/lambda[4])*exp(-0.5d0*((X-X335)/lambda[4])^2)
 FDEMc=(lambda[5]^2*(DEMc3/DEMc1)/sqrt(2.d0*!pi)/lambda[6])*exp(-0.5d0*((X-X131)/lambda[6])^2)
 Sk   = (DEMc*dX* FDEM(X,lambda))                                          ##Qkl
 SkA  = (DEMc*dX*FDEMa)                                                    ##Qkl
 SkB  = (DEMc*dX*FDEMb)                                                    ##Qkl
 SkC  = (DEMc*dX*FDEMc)                                                    ##Qkl
 Sk_1 =  SkA/lambda[0]
 Sk_2 = (DEMc*dX*FDEMa *  (X-lambda[1])/lambda[2]      /lambda[2])         ##Qkl
 Sk_3 = (DEMc*dX*FDEMa *(((X-lambda[1])/lambda[2])^2-1)/lambda[2])         ##Qkl                                            
 Sk_4 = 2.*SkB/lambda[3]
 Sk_5 = (DEMc*dX*FDEMb*(((X-X335)/lambda[4])^2-1)/lambda[4])               ##Qkl
 Sk_6 = 2.*SkC/lambda[5] 
 Sk_7 = (DEMc*dX*FDEMb*(((X-X131)/lambda[6])^2-1)/lambda[6])               ##Qkl                                                 
return,$
-2.*[total((FBE-Sk)*Sk_1),total((FBE-Sk)*Sk_2),total((FBE-Sk)*Sk_3),total((FBE-Sk)*Sk_4),total((FBE-Sk)*Sk_5),total((FBE-Sk)*Sk_6),total((FBE-Sk)*Sk_7)]
endif
;-
if parametrization eq 'equalizer' then begin

Sk          = (DEMc*dX* FDEM(X,lambda))                                    ##Qkl
FDEMv       = dblarr(ntbins,nband)
Sk_A        = dblarr(nband,nband)
gradient_phi= dblarr(nband)
for i=0,nband-1 do begin
FDEMv(*,i)     = (lambda[i]^2/sqrt(2.d0*!pi)/sigma_v[i])*exp(-0.5d0*((X-xmax[i])/sigma_v[i])^2)*(demcv[i]/demc)
Sk_A (*,i )    = 2.*(((DEMc*dX* FDEMv(*,i)))                                    ##Qkl)/lambda[i]
gradient_phi(i)= total((FBE-Sk)*Sk_A(*,i)) 
endfor
return,-2.*gradient_phi
endif

end


;--------------------------------------------------------------------------------


;--------------------------------------------------------------------------------
pro singletempsolution,sca,Tmin,FBE,logaT,Qkl,DEMc,Tc,Xc,nband,nTbins,X,dX
 PHI1 = fltarr(nTbins)

 DEM1 = Qkl##FBE/total(Qkl^2,1)
 For il=0,nTbins-1 do PHI1(il) = total((Qkl(*,il)*DEM1(il)-FBE)^2)

 imin = (where(PHI1 eq min(PHI1)))(0)
   Xc = X(imin)
   Tc = Xc*Tmin
 DEMc = DEM1(imin)

goto,skipgraph
window,0,xs=400,ys=400
!p.charsize=1.5
!p.multi=[0,1,2]
 plot,X*Tmin/1.e6,DEM1,/ylog
 plot,X*Tmin/1.e6,PHI1
oplot,[1,1]*Tc/1.e6,[0,max(PHI1)]
!p.multi=0
print,DEMc*dX(imin),Tc/1.e6
print,Xc,Tc,DEMc,PHI1(imin)
print,'FBE_284_synth/FBE_284_tom=',Qkl(0,imin)*DEMc/FBE(0)
print,'FBE_335_synth/FBE_335_tom=',Qkl(1,imin)*DEMc/FBE(1)
stop
skipgraph:
return
end

;--------------------------------------------------------------------------------
pro doubletempsolution,sca,Tmin,FBE,logaT,Qkl,DEMa,DEMb,DEMc,Ta,Tb,Xa,Xb,nband,nTbins,X,dX

PHI1  = fltarr(nTbins,nTbins)
DEMa_A= fltarr(nTbins,nTbins)
DEMb_A= fltarr(nTbins,nTbins)
xa_A  = fltarr(nTbins,nTbins)
xb_A  = fltarr(nTbins,nTbins)

coef_b= fltarr(nTbins,nTbins)
coef_d= total(Qkl^2,1)
coef_I= Qkl##FBE


for ia=0,nTbins-1 do begin
for ib=0,nTbins-1 do begin

coef_b(ia,ib)= total(qkl(*,ia)*qkl(*,ib))

xa_A(ia,ib)=x(ia)

xb_A(ia,ib)=x(ib)

DEMa_A(ia,ib)=(coef_d(ib)*coef_I(ia) - coef_b(ia,ib)*coef_I(ib)) / (coef_d(ia)*coef_d(ib)-coef_b(ia,ib)^2)

DEMb_A(ia,ib)=(coef_d(ia)*coef_I(ib) - coef_b(ia,ib)*coef_I(ia)) / (coef_d(ia)*coef_d(ib)-coef_b(ia,ib)^2)


PHI1  (ia,ib)=total((Qkl(*,ia)*DEMa_A(ia,ib)+Qkl(*,ib)*DEMb_A(ia,ib)-FBE)^2)

imin = (where(PHI1 eq min(PHI1)))(0)

DEMa = DEMa_A(imin)
DEMb = DEMb_A(imin)
DEMc = DEMa + DEMb
xa   = xa_A(imin) 
xb   = xb_A(imin)
Ta   = Tmin*Xa
Tb   = Tmin*Xb

endfor
endfor


return
end




;--------------------------------------------------------------------------------
pro custom_grid,instrument,logaT,Tmin,minX,maxX,Xold,skernel,nband
common customX,Xcustom;,kernelcustom

  dlogaT = logaT(1)-logaT(0)
  dlT    = 0.05
  dI     = fix(dlT/dlogaT)
  nI     = 5
;if instrument eq 'euvi' then begin
  Xcustom=[minX,maxX]

  for ib=0,nband-1 do begin
     q=reform(skernel(ib,*)) & Imq=(where(q eq max(q)))(0)
     Xcustom=[Xcustom,Xold(Imq+(indgen(nI)-(nI-1)/2)*dI)] 
  endfor
;endif

    Xcustom=Xcustom(sort(Xcustom))
    Lcustom=n_elements(Xcustom)
    Xcustom(0)=Xcustom(0)+(Xcustom(1)-Xcustom(0))/2.    
    Xcustom(Lcustom-1)=Xcustom(Lcustom-1)-(Xcustom(Lcustom-1)-Xcustom(Lcustom-2))/2.

    kernelcustom=fltarr(nband,Lcustom)
    for ib=0,nband-1 do $
    kernelcustom(ib,*)=interpol(skernel(ib,*),Xold,Xcustom,/spline)

    window
    plot,Tmin*Xold,skernel(2,*),/ylog,/nodata,xstyle=1,/xlog
    for ib=0,nband-1 do oplot,Tmin*Xold,skernel(ib,*)
    for ib=0,nband-1 do oplot,Tmin*Xcustom,kernelcustom(ib,*),psym=4,th=4
stop

return
end

;--------------------------------------------------------------------------------
pro Interpol_Qkl,logT,kernel,nband,nTbins,sca,Tmin,Tmax,Qkl,X,dX,L,logaT,instrument
common customX,Xcustom

;Select range of interest
 range=where( (logT ge alog10(Tmin)) AND (logT le alog10(Tmax)) )
 range=[min(range)-1,range,max(range)+1]

;Select the temperature range to be used:
;and correct Tmin and Tmax
 logaT  = [logT(   range)]
skernel = kernel(*,range)
;Tmin   = 10.^min(logaT) 
;Tmax   = 10.^max(logaT) 

;Make Xold = Told / Tmin
 Xold = 10.^logaT/Tmin
;minX=min(Xold)
;maxX=max(Xold)
 minX=Tmin/Tmin
 maxX=Tmax/Tmin

;Decide NEW grid
;make  X = T/Tmin
;make dX the array of bin widths
 if sca eq 0 then begin 
 DlogaT_L = (max(logaT)-min(logaT))  / float(nTbins)
  logaT_L = min(logaT) + DlogaT_L/2. + DlogaT_L * findgen(nTbins)
        X = 10.^logaT_L /Tmin
 logaT_La = min(logaT) +               DlogaT_L * findgen(nTbins)
 logaT_Lb = min(logaT) +  DlogaT_L   + DlogaT_L * findgen(nTbins)
       Xa = 10.^logaT_La/Tmin
       Xb = 10.^logaT_Lb/Tmin
       dX = Xb-Xa
    endif

  if sca eq 1 then begin
     DT_L = (Tmax-Tmin) / float(nTbins)
    
  T_L = Tmin + DT_L /2. + DT_L *findgen(nTbins)
        X = T_L/Tmin
       dX = (maxX-minX)/float(nTbins) + fltarr(nTbins)
  endif
  if sca eq 2 then begin
      custom_grid,instrument,logaT,Tmin,minX,maxX,Xold,skernel,nband
  endif

;Interpolate logaT into the X-grid
 logaT = interpol(logaT,Xold,X)

;Make the Qkl array, interpolated into the X-grid
 Qkl=dblarr(nband,nTbins)
 for ib=0,nband-1 do Qkl(ib,*) = interpol(skernel(ib,*),Xold,X);,/spline)

;Check interpolation:
goto,skip
    window

    if sca eq 0 then begin
    plot,logT,kernel(2,*),/ylog,/nodata,xstyle=1
    for ib=0,nband-1 do oplot,logT          ,kernel(ib,*)
    for ib=0,nband-1 do oplot,alog10(X*Tmin),   Qkl(ib,*),psym=4,th=4
    oplot,[1,1]*alog10(Tmin),[0,1]
    for il=0,ntbins-1 do $
    oplot,[1,1]*logaT_La(il),[.00001,1.]*max(kernel),th=0.5
    for il=0,ntbins-1 do $
    oplot,[1,1]*logaT_Lb(il),[.00001,1.]*max(kernel),linestyle=2,th=3
    endif

    if sca eq 1 then begin
    plot,10.^logT/1.e6,kernel(2,*),/ylog,xr=[0.5,5]
    for ib=0,nband-1 do oplot,10.^logT/1.e6   ,kernel(ib,*)
    for ib=0,nband-1 do oplot,X*Tmin/1.e6,        Qkl(ib,*),psym=4,th=4
    endif

    stop

skip:
return
end

;-------------------------------------------------
pro start_point_newpar4,new_guess,Xc,nn,new_guess_a,npar,multiple=multiple
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
common parametrizacion,parametrization

n1= 5
n2= 4
n3=11
n4=11

l1min =  0.2d
l1max = 10.0d
l2min = min(x)
l2max = max(x)
l3min = 3.*min(dX) ; 0.05e6/Tmin    ; 0.05MK
l3max = (max(x)-min(x))/4 ; All range
l4min = .5
l4max = 3

 l1= l1min + (l1max-l1min) *findgen(n1)/float(n1-1)
;l2= l2min + (l2max-l2min) *findgen(n2)/float(n2-1)
;l2= [l2min , Xc-(Xc-l2min)/2. , Xc , Xc+(l2max-Xc)/2. , l2max]
 l2 = [1.72000   ,   3.00000   ,   3.74000   ,   4.94000]
 l3= l3min + (l3max-l3min) *findgen(n3)/float(n3-1)
 l4= l4min + (l4max-l4min) *findgen(n4)/float(n4-1)

phiA = dblarr(n1,n2,n3,n4)
 l1A = dblarr(n1,n2,n3,n4)
 l2A = dblarr(n1,n2,n3,n4)
 l3A = dblarr(n1,n2,n3,n4)
 l4A = dblarr(n1,n2,n3,n4)

for i1=0,n1-1 do begin
for i2=0,n2-1 do begin
for i3=0,n3-1 do begin
for i4=0,n4-1 do begin
phiA(i1,i2,i3,i4)=PHI([l1(i1),l2(i2),l3(i3),l4(i4)])
 l1A(i1,i2,i3,i4)=l1(i1)
 l2A(i1,i2,i3,i4)=l2(i2)
 l3A(i1,i2,i3,i4)=l3(i3)
 l4A(i1,i2,i3,i4)=l4(i4)
endfor
endfor
endfor
endfor

ii= median( where(phiA eq min(phiA)) )
new_guess=[l1A(ii),l2A(ii),l3A(ii),l4A(ii)]

nn=1
new_guess_a=dblarr(nn,npar)
new_guess_a(0,*)=new_guess

return
end
;-------------------------------------------------

pro start_point_newpar5,new_guess,Xc,nn,new_guess_a,npar,multiple=multiple
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
common parametrizacion,parametrization
common intent,intento

if new_guess(0) ne -666. AND keyword_set(multiple) then goto,skip_min5D

if intento eq 0 then begin
n1= 5
n2= 4
n3= 11
n4= 11
n5= 11

l1min =  0.2d
l1max = 10.0d
l2min = min(x)
l2max = max(x)
l3min = 3.*min(dX) ; 0.05e6/Tmin    ; 0.05MK
l3max = (max(x)-min(x))/2 ; All range
l4min = .5
l4max = 3
l5min =-2
l5max = 2 
endif

if intento NE 0 then begin
n1= 5
n2= 4
n3= 7
n4= 21
n5= 11
l1min =  0.2d
l1max = 10.0d
l2min = min(x)
l2max = max(x)
l3min = 3.*min(dX) ; 0.05e6/Tmin    ; 0.05MK
l3max = (max(x)-min(x))/4 ; All range
l4min = .5
l4max = 3
l5min =-2
l5max = 2 
endif

 l1= l1min + (l1max-l1min) *findgen(n1)/float(n1-1)
;l2= l2min + (l2max-l2min) *findgen(n2)/float(n2-1)
;l2= [l2min , Xc-(Xc-l2min)/2. , Xc , Xc+(l2max-Xc)/2. , l2max]
 l2 = [1.72000   ,   3.00000   ,   3.74000   ,   4.94000]
 l3= l3min + (l3max-l3min) *findgen(n3)/float(n3-1)
 l4= l4min + (l4max-l4min) *findgen(n4)/float(n4-1)
 l5= l5min + (l5max-l5min) *findgen(n5)/float(n5-1)

phiA = dblarr(n1,n2,n3,n4,n5)
 l1A = dblarr(n1,n2,n3,n4,n5)
 l2A = dblarr(n1,n2,n3,n4,n5)
 l3A = dblarr(n1,n2,n3,n4,n5)
 l4A = dblarr(n1,n2,n3,n4,n5)
 l5A = dblarr(n1,n2,n3,n4,n5)

for i1=0,n1-1 do begin
for i2=0,n2-1 do begin
for i3=0,n3-1 do begin
for i4=0,n4-1 do begin
for i5=0,n5-1 do begin

phiA(i1,i2,i3,i4,i5)=PHI([l1(i1),l2(i2),l3(i3),l4(i4),l5(i5)])
 l1A(i1,i2,i3,i4,i5)=l1(i1)
 l2A(i1,i2,i3,i4,i5)=l2(i2)
 l3A(i1,i2,i3,i4,i5)=l3(i3)
 l4A(i1,i2,i3,i4,i5)=l4(i4)
 l5A(i1,i2,i3,i4,i5)=l5(i5)

endfor
endfor
endfor
endfor
endfor

ii= median( where(phiA eq min(phiA)) )
new_guess=[l1A(ii),l2A(ii),l3A(ii),l4A(ii),l5A(ii)]

nn=1
new_guess_a=dblarr(nn,npar)
new_guess_a(0,*)=new_guess

skip_min5D:

if keyword_set(multiple) then begin
print,'Explore from multiple starting points'
nn=9
new_guess_a=dblarr(nn,5)

new_guess_a(1,*)=[1,Xc,min(l3),min(l4),min(l5)]
new_guess_a(2,*)=[1,Xc,min(l3),min(l4),max(l5)]
new_guess_a(3,*)=[1,Xc,min(l3),max(l4),min(l5)]
new_guess_a(4,*)=[1,Xc,min(l3),max(l4),max(l5)]

new_guess_a(5,*)=[1,Xc,max(l3),min(l4),min(l5)]
new_guess_a(6,*)=[1,Xc,max(l3),min(l4),max(l5)]
new_guess_a(7,*)=[1,Xc,max(l3),max(l4),min(l5)]
new_guess_a(8,*)=[1,Xc,max(l3),max(l4),max(l5)]
endif

return
end

;------------------------------------------------------------------------------------------------
pro start_point,new_guess,Xc,nn,new_guess_a,npar,multiple=multiple
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
common parametrizacion,parametrization

;if new_guess(0) ne -666. AND keyword_set(multiple) then goto,skip_min3D

; Changed new guess!
 new_guess=[1.,Xc ,2. ]
;new_guess=[1.,1.5,1.5]

;goto,skip_min3D

n1= 10
n2= 7
n3= 10


phiA = dblarr(n1,n2,n3)


hmin =  0.1d 
hmax = 10.0d 
cmin = min(x)
cmax = max(x)
dmin = 3.*min(dX); 0.05e6/Tmin        ; 0.05MK
dmax = (max(x)-min(x))/2. ; All range

l1= hmin + (hmax-hmin) *findgen(n1)/float(n1-1)
;l2= cmin + (cmax-cmin) *findgen(n2)/float(n2-1)
l2=[cmin,xc - (xc-cmin)*2./3, xc - (xc-cmin)*1./3,xc,xc + (cmax -xc)*1./3,xc + (cmax -xc)*2./3 ,cmax]
l3= dmin + (dmax-dmin) *findgen(n3)/float(n3-1)


l1A = dblarr(n1,n2,n3)
l2A = dblarr(n1,n2,n3)
l3A = dblarr(n1,n2,n3)

for i1=0,n1-1 do begin
for i2=0,n2-1 do begin
for i3=0,n3-1 do begin
phiA(i1,i2,i3)=PHI([l1(i1),l2(i2),l3(i3)])
 l1A(i1,i2,i3)=l1(i1)
 l2A(i1,i2,i3)=l2(i2)
 l3A(i1,i2,i3)=l3(i3)
endfor
endfor
endfor

ii= median( where(phiA eq min(phiA)) )
new_guess=[l1A(ii),l2A(ii),l3A(ii)]

f1=abs(l1-new_guess(0)) & ii1=(where(f1 eq min(f1)))(0)
f2=abs(l2-new_guess(1)) & ii2=(where(f2 eq min(f2)))(0)
f3=abs(l3-new_guess(2)) & ii3=(where(f3 eq min(f3)))(0)

nn=1
new_guess_a=dblarr(nn,3)
new_guess_a(0,*)=new_guess

skip_min3D:

if keyword_set(multiple) then begin
print,'Explore from multiple starting points'
nn=9
new_guess_a=dblarr(nn,3)

new_guess_a(1,*)=[min(l1),min(l2),min(l3)]
new_guess_a(2,*)=[min(l1),min(l2),max(l3)]
new_guess_a(3,*)=[min(l1),max(l2),min(l3)]
new_guess_a(4,*)=[min(l1),max(l2),max(l3)]

new_guess_a(5,*)=[max(l1),min(l2),min(l3)]
new_guess_a(6,*)=[max(l1),min(l2),max(l3)]
new_guess_a(7,*)=[max(l1),max(l2),min(l3)]
new_guess_a(8,*)=[max(l1),max(l2),max(l3)]
endif

goto,skipprints
print,'indexes for minimum:',ii1,ii2,ii3
print,'lambda:',new_guess
stop
 for i=0,n1-1 do print,min(phiA(i,*,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,i,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,*,i))/min(phia)
 print,'-'
 stop
skipprints:

return
end
;-------------------------------------------------------------------

pro start_point_eq,new_guess,Xc,nn,new_guess_a,npar,multiple=multiple
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
common parametrizacion,parametrization

n=5
hmin =  0.1d^.5
hmax =  10.d^.5 
l    = hmin + (hmax-hmin) *findgen(n)/float(n-1)

if nband eq 3 then phiA = dblarr(n,n,n)
if nband eq 4 then phiA = dblarr(n,n,n,n)
if nband eq 5 then phiA = dblarr(n,n,n,n,n)
if nband eq 6 then phiA = dblarr(n,n,n,n,n,n)

    l1A = phiA
    l2A = phiA
    l3A = phiA
if nband ge 4 then begin
    l4A = phiA
    l5A = phiA
    l6A = phiA
endif

for i1=0,n-1 do begin
for i2=0,n-1 do begin
for i3=0,n-1 do begin
if nband eq 3 then begin
phiA(i1,i2,i3)=PHI([l(i1),l(i2),l(i3)])
 l1A(i1,i2,i3)=l(i1)
 l2A(i1,i2,i3)=l(i2)
 l3A(i1,i2,i3)=l(i3)
 goto,next3 
endif
for i4=0,n-1 do begin
if nband eq 4 then begin
phiA(i1,i2,i3,i4)=PHI([l(i1),l(i2),l(i3),l(i4)])
 l1A(i1,i2,i3,i4)=l(i1)
 l2A(i1,i2,i3,i4)=l(i2)
 l3A(i1,i2,i3,i4)=l(i3)
 l4A(i1,i2,i3,i4)=l(i4)
 goto,next4
endif
for i5=0,n-1 do begin
if nband eq 5 then begin
phiA(i1,i2,i3,i4,i5)=PHI([l(i1),l(i2),l(i3),l(i4),l(i5)])
 l1A(i1,i2,i3,i4,i5)=l(i1)
 l2A(i1,i2,i3,i4,i5)=l(i2)
 l3A(i1,i2,i3,i4,i5)=l(i3)
 l4A(i1,i2,i3,i4,i5)=l(i4)
 l5A(i1,i2,i3,i4,i5)=l(i5)
 goto,next4
endif
for i6=0,n-1 do begin
phiA(i1,i2,i3,i4,i5,i6)=PHI([l(i1),l(i2),l(i3),l(i4),l(i5),l(i6)])
 l1A(i1,i2,i3,i4,i5,i6)=l(i1)
 l2A(i1,i2,i3,i4,i5,i6)=l(i2)
 l3A(i1,i2,i3,i4,i5,i6)=l(i3)
 l4A(i1,i2,i3,i4,i5,i6)=l(i4)
 l5A(i1,i2,i3,i4,i5,i6)=l(i5)
 l6A(i1,i2,i3,i4,i5,i6)=l(i6)
endfor
next5:
endfor
next4:
endfor
next3:
endfor
endfor
endfor

ii= median( where(phiA eq min(phiA)) )

if nband eq 3 then new_guess=[l1A(ii),l2A(ii),l3A(ii)]
if nband eq 4 then new_guess=[l1A(ii),l2A(ii),l3A(ii),l4A(ii)]
if nband eq 5 then new_guess=[l1A(ii),l2A(ii),l3A(ii),l4A(ii),l5A(ii)]
if nband eq 6 then new_guess=[l1A(ii),l2A(ii),l3A(ii),l4A(ii),l5A(ii),l6A(ii)]

nn=1
new_guess_a=dblarr(nn,nband)
new_guess_a(0,*)=new_guess

skip_min3D:

if keyword_set(multiple) then begin
print,'Explore from multiple starting points'
nn=9
new_guess_a=dblarr(nn,3)

new_guess_a(1,*)=[min(l1),min(l2),min(l3)]
new_guess_a(2,*)=[min(l1),min(l2),max(l3)]
new_guess_a(3,*)=[min(l1),max(l2),min(l3)]
new_guess_a(4,*)=[min(l1),max(l2),max(l3)]

new_guess_a(5,*)=[max(l1),min(l2),min(l3)]
new_guess_a(6,*)=[max(l1),min(l2),max(l3)]
new_guess_a(7,*)=[max(l1),max(l2),min(l3)]
new_guess_a(8,*)=[max(l1),max(l2),max(l3)]
endif

goto,skipprints
print,'indexes for minimum:',ii1,ii2,ii3
print,'lambda:',new_guess
stop
 for i=0,n1-1 do print,min(phiA(i,*,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,i,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,*,i))/min(phia)
 print,'-'
 stop
skipprints:

return
end
;-------------------------------------------------

pro start_point_lorentz,new_guess,Xc,nn,new_guess_a,npar,multiple=multiple
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
common parametrizacion,parametrization


n1= 15
n2=  5
n3= 15

phiA = dblarr(n1,n2,n3)

hmin =  0.2d^.5 ;.3
hmax =  5.0d^.5 ;.5 
cmin = min(x)
cmax = max(x)
dmin = 3.*min(dX); 0.05e6/Tmin        ; 0.05MK
dmax = (max(x)-min(x))/2. ; All range

l1= hmin + (hmax-hmin) *findgen(n1)/float(n1-1)
l2= cmin + (cmax-cmin) *findgen(n2)/float(n2-1)
l3= dmin + (dmax-dmin) *findgen(n3)/float(n3-1)

l1A = dblarr(n1,n2,n3)
l2A = dblarr(n1,n2,n3)
l3A = dblarr(n1,n2,n3)

for i1=0,n1-1 do begin
for i2=0,n2-1 do begin
for i3=0,n3-1 do begin

phiA(i1,i2,i3)=PHI([l1(i1),l2(i2),l3(i3)])
 l1A(i1,i2,i3)=l1(i1)
 l2A(i1,i2,i3)=l2(i2)
 l3A(i1,i2,i3)=l3(i3)
endfor
endfor
endfor

ii= median( where(phiA eq min(phiA)) )
new_guess=[l1A(ii),l2A(ii),l3A(ii)]

f1=abs(l1-new_guess(0)) & ii1=(where(f1 eq min(f1)))(0)
f2=abs(l2-new_guess(1)) & ii2=(where(f2 eq min(f2)))(0)
f3=abs(l3-new_guess(2)) & ii3=(where(f3 eq min(f3)))(0)

nn=1
new_guess_a=dblarr(nn,3)
new_guess_a(0,*)=new_guess

skip_min3D:

if keyword_set(multiple) then begin
print,'Explore from multiple starting points'
nn=9
new_guess_a=dblarr(nn,3)

new_guess_a(1,*)=[min(l1),min(l2),min(l3)]
new_guess_a(2,*)=[min(l1),min(l2),max(l3)]
new_guess_a(3,*)=[min(l1),max(l2),min(l3)]
new_guess_a(4,*)=[min(l1),max(l2),max(l3)]

new_guess_a(5,*)=[max(l1),min(l2),min(l3)]
new_guess_a(6,*)=[max(l1),min(l2),max(l3)]
new_guess_a(7,*)=[max(l1),max(l2),min(l3)]
new_guess_a(8,*)=[max(l1),max(l2),max(l3)]
endif

goto,skipprints
print,'indexes for minimum:',ii1,ii2,ii3
print,'lambda:',new_guess
stop
 for i=0,n1-1 do print,min(phiA(i,*,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,i,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,*,i))/min(phia)
 print,'-'
 stop
skipprints:

return
end

;-------------------------------------------------


pro start_point_2,new_guess,Xc,nn,new_guess_a,multiple=multiple
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
common parametrizacion,parametrization

;goto,skip_min3D

na1= 5
na2= 5
na3= 5

nb1= 5
nb2= 5
nb3= 5

phiA = dblarr(na1,na2,na3,nb1,nb2,nb3)

hmin=0.1d
hmax=10.0d
cmin=min(x) / 3.
cmax=max(x) * 2.
dmin=   0.5*dx
dmax=nTbins*dx * 2.

la1= hmin + (hmax-hmin) *findgen(na1)/float(na1-1)
la2= cmin + (cmax-cmin) *findgen(na2)/float(na2-1)
la3= dmin + (dmax-dmin) *findgen(na3)/float(na3-1)

lb1=la1
lb2=la2
lb3=la3

la1A = dblarr(na1,na2,na3,nb1,nb2,nb3)
la2A = dblarr(na1,na2,na3,nb1,nb2,nb3)
la3A = dblarr(na1,na2,na3,nb1,nb2,nb3)

lb1A = dblarr(na1,na2,na3,nb1,nb2,nb3)
lb2A = dblarr(na1,na2,na3,nb1,nb2,nb3)
lb3A = dblarr(na1,na2,na3,nb1,nb2,nb3)

for ia1=0,na1-1 do begin
for ia2=0,na2-1 do begin
for ia3=0,na3-1 do begin
for ib1=0,nb1-1 do begin
for ib2=0,nb2-1 do begin
for ib3=0,nb3-1 do begin
phiA(ia1,ia2,ia3,ib1,ib2,ib3)=PHI([la1(ia1),la2(ia2),la3(ia3),lb1(ib1),lb2(ib2),lb3(ib3)])
 la1A(ia1,ia2,ia3,ib1,ib2,ib3)=la1(ia1)
 la2A(ia1,ia2,ia3,ib1,ib2,ib3)=la2(ia2)
 la3A(ia1,ia2,ia3,ib1,ib2,ib3)=la3(ia3)
 lb1A(ia1,ia2,ia3,ib1,ib2,ib3)=lb1(ib1)
 lb2A(ia1,ia2,ia3,ib1,ib2,ib3)=lb2(ib2)
 lb3A(ia1,ia2,ia3,ib1,ib2,ib3)=lb3(ib3)
endfor
endfor
endfor
endfor
endfor
endfor

ii=median( where(phiA eq min(phiA)) )
new_guess=[la1A(ii),la2A(ii),la3A(ii),lb1A(ii),lb2A(ii),lb3A(ii)]
stop

;f1=abs(l1-new_guess(0)) & ii1=(where(f1 eq min(f1)))(0)
;f2=abs(l2-new_guess(1)) & ii2=(where(f2 eq min(f2)))(0)
;f3=abs(l3-new_guess(2)) & ii3=(where(f3 eq min(f3)))(0)

nn=1
new_guess_a=dblarr(nn,6)
new_guess_a(0,*)=new_guess

;print,'My Guess:',reform([new_guess,min(phiA)])

if keyword_set(multiple) then begin
print,'Explore from multiple starting points'
nn=9
new_guess_a=dblarr(nn,3)
new_guess_a(0,*)=new_guess

new_guess_a(1,*)=[min(l1),min(l2),min(l3)]
new_guess_a(2,*)=[min(l1),min(l2),max(l3)]
new_guess_a(3,*)=[min(l1),max(l2),min(l3)]
new_guess_a(4,*)=[min(l1),max(l2),max(l3)]

new_guess_a(5,*)=[max(l1),min(l2),min(l3)]
new_guess_a(6,*)=[max(l1),min(l2),max(l3)]
new_guess_a(7,*)=[max(l1),max(l2),min(l3)]
new_guess_a(8,*)=[max(l1),max(l2),max(l3)]
endif

goto,skipprints
print,'indexes for minimum:',ii1,ii2,ii3
 for i=0,n1-1 do print,min(phiA(i,*,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,i,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,*,i))/min(phia)
 print,'-'
skipprints:

return
end


;+---------------------------------------------------------------
pro start_point_2c,new_guess,Xc,nn,new_guess_a,npar,multiple=multiple
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
common parametrizacion,parametrization
common gauss_cond,T335,X335,Q335_max,lambda_inicial,DEMc1,DEMc2,DEMc3,T131,X131,Q131_max,SIGMA,X284,T284

if new_guess(0) ne -666. AND keyword_set(multiple) then goto,skip_min3D


n1= 1
n2= 1
n3= 1
n4= 11
;n5= 11
n5= 1
phiA = dblarr(n1,n2,n3,n4,n5)

hmin =  0.10d^.5
hmax =  10.0d^.5
dmin =            3*dx(0)        
dmax = (max(x)-min(x))/2. ; All range


l1= [lambda_inicial(0)]
l2= [lambda_inicial(1)]
l3= [lambda_inicial(2)]
l4= hmin + (hmax-hmin) *findgen(n4)/float(n4-1)
;l5= dmin + (dmax-dmin) *findgen(n5)/float(n5-1)
 l5= [0.25*1.e6/Tmin] 

l1A = dblarr(n1,n2,n3,n4,n5)
l2A = dblarr(n1,n2,n3,n4,n5)
l3A = dblarr(n1,n2,n3,n4,n5)
l4A = dblarr(n1,n2,n3,n4,n5)
l5A = dblarr(n1,n2,n3,n4,n5)


for i1=0,n1-1 do begin
for i2=0,n2-1 do begin
for i3=0,n3-1 do begin
for i4=0,n4-1 do begin
for i5=0,n5-1 do begin
phiA(i1,i2,i3,i4,i5)=PHI([l1(i1),l2(i2),l3(i3),l4(i4),l5(i5)])
 l1A(i1,i2,i3,i4,i5)=l1(i1)
 l2A(i1,i2,i3,i4,i5)=l2(i2)
 l3A(i1,i2,i3,i4,i5)=l3(i3)
 l4A(i1,i2,i3,i4,i5)=l4(i4)
 l5A(i1,i2,i3,i4,i5)=l5(i5)
endfor
endfor
endfor
endfor
endfor

ii= median( where(phiA eq min(phiA)) )
new_guess=[l1A(ii),l2A(ii),l3A(ii),l4A(ii),l5A(ii)]

nn=1
new_guess_a=dblarr(nn,5)
new_guess_a(0,*)=new_guess

skip_min3D:

if keyword_set(multiple) then begin
print,'Explore from multiple starting points'
nn=9
new_guess_a=dblarr(nn,3)

new_guess_a(1,*)=[min(l1),min(l2),min(l3)]
new_guess_a(2,*)=[min(l1),min(l2),max(l3)]
new_guess_a(3,*)=[min(l1),max(l2),min(l3)]
new_guess_a(4,*)=[min(l1),max(l2),max(l3)]

new_guess_a(5,*)=[max(l1),min(l2),min(l3)]
new_guess_a(6,*)=[max(l1),min(l2),max(l3)]
new_guess_a(7,*)=[max(l1),max(l2),min(l3)]
new_guess_a(8,*)=[max(l1),max(l2),max(l3)]
endif

goto,skipprints
print,'indexes for minimum:',ii1,ii2,ii3
print,'lambda:',new_guess
stop
 for i=0,n1-1 do print,min(phiA(i,*,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,i,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,*,i))/min(phia)
 print,'-'
 stop
skipprints:

return
end
;+---------------------------------------------------------------
pro start_point_double,new_guess,Xc,nn,new_guess_a,npar,multiple=multiple
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
common parametrizacion,parametrization
common gauss_cond,T335,X335,Q335_max,lambda_inicial,DEMc1,DEMc2,DEMc3,T131,X131,Q131_max,SIGMA

if new_guess(0) ne -666. AND keyword_set(multiple) then goto,skip_min3D


n1= 1
n2= 1
n3= 1
n4= 11
n5= 1
n6= 11
;n6=1

phiA = dblarr(n1,n2,n3,n4,n5,n6)

hmin =  0.10d^.5
hmax =  10.0d^.5
dmin =            3*dx(0)        
dmax = (max(x)-min(x))/2. ; All range
;experimento 1
dmax = 0.5*1.e6/Tmin
;experimento 2
dmin= 0.2*1.e6/Tmin
;experimento 3
dmin=0.25*1.e6/Tmin

l1= [lambda_inicial(0)]
l2= [lambda_inicial(1)]
l3= [lambda_inicial(2)]
l4= hmin + (hmax-hmin) *findgen(n4)/float(n4-1)
l5= [X335]
l6= dmin + (dmax-dmin) *findgen(n6)/float(n6-1)
;l6= [0.25*1.e6/Tmin]

l1A = dblarr(n1,n2,n3,n4,n5,n6)
l2A = dblarr(n1,n2,n3,n4,n5,n6)
l3A = dblarr(n1,n2,n3,n4,n5,n6)
l4A = dblarr(n1,n2,n3,n4,n5,n6)
l5A = dblarr(n1,n2,n3,n4,n5,n6)
l6A = dblarr(n1,n2,n3,n4,n5,n6)

for i1=0,n1-1 do begin
for i2=0,n2-1 do begin
for i3=0,n3-1 do begin
for i4=0,n4-1 do begin
for i5=0,n5-1 do begin
for i6=0,n6-1 do begin
phiA(i1,i2,i3,i4,i5,i6)=PHI([l1(i1),l2(i2),l3(i3),l4(i4),l5(i5),l6(i6)])
 l1A(i1,i2,i3,i4,i5,i6)=l1(i1)
 l2A(i1,i2,i3,i4,i5,i6)=l2(i2)
 l3A(i1,i2,i3,i4,i5,i6)=l3(i3)
 l4A(i1,i2,i3,i4,i5,i6)=l4(i4)
 l5A(i1,i2,i3,i4,i5,i6)=l5(i5)
 l6A(i1,i2,i3,i4,i5,i6)=l6(i6)
endfor
endfor
endfor
endfor
endfor
endfor

ii= median( where(phiA eq min(phiA)) )
new_guess=[l1A(ii),l2A(ii),l3A(ii),l4A(ii),l5A(ii),l6A(ii)]


nn=1
new_guess_a=dblarr(nn,6)
new_guess_a(0,*)=new_guess

skip_min3D:

if keyword_set(multiple) then begin
print,'Explore from multiple starting points'
nn=9
new_guess_a=dblarr(nn,3)

new_guess_a(1,*)=[min(l1),min(l2),min(l3)]
new_guess_a(2,*)=[min(l1),min(l2),max(l3)]
new_guess_a(3,*)=[min(l1),max(l2),min(l3)]
new_guess_a(4,*)=[min(l1),max(l2),max(l3)]

new_guess_a(5,*)=[max(l1),min(l2),min(l3)]
new_guess_a(6,*)=[max(l1),min(l2),max(l3)]
new_guess_a(7,*)=[max(l1),max(l2),min(l3)]
new_guess_a(8,*)=[max(l1),max(l2),max(l3)]
endif

goto,skipprints
print,'indexes for minimum:',ii1,ii2,ii3
print,'lambda:',new_guess
stop
 for i=0,n1-1 do print,min(phiA(i,*,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,i,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,*,i))/min(phia)
 print,'-'
 stop
skipprints:

return
end

pro start_point_dgfw,new_guess,Xc,nn,new_guess_a,npar,multiple=multiple
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
common parametrizacion,parametrization
common gauss_cond,T335,X335,Q335_max,lambda_inicial,DEMc1,DEMc2,DEMc3,T131,X131,Q131_max,SIGMA,X284,T284

if new_guess(0) ne -666. AND keyword_set(multiple) then goto,skip_min3D


n1= 1
n2= 1
n3= 1
n4= 11
n5= 1
if nband eq 5 then n5=11

phiA = dblarr(n1,n2,n3,n4,n5)

hmin =  0.10d^.5
hmax =  10.0d^.5
xcmin= 2.e6/Tmin
xcmax= 3.e6/Tmin
dmin = 3*dx(0)        
dmax = (max(x)-min(x))/2. ; All range


l1= [lambda_inicial(0)]
l2= [lambda_inicial(1)]
l3= [lambda_inicial(2)]
l4= hmin + (hmax-hmin) *findgen(n4)/float(n4-1)
if nband eq 4 then l5= [X335]
if nband eq 5 then l5= xcmin + (xcmax-xcmin) *findgen(n5)/float(n5-1)

l1A = dblarr(n1,n2,n3,n4,n5)
l2A = dblarr(n1,n2,n3,n4,n5)
l3A = dblarr(n1,n2,n3,n4,n5)
l4A = dblarr(n1,n2,n3,n4,n5)
l5A = dblarr(n1,n2,n3,n4,n5)

for i1=0,n1-1 do begin
for i2=0,n2-1 do begin
for i3=0,n3-1 do begin
for i4=0,n4-1 do begin
for i5=0,n5-1 do begin
phiA(i1,i2,i3,i4,i5)=PHI([l1(i1),l2(i2),l3(i3),l4(i4),l5(i5)])
 l1A(i1,i2,i3,i4,i5)=l1(i1)
 l2A(i1,i2,i3,i4,i5)=l2(i2)
 l3A(i1,i2,i3,i4,i5)=l3(i3)
 l4A(i1,i2,i3,i4,i5)=l4(i4)
 l5A(i1,i2,i3,i4,i5)=l5(i5)
endfor
endfor
endfor
endfor
endfor

ii= median( where(phiA eq min(phiA)) )
new_guess=[l1A(ii),l2A(ii),l3A(ii),l4A(ii),l5A(ii)]

nn=1
new_guess_a=dblarr(nn,5)
new_guess_a(0,*)=new_guess

skip_min3D:

if keyword_set(multiple) then begin
print,'Explore from multiple starting points'
nn=9
new_guess_a=dblarr(nn,3)

new_guess_a(1,*)=[min(l1),min(l2),min(l3)]
new_guess_a(2,*)=[min(l1),min(l2),max(l3)]
new_guess_a(3,*)=[min(l1),max(l2),min(l3)]
new_guess_a(4,*)=[min(l1),max(l2),max(l3)]

new_guess_a(5,*)=[max(l1),min(l2),min(l3)]
new_guess_a(6,*)=[max(l1),min(l2),max(l3)]
new_guess_a(7,*)=[max(l1),max(l2),min(l3)]
new_guess_a(8,*)=[max(l1),max(l2),max(l3)]
endif

goto,skipprints
print,'indexes for minimum:',ii1,ii2,ii3
print,'lambda:',new_guess
stop
 for i=0,n1-1 do print,min(phiA(i,*,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,i,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,*,i))/min(phia)
 print,'-'
 stop
skipprints:

return
end
pro start_point_g2cfw,new_guess,Xc,nn,new_guess_a,npar,multiple=multiple
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
common parametrizacion,parametrization
common gauss_cond,T335,X335,Q335_max,lambda_inicial,DEMc1,DEMc2,DEMc3,T131,X131,Q131_max,SIGMA

if new_guess(0) ne -666. AND keyword_set(multiple) then goto,skip_min3D


n1= 1
n2= 1
n3= 1
n4= 11


phiA = dblarr(n1,n2,n3,n4)

hmin =  0.10d^.5
hmax =  10.0d^.5
dmin =            3*dx(0)        
dmax = (max(x)-min(x))/2. ; All range


l1= [lambda_inicial(0)]
l2= [lambda_inicial(1)]
l3= [lambda_inicial(2)]
l4= hmin + (hmax-hmin) *findgen(n4)/float(n4-1)




l1A = dblarr(n1,n2,n3,n4)
l2A = dblarr(n1,n2,n3,n4)
l3A = dblarr(n1,n2,n3,n4)
l4A = dblarr(n1,n2,n3,n4)



for i1=0,n1-1 do begin
for i2=0,n2-1 do begin
for i3=0,n3-1 do begin
for i4=0,n4-1 do begin
phiA(i1,i2,i3,i4)=PHI([l1(i1),l2(i2),l3(i3),l4(i4)])
 l1A(i1,i2,i3,i4)=l1(i1)
 l2A(i1,i2,i3,i4)=l2(i2)
 l3A(i1,i2,i3,i4)=l3(i3)
 l4A(i1,i2,i3,i4)=l4(i4)
endfor
endfor
endfor
endfor

ii= median( where(phiA eq min(phiA)) )
new_guess=[l1A(ii),l2A(ii),l3A(ii),l4A(ii)]


nn=1
new_guess_a=dblarr(nn,4)
new_guess_a(0,*)=new_guess

skip_min3D:

if keyword_set(multiple) then begin
print,'Explore from multiple starting points'
nn=9
new_guess_a=dblarr(nn,3)

new_guess_a(1,*)=[min(l1),min(l2),min(l3)]
new_guess_a(2,*)=[min(l1),min(l2),max(l3)]
new_guess_a(3,*)=[min(l1),max(l2),min(l3)]
new_guess_a(4,*)=[min(l1),max(l2),max(l3)]

new_guess_a(5,*)=[max(l1),min(l2),min(l3)]
new_guess_a(6,*)=[max(l1),min(l2),max(l3)]
new_guess_a(7,*)=[max(l1),max(l2),min(l3)]
new_guess_a(8,*)=[max(l1),max(l2),max(l3)]
endif

goto,skipprints
print,'indexes for minimum:',ii1,ii2,ii3
print,'lambda:',new_guess
stop
 for i=0,n1-1 do print,min(phiA(i,*,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,i,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,*,i))/min(phia)
 print,'-'
 stop
skipprints:

return
end
;+---------------------------------------------------------------
pro start_point_3c,new_guess,Xc,nn,new_guess_a,npar,multiple=multiple
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
common parametrizacion,parametrization
common gauss_cond,T335,X335,Q335_max,lambda_inicial,DEMc1,DEMc2,DEMc3,T131,X131,Q131_max,SIGMA

if new_guess(0) ne -666. AND keyword_set(multiple) then goto,skip_min3D


n1= 1
n2= 1
n3= 1
n4= 1
n5= 1
n6= 21
phiA = dblarr(n1,n2,n3,n4,n5,n6)

hmin =  0.10d^.5
hmax =  10.0d^.5

l1= [lambda_inicial(0)]
l2= [lambda_inicial(1)]
l3= [lambda_inicial(2)]
l4= [lambda_inicial(3)]
l5= [lambda_inicial(4)]
l6= hmin + (hmax-hmin) *findgen(n6)/float(n6-1)

l1A = dblarr(n1,n2,n3,n4,n5,n6)
l2A = dblarr(n1,n2,n3,n4,n5,n6)
l3A = dblarr(n1,n2,n3,n4,n5,n6)
l4A = dblarr(n1,n2,n3,n4,n5,n6)
l5A = dblarr(n1,n2,n3,n4,n5,n6)
l6A = dblarr(n1,n2,n3,n4,n5,n6)

for i1=0,n1-1 do begin
for i2=0,n2-1 do begin
for i3=0,n3-1 do begin
for i4=0,n4-1 do begin
for i5=0,n5-1 do begin
for i6=0,n6-1 do begin
phiA(i1,i2,i3,i4,i5,i6)=PHI([l1(i1),l2(i2),l3(i3),l4(i4),l5(i5),l6(i6)])
 l1A(i1,i2,i3,i4,i5,i6)=l1(i1)
 l2A(i1,i2,i3,i4,i5,i6)=l2(i2)
 l3A(i1,i2,i3,i4,i5,i6)=l3(i3)
 l4A(i1,i2,i3,i4,i5,i6)=l4(i4)
 l5A(i1,i2,i3,i4,i5,i6)=l5(i5)
 l6A(i1,i2,i3,i4,i5,i6)=l6(i6)
endfor
endfor
endfor
endfor
endfor
endfor

ii= median( where(phiA eq min(phiA)) )
new_guess=[l1A(ii),l2A(ii),l3A(ii),l4A(ii),l5A(ii),l6A(ii)]

nn=1
new_guess_a=dblarr(nn,6)
new_guess_a(0,*)=new_guess

skip_min3D:

if keyword_set(multiple) then begin
print,'Explore from multiple starting points'
nn=9
new_guess_a=dblarr(nn,3)

new_guess_a(1,*)=[min(l1),min(l2),min(l3)]
new_guess_a(2,*)=[min(l1),min(l2),max(l3)]
new_guess_a(3,*)=[min(l1),max(l2),min(l3)]
new_guess_a(4,*)=[min(l1),max(l2),max(l3)]

new_guess_a(5,*)=[max(l1),min(l2),min(l3)]
new_guess_a(6,*)=[max(l1),min(l2),max(l3)]
new_guess_a(7,*)=[max(l1),max(l2),min(l3)]
new_guess_a(8,*)=[max(l1),max(l2),max(l3)]
endif

goto,skipprints
print,'indexes for minimum:',ii1,ii2,ii3
print,'lambda:',new_guess
stop
 for i=0,n1-1 do print,min(phiA(i,*,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,i,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,*,i))/min(phia)
 print,'-'
 stop
skipprints:

return
end

;+---------------------------------------------------------------
pro start_point_7par,new_guess,Xc,nn,new_guess_a,npar,multiple=multiple
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
common parametrizacion,parametrization
common gauss_cond,T335,X335,Q335_max,lambda_inicial,DEMc1,DEMc2,DEMc3,T131,X131,Q131_max,SIGMA

if new_guess(0) ne -666. AND keyword_set(multiple) then goto,skip_min3D


n1= 1
n2= 1
n3= 1
n4= 1
n5= 1
n6= 5
n7= 5
phiA = dblarr(n1,n2,n3,n4,n5,n6,n7)

hmin =  0.1d^.5
hmax =  10d^.5
dmin =  X131*abs(1.-10.^.05)        ;.005
dmax =  X131*abs(1.-10.^.1)         ;.005

l1= [lambda_inicial(0)]
l2= [lambda_inicial(1)]
l3= [lambda_inicial(2)]
l4= [lambda_inicial(3)]
l5= [lambda_inicial(4)]
l6= hmin + (hmax-hmin) *findgen(n6)/float(n6-1)
l7= dmin + (dmax-dmin) *findgen(n7)/float(n7-1)

l1A = dblarr(n1,n2,n3,n4,n5,n6,n7)
l2A = dblarr(n1,n2,n3,n4,n5,n6,n7)
l3A = dblarr(n1,n2,n3,n4,n5,n6,n7)
l4A = dblarr(n1,n2,n3,n4,n5,n6,n7)
l5A = dblarr(n1,n2,n3,n4,n5,n6,n7)
l6A = dblarr(n1,n2,n3,n4,n5,n6,n7)
l7A = dblarr(n1,n2,n3,n4,n5,n6,n7)

for i1=0,n1-1 do begin
for i2=0,n2-1 do begin
for i3=0,n3-1 do begin
for i4=0,n4-1 do begin
for i5=0,n5-1 do begin
for i6=0,n6-1 do begin
for i7=0,n7-1 do begin
phiA(i1,i2,i3,i4,i5,i6,i7)=PHI([l1(i1),l2(i2),l3(i3),l4(i4),l5(i5),l6(i6),l7(i7)])
 l1A(i1,i2,i3,i4,i5,i6,i7)=l1(i1)
 l2A(i1,i2,i3,i4,i5,i6,i7)=l2(i2)
 l3A(i1,i2,i3,i4,i5,i6,i7)=l3(i3)
 l4A(i1,i2,i3,i4,i5,i6,i7)=l4(i4)
 l5A(i1,i2,i3,i4,i5,i6,i7)=l5(i5)
 l6A(i1,i2,i3,i4,i5,i6,i7)=l6(i6)
 l7A(i1,i2,i3,i4,i5,i6,i7)=l7(i7)
endfor
endfor
endfor
endfor
endfor
endfor
endfor

ii= median( where(phiA eq min(phiA)) )
new_guess=[l1A(ii),l2A(ii),l3A(ii),l4A(ii),l5A(ii),l6A(ii),l7A(ii)]

nn=1
new_guess_a=dblarr(nn,7)
new_guess_a(0,*)=new_guess

skip_min3D:

if keyword_set(multiple) then begin
print,'Explore from multiple starting points'
nn=9
new_guess_a=dblarr(nn,3)

new_guess_a(1,*)=[min(l1),min(l2),min(l3)]
new_guess_a(2,*)=[min(l1),min(l2),max(l3)]
new_guess_a(3,*)=[min(l1),max(l2),min(l3)]
new_guess_a(4,*)=[min(l1),max(l2),max(l3)]

new_guess_a(5,*)=[max(l1),min(l2),min(l3)]
new_guess_a(6,*)=[max(l1),min(l2),max(l3)]
new_guess_a(7,*)=[max(l1),max(l2),min(l3)]
new_guess_a(8,*)=[max(l1),max(l2),max(l3)]
endif

goto,skipprints
print,'indexes for minimum:',ii1,ii2,ii3
print,'lambda:',new_guess
stop
 for i=0,n1-1 do print,min(phiA(i,*,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,i,*))/min(phia)
 print,'-'
 for i=0,n2-1 do print,min(phiA(*,*,i))/min(phia)
 print,'-'
 stop
skipprints:

return
end

;-----------------------------------------------------------------
pro experiment,sca,Xc,logaT,Q171,Q195,Q284,noiselevel,box,exper,Ck
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny

 Nh = 10
 Nd = 10
 Nx = 10

 Mindx= dx/2.
 Maxdx= dx*nTbins

 Minh=0.1
 Maxh=10.

 Minx=min(x)
 Maxx=max(x)

 L0A = Minh +(Maxh -Minh )*findgen(Nh)/float(Nh-1)
 XcA = Minx +(Maxx -Minx )*findgen(Nx)/float(Nx-1)
 dxA = Mindx+(Maxdx-Mindx)*findgen(Nd)/float(Nd-1)

 DEMc = 1.e15

 openw,1,'./LDEM_expA-L'+strmid(string(nTbins),6,2)+'_NoiseLev-'+strmid(string(noiselevel),5,4)+'.txt' 
 printf,1,'Lambda0 Lambda1/Min-x  Lambda2/Dx      FBE_recovered/FBE_noisy           FBE_noisy/FBE_clean'
 print  , 'Lambda0 Lambda1/Min-x  Lambda2/Dx      FBE_recovered/FBE_noisy           FBE_noisy/FBE_clean'

 for ih=0,Nh-1 do begin
 L0=L0A(ih)

 for ic=0,nx-1 do begin
 xc=XcA(ic)
 for id=0,nd-1 do begin
 dd=dxa(id)

 lambda_input=[L0,xc,dd]
 
 FBE_clean = SynthFBE(lambda_input) ; This is the CLEAN forwardly computed FBE.

; Make a noisy FBE from FBE_clean:
  FBE_noisy = FBE_clean * (randomn(seed,nband)*noiselevel+1.)
; Assign FBE_noisy to FBE
  FBE=FBE_noisy

  normalize,FBE,Qkl,nband
  start_point,new_guess,Xc,nn,new_guess_a;,/multiple
  for i=0,nn-1 do $
  minPHI,smplx,sca,lambda_output,Tmin,dem_syn,reform(new_guess_a(i,*)),PHIv
  unnormalize,FBE,Qkl,nband

 FBE_recovered = SynthFBE(lambda_output) ; This is the FBE predicted by the obtained LDEM.

 printf,1,L0,Xc/MinX,dd/dx,[FBE_recovered/FBE_noisy],[FBE_noisy/FBE_clean],$
            format='(3(f8.3),"            ",3(f8.3),"          ",3(f8.3))'

   print, L0,Xc/MinX,dd/dx,[FBE_recovered/FBE_noisy],[FBE_noisy/FBE_clean],$
            format='(3(f8.3),"            ",3(f8.3),"          ",3(f8.3))'

 endfor; Sigma
 endfor; Centroid
 endfor; Area

 close,/all

return
end

;-----------------------------------
pro normalize,fbe,qkl,nband
common weights,Ck
fbe=fbe/Ck
for ik=0,nband-1 do qkl(ik,*)=qkl(ik,*)/Ck(ik)
return
end
pro unnormalize,fbe,qkl,nband
common weights,Ck
fbe=fbe*Ck
for ik=0,nband-1 do qkl(ik,*)=qkl(ik,*)*Ck(ik)
return
end

pro rotate,FBE_A,hdr,nband,ANGLE
for ib=0,nband-1 do begin
   tmp = reform(FBE_A(*,*,ib))
   tmp = rot(tmp,ANGLE,1,hdr.crpix1-1,hdr.crpix2-1,/pivot,missing=-999.)
   FBE_A(*,*,ib) = tmp
endfor
return
end
;------------------------------------------------------------------
pro compute_grid,hdr,ra,pa,ya,za,euvi=euvi,aia=aia

if keyword_set(aia ) then Rs=hdr.rsun_obs          ; Sun radius in arcsec
if keyword_set(euvi) then Rs=hdr.rsun              ; Sun radius in arcsec
px=hdr.cdelt1                                      ; Pixel size in arcsec

Rs=Rs/px                                           ; Sun radius in pixels
px=1./Rs                                           ; Pixel size in Rsun units

iy0=hdr.crpix1-1
iz0=hdr.crpix2-1

y = px * (findgen(hdr.naxis1)-iy0)
z = px * (findgen(hdr.naxis1)-iz0)

 u=1.+fltarr(hdr.naxis1)
ya=y#u
za=u#z

ra = sqrt(ya^2+za^2)
ta = fltarr(hdr.naxis1,hdr.naxis1)

p=where(ya gt 0.)
ta(p) = Acos( za(p) / ra(p) )
p=where(ya lt 0.)
ta(p) = 2.*!pi-Acos( za(p) / ra(p) )
p=where(ya eq 0. AND za gt 0.)
if p(0) ne -1 then ta(p)=0.
p=where(ya eq 0. AND za lt 0.)
if p(0) ne -1 then ta(p)=!pi
ta=2.*!pi-ta

PA=ta

return
end

;------------------------------------------------------------------
pro minPHI,smplx,sca,out,dem_syn,new_guess,PHIv,ftol
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
P=new_guess
ftol=1.d-4 ; convergence tolerance on the cost function value 
retry:
frprmn,P,ftol,iter,fret

;NEW CODE NEXT FOUR LINES
if fret eq -666. then begin
   ftol=ftol*5.
  ;print,ftol
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













;--------------------------------------------------------------------------------
FUNCTION SynthFBE_loops,lambda
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
synth_fbe=fltarr(nband)
for k=0,nband-1 do begin
  tempsum=0.d0
  for l=0,nTbins-1  do begin
  tempsum=tempsum+Qkl(k,l)*DEMc*FDEM(X[l],lambda)
  endfor ; Temp loop
  synth_fbe(k)=tempsum
endfor ; Detector loop
return,synth_fbe
end
;--------------------------------------------------------------------------------
FUNCTION PHI_loops,lambda
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
sum=0.d0
for k=0,nband-1 do begin
    tempsum=0.d0
    for l=0,nTbins-1  do begin
    tempsum=tempsum+Qkl(k,l)*DEMc*FDEM(X(l),lambda)
    endfor ; Temp loop
    sum=sum+(FBE(k)-tempsum)^2
endfor ; Detector loop
return,sum
end
;--------------------------------------------------------------------------------
FUNCTION gradPHI_loops,lambda
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny
S1=0.d0
S2=0.d0
S3=0.d0
for k=0,nband-1 do begin
Sk  =0.d0
Sk_1=0.d0
Sk_2=0.d0
Sk_3=0.d0
for l=0,nTbins-1   do begin
Sk  =Sk  +Qkl(k,l)* DEMc*FDEM(X(l),lambda)
Sk_1=Sk_1+Qkl(k,l)*(DEMc*FDEM(X(l),lambda)/lambda[0])
Sk_2=Sk_2+Qkl(k,l)*(DEMc*FDEM(X(l),lambda)/lambda[2]) * ((X(l)-lambda[1])/lambda[2])
Sk_3=Sk_3+Qkl(k,l)*(DEMc*FDEM(X(l),lambda)/lambda[2]) *(((X(l)-lambda[1])/lambda[2])^2-1)
endfor ; Temp loop
S1=S1+(FBE(k)-Sk)*Sk_1
S2=S2+(FBE(k)-Sk)*Sk_2
S3=S3+(FBE(k)-Sk)*Sk_3
endfor ; Detector loop
return,-2.*[S1,S2,S3]
end
;--------------------------------------------------------------------------------



;------------------------------------------------------------------
pro minPHI_OLD,smplx,sca,out,dem_syn,new_guess,in,exper,PHIv
common parameters,FBE,Qkl,X,dX,nTbins,DEMc,Tc,nDetect,nband,Tmin,tiny

;goto,skip_frprmn
P=new_guess
ftol=1.d-6 ; convergence tolerance on the cost function value 
frprmn,P,ftol,iter,fret
out=P
PHIv=fret
return

if exper eq 0 then return

print,'FRPRMN Parameters: ',out
print,'COST FUNCTION    : ',fret
print,'Number of iters  : ',iter
print,'Set-up Parameters: ',in
print,'Search Parameters: ',new_guess
print,'FRPRMN / Setup   : ',out/in
skip_frprmn:

goto,skip_dfp
P=new_guess
Gtol = 1.0d-2
DFPMIN, P, Gtol, Fmin, 'PHI', 'gradPHI', tolx=5.e-8 , eps=3.e-8, iter=iter, itmax=200,stepmax=10000.;,/double
out=P
print,'Set-up Parameters: ',in
print,'Search Parameters: ',new_guess
print,'DFP    Parameters: ',out
print,'DFP    / Setup   : ',out/in
print,iter
skip_dfp:

goto,skip_powell
 P = new_guess
xi = TRANSPOSE([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]])
POWELL, P, xi, 1.e-3, fmin, 'PHI'
out=P
print,'Set-up Parameters: ',in
print,'Search Parameters: ',new_guess
print,'Powell Parameters: ',out
print,'Powell / Setup   : ',out/in       
Skip_powell:

goto,skip_amoeba
;  print,smplx
  result = AMOEBA(1.e-3,FUNCTION_NAME='PHI',simplex=smplx,$
                  ncalls=iter,FUNCTION_VALUE = fval) 
  IF N_ELEMENTS(Result) EQ 1 THEN MESSAGE, 'AMOEBA failed to converge'
  out = smplx[*,0]
  print,'Set-up Parameters: ',in
  print,'Search Parameters: ',new_guess
  print,'Amoeba Parameters: ',out
  print,'Amoeba / Setup   : ',out/in         
  print,'----------------------------------------------------------'
  goto,fin
skip_amoeba:

 SigT_DT = in(2)/(x(1)-x(0))

 XX=min(x)+(max(x)-min(x))*findgen(1000)/999.
 if sca eq 1 then $
 plot,XX*Tmin/1.e6,FDEM(XX,out),xtitle='T [MK]',title='SigT/DT='+strmid(string(SigT_DT),0,10),$
 xr=[0.,3.5],xstyle=1,$
 yr=[-.1,1.1]*MAX([FDEM(XX,out),dem_syn/DEMc,FDEM(XX,new_guess),FDEM(XX,in)]),ystyle=1,th=1
 loadct,27
 oplot,XX*Tmin/1.e6,FDEM(XX,out),color=200
 loadct,0
 if sca eq 0 then $
 plot,10.^XX,FDEM(XX,out),xtitle='T/T!DMIN!N',title='DEM(T)/DEM!DC!N',$
 xr=[0.,max(10.^x)],xstyle=1,$
 yr=[0,MAX([FDEM(XX,out),dem_syn/DEMc,FDEM(XX,new_guess)])]*1.2,ystyle=1,th=1

 if sca eq 1 then oplot,x*Tmin/1.e6,    dem_syn/DEMc,th=2,psym=4 
 if sca eq 0 then oplot,10.^x,dem_syn/DEMc,th=2,psym=4 

  oplot,XX*Tmin/1.e6,FDEM(XX,in ),linestyle=2

; if sca eq 1 then oplot,[1,1]*in(1)    ,[0,10],th=2
; if sca eq 1 then oplot,[1,1]*10.^in(1),[0,10],th=2

 loadct,27
 oplot,xx*Tmin/1.e6,FDEM(XX,new_guess),linestyle=3,color=100
 loadct,0

 wait,0.1
 
; stop
fin:

return
end

;-------------------------------------------------
pro qprime, Q171, Q195, Q284, logaT, M, Qp, lT_M

nt=fix((size(logat))(1))
Q=fltarr(3,nt)
Q(0,*)=Q171
Q(1,*)=Q195
Q(2,*)=Q284

ind=intarr(M+1)
ind(0)=0
ind(M)=nt

DlT=(Max(logaT)-Min(logaT))/M
for i=1,M-1 do begin
lT0=float(i)*DlT+min(logaT)
f  =abs(logaT-lT0)
ind(i)=where(f eq min(f))
endfor

Qp = fltarr(3,M)

for r=0,2   do begin
for c=0,M-1 do begin
Qp(r,c)=total(Q(r,ind(c):ind(c+1)-1));/float(ind(c+1)-ind(c))
endfor
endfor

lT_M=DlT*findgen(M)+DlT/2.+min(logat)

stop
return
end

;-------------------------------------------------
pro qprime_explore2, Q171, Q195, Q284, logaT, M, Qp, lT_M


nt=fix((size(logat))(1))
Q=fltarr(2,nt)
Q(0,*)=Q171
Q(1,*)=Q195
;Q(2,*)=Q284

ind=intarr(M+1)
ind(0)=0
ind(M)=nt

;DlT=(Max(logaT)-Min(logaT))/M
;for i=1,M-1 do begin
;lT0=float(i)*DlT+min(logaT)
;f  =abs(logaT-lT0)
;ind(i)=where(f eq min(f))
;endfor
DlT=(Max(logaT)-Min(logaT))/M
Qp = fltarr(2,M)

n=(nt-3-1+1)*(nt-2-1+1)/2
results=fltarr(n,4)
z=0

;openw,1,'Tbins.txt'
for i=1  ,nt-2 do begin
ind(1)=i
for r=0,1   do begin
for c=0,M-1 do begin
Qp(r,c)=total(Q(r,ind(c):ind(c+1)-1));/float(ind(c+1)-ind(c))
endfor
endfor

lT_M=DlT*findgen(M)+DlT/2.+min(logat)

Qkl = Qp*1.d0
meanQkl=mean(Qkl)
Mkl=Qkl/meanQkl

SVDC, Mkl, W, U, V, /double   
;print, reform([z,10.^logaT(ind)/1.e6,W(0)/W(1)])
results(z,*)=reform([10.^logaT(ind)/1.e6,max(W)/min(W)])
if z eq 45 then goto,out
z=z+1
endfor

out:
Lt_M=[results(z,0)+0.5*(results(z,1)-results(z,0)),$
      results(z,1)+0.5*(results(z,2)-results(z,1))  ]


return
end

;-------------------------------------------------
pro qprime_explore, Q171, Q195, Q284, logaT, M, Qp, lT_M


nt=fix((size(logat))(1))
Q=fltarr(3,nt)
Q(0,*)=Q171
Q(1,*)=Q195
Q(2,*)=Q284

ind=intarr(M+1)
ind(0)=0
ind(M)=nt

;DlT=(Max(logaT)-Min(logaT))/M
;for i=1,M-1 do begin
;lT0=float(i)*DlT+min(logaT)
;f  =abs(logaT-lT0)
;ind(i)=where(f eq min(f))
;endfor
DlT=(Max(logaT)-Min(logaT))/M
Qp = fltarr(3,M)

n=(nt-3-1+1)*(nt-2-1+1)/2
results=fltarr(n,5)
z=0

;openw,1,'Tbins.txt'
for i=1  ,nt-3 do begin
for j=i+1,nt-2 do begin

ind(1)=i
ind(2)=j

for r=0,2   do begin
for c=0,M-1 do begin
Qp(r,c)=total(Q(r,ind(c):ind(c+1)-1));/float(ind(c+1)-ind(c))
endfor
endfor

lT_M=DlT*findgen(M)+DlT/2.+min(logat)

Qkl = Qp*1.d0
meanQkl=mean(Qkl)
Mkl=Qkl/meanQkl

SVDC, Mkl, W, U, V, /double   
;printf,1, reform([10.^logaT(ind)/1.e6,W(0)/W(2)])
;print,ind(1:M-1)-ind(0:M-2)
results(z,*)=reform([10.^logaT(ind)/1.e6,max(W)/min(W)])
if z eq 4344 then goto,out
z=z+1
endfor
endfor
;close,1
out:
Lt_M=[results(z,0)+0.5*(results(z,1)-results(z,0)),$
      results(z,1)+0.5*(results(z,2)-results(z,1)),$
      results(z,2)+0.5*(results(z,3)-results(z,2))  ]
;stop
return
end

;-------------------------------------------------
pro setsimplex,DX,L,sca,Tc,Tmin,X,smplx,new_guess
if sca eq 0 then Xc = alog10(Tc/Tmin)
if sca eq 1 then Xc =        Tc/Tmin

;Proposed sigmaX for simplex is "M" X-bins, enter "M"
M=10.
  smplx=[ [ 0.1,    Xc,DX      ], $
          [10.0,    Xc,DX      ], $
          [ 1.0,max(X),DX      ], $
          [ 1.0,    Xc,DX*M    ]   ] * 1.d0

a=new_guess(0)
b=new_guess(1)
c=new_guess(2)

  smplx=[ [    a ,    b ,    c ], $
          [1.1*a ,1.1*b ,    c ], $
          [    a ,1.1*b ,1.1*c ], $
          [1.1*a ,    b ,1.1*c ]  ]  * 1.d0
return
end

;-------------------------------------------------
pro setsimplex2,new_guess,smplx

a=new_guess(0)
b=new_guess(1)
c=new_guess(2)

  smplx=[ [    a ,    b ,    c ], $
          [1.1*a ,1.1*b ,    c ], $
          [    a ,1.1*b ,1.1*c ], $
          [1.1*a ,    b ,1.1*c ]  ]  * 1.d0
return
end


;--------------------EJEMPLO DE MINIMIZACION-----------
;-----------------------------------------------------
pro test
common parameter,a

a=[1.23e2,2.,5.e-2,1.e1]

lambda_guess=[1.,1.,1.,1.]

print,reform(lambda_guess/a),$
       format='(("lambda_guess / a:    "),4(e15.3))'
print

P     = lambda_guess
ftol  = 1.d-6 ; convergence tolerance on the cost function value 
frprmn,P,ftol,iter,fret
lambda_out=P      ; out: parameter array that minimized PHI
PHIv      =fret   ; the value PHI(out)

print
print,reform(lambda_out/a),$
       format='(("lambda_out   / a:    "),4(e15.3))'
print
print,'Iterations of FRPRMN:',iter
print,'Final value of PHI:  ',fret
print

return
end

;------------------------------------------------------
function PHI_test,lambda
common parameter,a
C=a;*0.+1.
return,total((lambda-a)^2/c^2)
end

function gradPHI_test,lambda
common parameter,a
c=a;*0.+1.
return,2.*(lambda-a)/c^2
end
;------------------------------------------------------

pro grilla,nr,nth,np,rad,lat,lon
drad=0.01
dt  =2.
rad =   1. + drad/2. + drad * findgen(Nr)
lat = -90. + dt  /2. + dt   * findgen(nth)
lon =   0. + dt  /2. + dt   * findgen(np)
return
end



pro second_gaussian_guess,sca,Tmin,FBE,logaT,Qkl,DEMc,Tc,Xc,nband,nTbins,X,dX,sigma

n1= 11
n2= 11

sigma=(0.37/2)*1.e6/Tmin

phiA = dblarr(n1,n2)


hmin =  0.1d 
hmax = 10.0d 
cmin = 1.2e6/Tmin
cmax = 2.3e6/Tmin


l1= hmin + (hmax-hmin) *findgen(n1)/float(n1-1)
l2= cmin + (cmax-cmin) *findgen(n2)/float(n2-1)


l1A = dblarr(n1,n2)
l2A = dblarr(n1,n2)


for i1=0,n1-1 do begin
   for i2=0,n2-1 do begin
      f=(l1(i1)/sqrt(2.d0*!pi)/sigma)*exp(-0.5d0*((X-l2(i2))/sigma)^2)
      phiA(i1,i2)=total((FBE-(DEMc*dX*f)##Qkl)^2)
      l1A(i1,i2)=l1(i1)
      l2A(i1,i2)=l2(i2)
   endfor
endfor


ii= median( where(phiA eq min(phiA)) )
new_guess=[l1A(ii),l2A(ii)]

DEMc =DEMc*l1A(ii)
Xc   =l2A(ii)
Tc   =Tmin*Xc

;goto,skipgraph
;window,0,xs=400,ys=400
;!p.charsize=1.5
;!p.multi=[0,1,2]
; plot,X*Tmin/1.e6,DEM1,/ylog
; plot,X*Tmin/1.e6,PHI1
;oplot,[1,1]*Tc/1.e6,[0,max(PHI1)]
;!p.multi=0
print,Tc/1.e6,DEMc,PHIA(ii)
f=(l1A(ii)/sqrt(2.d0*!pi)/sigma)*exp(-0.5d0*((X-l2A(ii))/sigma)^2)
print,'FBE_284_synth/FBE_284_tom=',((DEMc*dX*f)##Qkl)(0)/FBE(0)
print,'FBE_335_synth/FBE_335_tom=',((DEMc*dX*f)##Qkl)(1)/FBE(1)
stop
skipgraph:

return
end

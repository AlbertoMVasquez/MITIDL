

pro carrington_maps

  common comunes,tm,wt,nband,demc,PHI,parametrizacion,Tmin,Tmax,nr,nth,np,rad,lat,lon,lambda,WTc
  common results_tomo,tfbe,sfbe,N_e

  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH

  rad_range=[ 1.1 , 1.12]       ; Rsun
  lat_range=[- 40 ,   0.]       ; deg
  lon_range=[ 250 , 350.]       ; deg
  file='LDEM.v3_CR2198_l.70.90.90_h1_reduced_Rmin1.00_Rmax1.26_Nr26_InstRmax1.26_bf4_r3d_B_chianti.ioneq_sun_coronal_1992_feldman_ext.abundaia3_171_gauss1_lin_Norm-median_singlStart'

  read_ldem,file,/ldem,/gauss1
  ratio=sfbe/tfbe
  R=total( abs(1.-ratio), 4 ) / float(nband)
  ZDA = where(demc eq -999.)
  CNS = where(demc ne -999. AND (R      gt 0.25))


  superhigh=0.25
  superlow=0.01
  p=where(demc ne -999. AND R ge superhigh) & if p(0) ne -1 then R(p)=superhigh
  p=where(demc ne -999. AND R le superlow)  & if p(0) ne -1 then R(p)=superlow
  if ZDA(0) ne -1 then R(ZDA)=0.
  

  Nesat = N_e
  Tmsat = Tm
  if ZDA(0) ne -1 then begin
     Nesat(ZDA)=0.
     Tmsat(ZDA)=0.
  endif

  if CNS(0) ne -1 then begin
     Nesat(CNS)=0.
     Tmsat(CNS)=0.
  endif

  
  r0A=[1.025,1.065,1.105]


  xdisplay,map=Nesat,file='Ne_CR2198',nr=26,nt=90,rmin=1.0,rmax=1.26,r0A=r0A,win=0,minA=[0.,0.,0.],maxA=[2.2, 1.6, 1.0], clrtbl= 4,titulo='Ne [10!U8!Ncm!U-3!N]',units=1.e8,$
           rad_range=rad_range,lat_range=lat_range,lon_range=lon_range
  
  xdisplay,file='Tm_CR2198',nr= 26,nt=90,rmin=1.0,rmax=1.26,$
           r0A=r0A,win=0,titulo='CR-2198 DEMT: Tm [MK]',minA=[0,0,0],maxA=[2.5,2.5,2.5],$
           clrtbl=5,units=1.e6,map=tmsat,rad_range=rad_range,lat_range=lat_range,lon_range=lon_range


  xdisplay,map=R,file='R_CR2198',nr=26,nt=90,rmin=1.0,rmax=1.26,r0A=r0A,win=0,minA=[0.,0.,0.],maxA=[0.25, 0.25, 0.25], clrtbl= 12,titulo='R',units=1.,$
           rad_range=rad_range,lat_range=lat_range,lon_range=lon_range
  
  


  return
end

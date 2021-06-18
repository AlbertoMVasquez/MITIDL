

pro wrapper

  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH
  rad_range=[ 1.2 , 1.22]       ; Rsun
  lat_range=[- 40 ,   0.]       ; deg
  lon_range=[ 250 , 350.]       ; deg
  file='LDEM.v3_CR2198_l.70.90.90_h1_reduced_Rmin1.00_Rmax1.26_Nr26_InstRmax1.26_bf4_r3d_B_chianti.ioneq_sun_coronal_1992_feldman_ext.abundaia3_171_gauss1_lin_Norm-median_singlStart'
  median_values_box,rad_range,lat_range,lon_range,file
  return
end


pro median_values_box,rad_range,lat_range,lon_range,file
  common comunes,tm,wt,nband,demc,PHI,parametrizacion,Tmin,Tmax,nr,nth,np,rad,lat,lon,lambda,WTc
  common results_tomo,tfbe,sfbe,N_e


  read_ldem,file,/ldem,/gauss1

  nband=3
  x171 = tfbe(*,*,*,0)
  x193 = tfbe(*,*,*,1)
  x211 = tfbe(*,*,*,2)
  ratio=sfbe/tfbe
  R=total( abs(1.-ratio), 4 ) / float(nband)
  ZDA = where(demc eq -999.)
  CNS = where(demc ne -999. AND (R      gt 0.25))
  grilla,nr,nth,np,radA,latA,lonA
  
  OKvoxel = demc*0. + 1.
  OKvoxel(ZDA) = 0.
  OKvoxel(CNS) = 0.


   index = where (OKvoxel eq 1. and $
                   radA gt rad_range(0) and radA lt rad_range(1) and $
                   latA gt lat_range(0) and latA lt lat_range(1) and $
                   lonA gt lon_range(0) and lonA lt lon_range(1)      )  

  print,'median value of Nm in the box: [10^8 cm-3]',median(n_e(index))/1.e8
  print,'median value of Tm in the box: [MK]       ',median(tm (index))/1.e6
  print,'median value of WT in the box: [MK]       ',median(wt (index))/1.e6
  print
  print,'median value of FBE 171 A in the box:',median(x171(index))
  print,'median value of FBE 193 A in the box:',median(x193(index))
  print,'median value of FBE 211 A in the box:',median(x211(index))
  print

  dir = '/data1/tomography/bindata/'
  ;xread,dir=dir,file='x.comp1074.dynamics.Dt2_CR2198.bf2.ri1.00.ro1.50_50_90_180_r3D_1.7_IRMIN_1.09',nr=50,nt=90,np=180,map=comp1074
  xread,dir=dir,file='x.comp1074_Rmin1.0_Rmax1.5_IRmin1.09_IRmax1.3_50x90x180_BF2_L1.7',nr=50,nt=90,np=180,map=comp1074
  ;xread,dir=dir,file='x.comp1079.dynamics.Dt2_CR2198.bf2.ri1.00.ro1.50_50_90_180_r3D_L3.0_IRMIN_1.09',nr=50,nt=90,np=180,map=comp1079
  xread,dir=dir,file='x.comp1079_Rmin1.0_Rmax1.5_IRmin1.09_IRmax1.3_50x90x180_BF2_L3.0',nr=50,nt=90,np=180,map=comp1079
  nr = 50 & nth = 90 & np = 180
  grilla,nr,nth,np,radA,latA,lonA  
  index = where (comp1074 gt 0. and $
                 radA gt rad_range(0) and radA lt rad_range(1) and $
                 latA gt lat_range(0) and latA lt lat_range(1) and $
                 lonA gt lon_range(0) and lonA lt lon_range(1)      )
    print,'median value of CoMP 1074 in the box:',median(comp1074(index))
    
    index = where (comp1079 gt 0. and $
                   radA gt rad_range(0) and radA lt rad_range(1) and $
                   latA gt lat_range(0) and latA lt lat_range(1) and $
                   lonA gt lon_range(0) and lonA lt lon_range(1)      )
    print,'median value of CoMP 1079 in the box:',median(comp1079(index))
    print

  return
end


pro grilla,nr,nth,np,radA,latA,lonA
  drad=0.01
  dt  =2.
  rad =   1. + drad/2. + drad * findgen(Nr)
  lat = -90. + dt  /2. + dt   * findgen(nth)
  lon =   0. + dt  /2. + dt   * findgen(np)
  radA=fltarr(nr,nth,np)
  latA=fltarr(nr,nth,np)
  lonA=fltarr(nr,nth,np)
  for irad=0,nr -1 do radA(irad,*   ,*   )=rad(irad)
  for ilat=0,nth-1 do latA(*   ,ilat,*   )=lat(ilat)
  for ilon=0,np -1 do lonA(*   ,*   ,ilon)=lon(ilon)
  return
end

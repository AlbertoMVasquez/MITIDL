

pro wrapper

  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH
  rad_range=[ 1.1 , 1.12]       ; Rsun
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

  print,'median value of Ne in the box: [10^8 cm-3]',median(n_e(index))/1.e8
  print,'median value of Te in the box: [MK]       ',median(tm (index))/1.e6
  print
  print,'median value of FBE 171 A in the box:',median(x171(index))
  print,'median value of FBE 193 A in the box:',median(x193(index))
  print,'median value of FBE 211 A in the box:',median(x211(index))

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

pro analyze_xfiles

  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH

; Load tomographies into memory
  dir = '/data1/tomography/bindata/'
  xread,dir=dir,file='x_KCOR.CR2198.13imgs-reduced.bf2.ri1.05.ro2.25_Inst_1.09_2.00_120_90_180_dropneg_r3D_l1e-4',nr=120,nt=90,np=180,map=y0
  xread,dir=dir,file='x.comp1074.dynamics.Dt2_CR2198.bf2.ri1.00.ro1.50_50_90_180_r3D_1.7_IRMIN_1.09',nr=50,nt=90,np=180,map=comp1074
  xread,dir=dir,file='x.comp1079.dynamics.Dt2_CR2198.bf2.ri1.00.ro1.50_50_90_180_r3D_L3.0_IRMIN_1.09',nr=50,nt=90,np=180,map=comp1079
  xread,dir=dir,file='x_aia.171.cr2198.26x90_bf4_ri.00_ro1.09_h1_Oldset_r3d_reduced_L0.70',nr=26,nt=90,np=180,map=FBE171
  xread,dir=dir,file='x_aia.193.cr2198.26x90_bf4_ri.00_ro1.09_h1_Oldset_r3d_reduced_L0.90',nr=26,nt=90,np=180,map=FBE193
  xread,dir=dir,file='x_aia.211.cr2198.26x90_bf4_ri.00_ro1.09_h1_Oldset_r3d_reduced_L0.90',nr=26,nt=90,np=180,map=FBE211

; Bring FBEs into CGS units:
  FBE171 = FBE171/0.0696
  FBE193 = FBE193/0.0696
  FBE211 = FBE211/0.0696

; Heights at which to make Carrmaps:
  r0A = [1.1,1.2]
  
; Box to analyze:
  rad_range=[ 1.1 , 1.12]       ; Rsun
  rad_range=[ 1.2 , 1.22]       ; Rsun
  lat_range=[- 40 ,   0.]       ; deg
  lon_range=[ 250 , 350.]       ; deg

; Common settings:
   nt   = 90
   np   = 2*nt
   scalefactor=4
   win=0

; CoMP settings:
   file = 'x.comp1074.dynamics.Dt2_CR2198.bf2.ri1.00.ro1.50_50_90_180_r3D_1.7_IRMIN_1.09'
   Rmin   = 1.0
   Rmax   = 1.5
   Nr     =  50
  instrument = 'COMP'
  titulo='CR-2198 CoMP 1074'
   xdisplay,dir=dir,map=comp1074,nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,r0A=r0A,units=units,titulo=titulo,$
            rad_range=rad_range,lat_range=lat_range,lon_range=lon_range,file=file

   file = 'x.comp1079.dynamics.Dt2_CR2198.bf2.ri1.00.ro1.50_50_90_180_r3D_L3.0_IRMIN_1.09'
   Rmin   = 1.0
   Rmax   = 1.5
   Nr     =  50
  instrument = 'COMP'
  titulo='CR-2198 CoMP 1079'
   xdisplay,dir=dir,map=comp1079,nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,r0A=r0A,units=units,titulo=titulo,$
            rad_range=rad_range,lat_range=lat_range,lon_range=lon_range,file=file

   stop

   
; KCOR settings:

   file='x_KCOR.CR2198.13imgs-reduced.bf2.ri1.05.ro2.25_Inst_1.09_2.00_120_90_180_dropneg_r3D_l1e-4'
   rmin = 1.05
   rmax = 2.25
   nr   = 120
   instrument = 'kcor'
; KCOR analysis:   
   units=1.e8
   titulo='CR-2198 VL Ne [10!U8!N cm!U-3!N]'
   xdisplay,dir=dir,map=y0,nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,r0A=r0A,units=units,titulo=titulo,$
            rad_range=rad_range,lat_range=lat_range,lon_range=lon_range,file=file



; AIA settings:
   rmin = 1.00
   rmax = 1.26
   nr   = 26
   instrument = 'aia'
   units=1.   

; AIA analysis:

   file  ='x_aia.171.cr2198.26x90_bf4_ri.00_ro1.09_h1_Oldset_r3d_reduced_L0.70'
   titulo='CR-2198 FBE-171 [PH cm!U-3!N s!U-1!N sr!U-1!N]'
   xdisplay,dir=dir,map=FBE171,nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,r0A=r0A,units=units,titulo=titulo,$
            rad_range=rad_range,lat_range=lat_range,lon_range=lon_range,file=file
   
   file  ='x_aia.193.cr2198.26x90_bf4_ri.00_ro1.09_h1_Oldset_r3d_reduced_L0.90'   
   titulo='CR-2198 FBE-193 [PH cm!U-3!N s!U-1!N sr!U-1!N]'
   xdisplay,dir=dir,map=FBE193,nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,r0A=r0A,units=units,titulo=titulo,$
            rad_range=rad_range,lat_range=lat_range,lon_range=lon_range,file=file
   
   file  ='x_aia.211.cr2198.26x90_bf4_ri.00_ro1.09_h1_Oldset_r3d_reduced_L0.90'
   titulo='CR-2198 FBE-211 [PH cm!U-3!N s!U-1!N sr!U-1!N]'
   xdisplay,dir=dir,map=FBE211,nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,r0A=r0A,units=units,titulo=titulo,$
            rad_range=rad_range,lat_range=lat_range,lon_range=lon_range,file=file

   stop
  return
end

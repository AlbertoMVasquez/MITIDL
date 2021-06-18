pro xdisplay_comp


; Heights at which to make Carrmaps:
  r0A = [1.1,1.2]

; expand path to use the x-tools in SolarTom_idl directory
  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH

; lat and lon grid parameters:
  nt     = 90
  np     = 2*nt  
; radial grid parameters:
  Rmin   = 1.0
  Rmax   = 1.5
  Nr     =  50

; CoMP Maps;
  xread,file='x.comp1074.dynamics.Dt2_CR2198.bf2.ri1.00.ro1.50_50_90_180_r3D_1.7_IRMIN_1.09',nr=nr,nt=nt,np=np,map=comp1074_old
  file = 'x.comp1074.old'
  titulo='CoMP 1074 (OLD)'
  xdisplay,map=comp1074_old,nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,r0A=r0A,units=units,titulo=titulo,file=file,minS=minS,maxS=maxS

  xread,file='x.comp1074_Rmin1.0_Rmax1.5_IRmin1.09_IRmax1.3_50x90x180_BF2_L1.7',nr=nr,nt=nt,np=np,map=comp1074_new
  file  = 'x.comp1074.new'
  titulo='CoMP 1074 (NEW)'
  xdisplay,map=comp1074_new,nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,r0A=r0A,units=units,titulo=titulo,file=file,minA=minS,maxA=maxS
   
  xread,file='x.comp1079.dynamics.Dt2_CR2198.bf2.ri1.00.ro1.50_50_90_180_r3D_L3.0_IRMIN_1.09',nr=nr,nt=nt,np=np,map=comp1079_old
  file  = 'x.comp1079.old'
  titulo='CoMP 1079 (OLD)'
  xdisplay,map=comp1079_old,nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,r0A=r0A,units=units,titulo=titulo,file=file,minA=minS,maxA=maxS

  xread,file='x.comp1079_Rmin1.0_Rmax1.5_IRmin1.09_IRmax1.3_50x90x180_BF2_L3.0',nr=nr,nt=nt,np=np,map=comp1079_new
  file  = 'x.comp1079.new'
  titulo='CoMP 1079 (NEW)'
  xdisplay,map=comp1079_new,nr=nr,nt=nt,np=np,rmin=rmin,rmax=rmax,r0A=r0A,units=units,titulo=titulo,file=file,minA=minS,maxA=maxS


  rmin = 1.07 
  rmax = 1.32
  Nr   = 25


  rmin = 1.07 
  rmax = 1.22
  Nr   = 15



  
  return
end

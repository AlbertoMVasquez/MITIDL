
pro wrapper
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common NT_limits, Ne0_Limits, Te0_Limits

  ;COMMON FOV of all instruments 
  rmin= 1.05
  rmax= 1.25
  ;This vectors define the order of the measurement
  i_mea_vec           =[0       ,0       ,1    ,1    ,1    ]
  ion_label_vec       =['fexiii','fexiii',''   ,''   ,''   ]
  line_wavelength_vec =['10747' ,'10801' ,''   ,''   ,''   ]
  instrument_label_vec=[''      ,''      ,'aia','aia','aia']
  band_label_vec      =[''      ,''      ,'171','193','211']
  ;File name of x-tomographic products of all instruments 
  ;(need to be consistent with above)
  xfiles   = ['x_KCOR.CR2198.13imgs-reduced.bf2.ri1.05.ro2.25_Inst_1.09_2.00_120_90_180_dropneg_r3D_l1e-4',$
              'x.comp1074.dynamics.Dt2_CR2198.bf2.ri1.00.ro1.50_50_90_180_r3D_1.7_IRMIN_1.09',$
              'x.comp1079.dynamics.Dt2_CR2198.bf2.ri1.00.ro1.50_50_90_180_r3D_L3.0_IRMIN_1.09',$
              'x_aia.171.cr2198.26x90_bf4_ri.00_ro1.09_h1_Oldset_r3d_reduced_L0.70',$
              'x_aia.193.cr2198.26x90_bf4_ri.00_ro1.09_h1_Oldset_r3d_reduced_L0.90',$
              'x_aia.211.cr2198.26x90_bf4_ri.00_ro1.09_h1_Oldset_r3d_reduced_L0.90']

 ; Limits for the grid. These are DYNAMICAL as to adjust to future changes in G-tables:
 ; Ne0_Limits = [max([min(Ne1),min(Ne2),min(Ne3),min(Ne4),min(Ne5)]),min([max(Ne1),max(Ne2),max(Ne3),max(Ne4),max(Ne5)])]
 ; Te0_Limits = [max([min(Te1),min(Te2),min(Te3),min(Te4),min(Te5)]),min([max(Te1),max(Te2),max(Te3),max(Te4),max(Te5)])]
 ; Restricted Ne and Te ranges 
  Ne0_Limits = [1.0e6,5.0e9]
  Te0_Limits = [0.5e6,5.0e6]

  method=3 ; BFGS method

  
  MIT,rmin,rmax,xfiles,/Riemann,min_method=method

  return
end


;================================================================
; Master code to Multi-Instrument Tomography (MIT)

;===============================================================

pro MIT,rmin,rmax,xfiles,min_method=min_method,riemann=riemann,$
        NNe_provided=NNe_provided,NTe_provided=NTe_provided 
  common measurement_vectors,i_mea_vec,ion_label_vec,line_wavelength_vec,instrument_label_vec,band_label_vec
  common NT_limits, Ne0_Limits, Te0_Limits

  if not keyword_set(min_method) then begin
     print,'minimization method (min_method keyword) not selected:'
     print,'1: Downhill Simplex'
     print,'2: Powell          '
     print,'3: BFGS            '
     print,'4: Polak-Ribiere   '
     print,'5: minimization with constrains'
     return
  endif
  
  if keyword_set(Riemann) then begin 
     Phi_name     ='cost_function_cs'
     grad_Phi_name='grad_cost_function_cs'
  endif else begin
     Phi_name     ='cost_function'
     grad_Phi_name='grad_cost_function'
  endelse
  
  print,'Minimization of ',Phi_name
  print

  if min_method eq 1 then print,'Downhill simplex Method'
  if min_method eq 2 then print,'Powell method'
  if min_method eq 3 then print,'BFGS method'
  IF min_method eq 4 then begin
     PRINT,'Polak-Ribiere method'
     if not keyword_set(riemann) then begin
        command1='cp phi1.pro phi.pro'
        command2='cp gradphi1.pro gradphi.pro'
     endif else begin
        command1='cp phi2.pro phi.pro'
        command2='cp gradphi2.pro gradphi.pro'
     endelse
     print,command1 & spawn,command1
     print,command2 & spawn,command2
     print
  endif
  IF min_method eq 5 then print,"Constrained Minimization"

  
  if keyword_set(Riemann) then begin
     if not keyword_set(NNe_provided) then NNe_provided = 100
     if not keyword_set(NTe_provided) then NTe_provided = 100
     
     if keyword_set(   uniform) then make_grid,   /uniform,NNe_provided=NNe_provided,NTe_provided=NTe_provided
     if keyword_set(loguniform) then make_grid,/loguniform,NNe_provided=NNe_provided,NTe_provided=NTe_provided
     if keyword_set( lnuniform) then make_grid, /lnuniform,NNe_provided=NNe_provided,NTe_provided=NTe_provided
     if not keyword_set(uniform) and not keyword_set(loguniform) and not keyword_set(lnuniform) then begin
        print,'choose a grid (uniform, loguniform, or lnuniform)'
        return
     endif
  endif


  


; Load tomographies into memory
  xread,file=xfiles(0),nr=120,nt=90,np=180,map=KCOR
  xread,file=xfiles(1),nr=50 ,nt=90,np=180,map=comp1074
  xread,file=xfiles(2),nr=50 ,nt=90,np=180,map=comp1079
  xread,file=xfiles(3),nr=26 ,nt=90,np=180,map=FBE171
  xread,file=xfiles(4),nr=26 ,nt=90,np=180,map=FBE193
  xread,file=xfiles(5),nr=26 ,nt=90,np=180,map=FBE211
  
; Bring FBEs into CGS units:
  FBE171 = FBE171/0.0696
  FBE193 = FBE193/0.0696
  FBE211 = FBE211/0.0696

; This is correct ??? something more??
  comp1074=comp1074/6.96e10
  comp1079=comp1079/6.96e10

;  KCOR:
   rmin_kcor = 1.05
   rmax_kcor = 2.25
   Nr_kcor   = 120
   dr_kcor   = (Rmax_kcor - Rmin_kcor)/Nr_kcor 
   r_kcor    = Rmin_kcor + dr_kcor * indgen(Nr_kcor) + dr_kcor/2
   i_kcor = where ( r_kcor ge rmin and r_kcor le rmax ) 
   kcor = reform ( kcor(i_kcor,*,*))


;  CoMP:
   Rmin_comp   = 1.0
   Rmax_comp   = 1.5
   Nr_comp     =  50
   dr_comp     = (Rmax_comp - Rmin_comp)/Nr_comp 
   r_comp      = Rmin_comp + dr_comp * indgen(Nr_comp) + dr_comp/2
   i_comp = where ( r_comp ge rmin and r_comp le rmax ) 
   comp1074 = reform(comp1074 ( i_comp,*,*))
   comp1079 = reform(comp1079 ( i_comp,*,*))

    
;  AIA settings:
   rmin_aia = 1.00
   rmax_aia = 1.26
   nr_aia   = 26
   dr_aia   = (Rmax_aia - Rmin_aia)/Nr_aia 
   r_aia    = Rmin_aia + dr_aia * indgen(Nr_aia) + dr_aia/2
   i_aia  = where ( r_aia ge rmin and r_aia le rmax ) 
   FBE171 = reform ( FBE171(i_aia,*,*) )
   FBE193 = reform ( FBE193(i_aia,*,*) )
   FBE211 = reform ( FBE211(i_aia,*,*) )
   

;  tomographic grid common to all instrument:
   Nr  = n_elements(i_aia)
   nt  = 90
   np  = 180
   dt  = 2.
   dp  = dt 
   rad =  r_aia(i_aia)
   lat = -90. + dt/2. + dt * findgen(Nt)
   lon =   0. + dp/2. + dp * findgen(Np)


   Npar      = 6
   Ndat      = 6
   parA      = dblarr(Nr,Nt,Np,npar) -666.
   yA        = dblarr(Nr,Nt,Np,ndat) 
   ysynthA   = dblarr(Nr,Nt,Np,ndat) 
;--------------------------------

   r_array = rad
   load_sk_array,Ne_array,Te_array,r_array,sk_A
   change_units_grid

   irad1=0                      
   irad2=Nr-1                   
   ;if keyword_set(oneheight) then begin
   ;   irad1=oneheight
   ;   irad2=irad1
   ;endif
   Dirad=1
   ilat1=0
   ilat2=nt-1
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
      sk_over_fip_factor = sk_A(*,*,*,ir)
      for ith=ilat1,ilat2,Dilat do begin
         for ip =ilon1,ilon2,Dilon do begin

           ;Load y0 and y in the voxel
            y0 = kcor(ir,ith,ip) 
            y  =[comp1074(ir,ith,ip),comp1079(ir,ith,ip),$
                 FBE171(ir,ith,ip),FBE193(ir,ith,ip),FBE211(ir,ith,ip)]
            if min([y0,y]) le 0. then goto,skipmin ; SKIP ZDAs

           ;INITIAL GUESS:
            make_guess_ini_new_units,guess_ini,PHIguess
           
           ;MINIMIZATION BLOCK
            minimizador,phi_name,grad_phi_name,guess_ini,P,min_method=min_method
            parA(ir,ith,ip,*) = P ; save the parameters vector in a 3D array 
            skipmin:
 
         endfor                 ; IP  loop
      endfor                    ; ITH loop
   endfor                       ; IR loop

;--------------------------------

   

   stop
   RETURN
END

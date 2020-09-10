pro show_exp_contrl
 
; Experiment suffix (con o sin ruido)
  suffix = 'con_ruido'
  suffix_file = '_'+suffix
  suffix_title =' ('+suffix+')'

  
  print,'exp. controlado ',suffix
  print
  
; Directory where is exp_contrl.out and 
; where the plots are saved.
  dir = '~/Downloads/exp_contrl/'+suffix+'/'
  
; Store the data in memory 
  restore,dir+'exp_contr.out'

; IN and OUT arrays of Nm, Tm, FIP, sigN, sigT, q  
 ;npar=(size(par_in))(4)
  nm_in   =reform(par_in  (*,*,*,0))/1.e8
  nm_out  =reform(par_out (*,*,*,0))
  fip_in  =reform(par_in  (*,*,*,1))
  fip_out =reform(par_out (*,*,*,1))
  Tm_in   =reform(par_in  (*,*,*,2))/1.e6
  Tm_out  =reform(par_out (*,*,*,2))
  sigT_in =reform(par_in  (*,*,*,3))/1.e6
  sigT_out=reform(par_out (*,*,*,3))
  sigN_in =reform(par_in  (*,*,*,4))/1.e8
  sigN_out=reform(par_out (*,*,*,4))
  q_in    =reform(par_in  (*,*,*,5))
  q_out   =reform(par_out (*,*,*,5))


; Plots settings
  !PATH = Expand_Path('+~/idlfiles/coyote/') + ':' + !PATH
  !p.background=255
  !p.color=0

; =============== Nm ================== 
  print,'Nm comparison'
  rel_diff = (nm_out-nm_in)/nm_in 
  print_rel_diff,rel_diff,0.1
  print_rel_diff,rel_diff,0.2
  print_rel_diff,rel_diff,0.3
  window,xs=500,ys=500
  cghistoplot,rel_diff,title='Nm rel diff.'+suffix_title,xrange=[-1,1],binsize=1.e-1
  record_jpg,dir,'Nm'+suffix_file+'.jpg'
  print,'--------------------------------------------------------------'
  print
  
; =============== Tm ================== 
  print,'Tm comparison'
  rel_diff = (tm_out-tm_in)/tm_in
  print_rel_diff,rel_diff,0.1
  print_rel_diff,rel_diff,0.2
  print_rel_diff,rel_diff,0.3
  window,xs=500,ys=500
  cghistoplot,rel_diff,title='Tm rel diff.'+suffix_title,xrange=[-1,1],binsize=1.e-1
  record_jpg,dir,'Tm'+suffix_file+'.jpg'
  print,'--------------------------------------------------------------'
  print

; =============== FIP ================== 
  print,'fip comparison'
  rel_diff = (fip_out-fip_in)/fip_in
  print_rel_diff,rel_diff,0.1
  print_rel_diff,rel_diff,0.2
  print_rel_diff,rel_diff,0.3
  window,xs=500,ys=500
  cghistoplot,rel_diff,title='FIP rel diff.'+suffix_title,xrange=[-2,2],binsize=1.e-1
  record_jpg,dir,'fip_factor'+suffix_file+'.jpg'
  print,'--------------------------------------------------------------'
  print
  
; =============== sigN ================== 
   print,'sigN comparison'
   rel_diff = (sigN_out-sigN_in)/sigN_in
   print_rel_diff,rel_diff,0.1
   print_rel_diff,rel_diff,0.2
   print_rel_diff,rel_diff,0.3
   window,xs=500,ys=500
   cghistoplot,rel_diff,title='sigN rel diff'+suffix_title,xrange=[-2,2],binsize=1.e-1
   record_jpg,dir,'sigN'+suffix_file+'.jpg'
   print,'--------------------------------------------------------------'
   print
   
; =============== sigT ================== 
   print,'sigT comparison'
   rel_diff = (sigT_out-sigT_in)/sigT_in
   print_rel_diff,rel_diff,0.1
   print_rel_diff,rel_diff,0.2
   print_rel_diff,rel_diff,0.3
   window,xs=500,ys=500
   cghistoplot,rel_diff,title='sigT rel diff'+suffix_title,xrange=[-2,2],binsize=1.e-1
   record_jpg,dir,'sigT'+suffix_file+'.jpg'
   print,'--------------------------------------------------------------'
   print


; =============== q ================== 
   print,'q comparison'
   rel_diff = (q_out-q_in)/q_in
   print_rel_diff,rel_diff,0.1
   print_rel_diff,rel_diff,0.2
   print_rel_diff,rel_diff,0.3
   window,xs=500,ys=500
   cghistoplot,rel_diff,title='q rel diff'+suffix_title,xrange=[-2,2],binsize=1.e-1
   record_jpg,dir,'q'+suffix_file+'.jpg'
   print,'--------------------------------------------------------------'
   print

  return
end


pro print_rel_diff,rel_diff,frac


  print,'--------------------------------------------------------------'
  index = where (abs(rel_diff) lt frac)
  print,'fraction of voxel with rel. diff. within ',frac*100,' %:',$
     float(n_elements(index))/n_elements(rel_diff)
  

  return
end

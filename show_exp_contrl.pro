; show_exp_contrl,dir='exp_contrl_old'
; show_exp_contrl,dir='exp_contrl_old',/ruido
; show_exp_contrl,dir='exp_contrl'
; show_exp_contrl,dir='exp_contrl',/ruido
; show_exp_contrl,dir='exp_contrl2'
; show_exp_contrl,dir='exp_contrl_amoeba'
; ::::::::::::::::::::V2::::::::::::::::
; show_exp_contrl,dir='exp_contrl_v2'
; show_exp_contrl,dir='exp_contrl_v2',/ruido
; show_exp_contrl,dir='exp_contrl_v2_exp_B'
; show_exp_contrl,dir='exp_contrl_v2_exp_B',/ruido
; show_exp_contrl,dir='exp_contrl_v2_100x100'
; show_exp_contrl,dir='exp_contrl_v2_amoeba'
; show_exp_contrl,dir='exp_contrl_v2_amoeba',/ruido

pro show_exp_contrl,ruido=ruido,dir=dir


 
; Experiment suffix (con o sin ruido)
  if     keyword_set(ruido) then  suffix = 'con_ruido'
  if NOT keyword_set(ruido) then  suffix = 'sin_ruido'

  suffix_file = '_'+suffix
  suffix_title =' ('+suffix+')'

  
  print,'exp. controlado ',suffix
  print
  
; Directory where is exp_contrl.out and 
; where the plots are saved.
  if not keyword_set(dir) then dir  = '~/Downloads/exp_contrl/'+suffix+'/'
  if     keyword_set(dir) then dir  = '~/Downloads/'+dir+'/'+suffix+'/'

; Store the data in memory 
  ;restore,dir+'exp_contr.out'
  restore,dir+'mit_exp_contrl.out'

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



  ;plot_scatter_and_historatio,nm_in,nm_out,xsuffix='modeled',ysuffix='reconstructed',filename='Nm_'+suffix,titulo='Nm ('+suffix+')'
  ;plot_scatter_and_historatio,tm_in,tm_out,xsuffix='modeled',ysuffix='reconstructed',filename='Tm_'+suffix,titulo='Tm ('+suffix+')'
  ;plot_scatter_and_historatio,fip_in,fip_out,xsuffix='modeled',ysuffix='reconstructed',filename='FIP_'+suffix,titulo='FIP ('+suffix+')'
  ;plot_scatter_and_historatio,sigN_in,sigN_out,xsuffix='modeled',ysuffix='reconstructed',filename='sigN_'+suffix,titulo='sigN ('+suffix+')'
  ;plot_scatter_and_historatio,sigT_in,sigT_out,xsuffix='modeled',ysuffix='reconstructed',filename='sigT_'+suffix,titulo='sigT ('+suffix+')'
  ;plot_scatter_and_historatio,q_in,q_out,xsuffix='modeled',ysuffix='reconstructed',filename='q_'+suffix,titulo='q ('+suffix+')'
  ;return
  
  
  ;index = where (abs((sigN_in -sigN_out)/sigN_in) gt 0.5)
  ;index = where( abs(sigT_in/Tm_in - sigN_in/Nm_in) lt 1.e-1  )
  ;stop


; Plots settings
  !PATH = Expand_Path('+~/idlfiles/coyote/') + ':' + !PATH
  !p.background=255
  !p.color=0

; =============== Nm ================== 
  print,'Nm comparison'
  rel_diff = ((nm_out-nm_in)/nm_in);(index) 
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
  rel_diff = ((tm_out-tm_in)/tm_in);(index)
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
  rel_diff = ((fip_out-fip_in)/fip_in);(index)
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
   rel_diff = ((sigN_out-sigN_in)/sigN_in);(index)
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
   rel_diff = ((sigT_out-sigT_in)/sigT_in);(index)
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
   rel_diff = ((q_out-q_in)/q_in);(index)
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

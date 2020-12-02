; Examples of call sequence:
; show_exp_contrl,dir='exp_contrl_v2',suffix_exp='(SR) G.Conj ftol1e-4'
; show_exp_contrl,dir='exp_contrl_v2',/ruido,suffix_exp='(CR) G.Conj ftol1e-4'
; show_exp_contrl,dir='exp_contrl_v2_amoeba',suffix_exp='(SR) AMOEBA ftol1e-4'
; show_exp_contrl,dir='exp_contrl_v2_amoeba',suffix_exp='(CR) AMOEBA ftol1e-4',/ruido
; show_exp_contrl,dir='exp_contrl_v2_ftol1e-6',suffix_exp='(SR) G.Conj ftol1e-6'
; show_exp_contrl,dir='exp_contrl_v2_ftol1e-6',/ruido,suffix_exp='(CR) G.Conj ftol1e-6'
; show_exp_contrl,dir='exp_contrl_v2_amoeba_ftol1e-6',suffix_exp='(SR) AMOEBA ftol1e-6'
; show_exp_contrl,dir='exp_contrl_v2_amoeba_ftol1e-6',/ruido,suffix_exp='(CR) AMOEBA ftol1e-6'
; show_exp_contrl,dir='exp_contrl_v2_ftol1e-1e-4_TOL1e-7',suffix_exp='(SR) G.Conj ftol1e-4 tol1e-7'
; show_exp_contrl,dir='exp_contrl_v2_ftol1e-1e-6_TOL1e-4',suffix_exp='(SR) G.Conj ftol1e-6 tol1e-4'
; show_exp_contrl,dir='exp_contrl_v2_ftol1e-1e-5_TOL1e-4',suffix_exp='(SR) G.Conj ftol1e-5 tol1e-4'
; show_exp_contrl,dir='exp_contrl_v2_amoeba_ftol1e-8',suffix_exp='(SR) AMOEBA ftol1e-8'
; show_exp_contrl,dir='exp_contrl_v2_amoeba_ftol1e-10',suffix_exp='(SR) AMOEBA ftol1e-10'
; show_exp_contrl,dir='exp_contrl_v2_ftol1e-8TOL1e-4',suffix_exp='(SR) G.Conj ftol1e-8 tol=1e-4'
; show_exp_contrl,dir='exp_contrl_v2_ftol1e-8TOL1e-6',suffix_exp='(SR) G.Conj ftol1e-8 tol=1e-6'
; show_exp_contrl,dir='exp_contrl_v2_ftol1e-8TOL1e-2',suffix_exp='(SR) G.Conj ftol1e-8 tol=1e-2'
; show_exp_contrl,dir='exp_contrl_v2_ftol1e-10TOL1e-4',suffix_exp='(SR) G.Conj ftol1e-10 tol=1e-4'

pro wrapper

  show_exp_contrl,dir='exp_contrl_v2_ftol1e-8TOL1e-6',suffix_exp='(CR 0.05) G.Conj ftol1e-8 tol=1e-6',file='mit_exp_contrl_con_ruido_0.05.out',fnoise_suffix='0.05',/ruido

  return
end

pro show_exp_contrl,ruido=ruido,dir=dir,file=file,suffix_exp=suffix_exp,fnoise_suffix=fnoise_suffix
  common units,ne_unit,te_unit

IF NOT keyword_set(suffix_exp) then suffix_exp=''
 
; (con o sin ruido)
  if     keyword_set(ruido) then  suffix = 'con_ruido'
  if NOT keyword_set(ruido) then  suffix = 'sin_ruido'

  if     keyword_set(fnoise_suffix) then noise_suffix = suffix+'_'+fnoise_suffix
  if not keyword_set(fnoise_suffix) then noise_suffix = suffix
  
  suffix_file = '_'+suffix
  suffix_title =' ('+suffix+')'

  
  print,'exp. controlado ',suffix
  print
  
; Directory where is exp_contrl.out and 
; where the plots are saved.
  if not keyword_set(dir) then STOP
  if     keyword_set(dir) then dir  = '/data1/DATA/MIT/'+dir+'/'+suffix+'/'

; Store the data in memory 
  ;restore,dir+'exp_contr.out'
  if not keyword_set(file) then file = 'mit_exp_contrl.out'
  restore,dir+file

  load_units

; IN and OUT arrays of Nm, Tm, FIP, sigN, sigT, q  
 ;npar=(size(par_in))(4)
  nm_in   =reform(par_in  (*,*,*,0))/ne_unit
  nm_out  =reform(par_out (*,*,*,0))
  fip_in  =reform(par_in  (*,*,*,1))
  fip_out =reform(par_out (*,*,*,1))
  Tm_in   =reform(par_in  (*,*,*,2))/te_unit
  Tm_out  =reform(par_out (*,*,*,2))
  sigT_in =reform(par_in  (*,*,*,3))/te_unit
  sigT_out=reform(par_out (*,*,*,3))
  sigN_in =reform(par_in  (*,*,*,4))/ne_unit 
  sigN_out=reform(par_out (*,*,*,4))
  q_in    =reform(par_in  (*,*,*,5))
  q_out   =reform(par_out (*,*,*,5))

  rel_diff_Nm   = (Nm_out-Nm_in)/Nm_in
  rel_diff_fip  = (fip_out-fip_in)/fip_in
  rel_diff_Tm   = (Tm_out-Tm_in)/Tm_in
  rel_diff_sigT = (sigT_out-sigT_in)/sigT_in
  rel_diff_sigN = (sigN_out-sigN_in)/sigN_in
  rel_diff_q    = (q_out-q_in)/q_in

  ;p = where(abs(rel_diff_q) gt 0.5)
  ;p2 = where(abs(rel_diff_q) gt 0.5 and abs(q_out) gt 0.95)

  ;index = where (scoreR lt 0.05)
  index = where (scoreR gt 0. and Nm_out gt 0.) ; this index select all the experiments
  print,'N=',n_elements(index)
  print

  plot_scatter_and_historatio,nm_in(index),nm_out(index),xsuffix='modeled',ysuffix='reconstructed',filename='Nm_'+noise_suffix,titulo='Nm '+suffix_exp,dir=dir
  plot_scatter_and_historatio,tm_in(index),tm_out(index),xsuffix='modeled',ysuffix='reconstructed',filename='Tm_'+noise_suffix,titulo='Tm '+suffix_exp,dir=dir
  plot_scatter_and_historatio,fip_in(index),fip_out(index),xsuffix='modeled',ysuffix='reconstructed',filename='FIP_'+noise_suffix,titulo='FIP '+suffix_exp,dir=dir
  plot_scatter_and_historatio,sigN_in(index),sigN_out(index),xsuffix='modeled',ysuffix='reconstructed',filename='sigN_'+noise_suffix,titulo='sigN '+suffix_exp,dir=dir
  plot_scatter_and_historatio,sigT_in(index),sigT_out(index),xsuffix='modeled',ysuffix='reconstructed',filename='sigT_'+noise_suffix,titulo='sigT '+suffix_exp,dir=dir
  plot_scatter_and_historatio,q_in(index),q_out(index),xsuffix='modeled',ysuffix='reconstructed',filename='q_'+noise_suffix,titulo='q '+suffix_exp,dir=dir
  ;return
  ;index = where (abs((sigN_in -sigN_out)/sigN_in) gt 0.5)
  ;index = where( abs(sigT_in/Tm_in - sigN_in/Nm_in) lt 1.e-1  )
  ;stop


; Plots settings
  !PATH = Expand_Path('+~/idlfiles/coyote/') + ':' + !PATH
  !p.background=255
  !p.color=0

  goto,skip_scatter_R
  ; SCATTER-PLOTS: Rel Diff. VS score R
  plot,abs(rel_diff_Nm),scoreR,psym=4,xtitle='rel diff.',ytitle='score R',title='Nm '+suffix_exp,/xlog,/ylog,charsize=1.5
  record_jpg,dir,'Nm_scatter'+suffix_title+'.jpg'
  plot,abs(rel_diff_Tm),scoreR,psym=4,xtitle='rel diff.',ytitle='score R',title='Tm '+suffix_exp,/xlog,/ylog,charsize=1.5
  record_jpg,dir,'Tm_scatter'+suffix_file+'.jpg'
  plot,abs(rel_diff_fip),scoreR,psym=4,xtitle='rel diff.',ytitle='score R',title='FIP '+suffix_exp,/xlog,/ylog,charsize=1.5
  record_jpg,dir,'FIP_scatter'+suffix_file+'.jpg'
  plot,abs(rel_diff_sigT),scoreR,psym=4,xtitle='rel diff.',ytitle='score R',title='sigT '+suffix_exp,/xlog,/ylog,charsize=1.5
  record_jpg,dir,'sigT_scatter'+suffix_file+'.jpg'
  plot,abs(rel_diff_sigN),scoreR,psym=4,xtitle='rel diff.',ytitle='score R',title='sigN '+suffix_exp,/xlog,/ylog,charsize=1.5
  record_jpg,dir,'sigN_scatter'+suffix_file+'.jpg'
  plot,abs(rel_diff_q),scoreR,psym=4,xtitle='rel diff.',ytitle='score R',title='q '+suffix_exp,/xlog,/ylog,charsize=1.5
  record_jpg,dir,'q_scatter'+suffix_file+'.jpg'
  skip_scatter_R:

 openw,1,dir+'tabla.txt'
 printf,1,'======================'
 printf,1,suffix_exp
 printf,1,'======================'
 printf,1,''
; =============== Nm ================== 
  print,'Nm comparison'
  printf,1,'Nm comparison'
  rel_diff = rel_diff_Nm(index) 
  print_rel_diff,rel_diff,0.1
  print_rel_diff,rel_diff,0.2
  print_rel_diff,rel_diff,0.3
  ;window,xs=500,ys=500
  ;cghistoplot,rel_diff,title='d_Nm '+suffix_exp,xrange=[-1,1],binsize=1.e-1
  ;record_jpg,dir,'Nm'+suffix_file+'.jpg'
  print,'--------------------------------------------------------------'
  printf,1,'--------------------------------------------------------------'
  print
  printf,1,''
 
 

; =============== Tm ================== 
  print,'Tm comparison'
  printf,1,'Tm comparison'
  rel_diff = rel_diff_Tm(index)
  print_rel_diff,rel_diff,0.1
  print_rel_diff,rel_diff,0.2
  print_rel_diff,rel_diff,0.3
  ;window,xs=500,ys=500
  ;cghistoplot,rel_diff,title='d_Tm '+suffix_exp,xrange=[-1,1],binsize=1.e-1
  ;record_jpg,dir,'Tm'+suffix_file+'.jpg'
  print,'--------------------------------------------------------------'
  printf,1,'--------------------------------------------------------------'
  print
  printf,1,''

; =============== FIP ================== 
  print,'fip comparison'
  printf,1,'fip comparison'
  rel_diff = rel_diff_fip(index)
  print_rel_diff,rel_diff,0.1
  print_rel_diff,rel_diff,0.2
  print_rel_diff,rel_diff,0.3
  ;window,xs=500,ys=500
  ;cghistoplot,rel_diff,title='d_FIP '+suffix_exp,xrange=[-2,2],binsize=1.e-1
  ;record_jpg,dir,'fip_factor'+suffix_file+'.jpg'
  print,'--------------------------------------------------------------'
  printf,1,'--------------------------------------------------------------'
  print
  printf,1,''

; =============== sigN ================== 
   print,'sigN comparison'
   printf,1,'sigN comparison'
   rel_diff = rel_diff_sigN(index)
   print_rel_diff,rel_diff,0.1
   print_rel_diff,rel_diff,0.2
   print_rel_diff,rel_diff,0.3
   ;window,xs=500,ys=500
   ;cghistoplot,rel_diff,title='d_sigN '+suffix_exp,xrange=[-2,2],binsize=1.e-1
   ;record_jpg,dir,'sigN'+suffix_file+'.jpg'
   print   ,'--------------------------------------------------------------'
   printf,1,'--------------------------------------------------------------'
   print
   printf,1,''
   
; =============== sigT ================== 
   print,'sigT comparison'
   printf,1,'sigT comparison'
   rel_diff = rel_diff_sigT(index)
   print_rel_diff,rel_diff,0.1
   print_rel_diff,rel_diff,0.2
   print_rel_diff,rel_diff,0.3
   ;window,xs=500,ys=500
   ;cghistoplot,rel_diff,title='d_sigT '+suffix_exp,xrange=[-2,2],binsize=1.e-1
   ;record_jpg,dir,'sigT'+suffix_file+'.jpg'
   print,'--------------------------------------------------------------'
   printf,1,'--------------------------------------------------------------'
   print
   printf,1,''
   

; =============== q ================== 
   print,'q comparison'
   printf,1,'q comparison'
   rel_diff = rel_diff_q(index)
   print_rel_diff,rel_diff,0.1
   print_rel_diff,rel_diff,0.2
   print_rel_diff,rel_diff,0.3
   ;window,xs=500,ys=500
   ;cghistoplot,rel_diff,title='d_q '+suffix_exp,xrange=[-2,2],binsize=1.e-1
   ;record_jpg,dir,'q'+suffix_file+'.jpg'
   print,'--------------------------------------------------------------'
   printf,1,'--------------------------------------------------------------'
   print
   printf,1,''

   close,1

  return
end


pro print_rel_diff,rel_diff,frac


  print,'--------------------------------------------------------------'
  printf,1,'--------------------------------------------------------------'
  index = where (abs(rel_diff) lt frac)
  print,'fraction of voxel with rel. diff. within ',frac*100,' %:',$
     float(n_elements(index))/n_elements(rel_diff)

  printf,1,'fraction of voxel with rel. diff. within ',frac*100,' %:',$
         float(n_elements(index))/n_elements(rel_diff)


  return
end

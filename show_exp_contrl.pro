pro show_exp_contrl


  suffix = 'sin_ruido'

  dir = '~/Downloads/exp_contrl/'+suffix+'/'
  suffix_file = '_'+suffix
  restore,dir+'exp_contr.out'
  suffix_title =' ('+suffix+')'

  npar=(size(par_in))(4)
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
 

  !PATH = Expand_Path('+~/idlfiles/coyote/') + ':' + !PATH
  !p.background=255
  !p.color=0
  window,xs=500,ys=500
  cghistoplot,(nm_out-nm_in)/nm_in,title='Nm rel diff.'+suffix_title,xrange=[-1,1],binsize=1.e-1
  record_jpg,dir,'Nm'+suffix_file+'.jpg'
  window,xs=500,ys=500
  cghistoplot,(tm_out-tm_in)/tm_in,title='Tm rel diff.'+suffix_title,xrange=[-1,1],binsize=1.e-1
  record_jpg,dir,'Tm'+suffix_file+'.jpg'
  window,xs=500,ys=500
  cghistoplot,(fip_out-fip_in)/fip_in,title='FIP rel diff.'+suffix_title,xrange=[-2,2],binsize=1.e-1
  record_jpg,dir,'fip_factor'+suffix_file+'.jpg'
  window,xs=500,ys=500
  cghistoplot,(sigT_out-sigT_in)/sigT_in,title='sigT rel diff.'+suffix_title,xrange=[-2,2],binsize=1.e-1
  record_jpg,dir,'sigT'+suffix_file+'.jpg'
  window,xs=500,ys=500
  cghistoplot,(sigN_out-sigN_in)/sigN_in,title='sigN rel diff'+suffix_title,xrange=[-2,2],binsize=1.e-1
  record_jpg,dir,'sigN'+suffix_file+'.jpg'
  window,xs=500,ys=500
  cghistoplot,(q_out-q_in)/q_in,title='q rel diff.'+suffix_title,xrange=[-2,2],binsize=1.e-1
  record_jpg,dir,'q'+suffix_file+'.jpg'
 

  return
end

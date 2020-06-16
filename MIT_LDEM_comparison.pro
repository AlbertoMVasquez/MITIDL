pro MIT_LDEM_comparison
  
  restore,'~/Downloads/param.out'
  restore,'~/Downloads/ldem_mit.out'
  ndat=(size(demt_A))(4)
  nmdem=reform(demt_A(*,*,*,ndat-4))
  tmdem=reform(demt_A(*,*,*,ndat-3))
  wtdem=reform(demt_A(*,*,*,ndat-2))
  !PATH = Expand_Path('+~/idlfiles/coyote/') + ':' + !PATH
  !p.background=255
  !p.color=0
  window,xs=500,ys=500
  cghistoplot,nmdem/nma,title='Nm ratio DEMT/MIT'
  record_jpg,'~/Downloads/','Nm.jpg'
  window,xs=500,ys=500
  cghistoplot,tmdem/tma,title='Tm ratio DEMT/MIT'
  record_jpg,'~/Downloads/','Tm.jpg'
  window,xs=500,ys=500
  cghistoplot,wtdem/sta,title='WT/sigma_T'
  record_jpg,'~/Downloads/','WT.jpg'
  return
end

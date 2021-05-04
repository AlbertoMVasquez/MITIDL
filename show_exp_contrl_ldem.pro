pro show_exp_contrl_ldem
common comunes,tm,wt,nband,demc,PHI,parametrizacion,Tmin,Tmax,nr,nth,np,rad,lat,lon,lambda,WTc


  dir_in = '/data1/tomography/bindata/'
  ;file_in= 'param_input_exp_contrl_ldem_sin_ruido.out'
  file_in = 'param_input_exp_contrl_ldem_con_ruido_0.1.out'
  restore,dir_in+file_in
  
  Ne0_in  = par_IN(*,*,*,0)
  Tc_in   = par_IN(*,*,*,1)
  sigT_in = par_IN(*,*,*,2)

  ;file_out = 'LDEM_AIA3_exp_contrl_sin_ruido'
  file_out = 'LDEM_AIA3_exp_contrl_con_ruido_0.1'
  read_ldem,file_out,/ldem,/gauss1

   Ne0_out  = sqrt(demc(*,*,*)*lambda(*,*,*,0))
   Tc_out   = lambda(*,*,*,1)*Tmin/1.e6
   sigT_out = lambda(*,*,*,2)*Tmin/1.e6

   print,'Ne0 ratio out/in:',transpose(Ne0_out/Ne0_in)
   print,'Tc ratio out/in:',transpose(Tc_out/Tc_in)
   print,'sigT ratio out/in:',transpose(sigT_out/sigT_in)

   dir = '/data1/Downloads/'
   
   plot_scatter_and_historatio,ne0_in,ne0_out,xsuffix='modeled',ysuffix='reconstructed',filename='Ne0_con_ruido',titulo='Ne0 [cm-3] (CR 0.1)',dir=dir
   plot_scatter_and_historatio,Tc_in,Tc_out,xsuffix='modeled',ysuffix='reconstructed',filename='Tc_con_ruido',titulo='Tc [MK] (CR 0.1)',dir=dir
   plot_scatter_and_historatio,sigT_in,sigT_out,xsuffix='modeled',ysuffix='reconstructed',filename='sigT_con_ruido',titulo='SigT [MK] (CR 0.1)',dir=dir
  return
end

pro load_demt,file_demt,nm_demt_array,tm_demt_array,wt_demt_array
  common comunes,tm,wt,nband,demc,PHI,parametrizacion,Tmin,Tmax,nr,nth,np,rad,lat,lon,lambda,WTc
  common results_tomo,tfbe,sfbe,N_e
  
    read_ldem,file_demt,/ldem,/gauss1
    nm_demt_array = N_e;/1.e8
    tm_demt_array = Tm ;/1.e6
    wt_demt_array = WT ;/1.e6

  return
end

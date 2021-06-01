pro read_two_images
  dir1='/data1/tomography/DATA/comp/1074/CR2198/Full_Data/20171203.comp.1074.daily_dynamics/'
  dir2='/data1/tomography/DATA/comp/1074/CR2198/Full_Data/20171203.comp.1074.daily_dynamics.ALBERT/'
  file1='20171203.184938.comp.1074.dynamics_avg_total_intensity_t017.8-Dt2.00.fts'
  file2='20171203.184938.comp.1074.dynamics_avg_total_intensity_t017.8-Dt2.00.fts'
  
  mreadfits,dir1+file1,hdr1,img1
  mreadfits,dir2+file2,hdr2,img2

  p1 = where(img1 gt 0)
  p2 = where(img2 gt 0)

  print,'ima 1:',mean(img1(p1))
  print,'ima 2:',mean(img2(p2))
  STOP
  return
end

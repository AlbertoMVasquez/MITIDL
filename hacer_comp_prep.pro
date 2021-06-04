
   ; este script ejecuta comp_prep > Necesita SolSoft
   ; IDL> 
   ;  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH
   ;  .r comp_prep
   ;  .r hacer_comp_prep

pro hacer_comp_prep
  r0=[1.1,1.3]
  ;comp_line = '1074' 
  comp_line = '1079' 
  year = '2017'
  month= '12'
  day_v = ['03','04','05','07','08','09','10','11','12','13','14','15','16']
  ndays = n_elements(day_v)
  dir_out   = '/data1/tomography/DATA/comp/'+comp_line+'/CR2198/PREP/'
  ; loop para aplicar avg_compute_dynamics 
  ; 03-12-2017 > 16-12-2017
  for i=0,ndays-1 do begin
     day  = day_v[i]
     date      = year+month+day
     data_dir  = '/data1/tomography/DATA/comp/'+comp_line+'/CR2198/Full_Data/'+date+'.comp.'+comp_line+'.daily_dynamics/'  
     compute_avg_dynamics,data_dir=data_dir,file_list='list.txt',window_lapse=2.,init_hour=17.8,/dynamics,dir_out=dir_out,r0=r0
  endfor
  return
 end

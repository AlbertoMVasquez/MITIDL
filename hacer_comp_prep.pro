

   ; este script ejecuta comp_prep > Necesita SolSoft
   ; IDL> @hacer_comp_prep.pro

    !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH
    .r comp_prep
    
    ;comp_prep,data_dir='/data1/tomography/DATA/comp/1074/CR2198/',file_list='list_test.txt',r0=[1.05,1.07],/dynamics
 
    ; En realidad, esta es la rutina de procesamiento principal...
    compute_avg_dynamics,data_dir='/data1/tomography/DATA/comp/1074/CR2198/Full_Data/20171203.comp.1074.daily_dynamics/'$
                         ,file_list='list.txt',window_lapse=2.,init_hour=17.8,/dynamics,$
                         dir_out='/data1/tomography/DATA/comp/1074/CR2198/test/'

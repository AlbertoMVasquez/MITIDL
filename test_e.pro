pro test_e
  
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  
  r0=1.2
  fip_factor=1.
  Tem=1.35e6
  Nem=1.3e8
  SiTe=0.25e6
  SigNe=0.25e8
  q=0.5 

  Ne0_Limits = [1.e7, 1.e9]
  Te0_Limits = [1.e6, 3.e6]

  print, e_function( Ne0_Limits , Te0_Limits )

  return
end

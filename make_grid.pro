;-------------------------------------------------------
; Esta rutina crea una grilla fija y uniforme de Ne y
; Te para hacer las integrales dobles por CS.

;-------------------------------------------------------


pro make_grid

  common dimensions, NTe, NNe
  common NT_limits, Ne0_Limits, Te0_Limits
  common NT_arrays,Ne_array,Te_array

  nmin=Ne0_Limits(0)
  nmax=Ne0_Limits(1)
  tmin=te0_Limits(0)
  tmax=te0_Limits(1)
  
  dNe=(nmax -nmin)/NNe
  dTe=(Tmax -Tmin)/NTe

  Ne_array = Nmin + dNe* findgen(NNe) + dNe/2 
  Te_array = tmin + dTe* findgen(NTe) + dTe/2 


  return
end

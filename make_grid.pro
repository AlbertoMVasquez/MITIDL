;-------------------------------------------------------
; Esta rutina crea una grilla fija y uniforme de Ne y
; Te para hacer las integrales dobles por CS.

; HISTORY
; V1.0 F.A. Nuevo, IAFE, April-2020
; V1.1 A.M. Vasquez, IAFE, April-2020
;
;-------------------------------------------------------
pro make_grid

  common dimensions, NTe, NNe
  common NT_limits, Ne0_Limits, Te0_Limits
  common NT_arrays,Ne_array,Te_array,dNe_array,dTe_array
  common tables,Te1,Te2,Te3,Te4,Te5,Ne1,Ne2,Ne3,Ne4,Ne5,G1,G2,G3,G4,G5,r1,r2
 

 ; Create a fixed grid of Ne and Te
  NNe=100;200
  NTe=80;200

  print,'Created Grid:',NNe,'    x',NTe
  print
  
; Valores originales de Fede
; USAR ESTOS VALORES, VER TEST_INTEGRAL_LIMITS
  Ne0_Limits = [1.e6, 5.e9]
  Te0_Limits = [0.5e6,5.0e6]

; Valores dinámicos razonables:
;  Ne0_Limits = [max([min(Ne1),min(Ne2),min(Ne3),min(Ne4),min(Ne5)]),min([max(Ne1),max(Ne2),max(Ne3),max(Ne4),max(Ne5)])]
;  Te0_Limits = [max([min(Te1),min(Te2),min(Te3),min(Te4),min(Te5)]),min([max(Te1),max(Te2),max(Te3),max(Te4),max(Te5)])]

  Ne_min=Ne0_Limits(0)
  Ne_max=Ne0_Limits(1)
  Te_min=Te0_Limits(0)
  Te_max=Te0_Limits(1)

  Negrid = Ne_min + (Ne_max-Ne_min) * findgen(NNe+1)/float(NNe)
  Tegrid = Te_min + (Te_max-Te_min) * findgen(NTe+1)/float(NTe)

  Ne_array  = (Negrid(0:NNe-1)+Negrid(1:NNe))/2.
  Te_array  = (Tegrid(0:NTe-1)+Tegrid(1:NTe))/2.
  dNe_array = Negrid(1:NNe)-Negrid(0:NNe-1)
  dTe_array = Tegrid(1:NTe)-Tegrid(0:NTe-1)

  ;dNe=(Ne_max - Ne_min)/NNe
  ;dTe=(Te_max - Te_min)/NTe
  ;Ne_array = Ne_min + dNe/2. + dNe*findgen(NNe)
  ;Te_array = Te_min + dTe/2. + dTe*findgen(NTe)
  return
end

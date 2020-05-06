;-------------------------------------------------------
; Create a fixed grid of NNe and NTe values of Ne and Te. 
; To be used for Riemann cuadratures approach.
;
; HISTORY
; V1.0 F.A. Nuevo, IAFE, April-2020
; V1.1 A.M. Vasquez, IAFE, April-2020
;           Set dynamical grid limits.
; V1.2 A.M. Vasquez, IAFE, April-2020
;           Allowed for uniform or log-uniform grids.
;
;-------------------------------------------------------
pro make_grid,uniform=uniform,loguniform=loguniform,NNe_provided=NNe_provided,NTe_provided=NTe_provided

  common dimensions, NTe, NNe
  common NT_limits, Ne0_Limits, Te0_Limits
  common NT_arrays,Ne_array,Te_array,dNe_array,dTe_array,dTN
  common tables,Te1,Te2,Te3,Te4,Te5,Ne1,Ne2,Ne3,Ne4,Ne5,G1,G2,G3,G4,G5,r1,r2

  if not keyword_set(NNe_provided) then NNe_provided = 100
  if not keyword_set(NTe_provided) then NTe_provided = 100

  NNe=NNe_provided
  NTe=NTe_provided

  if keyword_set(   uniform) then print,'Created     Uniform Grid NTe x NNe = :',NTe,'    x',NNe
  if keyword_set(loguniform) then print,'Created Log-Uniform Grid NTe x NNe = :',NTe,'    x',NNe
  print
  
; Limits for the grid. These are DYNAMICAL as to adjust to future changes in G-tables:
  Ne0_Limits = [max([min(Ne1),min(Ne2),min(Ne3),min(Ne4),min(Ne5)]),min([max(Ne1),max(Ne2),max(Ne3),max(Ne4),max(Ne5)])]
  Te0_Limits = [max([min(Te1),min(Te2),min(Te3),min(Te4),min(Te5)]),min([max(Te1),max(Te2),max(Te3),max(Te4),max(Te5)])]
  
  Ne_min=Ne0_Limits(0)
  Ne_max=Ne0_Limits(1)
  Te_min=Te0_Limits(0)
  Te_max=Te0_Limits(1)

; Arrays Ne_grid and Te_grid of Ne and Te values at NNe+1 and NTe+1 CELL-BOUNDARIES (the grid):
; A) Uniform grid:
  if keyword_set(uniform) then begin
     Ne_grid = Ne_min + (Ne_max-Ne_min) * findgen(NNe+1)/float(NNe)
     Te_grid = Te_min + (Te_max-Te_min) * findgen(NTe+1)/float(NTe)
  endif
; B) Log-Uniform grid:
  if keyword_set(loguniform) then begin
     log10_Ne_grid = alog10(Ne_min) + (alog10(Ne_max)-alog10(Ne_min)) * findgen(NNe+1)/float(NNe)
     log10_Te_grid = alog10(Te_min) + (alog10(Te_max)-alog10(Te_min)) * findgen(NTe+1)/float(NTe)
     Ne_grid       = 10.^log10_Ne_grid
     Te_grid       = 10.^log10_Te_grid
  endif

; 1D-Arrays Ne_array and Te_array of Ne and Te values at NNe and NTe CELL-CENTERS:
  Ne_array  = (Ne_grid(0:NNe-1)+Ne_grid(1:NNe))/2.
  Te_array  = (Te_grid(0:NTe-1)+Te_grid(1:NTe))/2.

; 1D-Arrays dNe_array and dTe_array of NNe and NTe CELL-WIDTHS:
  dNe_array = Ne_grid(1:NNe)-Ne_grid(0:NNe-1)
  dTe_array = Te_grid(1:NTe)-Te_grid(0:NTe-1)

; 2D-Array dTN of NTe x NNe dimensions with elements
; dTN(i,j) = dTe_array(i)*dNe_array(j):
  dTN = dTe_array#dNe_array     ; array NTe x NNe (col x row)
  
  return
end

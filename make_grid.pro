;-------------------------------------------------------
; Create a fixed grid of NNe and NTe values of Ne and Te. 
; To be used for Riemann cuadratures approach.
;
; KEYWORDS;
; uniform : if keyword set the grid is uniform in Ne and Te
; loguniform: if keyword set the grid is uniform in log10Ne and
; log10Te
; lnuniform: if keyword set the grid is uniform in lnNe and lnTe.
; Also, use the Jacobian of the transformation to calculate dNe_array
; and dTe_array.
; NNe_provided: number of points in the Ne grid
; NTe_provided: number of points in the Te grid

; HISTORY
; V1.0 F.A. Nuevo, IAFE, April-2020
; V1.1 A.M. Vasquez, IAFE, April-2020
;           Set dynamical grid limits.
; V1.2 A.M. Vasquez, IAFE, May-2020
;           Allowed for uniform or log-uniform grids.
; V1.3 F.A. Nuevo, IAFE, May-2020
;           Add ln-uniform grid and use the Jacobian of
;           transformation to calculate dNe_array and dTe_array
;-------------------------------------------------------
pro make_grid,uniform=uniform,$
              lnuniform=lnuniform,$
              loguniform=loguniform,$
              NNe_provided=NNe_provided,NTe_provided=NTe_provided

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
  if keyword_set( lnuniform) then print,'Created  Ln-Uniform Grid NTe x NNe = :',NTe,'    x',NNe
  print
    
  Ne_min=Ne0_Limits(0)
  Ne_max=Ne0_Limits(1)
  Te_min=Te0_Limits(0)
  Te_max=Te0_Limits(1)

; Arrays Ne_grid and Te_grid of Ne and Te values at NNe+1 and NTe+1 CELL-BOUNDARIES (the grid):
; A) Uniform grid:
  if keyword_set(uniform) then begin
     dNe     = (Ne_max-Ne_min)/float(NNe) 
     dTe     = (Te_max-Te_min)/float(NTe) 
     Ne_grid = Ne_min + dNe * findgen(NNe+1)
     Te_grid = Te_min + dTe * findgen(NTe+1)
  endif
; B) Log-Uniform grid:
  if keyword_set(loguniform) then begin
     dlog10Ne      = (alog10(Ne_max)-alog10(Ne_min))/float(NNe)
     dlog10Te      = (alog10(Te_max)-alog10(Te_min))/float(NTe)
     log10_Ne_grid = alog10(Ne_min) + dlog10Ne * findgen(NNe+1)
     log10_Te_grid = alog10(Te_min) + dlog10Te * findgen(NTe+1)
     Ne_grid       = 10.^log10_Ne_grid
     Te_grid       = 10.^log10_Te_grid
  endif
; C) Ln-Uniform grid:
  if keyword_set(lnuniform) then begin
     dlnNe      = (alog(Ne_max)-alog(Ne_min)) /float(NNe)
     dlnTe      = (alog(Te_max)-alog(Te_min)) /float(NTe)
     ln_Ne_grid = alog(Ne_min) + dlnNe * findgen(NNe+1)
     ln_Te_grid = alog(Te_min) + dlnTe * findgen(NTe+1)
     Ne_grid    = exp(ln_Ne_grid)
     Te_grid    = exp(ln_Te_grid)
  endif


; 1D-Arrays Ne_array and Te_array of Ne and Te values at NNe and NTe CELL-CENTERS:
  Ne_array  = (Ne_grid(0:NNe-1)+Ne_grid(1:NNe))/2.
  Te_array  = (Te_grid(0:NTe-1)+Te_grid(1:NTe))/2.

  if not keyword_set(lnuniform) then begin
; 1D-Arrays dNe_array and dTe_array of NNe and NTe CELL-WIDTHS:
     dNe_array = Ne_grid(1:NNe)-Ne_grid(0:NNe-1)
     dTe_array = Te_grid(1:NTe)-Te_grid(0:NTe-1)
; 2D-Array dTN of NTe x NNe dimensions with elements
; dTN(i,j) = dTe_array(i)*dNe_array(j):
     dTN = dTe_array#dNe_array  ; array NTe x NNe (col x row)
  endif

; dTN is calculated using the jacobian of transformation: 
; d(lnNe) = 1/Ne * dNe 
; d(lnTe) = 1/Te * dTe
  if keyword_set(lnuniform) then $
  dTN = (Te_array#Ne_array) * dlnNe * dlnTe
  
  return
end

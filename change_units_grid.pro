
pro change_units_grid,direction=direction
; IF direction= 1:  cm-3, K      > 10^8cm-3, MK
; IF direction=-1:  10^8cm-3, MK > cm-3, K
  common NT_arrays,Ne_array,Te_array,dNe_array,dTe_array,dTN
  common units,ne_unit,te_unit
  
  if not keyword_set(direction) then direction=1
; -----------------------------------------------------
  if direction eq 1 then begin
   ; Ne and Te grid 
     Ne_array = Ne_array /ne_unit
     Te_array = Te_array /te_unit
     dNe_array= dNe_array/ne_unit
     dTe_array= dTe_array/te_unit
     dTN      = dTN/te_unit/ne_unit
  endif

  if direction eq -1 then begin
   ; Ne and Te grid 
     Ne_array = Ne_array*ne_unit
     Te_array = Te_array*te_unit
     dNe_array= dNe_array*ne_unit
     dTe_array= dTe_array*te_unit
     dTN      = dTN*te_unit*ne_unit
  endif
; -----------------------------------------------------
  return
end


pro change_units_grid,direction=direction
; IF direction= 1:  cm-3, K      > 10^8cm-3, MK
; IF direction=-1:  10^8cm-3, MK > cm-3, K
  common NT_arrays,Ne_array,Te_array,dNe_array,dTe_array,dTN
  
  if not keyword_set(direction) then direction=1
; -----------------------------------------------------
  if direction eq 1 then begin
   ; Ne and Te grid 
     Ne_array = Ne_array /1.d8
     Te_array = Te_array /1.d6 
     dNe_array= dNe_array/1.d8
     dTe_array= dTe_array/1.d6
     dTN      = dTN/1.d6 /1.d8
  endif

  if direction eq -1 then begin
   ; Ne and Te grid 
     Ne_array = Ne_array*1.d8
     Te_array = Te_array*1.d6 
     dNe_array= dNe_array*1.d8
     dTe_array= dTe_array*1.d6
     dTN      = dTN*1.d6*1.d8
  endif
; -----------------------------------------------------
  return
end

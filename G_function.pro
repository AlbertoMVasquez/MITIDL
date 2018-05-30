function G_function, Te0, Ne0, r0 ; this is the ordering of the indexes in G table.
  common G_table,G,T_e,N_e,r,photT

  goto,tri_linear_interpolator
  
  ; Simply assign closest value, I will next implement a 3-linear interpolator
  fTe = abs(Te0-T_e) & iTe = median(where(fTe eq min(fTe)))
  fNe = abs(Ne0-N_e) & iNe = median(where(fNe eq min(fNe)))
  fr  = abs(r0 -r  ) & ir  = median(where(fr  eq min(fr )))
  if iTe eq -1 or iNe eq -1 or ir eq -1 then begin
     print,'Can not assign value to G_function'
     stop
  endif
  G_value = G(iTe,iNe,ir)
  goto,exit
  
  tri_linear_interpolator:
  ; Tri-linear interpolator
  ; Rename dimensions. X is Ne.
  xa = N_e & x0 = Ne0
  ya = T_e & y0 = Te0
  za = r   & z0 = r0
  ; Find the x-planes that surround x0
  ixA = max(where(xa le x0))
  if ixA eq -1 then begin
     print,'Ne value is out of range' ; Let's talk if we need something better here.
     stop
  endif
  ixB=ixA+1
  ; Define two 2D y-z arrays at x-planes that surround x0. 
  DATA_ARRAY_xA = reform(G(*,ixA,*))
  DATA_ARRAY_xB = reform(G(*,ixB,*))
  ; Bi-linearly interpolate in xA and xB planes.   
  G_xA = findval2D( DATA_ARRAY_xA ,ya ,za , y0, z0 )
  G_xB = findval2D( DATA_ARRAY_xB ,ya ,za , y0, z0 )
  ; Linearly interpolate the value og G along x, betwee xA and xB planes.
  G_value = G_xA + (G_xB-G_xA)*(x0-xa(ixA))/(xa(ixB)-xa(ixA))

  exit:
  return,G_value
end

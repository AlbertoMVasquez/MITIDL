;---------------------------------------------------------------------
; This routine make a grid in the variables-space and found the 
; the value in the grid the minimizes the cost function to use it
; as initial guess in the minimization. For Nm, Tm, sigN, and sigT
; uses the DEMT Nm, Tm and WT values. 
; Es significativamente mas rapida que make_guess_ini...
;
; OUTPUT:

; GUESS_INI: initial guess. 

; PHIGUESS: the value of the cost function in the initial guess.

; HISTORY:
; V1.0, Federico A. Nuevo, August, IAFE.

;---------------------------------------------------------------------


pro make_guess_ini_with_demt,nm_demt,tm_demt,wt_demt,guess_ini,PHIguess
  common tomographic_measurements, y0, y  

  print,'calculating initial guess...'
  tstart     = systime(/seconds)


  ; range of the variable-grid  
                                         n1=1
  fip_range         = [0.1, 5.0]       & n2=5 ; 0.5 > 5.
                                         n3=1
                                         n4=1
                                         n5=1
  q_range           = [0.1,0.8]        & n6=5

  ; array with the values of the cost function
  phiA = dblarr(n1,n2,n3,n4,n5,n6)
  
 ; variable-grid
  Nemv  = [y0];[Nm_demt]
  fipv  = fip_range (0) + (fip_range (1) -fip_range (0))*findgen(n2)/float(n2-1)
  Temv  = [Tm_demt]
  SigTv = [WT_demt]
  SigNv = [sqrt(abs(y0^2 - Nm_demt^2))]
     qv = q_range   (0) + (q_range   (1) -q_range   (0))*findgen(n6)/float(n6-1)
     
     
     for i1=0,n1-1 do begin
        for i2=0,n2-1 do begin
           for i3=0,n3-1 do begin
              for i4=0,n4-1 do begin
                 for i5=0,n5-1 do begin
                    for i6=0,n6-1 do begin
                       ; evaluates the cost function in the variable-grid
                       phiA(i1,i2,i3,i4,i5,i6) = cost_function_cs ([Nemv(i1), fipv(i2), Temv(i3), SigTv(i4), SigNv(i5), qv(i6)])
                    endfor
                 endfor
              endfor
           endfor
        endfor
     endfor
     
     ; found the minimum
     ii = median (where (phiA eq min(phiA)))
     PHIguess = phiA(ii)
     p  = ARRAY_INDICES(phiA, ii)
     ; initial guess
     guess_ini     = [Nemv(p(0)),fipv(p(1)),Temv(p(2)),sigTv(p(3)),sigNv(p(4)),qv(p(5))]

     t_elapsed  = systime(/seconds)-tstart
     print,'Elapsed time:',t_elapsed

     
     return
end

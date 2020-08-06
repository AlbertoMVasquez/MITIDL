;---------------------------------------------------------------------
; This routine make a grid in the variables-space and found the 
; the value in the grid the minimizes the cost function to use it
; as initial guess in the minimization.
;
; OUTPUT:

; GUESS_INI: initial guess. 

; PHIGUESS: the value of the cost function in the initial guess.

; HISTORY:
; V1.0, Federico A. Nuevo, June, IAFE.

;---------------------------------------------------------------------




pro make_guess_ini,guess_ini,PHIguess
  common tomographic_measurements, y0, y  

  print,'calculating initial guess...'
  tstart     = systime(/seconds)

  MK = 1.e6	
 
  ; de esta forma guess_ini tiene que coincidir
  ; con par_orig
  ;Nem_range         = [1.5 ,2.0]*1.e8  & n1=3
  ;fip_range         = [1.0, 1.2]       & n2=3
  ;Tem_range         = [1.0, 1.6]*MK    & n3=3
  ;SigT_range        = [0.1, 0.9]*MK    & n4=3
  ;sigN_range        = [0.1, 0.9]*1.e8  & n5=3
  ;q_range           = [0.1, 0.9]       & n6=3

  ; range of the variable-grid  
  Nem_range         = [0.25,4.0]*y0    & n1=5
  fip_range         = [0.1, 5.0]       & n2=5 ; 0.5 > 5.
  Tem_range         = [0.5, 2.5]*MK    & n3=5
  SigT_range        = [0.1, 1.0]*MK    & n4=5
  sigN_range        = [0.1,2.0]*y0     & n5=5
  q_range           = [0.1,0.8]        & n6=5

  ; array with the values of the cost function
  phiA = dblarr(n1,n2,n3,n4,n5,n6)
  
 ; variable-grid
  Nemv  = Nem_range (0) + (Nem_range (1) -Nem_range (0))*findgen(n1)/float(n1-1)
  fipv  = fip_range (0) + (fip_range (1) -fip_range (0))*findgen(n2)/float(n2-1)
  Temv  = Tem_range (0) + (Tem_range (1) -Tem_range (0))*findgen(n3)/float(n3-1)
  SigTv = SigT_range(0) + (sigT_range(1) -sigT_range(0))*findgen(n4)/float(n4-1)
  SigNv = SigN_range(0) + (sigN_range(1) -sigN_range(0))*findgen(n5)/float(n5-1)
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

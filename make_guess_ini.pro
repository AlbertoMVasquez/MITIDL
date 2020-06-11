pro make_guess_ini,y0,guess_ini
  
  
  MK = 1.e6	
  
  Nem_range         = [0.25,4.0]*y0    & n1=6
  fip_range         = [0.1, 10.]       & n2=6
  Tem_range         = [0.5, 2.5]*MK    & n3=6
  SigT_range        = [0.1, 1.0]*MK    & n4=6
  sigN_range        = [0.1,2.0]*y0     & n5=6
  q_range           = [0.1,0.9]        & n6=6
 
  phiA = dblarr(n1,n2,n3,n4,n5,n6)
  ;NemA = dblarr(n1,n2,n3,n4,n5,n6)
  ;fipA = dblarr(n1,n2,n3,n4,n5,n6)
  ;TemA = dblarr(n1,n2,n3,n4,n5,n6)
  ;SigTA= dblarr(n1,n2,n3,n4,n5,n6) 
  ;SigNA= dblarr(n1,n2,n3,n4,n5,n6)
  ;qA   = dblarr(n1,n2,n3,n4,n5,n6)

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
                       
                       phiA(i1,i2,i3,i4,i5,i6) = cost_function_cs ([Nemv(i1), fipv(i2), Temv(i3), SigTv(i4), SigNv(i5), qv(i6)])
                       ;NemA (i1,i2,i3,i4,i5,i6) = Nemv(i1)
                       ;fipA (i1,i2,i3,i4,i5,i6) = fipv(i2)
                       ;TemA (i1,i2,i3,i4,i5,i6) = Temv(i3)
                       ;SigTA(i1,i2,i3,i4,i5,i6) = SIGTv(i4)
                       ;SigNA(i1,i2,i3,i4,i5,i6) = SIGNv(i5)
                       ;qA   (i1,i2,i3,i4,i5,i6) = qv   (i6)
                    endfor
                 endfor
              endfor
           endfor
        endfor
     endfor

     ii = median (where (phiA eq min(phiA)))
     p  = ARRAY_INDICES(phiA, ii)
     ;guess_ini = [NemA(ii),fipA(ii),TemA(ii),sigTA(ii),sigNA(ii),qA(ii)]
     guess_ini = [Nemv(p(0)),fipv(p(1)),Temv(p(2)),sigTv(p(3)),sigNv(p(4)),qv(p(5))]

     return
end

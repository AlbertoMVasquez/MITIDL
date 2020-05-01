pro test_int2d
  
 Nx = 200
 Ny = 200

 xmin = 1.
 xmax = 5.

 ymin = 10.
 ymax = 30.

 ; Esta es la grilla de los BORDES de las celdas 
 ; son Nx y Ny puntos respectivamente
 ; Aquí lo hago UNIFORME, pero podría no serlo
 xgrid = xmin + (xmax-xmin) * findgen(Nx)/float(Nx-1)
 ygrid = ymin + (ymax-ymin) * findgen(Ny)/float(Ny-1)

 ; El resto del código NO ASUME uniformidad:

 ; Esta es la grilla de los CENTROS de las celdas
 ; Son Nx-1 y Ny-1 puntos respectivamente
 x_array = (xgrid(0:Nx-2)+xgrid(1:Nx-1))/2.
 y_array = (ygrid(0:Ny-2)+xgrid(1:Ny-1))/2.

 ; Estos son los arrays dx_array y dy_array
 ; Son Nx-1 y Ny-1 puntos respectivamente
 dx_array = xgrid(1:Nx-1)-xgrid(0:Nx-2)
 dy_array = ygrid(1:Ny-1)-ygrid(0:Ny-2)
 
 ; Notá como x_array, y_array, dx_array, dy_array tienen idénticas dimensiones

 ; Defino ahora una función 2D f_array, arbitraria, en los puntos x_array, y_array, 
 ; o sea que tiene dimensionalidad (Nx-1)x(Ny-1)

 f_array = dblarr(Nx-1,Ny-1)
 for ix=0,Nx-2 do $
    for iy=0,Ny-2 do $
       f_array(ix,iy) = cos((x_array(ix) + y_array(iy))^2)
 
 ; Y ahora calculo su Integral 2D sobre todo el dominio de x_array, y_array:

 tstart = systime(/seconds)
 suma2D = 0.d
 for ix=0,Nx-2 do $
 for iy=0,Ny-2 do $
    suma2D = suma2D + f_array(ix,iy) * dx_array(ix) * dy_array(iy)
 Tcalc_suma2D = systime(/seconds) - tstart
 
 ; Ahora lo calculo usando operadores de array de IDL:
 ; Referencia:   https://www.harrisgeospatial.com/docs/matrix_operators.html
 tstart = systime(/seconds)
 suma2D_cool = total(dy_array*(f_array#dx_array))
 Tcalc_suma2D_cool = systime(/seconds) - tstart

; verifico que dan igual y comparo tiempos de cálculo
 print,'suma2D     :',suma2D     ,'  calculada en',Tcalc_suma2D     ,'[seg]'    
 print,'suma2D_cool:',suma2D_cool,'  calculada en',Tcalc_suma2D_cool,'[seg]'    
 print,'Cociente de Sumas (loops-sobre-#):',suma2D/suma2D_cool
 print,'Cociente de Tiempos de Cálculo de Sumas (loops-sobre-#):',Tcalc_suma2D/Tcalc_suma2D_cool

 return
 end


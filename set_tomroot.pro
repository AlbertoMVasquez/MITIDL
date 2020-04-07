;---------------------------------------------------------------------
;
; Brief description:
;
; Sets up tomroot directory for all routines.
;
; History:  V1.0, Alberto M. Vasquez, CLaSP, Spring-2018.
;
;---------------------------------------------------------------------

pro set_tomroot
  common directories,tomroot
  tomroot = '/data1/tomography_dev/'
 ;tomroot = '/data1/tomography/' ;04-02-20, modificado para usar en mi PC
  return
end

  

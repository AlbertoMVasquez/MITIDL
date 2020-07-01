;---------------------------------------------------------------
;  Evaluate function along a line (p. 413):
;---------------------------------------------------------------
function f1dim,x
  common f1com,pcom,xicom
;  if finite(PHI(pcom+x*xicom)) ne 1 then stop ; para detectar NaNs
  return,PHI(pcom+x*xicom)
end

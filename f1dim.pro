;---------------------------------------------------------------
;  Evaluate function along a line (p. 413):
;---------------------------------------------------------------
function f1dim,x
common f1com,pcom,xicom
return,PHI(pcom+x*xicom)
end

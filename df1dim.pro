;---------------------------------------------------------------
;  Evaluate derivative along a line (p. 417):
;---------------------------------------------------------------
function df1dim,x
common f1com,pcom,xicom
if finite(total( (gradPHI(pcom+x*xicom)) * xicom)) ne 1 then stop ; para detectar NaNs
return,total( (gradPHI(pcom+x*xicom)) * xicom)
end

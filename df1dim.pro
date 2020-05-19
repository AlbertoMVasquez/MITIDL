;---------------------------------------------------------------
;  Evaluate derivative along a line (p. 417):
;---------------------------------------------------------------
function df1dim,x
common f1com,pcom,xicom
return,total( (gradPHI(pcom+x*xicom)) * xicom)
end

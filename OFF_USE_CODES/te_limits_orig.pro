FUNCTION te_limits, Ne0
  common Ylimits,Y_Limits
  RETURN,Y_Limits+[Ne0^0,0.]
END

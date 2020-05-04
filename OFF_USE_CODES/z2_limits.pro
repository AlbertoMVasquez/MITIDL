function z2_limits, Z1
  common parameters, r0, fip_factor, Tem, Nem, SigTe, SigNe, q
  common Ylimits,Y_Limits
  return,(Y_Limits-q*Z1)/sqrt(1.-q^2)
end

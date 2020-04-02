
  dp  = 1.e-6 * p

   phi        = cost_function(p)
  dphi        = cost_function(p+dp) - phi
  dphi_grad   = total( grad_cost_function(p)*dp )

  print,' Delta(c) / c               = ', dp / p
  print,' Delta(Phi) / Phi(c)        = ', dphi / phi
  print,'gradPhi*Delta_c / Delta_Phi =' , dphi_grad / dphi




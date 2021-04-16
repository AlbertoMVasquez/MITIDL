function df_penalty_q,q
  
  big = 1000.d
  pow = 20.d

  RESULT = pow * big * q ^ (pow-1)

  return,RESULT
END

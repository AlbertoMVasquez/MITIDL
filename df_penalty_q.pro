function df_penalty_q,q
  
  big = 10.d
  pow = 20.d

  RESULT = pow * big * q ^ (pow-1)

  return,RESULT
END

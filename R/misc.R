Quan_Fry <- function(wv, Sal=35, Temp=26) {
  n0 = 1.31405
  n1 = 1.779E-4
  n2 = -1.05E-6
  n3 = 1.6E-8
  n4 = -2.02E-6
  n5 = 15.868
  n6 = 0.01155
  n7 = -0.00423
  n8 = -4382
  n9 = 1.1455E6
  n = n0 + (n1+n2*Temp+n3*Temp^2)*Sal + n4*Temp^2 +
    (n5 + n6*Sal + n7*Temp) / wv  + n8/wv^2 + n9/wv^3
  return(n)
}



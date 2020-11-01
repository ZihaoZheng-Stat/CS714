f = function(a,b){
  
  k = complex(real = a+2, imaginary = b)
  
  z = polyroot(c(1, -k, 1))
  
  m1 = Mod(z[1]); m2 = Mod(z[2])
  
  m = max(c(m1,m2))
  
  return(m<=1)
}

a = seq(-5,2, by = 0.01)

b = a

tt = matrix(NA, ncol = length(b), nrow = length(a))

for (i in 1:length(a)) {
  
  for (j in 1:length(b)) {
    
    tt[i,j] = f(a[i],b[j])
  }
}


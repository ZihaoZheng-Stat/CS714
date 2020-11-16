## I write a function for chebyshev grid and chebyshev matrix

cheb = function(N){
  
  ## first, the chebyshev grid
  
  x = cos(c(0:N)/N * pi)
  
  ## then, based on the chebyshev grid, I work on the chebyshev matrix
  
  w = matrix(NA, (N+1), (N+1))
  
  c = c(2,rep(1,N-1),2)
  
  ## the off-diagnal entry
  
  for (i in 1:(N+1)) {
    
    for (j in 1:(N+1)) {
      
      w[i,j] = (c[i]/c[j])*((-1)^(i+j))/(x[i]-x[j])
    }
    
  }
  
  ## modify the diagonal
  
  for (i in 2:N) {
    
    w[i,i] = -x[i]/(2*(1-x[i]^2))
    
  }
  ## the w_00 and w_NN
  
  w[1,1] = (2*N^2 + 1)/6; w[N+1,N+1] = -(2*N^2 + 1)/6
  
  
  return(list(cheb.grid = x, cheb.matrix = w))
  
  
}
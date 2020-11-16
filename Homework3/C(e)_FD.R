BB = seq(0.5,10,by = 0.5)

Error = NULL

for (kk in 1:length(BB)) {
  
  B = BB[kk]
  
  ## the new function
  
  f = function(x,y){return(sin(B*pi*x)*sin(B*pi*y))}
  
  Nx = 64; Nt = 2*Nx ## we have many choices for Nx and Nt later on
  
  dx = 1/Nx; dt = 1/Nt
  
  x <- y <- c(-Nx:Nx)/Nx
  
  ## I will always use U to save the answer, at each time (it is indeed U^n)
  
  U0 = matrix(0, ncol = Nx, nrow = Nx)
  
  U1 = U0
  
  for (j1 in 1:Nx) {
    
    for (j2 in 1:Nx) {
      
      U1[j1,j2] = f(x=dx*j1,y=dx*j2)*dt
    }
  }
  
  ## I will creat a function to update U
  
  updateU = function(U.pre, U.cur){
    
    ## U.pre: U^{n-1}
    ## U.cur: U^{n}
    ## want to output: U^{n+1}
    
    UU = U0
    
    for(j1 in 2:(Nx-1)){
      for(j2 in 2:(Nx-1)){
        UU[j1,j2] = (dt/dx)^2*(U.cur[j1-1,j2] + U.cur[j1+1,j2] + 
                                 U.cur[j1,j2-1] + U.cur[j1,j2+1] - 
                                 4*U.cur[j1,j2]) + 2*U.cur[j1,j2] - U.pre[j1,j2]
      }
    }
    
    return(list(U.pre = U.cur, U.cur = UU))
    
  }
  
  ## start iterate
  
  U.pre = U0; U.cur = U1
  
  u.update = updateU(U.pre, U.cur)
  
  for (i in 1:(round(1/dt))) {
    
    u.update = updateU(u.update$U.pre, u.update$U.cur)
  }
  
  uu = u.update$U.cur
  
  fstar = function(x,y,t){return(1/(sqrt(2)*B*pi)*sin(B*pi*x)*sin(B*pi*y)*sin(sqrt(2)*B*pi*t))}
  
  Utrue = U0
  
  for (j1 in 1:Nx) {
    
    for (j2 in 1:Nx) {
      
      Utrue[j1,j2] = fstar(x=dx*j1,y=dx*j2,t=1)
    }
  }
  
  Error[kk] =  max(abs(Utrue - uu))
  
  print(kk)
  
}


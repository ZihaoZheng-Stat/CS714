## first, give a function f for use later

f = function(x){return(exp(-400*(x-0.5)^2))}

Nx = 100; Nt = 2*Nx ## we have many choices for Nx and Nt later on

dx = 1/Nx; dt = 1/Nt

x <- y <- c(1:Nx)/Nx

## I will always use U to save the answer, at each time (it is indeed U^n)

U0 = matrix(0, ncol = Nx, nrow = Nx)

U1 = U0

for (j1 in 1:Nx) {
  
  for (j2 in 1:Nx) {
    
    U1[j1,j2] = f(dx*j1)*f(dx*j2)*dt
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

for (i in 1:1/Nt) {
  
  u.update = updateU(u.update$U.pre, u.update$U.cur)
}

#a = NULL
#a[[7]] = u.update$U.cur

#persp(x,y,u.update$U.cur)

## evaluate the error

x.true = a[[7]]

E = NULL

for (i in 1:6) {
  
  ok = seq(1,6400, by = 2^i)
  
  x.approx = x.true[ok,ok]
  
  e.matrix = x.approx - a[[7-i]]
  
  E[i] = max(e.matrix)
  
}

Dx = (1/3200)*(c(1,2,4,8,16,32))


plot(log(Dx), log(E), xlab = "log of dx", ylab = "log of error", type = "b")

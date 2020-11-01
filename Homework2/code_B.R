yfun = function(x){
  
  return(exp(-400*(x-0.5)^2))
  
}

Nstar = 10^4

xstar = c(0:Nstar)/Nstar

ystar = yfun(xstar)


l = 10; u = 200

N = l

notdone = TRUE

while (notdone) {
  
  x = c(0:N)/N
  
  b = approxfun(x, yfun(x))
  
  e = max(abs(ystar - b(xstar)))
  
  if(e > 0.01){
    l = N
    N = round(0.5*(N + u))
  }
  
  if(e < 0.01){
    u = N
    N = round(0.5*(N + l))
  }
  
  notdone = l < u - 1
  
  print(N)
  
}

plot(yfun, lwd = 2, ylab = "f(x)")

## demonstrate N = 6, 10, 20

N = 6

x = c(0:N)/N

b = approxfun(x, yfun(x))

curve(b, add = T, col = "red", lwd = 1)

N = 10

x = c(0:N)/N

b = approxfun(x, yfun(x))

curve(b, add = T, col = "blue", lwd = 1)

N = 20

x = c(0:N)/N

b = approxfun(x, yfun(x))

curve(b, add = T, col = "gold", lwd = 1)

legend("topleft", legend = c("True", "N = 6", "N = 10", "N = 20"),
       col = c("black", "red", "blue", "gold"), lwd = 2)

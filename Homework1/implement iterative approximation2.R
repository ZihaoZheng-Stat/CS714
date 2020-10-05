f = function(y){
  return(cos(2*pi*y))
}

f = function(y){
  
  if(cos(2*pi*y)>=0){
    return(1)
  }
  
  if(cos(2*pi*y)<0){
    return(-1)
  }
}


### set mm  and get the "correct" answer

mm = 4

M = 2^mm - 1; N = 2^mm - 1

hx = 1/(M+1); hy = 1/(N+1)

A = matrix(0, nrow = M, ncol = M)

A[1,1] = -2; A[1,2] = 1

A[M, (M-1)] = 1; A[M, M] = -2

for (i in 2:(M-1)) {
  
  A[i,i] = -2
  
  A[i,(i-1)] <- A[i, (i+1)] <- 1
}

### for the matrix B

B = matrix(0, nrow = N, ncol = N)

B[1,1] = -1; B[1,2] = 1

B[N, (N-1)] = 1; B[N, N] = -1

for (i in 2:(N-1)) {
  
  B[i,i] = -2
  
  B[i,(i-1)] <- B[i, (i+1)] <- 1
}

### two important identity matrix

Ix = diag(1, nrow = M)
Iy = diag(1, nrow = N)

aa = kronecker(A, Iy)

bb = kronecker(Ix, B)

### the linear system is the following

### Qx = c

Q = (1/hx^2)*aa + (1/hy^2)*bb

c = rep(0, M*N)

c[1:N] = -sapply(hy*c(1:N),f)

correctans = as.numeric(solve(Q) %*% c)

### Starting here, I evaluated the iteration property

Q1 = diag(diag(Q)); Q1_inv = diag(1/diag(Q))

Q2 = Q - Q1

x = rep(0, length(c)) ## an initial value

initial.error = sum((x - correctans)^2)

Error = NULL

It = 1000 ## a larger number

for (it in 1:It) {
  
  x = Q1_inv %*% (c - Q2 %*% x)
  
  Error[it] = sum((x - correctans)^2)
  
  print(it)
  
}

par(mfrow = c(2,2))
plot(c(1:It), Error,
     xlab = "Iteration step", ylab = "Error", 
     main = "hx = hy = 1/16")
plot(c(1:It), Error/initial.error,
     xlab = "Iteration step", ylab = "ratio between error and initial error",
     main = "hx = hy = 1/16", log = "y")
abline(h = 0.001, col = "red", lwd = 2)
abline(v = min(c(1:It)[(Error/initial.error)<=0.001]), col = "red", lwd = 2)


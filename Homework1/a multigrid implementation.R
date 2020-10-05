f = function(y){
  return(cos(2*pi*y))
}

### set mm  and get the "correct" answer

mm = 6

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

c[1:N] = -f(hy*c(1:N))

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

###########################################################################
###########################################################################

### multigrid approach

### step 1, take 3 iteration step on Qx = c

x = rep(0, length(c)) ## an initial value

Error2 = NULL

It = 3

for (it in 1:It) {
  
  x = Q1_inv %*% (c - Q2 %*% x)
  
  Error2[it] = sum((x - correctans)^2)
  
  print(it)
  
}

### step 2, compute the residual

r = c - Q %*% x

### step 3, Coarsen the residual

r.matrix = matrix(NA, ncol = M, nrow = N)

for (j in 1:M) {
  
  r.matrix[,j] = r[(N*(j-1)+1):(N*j)]
}

rownames(r.matrix) = c(1:N)*hy
colnames(r.matrix) = c(1:M)*hx

r2.matrix = r.matrix[seq(2,M, by = 2), seq(2,M,by = 2)]

r2 = NULL

for (j in 1:dim(r2.matrix)[1]) {
  
  r2 = c(r2, r2.matrix[,j])
  
}

### step 4, solve residual equation

mm = 5

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

Q.sub = (1/hx^2)*aa + (1/hy^2)*bb

### solve Q.sub x = -r2

x.sub = solve(Q.sub) %*% (-r2)

### step 5, Interpolate (x.sub) grid function back to the original grid

x.sub.matrix = matrix(NA, ncol = M, nrow = N)

for (j in 1:M) {
  
  x.sub.matrix[,j] = x.sub[(N*(j-1)+1):(N*j)]
}

rownames(x.sub.matrix) = c(1:N)*hy
colnames(x.sub.matrix) = c(1:M)*hx

mm = 6

M = 2^mm - 1; N = 2^mm - 1

hx = 1/(M+1); hy = 1/(N+1)

x.matrix = matrix(0, ncol = M, nrow = N)

rownames(x.matrix) = c(1:N)*hy
colnames(x.matrix) = c(1:M)*hx

x.matrix[colnames(x.matrix)%in%colnames(x.sub.matrix),
         colnames(x.matrix)%in%colnames(x.sub.matrix)] = x.sub.matrix

xx = NULL

for (j in 1:dim(x.matrix)[1]) {
  
  xx = c(xx, x.matrix[,j])
  
}

### step 6, get a better approximation

x.new = x - xx

sum((x - correctans)^2)
sum((x.new - correctans)^2)
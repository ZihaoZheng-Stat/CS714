mm = 7

f = function(y){
  return(cos(2*pi*y))
}


M = 2^mm - 1; N = 2^mm - 1

hx = 1/(M+1); hy = 1/(N+1)

## U is a matrix, where each row is u_i1, ..., u_iN

U  = matrix(0, nrow = M, ncol = N)

## The huge matrix Q of size MN * MN (I do not form this matrix)

## But I would like to form some pieces

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

c = rep(0, M*N)

c[1:N] = -f(hy*c(1:N))

### The diagnal term of matrix Q (diagnal block of size N*N)

diag_Q = (-2/hx^2)*diag(1, nrow = N, ncol = N) + B/hy^2

offdiag_Q = diag(1, nrow = N, ncol = N)/hx^2

### When using Jacobi algorithm, we need to consider the diagnal term of diagnal matrix diag_Q

Q1 = diag(diag_Q)

### Then I update the matrix U, using Jacobi algorithm

for (it in 1:20) {
  
  U.tmp = U
  
  U.tmp[1,] = (1/Q1)*(c[1:N] - ((diag_Q - Q1) %*% U[1,] + 
                                  offdiag_Q %*% U[2,]))
  
  for (i in 2:(M-1)) {
    
    U.tmp[i,] = (1/Q1)*(rep(0,N) - (offdiag_Q %*% U[(i-1),] + 
                                      (diag_Q - Q1) %*% U[i,] + 
                                      offdiag_Q %*% U[(i+1),]))
  }
  
  U.tmp[M,] = (1/Q1)*(rep(0,N) - ((diag_Q - Q1) %*% U[M,] + 
                                    offdiag_Q %*% U[(M-1),]))
  
  U = U.tmp
  
  print(it)
}



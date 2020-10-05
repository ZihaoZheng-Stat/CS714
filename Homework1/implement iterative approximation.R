
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

### set the list to store all the answer (with different grid)

Ans = NULL

H = NULL

for (mm in 2:6) {
  
  M = 2^mm - 1; N = 2^mm - 1
  
  hx = 1/(M+1); hy = 1/(N+1)
  
  H[mm] = hx
  
  ### for the matrix A, the diagnal is -2, upper diagnal is 1 and lower diagnal is also 1
  
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
  
  c[1:N] = -sapply(hy*c(1:N), f)
  
  ### this is by using the direct approach, 
  
  ans = as.numeric(solve(Q) %*% c)
  
  ### but we can still use Jacobbi method, Q = Q1 + Q2 where Q1 is the diagnal matrix
  
  #Q1 = diag(diag(Q)); Q1_inv = diag(1/diag(Q))
  
  #Q2 = Q - Q1
  
  #x = rep(0, length(c)) ## an initial value
  
  #It = 2000 ## maybe a larger number, but experiments show 2000 is far larger to be converged(will be discussed later)
  
  #for (it in 1:It) {
    
  #  x = Q1_inv %*% (c - Q2 %*% x)
    
  #  #print(it)
    
  #}
  
  #ans = x
  
  ### in order to compare the solution on a nested grid, get a matrix of ans
  
  ans.matrix = matrix(NA, ncol = M, nrow = N)
  
  for (j in 1:M) {
    
    ans.matrix[,j] = ans[(N*(j-1)+1):(N*j)]
  }
  
  rownames(ans.matrix) = c(1:N)*hy
  colnames(ans.matrix) = c(1:M)*hx
  
  Ans[[mm]] = ans.matrix
  
  print(mm)
  
}
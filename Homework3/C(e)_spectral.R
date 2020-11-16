B = 1

f = function(x,y){return(sin(B*pi*x)*sin(B*pi*y))}

Nx = 64; Nt = Nx^2 ## we have many choices for Nx and Nt later on

dx = 1/Nx; dt = 1/Nt

cheb.tmp = cheb(Nx)

x <- y <- cheb.tmp$cheb.grid

w = cheb.tmp$cheb.matrix

w = w %*% w ## for second derivative

## for laplacian, we need the outer product

L = kronecker(w,diag(1, nrow = Nx + 1)) + kronecker(diag(1, nrow = Nx + 1),w)

## I will always use U to save the answer, at each time (it is indeed U^n)

U0 = matrix(0, ncol = Nx+1, nrow = Nx+1)

U1 = U0

for (j1 in 1:Nx+1) {
  
  for (j2 in 1:Nx+1) {
    
    U1[j1,j2] = f(x=dx*j1,y=dx*j2)*dt
    
    ## will later adjust to fourth order accuracy
  }
}

## for further update, I need to work on a vectorized 

u.old = c(U0)
u = c(U1)

t.step = round(1/dt)

for (t in 1:20) {
  
  #u.new = 2*u - u.old + dt^2*dx^2*L%*%u
  
  u.new = 2*u - u.old + dt^2*L%*%u
  
  #u.new = 2*u - u.old + dt^2*L%*%u + dt^4*L%*%L%*%u/12
  
  u.old = u
  
  u = u.new
  
  print(t)
}

U.matrix = U0

for (i in 1:(Nx+1)) {
  
  U.matrix[i,] = u[((Nx+1)*(i-1)+1) : ((Nx+1)*i)]
  
}

fstar = function(x,y,t){return(1/(sqrt(2)*B*pi)*sin(B*pi*x)*sin(B*pi*y)*sin(sqrt(2)*B*pi*t))}

Utrue = U0

for (j1 in 1:Nx) {
  
  for (j2 in 1:Nx) {
    
    Utrue[j1,j2] = fstar(x=dx*j1,y=dx*j2,t=dt*20)
  }
}

a = abs(Utrue)
b = abs(U.matrix)
max(abs(a-b))
max(abs(Utrue - U.matrix))
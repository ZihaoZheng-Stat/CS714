N = 500

x = matrix(0, nrow = 100, ncol = 100)

x[1,1:3] = c(6,-4,1)
x[2,1:4] = c(-4,6,-4,1)
x[100,98:100] = c(1,-4,6)
x[99,97:100] = c(1,-4,6,-4)

for(i in 3:97){
  
  x[i,((i-2):(i+2))] = c(1,-4,6,-4,1)
}

a = eigen(x)$values; max(a);min(a)
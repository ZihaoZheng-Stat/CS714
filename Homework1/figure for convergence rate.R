### analyze the answer

load("E://Fall 2020//CS 714//homework 1//answer.RData")

### choose the finest one as the "correct answer", and check the error

### note that we are constructing the "nested grid"

#a1 = ans.matrix1

a1 = Ans[[2]]

#a2 = ans.matrix2

a2 = Ans[[3]]

#a3 = ans.matrix3

a3 = Ans[[4]]

#a4 = ans.matrix4

a4 = Ans[[5]]

#a5 = ans.matrix5

a5 = Ans[[6]]

#a6 = ans.matrix6

#a6 = Ans[[7]]

a11 = a6[rownames(a6) %in% rownames(a1), colnames(a6) %in% colnames(a1)]
t1 = max(abs(a1 - a11))

a21 = a6[rownames(a6) %in% rownames(a2), colnames(a6) %in% colnames(a2)]
t2 = max(abs(a2 - a21))

a31 = a6[rownames(a6) %in% rownames(a3), colnames(a6) %in% colnames(a3)]
t3 = max(abs(a3 - a31))

a41 = a6[rownames(a6) %in% rownames(a4), colnames(a6) %in% colnames(a4)]
t4 = max(abs(a4 - a41))

a51 = a6[rownames(a6) %in% rownames(a5), colnames(a6) %in% colnames(a5)]
t5 = max(abs(a5 - a51))

x = c(1/2, 1/4, 1/8, 1/16, 1/32) ## this is h
x = H[2:6]
y = c(t1, t2, t3, t4, t5) ## this is error

plot(log(x), log(y), xlab = "log of h", ylab = "log of max-norm error",
     type = "b")

### check the slope

xx = log(x)

yy = log(y)

lm(yy[1:4]~xx[1:4])
(log(y[5]) - log(y[4]))/(log(x[5]) - log(x[4]))

### plot the figure

plot(x,y, log = "xy", type = "b",
     xlab = "h", ylab = "error")
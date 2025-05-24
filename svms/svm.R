rm(list = ls())
if (length(dev.list())) {
  dev.off()
}

library("kernlab")

set.seed(203)

C1_LABEL = -1
C2_LABEL = 1

s = 0.7
# h = 0.1
# C = 10

# aberturas das gaussianas: quanto maior s, mais ela espalha
nc <- 50 # numero de pontos de cada classe
k_list <- c(10)

# gera os dados - fora do loop, para que sejam os mesmo dados para todos os k
# centro dessa gaussiana esta na posicao (2,2)
xc1 <- matrix(rnorm(nc * 2) , ncol = 2) * s + t (matrix(c (2 , 2) , nrow =
                                                          2, ncol = nc))
# centro dessa outra gaussiana esta na posicao (4,4)
xc2 <- matrix(rnorm(nc * 2) , ncol = 2) * s + t (matrix(c (4 , 4) , nrow =
                                                          2, ncol = nc))

# matriz de entrada
X <- rbind(xc1, xc2)

# rotulos de entrada
yc1 <- rep(C1_LABEL, nc)
yc2 <- rep(C2_LABEL, nc)
Y <- c(yc1, yc2)

h_list <- c(0.001, 0.01, 0.1, 1)
C_list <- c(10)

# par(mfrow=c(2,2))

for (h in h_list) {
  for (C in C_list) {
    plot(
      NULL,
      main = paste("SVM: h = ", h, " C = ", C),
      xlab = "x1",
      ylab = "x2",
      ylim = c(0, 6),
      xlim = c(0, 6),
      cex.main = 2,
      cex.axis = 2,
      cex.lab = 2
    )
    
    points(xc1, col = "red", lwd = 2)
    points(xc2, col = "blue", lwd = 2)
    
    svm_train <- ksvm(X, Y, type='C-bsvc', kenel='rbfdot',
                      kpar=list(sigma=h), C=C)
    
    yhat <- predict(svm_train, X, type="response")
    
    a <- alpha(svm_train)
    ai <- SVindex(svm_train)
    nsvec <- nSV(svm_train)
    points(X[ai, 1], X[ai, 2], col="black", lwd = 4)
  }
}

dall <- as.matrix(dist(X, diag = T, upper = T))
kall <- exp(-(dall * dall) / (2 * (h)^2)) # essa é a matriz de kernel
heatmap(kall)

k11 <- kall[(1:50),(1:50)]
k12 <- kall[(1:50),(51:100)]
k21 <- kall[(51:100),(1:50)]
k22 <- kall[(51:100),(51:100)]

p11 <- rowSums(k11)
p12 <- rowSums(k12)
p21 <- rowSums(k21)
p22 <- rowSums(k22)

p1<-cbind(p11, p12) / 50
p2<-cbind(p21, p22) / 50

pall <- cbind(p1, p2)

# normaliza

# obs: valor de h que maximiza a soma dos elementos da matriz de covariancia
# (soma dos rows) -> interessante par ao trabalho
# plota essa função

plot(p1[,1], p1[,2], col = "red")
par(new = T)
plot(p2[,1], p2[,2], col = "blue")




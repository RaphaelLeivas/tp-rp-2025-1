rm(list = ls())
if (length(dev.list())) {
  dev.off()
}

library("kernlab")

set.seed(203)

distance_two_points <- function(x1, x2) {
  return (sqrt(sum((x1 - x2)^2)))
}

C1_LABEL = -1
C2_LABEL = 1

s = 0.75

# aberturas das gaussianas: quanto maior s, mais ela espalha
nc <- 250 # numero de pontos de cada classe
k_list <- c(10)

# gera os dados - fora do loop, para que sejam os mesmo dados para todos os k
# centro dessa gaussiana esta na posicao (2,2)
xc1 <- matrix(rnorm(nc * 2) , ncol = 2) * s + t (matrix(c (2 , 2) , nrow =
                                                          2, ncol = nc))
# centro dessa outra gaussiana esta na posicao (4,4)
xc2 <- matrix(rnorm(nc * 2) , ncol = 2) * s + t (matrix(c (4 , 4) , nrow =
                                                          2, ncol = nc))

plot(
  NULL,
  main = paste("Dados de entrada - s = ", s),
  xlab = "x1",
  ylab = "x2",
  ylim = c(0, 6),
  xlim = c(0, 6),
  cex.main = 2,
  cex.axis = 2,
  cex.lab = 2
)

points(xc1[,1], xc1[,2], col = "red", lwd = 2)
points(xc2[,1], xc2[,2], col = "blue", lwd = 2)

# matriz de entrada
X <- rbind(xc1, xc2)

# rotulos de entrada
yc1 <- rep(C1_LABEL, nc)
yc2 <- rep(C2_LABEL, nc)
Y <- c(yc1, yc2)

h_list <- seq(0.1, 5, 0.2)
C <- 10

dist_arr <- c()
area1_arr <- c()

for (h in h_list) {
  dall <- as.matrix(dist(X, diag = T, upper = T))
  kall <- exp(-(dall * dall) / (2 * (h)^2)) # essa é a matriz de kernel
  
  first_half <- nc
  second_hald <- nc + 1
  end <- 2 * nc
  
  k11 <- kall[(1:first_half),(1:first_half)]
  k12 <- kall[(1:first_half),(second_hald:end)]
  k21 <- kall[(second_hald:end),(1:first_half)]
  k22 <- kall[(second_hald:end),(second_hald:end)]
  
  p11 <- rowSums(k11)
  p12 <- rowSums(k12)
  p21 <- rowSums(k21)
  p22 <- rowSums(k22)
  
  p1<-cbind(p11, p12) / nc
  p2<-cbind(p21, p22) / nc
  
  pall <- cbind(p1, p2)
  
  plot(
    NULL,
    main = paste("SVM: h = ", h),
    xlab = "x1",
    ylab = "x2",
    ylim = c(0, 1),
    xlim = c(0, 1),
    cex.main = 2,
    cex.axis = 2,
    cex.lab = 2
  )

  points(p1[,1], p1[,2], col = "red", lwd = 2)
  points(p2[,1], p2[,2], col = "blue", lwd = 2)
  
  # calcula a distancia entre os centros para cada um
  m1 <- colMeans(p1)
  m2 <- colMeans(p2)
  
  # calcula os desvios padroes
  a1 <- sd(p1[which(p1[,1] >= p1[,2]),1]) # somente os abaixo da reta
  b1 <- sd(p1[which(p1[,1] >= p1[,2]),2])
  area1 <- pi * a1 * b1
  
  dist_arr <- c(dist_arr, distance_two_points(m1, m2))
  area1_arr <- c(area1_arr, area1)
}

plot(h_list, dist_arr, type = "b",
     main = "Distância entre centros x h",
     xlab = "h",
     ylab = "Distância",
     cex.main = 2,
     cex.axis = 2,
     cex.lab = 2,
     lwd = 2, lty = 1, col = "red")
par(new = T)
plot(h_list, area1_arr, type = "b",
     main = "Distância entre centros x h",
     xlab = "h",
     ylab = "Distância",
     cex.main = 2,
     cex.axis = 2,
     cex.lab = 2,
     lwd = 2, lty = 1, col = "orange")


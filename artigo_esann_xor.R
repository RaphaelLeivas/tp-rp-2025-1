rm(list = ls())
if (length(dev.list())) {
  dev.off()
}

library("kernlab")
library("mlbench")

set.seed(203)

C1_LABEL = -1
C2_LABEL = 1

distance_two_points <- function(x1, x2) {
  return (sqrt(sum((x1 - x2)^2)))
}

N <- 100
n <- 2

m1 <- c(2,2)
m2 <- c(4,4)
m3 <- c(2,4)
m4 <- c(4,2)

variancia = 0.3

g1 <- matrix(rnorm(2 * N), ncol = n, nrow = N)*variancia + matrix(m1, nrow = N, ncol = n, byrow = T)
g2 <- matrix(rnorm(2 * N), ncol = n, nrow = N)*variancia + matrix(m2, nrow = N, ncol = n, byrow = T)
g3 <- matrix(rnorm(2 * N), ncol = n, nrow = N)*variancia + matrix(m3, nrow = N, ncol = n, byrow = T)
g4 <- matrix(rnorm(2 * N), ncol = n, nrow = N)*variancia + matrix(m4, nrow = N, ncol = n, byrow = T)

xc1 <- rbind(g1, g2)
xc2 <- rbind(g3, g4)

nc1 <- nrow(xc1)
nc2 <- nrow(xc2)
nc_total <- nc1 + nc2

yc1 <- rep(C1_LABEL, nc1)
yc2 <- rep(C2_LABEL, nc2)

plot(
  NULL,
  main = "Espaço de Variáveis",
  xlab = "x1",
  ylab = "x2",
  ylim = c(0, 6),
  xlim = c(0, 6),
  cex.main = 2,
  cex.axis = 2,
  cex.lab = 1.5
)

points(xc1[,1], xc1[,2], col = "red", lwd = 2)
points(xc2[,1], xc2[,2], col = "blue", lwd = 2)

# matriz de entrada
X <- rbind(xc1, xc2)
Y <- c(yc1, yc2)

h_list <- seq(0.01, 4, 1)
C <- 10

dist_arr <- c()
area1_arr <- c()

for (h in h_list) {
  dall <- as.matrix(dist(X, diag = T, upper = T))
  kall <- exp(-(dall * dall) / (2 * (h)^2)) # essa é a matriz de kernel
  
  k11 <- kall[1:nc1, 1:nc1]
  k12 <- kall[(1:nc1),((nc1+1):nc_total)]
  k21 <- kall[((nc1+1):nc_total),(1:nc1)]
  k22 <- kall[((nc1+1):nc_total),((nc1+1):nc_total)]
  
  p11 <- rowSums(k11) / dim(k11)[2]
  p12 <- rowSums(k12) / dim(k12)[2]
  p21 <- rowSums(k21) / dim(k21)[2]
  p22 <- rowSums(k22) / dim(k22)[2]
  
  p1<-cbind(p11, p12) 
  p2<-cbind(p21, p22)
  
  plot(
    NULL,
    main = paste("SVM: h = ", h),
    xlab = "p1",
    ylab = "p2",
    ylim = c(0, 1),
    xlim = c(0, 1),
    cex.main = 2,
    cex.axis = 2,
    cex.lab = 1.5
  )

  points(p1[,1], p1[,2], col = "red", lwd = 2)
  points(p2[,1], p2[,2], col = "blue", lwd = 2)
  
  # calcula a distancia entre os centros para cada um
  m1 <- colMeans(p1)
  m2 <- colMeans(p2)
  
  # calcula os desvios padroes
  # a1 <- sd(p1[which(p1[,1] >= p1[,2]),1]) # somente os abaixo da reta
  # b1 <- sd(p1[which(p1[,1] >= p1[,2]),2])
  # 
  # a2 <- sd(p2[which(p2[,1] <= p2[,2]),1]) # somente os acima da reta
  # b2 <- sd(p2[which(p2[,1] <= p2[,2]),2])
  
  a1 <- sd(p1[,1])
  b1 <- sd(p1[,2])
  a2 <- sd(p2[,1])
  b2 <- sd(p2[,2])
  
  area1 <- pi * a1 * b1
  area2 <- pi * a2 * b2
  
  dist_arr <- c(dist_arr, distance_two_points(m1, m2))
  area1_arr <- c(area1_arr, area1 + area2)
  
  # treina o modelo svm
  svm_train <- ksvm(X, Y, type='C-bsvc', kernel='rbfdot',
                    kpar=list(sigma=h), C=C)
  
  # monta o grid
  seqi <- seq(0, 6, 0.1)
  seqj <- seq(0, 6, 0.1)
  M1 <- matrix(1, nrow =  length(seqi), ncol = length(seqj))

  ci <- 0

  for (i in seqi) {
    ci <- ci + 1
    cj <- 0
    
    for (j in seqj) {
      cj <- cj +1
      x <- matrix(c(i, j), byrow = T, ncol = 2)
      M1[ci, cj] <- predict(svm_train, x, type="response")
    }
  }
  
  filled.contour(seqi, seqj, M1, nlevels  = 1, lwd = 2, lty = 1, axes = T,
                 col = c("#FF000099", "#0000FF99"), 
                 main = paste("Espaço de Variáveis: h = ", h),
                 xlab = "x1",
                 ylab = "x2", 
                 plot.axes = {
                   axis(1)
                   axis(2)
                   contour(seqi, seqj, M1, nlevels  = 1, lwd = 2, add = TRUE)
                   points(xc1[,1], xc1[,2], col = "red", lwd = 2)
                   points(xc2[,1], xc2[,2], col = "blue", lwd = 2)
                 })
}

BUG()

DIST_COL <- "black"
AREA_COL <- "red"

par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
plot(h_list, dist_arr, ylab = "Distância entre centros", xlab = "h",
     type = "b", col = DIST_COL, lwd = 2, cex.main = 1.5,
     cex.axis = 1.5,pch = 16) # first plot
par(new = TRUE)
plot(h_list, area1_arr, axes = FALSE, bty = "n", xlab = "", ylab = "",
     main = "Métricas x h",
     type = "b", col = AREA_COL, lwd = 2, cex.main = 1.5,
     cex.axis = 1.5,
     cex.lab = 1.5, pch = 15)
axis(side=4, at = pretty(range(area1_arr)), col=AREA_COL, cex.axis = 1.5, col.axis=AREA_COL)
mtext("Área Elipse", side=4, line=3, col=AREA_COL)

## Add Legend
legend("topright",legend=c("Distância entre Centros","Área da Elipse"),
       text.col=c(DIST_COL,AREA_COL),pch=c(16,15),col=c(DIST_COL,AREA_COL))



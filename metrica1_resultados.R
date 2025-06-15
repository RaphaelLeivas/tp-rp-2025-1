rm(list = ls())
if (length(dev.list())) {
  dev.off()
}

library("kernlab")
library("mlbench")

set.seed(203)

C1_LABEL = 1
C2_LABEL = -1


distance_two_points <- function(x1, x2) {
  return (sqrt(sum((x1 - x2)^2)))
}

# pega os dados da package mlbench
data("BreastCancer")
data2 <- BreastCancer

# dados do Zoo
data("Zoo")
data3 <- Zoo

# dados do Vehicle
data("Vehicle")
data4 <- Vehicle

# dados do Glass
data("Glass")
data5 <- Glass

# ------- TRATAMENTO BREASTCANCER ---------

# Realiza o tratamento dos dados para remoção de NA
data2 <- data2[complete.cases(data2),]

start_variables_column <- 2
end_variables_column <- 10
label_column <- 11

X <- data.matrix(data2[, start_variables_column:end_variables_column])

N <- nrow(X) # numero de amostras
n <- ncol(X) # numero de variaveis

Y <- matrix(NA, nrow = N, ncol = 1)

for (i in 1:N) {
  # ultima coluna é a coluna com a classe
  if (data2[i, label_column] == "benign") {
    Y[i] <- C1_LABEL
  } else {
    Y[i] <- C2_LABEL
  }
}

# ------- TRATAMENTO BREASTCANCER ---------

# ------- TRATAMENTO ZOO ---------
# data3 <- data3[complete.cases(data3),]
# 
# label_column <- ncol(data3)
# legs_column <- 13
# for (j in 1:(ncol(data3) - 1)) {
#   if (j == legs_column) next
# 
#   data3[, j] <- as.integer(as.logical(data3[, j]))
# }
# 
# X <- data3[, 1:(ncol(data3) - 1)]
# Y <- c()
# 
# for (i in 1:nrow(data3)) {
#   if (data3[i, ncol(data3)] == "mammal") {
#     Y <- c(Y, C1_LABEL)
#   } else {
#     Y <- c(Y, C2_LABEL)
#   }
# }

# ------- TRATAMENTO ZOO ---------

# ------- TRATAMENTO VEHICLE ---------
# data4 <- data4[complete.cases(data4),]
# 
# label_column <- ncol(data4)
# 
# X <- data4[, 1:(ncol(data4) - 1)]
# Y <- c()
# 
# for (i in 1:nrow(data4)) {
#   if (data4[i, ncol(data4)] == "bus") {
#     Y <- c(Y, C1_LABEL)
#   } else {
#     Y <- c(Y, C2_LABEL)
#   }
# }

# ------- TRATAMENTO VEHICLE ---------

# ------- TRATAMENTO GLASS ---------
# data5 <- data5[complete.cases(data5),]
# 
# label_column <- ncol(data5)
# 
# X <- data5[, 1:(ncol(data5) - 1)]
# Y <- c()
# 
# for (i in 1:nrow(data5)) {
#   if (data5[i, ncol(data5)] == 1) {
#     Y <- c(Y, C1_LABEL)
#   } else {
#     Y <- c(Y, C2_LABEL)
#   }
# }

# ------- TRATAMENTO GLASS ---------

# ------- TRATAMENTO XOR ---------

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

# matriz de entrada
X <- rbind(xc1, xc2)
Y <- c(yc1, yc2)

# ------- TRATAMENTO XOR ---------

# junta tudo na matriz dos dados
all_data <- cbind(X, Y)

# embaralha a matriz dos dados de entrada - remove bias de coleta
all_data <- all_data[sample.int(nrow(all_data)), ]

xc1 <- all_data[which(all_data[,ncol(all_data)] == C1_LABEL), 1:(ncol(all_data) - 1)]
xc2 <- all_data[which(all_data[,ncol(all_data)] == C2_LABEL), 1:(ncol(all_data) - 1)]
nc1 <- nrow(xc1)
nc2 <- nrow(xc2)
nc_total <- nc1 + nc2

# matriz de entrada
X <- rbind(xc1, xc2)

h_list <- seq(0.1, 5, 0.2) # XOR
# h_list <- seq(0.1, 20, 1) # BreastCancer
# h_list <- seq(0.1, 5, 0.1) # Zoo
# h_list <- seq(0.1, 500, 25) # Vehicle
# h_list <- seq(0.1, 10, 0.2) # Glass
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

  number_of_sds <- 2
  
  a1 <- number_of_sds * sd(p1[,1])
  b1 <- number_of_sds * sd(p1[,2])
  a2 <- number_of_sds * sd(p2[,1])
  b2 <- number_of_sds * sd(p2[,2])
  
  area1 <- pi * a1 * b1
  area2 <- pi * a2 * b2
  
  dist_arr <- c(dist_arr, distance_two_points(m1, m2))
  area1_arr <- c(area1_arr, (area1 + area2) / 2)
  
  xc <- m1[1] # center x_c or h
  yc <- m1[2] # y_c or k
  a <- a1 # major axis length
  b <- b1 # minor axis length
  phi <- pi / 3 # angle of major axis with x axis phi or tau
  
  t <- seq(0, 2*pi, 0.01) 
  x <- xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi)
  y <- yc + a*cos(t)*cos(phi) + b*sin(t)*cos(phi)

  lines(x,y,pch=19, col='red', lwd = 4, lty = "dashed")
  
  xc <- m2[1] # center x_c or h
  yc <- m2[2] # y_c or k
  a <- a2 # major axis length
  b <- b2 # minor axis length
  phi <- pi / 3 # angle of major axis with x axis phi or tau
  
  t <- seq(0, 2*pi, 0.01) 
  x <- xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi)
  y <- yc + a*cos(t)*cos(phi) + b*sin(t)*cos(phi)
  
  lines(x,y,pch=19, col='blue', lwd = 4, lty = "dashed")
}

DIST_COL <- "black"
AREA_COL <- "red"

par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
plot(h_list, dist_arr, ylab = "Distância entre centros", xlab = "h",
     type = "b", col = DIST_COL, lwd = 2, cex.main = 1.5,
     cex.axis = 1.5,pch = 16, cex.lab = 1.5) # first plot
par(new = TRUE)
plot(h_list, area1_arr, axes = FALSE, bty = "n", xlab = "", ylab = "",
     main = "Métricas x h - XOR",
     type = "b", col = AREA_COL, lwd = 2, cex.main = 1.5,
     cex.axis = 1.5, cex.lab = 1.5,
     cex.lab = 1.5, pch = 15)
axis(side=4, at = pretty(range(area1_arr)), col=AREA_COL, cex.axis = 1.5, col.axis=AREA_COL)
mtext("Área Elipse", side=4, line=3, col=AREA_COL)

# Add Legend
# legend("bottom",legend=c("Distância entre Centros","Área da Elipse"),
#        text.col=c(DIST_COL,AREA_COL),pch=c(16,15),col=c(DIST_COL,AREA_COL))

print(paste("Max Referencia = ", h_list[which(dist_arr == max(dist_arr))]))
print(paste("Max Elipses = ", h_list[which(area1_arr == max(area1_arr))]))



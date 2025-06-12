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

# pega os dados da package mlbench
data("BreastCancer")
data2 <- BreastCancer

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


# junta tudo na matriz dos dados
all_data <- cbind(X, Y)

# embaralha a matriz dos dados de entrada - remove bias de coleta
all_data <- all_data[sample.int(nrow(all_data)), ]

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


xc1 <- all_data[which(all_data[,ncol(all_data)] == C1_LABEL), start_variables_column:end_variables_column]
xc2 <- all_data[which(all_data[,ncol(all_data)] == C2_LABEL), start_variables_column:end_variables_column]
nc1 <- nrow(xc1)
nc2 <- nrow(xc2)
nc_total <- nc1 + nc2

# plot(
#   NULL,
#   main = paste("Dados de entrada - s = ", s),
#   xlab = "x1",
#   ylab = "x2",
#   ylim = c(0, 6),
#   xlim = c(0, 6),
#   cex.main = 2,
#   cex.axis = 2,
#   cex.lab = 2
# )

# points(xc1[,1], xc1[,2], col = "red", lwd = 2)
# points(xc2[,1], xc2[,2], col = "blue", lwd = 2)

# matriz de entrada
X <- rbind(xc1, xc2)

# rotulos de entrada
yc1 <- rep(C1_LABEL, nc)
yc2 <- rep(C2_LABEL, nc)
Y <- c(yc1, yc2)

h_list <- seq(0.1, 40, 2)
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
}

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

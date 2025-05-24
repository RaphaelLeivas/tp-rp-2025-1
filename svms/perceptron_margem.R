rm(list = ls())
if (length(dev.list())) {
  dev.off()
}

library("plot3D")

set.seed(203)

# The margin of a linear classifier is the distance from the decision boundary 
# to the nearest training point

C1_LABEL = 1
C2_LABEL = 0

calc_margin_braga <- function(w, X_aug, Y) {
  
  margins <- (X_aug %*% w)
  
  imax_neg = which(margins == max(margins[margins < 0]))[1]
  
  imin_pos = which(margins == min(margins[margins > 0]))[1]
  
  return(list(imax_neg, imin_pos, margins))  # margem geométrica
}

s = 0.2

# aberturas das gaussianas: quanto maior s, mais ela espalha
nc <- 125 # numero de pontos de cada classe
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

# gera o grid
x1grid <- seq(0, 6, 0.1)
x2grid <- seq(0, 6, 0.1)
grid_matrix <- matrix(NA, nrow = length(x1grid), ncol = length(x2grid))

plot(
  NULL,
  main = paste("Perceptron: s = ", s),
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

w_list <- seq(5, 7, 0.1)
marg_arr <- c()

for (w in w_list) {
  Xaug <- cbind(1, X[,1], X[,2])
  w_otimo <- matrix(c(-w,1,1), ncol = 1, nrow = 3)
  
  b <- -w_otimo[1]/w_otimo[3]
  a <- -w_otimo[2]/w_otimo[3]
  
  lines(x1grid, a * x1grid + b)
  
  retlist <- calc_margin_braga(w_otimo, Xaug)
  min_pos <- retlist[[3]][retlist[[2]]]
  max_neg <- retlist[[3]][retlist[[1]]]
  
  marg_metrica <- abs(abs(min_pos) - abs(max_neg) / 2)
  
  marg_arr <- c(marg_arr, marg_metrica)
}

plot(w_list, marg_arr, type="b", lty="dashed", col = "black", lwd = 2,
     xlab = "W", ylab = "Diff Margem", main = "Diferença Margem x W")



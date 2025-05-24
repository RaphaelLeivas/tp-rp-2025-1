rm(list = ls())
if (length(dev.list())) {
  dev.off()
}

library("plot3D")

set.seed(203)

f <- function(x,m) {
  return ((x-m)^2)
}

xrange <- seq(0, 10, 0.1)
f1 <- f(xrange, 3)
f2 <- f(xrange, 6)

plot(xrange, f1, type = 'l', xlim=c(0,10), ylim=c(0, 100), col = "red", lwd = 2)
par(new = T)
plot(xrange, f2, type = 'l', xlim=c(0,10), ylim=c(0, 100), col = "blue", lwd = 2)


# plota no espaço de soluções
colvec <- c('black', 'green')
icol <- ((xrange >= 3) & (xrange <= 6))*1 + 1
plot(f1, f2, col = colvec[icol], lwd = 2)

rm(list = ls())
if (length(dev.list())) {
  dev.off()
}

source("C:\\dev\\padroes-2025-1\\utils\\trainperceptron.R")

s = 0.25

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
yc1 <- rep(1, nc)
yc2 <- rep(0, nc)
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

# treina o perceptron
eta <- 0.01
tol <- 0.01
maxepocas <- 1000
retlist <- trainperceptron(X, Y, eta, tol, maxepocas, 1)
w <- retlist[[1]] # w na posicao n+1 é o theta
evec <- retlist[[2]]

for (i in 1:length(x1grid)) {
  for (j in 1:length(x2grid)) {
    current_point_grid <- matrix(c(x1grid[i], x2grid[j]), ncol = 2)
    
    xtaug <- matrix(c(-1, current_point_grid[1,1], current_point_grid[1,2]), byrow = T, ncol = 1)
    
    # assume a reta otima
    # passa pelo 3,3
    w_otimo <- matrix(c(6,1,1), ncol = 1, nrow = 3)
    
    # calcula a saida do perceptron
    grid_matrix[i, j] = sign(t(w_otimo) %*% xtaug)
  }
}

contour2D(
  grid_matrix,
  x1grid,
  x2grid,
  nlevels = 1,
  xlim = c(0, 6),
  ylim = c(0, 6),
  add = T,
  col = "green",
  lwd = 2
)

w_list <- seq(5, 7, 0.1)
margem_vec <- c()
for (curr_w in w_list) {
  w_otimo <- matrix(c(curr_w,1,1), ncol = 1, nrow = 3)
  Xaug <- cbind(-1, X)
  distall <- Xaug %*% w_otimo
  plot(distall)
  
  # margem: media (soma e divide por 2) do max da neg e min da pos (das distancias)
  # varia o w_otimo (reta) para ver como a margem muda
  # varia o w de 5 a 7
  # mais info no cap 11 das notas de aula
  
  max_negs <- -Inf
  min_pos <- Inf
  
  for (i in 1:nrow(distall)) {
    curr_dist <- distall[i]
    
    if (curr_dist < 0 && curr_dist >= max_negs) {
      max_negs <- curr_dist
    }
    
    if (curr_dist > 0 && curr_dist <= min_pos) {
      min_pos <- curr_dist
    }
  }
  
  margem <- (min_pos - max_negs) / 2
  margem_vec <- c(margem_vec, margem)
}
plot(margem_vec)


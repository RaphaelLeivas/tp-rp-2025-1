# Limpeza inicial
rm(list = ls())
if (length(dev.list())) {
  dev.off()
}

# Bibliotecas
library(kernlab)
library(mlbench)
library(cluster)

set.seed(203)

# Função de plotagem
plot_similarity_space <- function(p1, p2, title) {
  plot(NULL, main = title,
       xlab = "Similaridade Média com Classe -1",
       ylab = "Similaridade Média com Classe +1",
       xlim = c(0, 1), ylim = c(0, 1),
       cex.main = 1.2, cex.axis = 1, cex.lab = 1)
  points(p1[,1], p1[,2], col = "red", pch = 16)
  points(p2[,1], p2[,2], col = "blue", pch = 16)
  legend("topright", legend = c("Classe -1", "Classe +1"),
         col = c("red", "blue"), pch = 16, bty = "n")
}

# Função principal do experimento
run_projection_experiment <- function(X, Y, dataset_name, h_ref = NULL) {
  cat(sprintf("\n>>> Iniciando análise para o dataset: %s\n", dataset_name))
  C1_LABEL <- -1; C2_LABEL <- 1
  xc1 <- X[Y == C1_LABEL, ]; xc2 <- X[Y == C2_LABEL, ]
  nc1 <- nrow(xc1); nc2 <- nrow(xc2)
  if (nc1 < 2 || nc2 < 2) stop("Uma das classes não tem amostras suficientes.")
  X_ord <- rbind(xc1, xc2)

  dall <- as.matrix(dist(X_ord))
  
  # Grade de busca híbrida para h
  base_h <- seq(0.1, 200, length.out = 200)
  if (!is.null(h_ref)) {
    around_ref <- seq(h_ref / 5, h_ref * 5, length.out = 200)
    h_list <- sort(unique(c(base_h, around_ref, h_ref)))
  } else {
    h_list <- sort(unique(base_h))
  }
  cat(sprintf("Analisando %d valores de h.\n", length(h_list)))

  sil_arr <- numeric(length(h_list))
  proj_list <- vector("list", length(h_list))
  
  for (i in seq_along(h_list)) {
    h <- h_list[i]
    K <- exp(- (dall * dall) / (2 * h^2))
    k11 <- K[1:nc1, 1:nc1]; k12 <- K[1:nc1, (nc1+1):(nc1+nc2)]; k22 <- K[(nc1+1):(nc1+nc2), (nc1+1):(nc1+nc2)]
    p11 <- rowSums(k11)/nc1; p12 <- rowSums(k12)/nc2
    p21 <- colSums(k12)/nc1; p22 <- rowSums(k22)/nc2
    p1 <- cbind(p11, p12); p2 <- cbind(p21, p22)
    proj_list[[i]] <- list(p1 = p1, p2 = p2) 

    pts <- rbind(p1, p2)
    lbls <- c(rep(C1_LABEL, nc1), rep(C2_LABEL, nc2))
    sil <- silhouette(as.integer(as.factor(lbls)), dist(pts))
    sil_arr[i] <- mean(sil[, 'sil_width'])
  }

  df_sil <- data.frame(h = h_list, sil = sil_arr)
  top10 <- head(df_sil[order(-df_sil$sil), ], 10)
  print("Top 10 valores de h por índice de silhueta:")
  print(top10)

  best_idx <- which.max(sil_arr)
  best_h <- h_list[best_idx]
  
  # Plot Silhueta vs h
  plot_limit <- min(length(h_list), 50) # Garante que não vai estourar o limite
  h_plot <- h_list[1:plot_limit]
  sil_plot <- sil_arr[1:plot_limit]
  plot(h_plot, sil_plot, type = 'l', lwd = 2, col = 'darkgreen',
       xlab = 'Parâmetro h', ylab = 'Silhueta Média',
       main = paste('Silhueta vs h (Primeiros', plot_limit, 'valores) -', dataset_name))
  points(h_plot, sil_plot, pch = 16)
  
  if (best_h %in% h_plot) abline(v = best_h, col = 'red', lty = 2, lwd = 2)
  if (!is.null(h_ref) && h_ref %in% h_plot) abline(v = h_ref, col = 'purple', lty = 3, lwd = 2)
  
  legend('topright', bty="n",
         legend = c(if(best_h %in% h_plot) sprintf('Melhor h = %.3f', best_h),
                    if(!is.null(h_ref) && h_ref %in% h_plot) sprintf('h_ref = %.3f', h_ref)),
         col = c('red', 'purple')[c(best_h %in% h_plot, !is.null(h_ref) && h_ref %in% h_plot)],
         lty = c(2, 3)[c(best_h %in% h_plot, !is.null(h_ref) && h_ref %in% h_plot)],
         lwd = 2)

  # Plot projeções para h selecionados
  select_h <- sort(unique(c(
    h_list[1], h_list[2],
    if (!is.null(h_ref)) h_ref,
    best_h,
    h_list[length(h_list)-1], h_list[length(h_list)]
  )))
  
  par(mfrow = c(2, ceiling(length(select_h)/2)), mar = c(4,4,3,1))
  
  for (h_val in select_h) {
    idx_to_plot <- which(h_list == h_val)[1] 
    pr <- proj_list[[idx_to_plot]]
    title <- sprintf('%s - h=%.3f%s', dataset_name, h_val,
                     if (!is.null(h_ref) && h_val == h_ref) ' (ref)' else if (h_val == best_h) ' (best)' else '')
    plot_similarity_space(pr$p1, pr$p2, title)
  }
  par(mfrow = c(1,1))
}

# --- Blocos de execução ---

# Dataset Zoo
data(Zoo, package = 'mlbench')
Xz <- data.matrix(Zoo[, -c(1,17)])
Yz <- ifelse(Zoo$type == 'mammal', 1, -1)
run_projection_experiment(Xz, Yz, 'Zoo', h_ref = 1.7)

# Dataset Vehicle
data(Vehicle, package = 'mlbench')
Xv <- data.matrix(Vehicle[, -ncol(Vehicle)])
Yv <- ifelse(Vehicle$Class == 'bus', 1, -1)
run_projection_experiment(Xv, Yv, 'Vehicle', h_ref = 50)

# Dataset Glass
data(Glass, package = 'mlbench')
Xg <- data.matrix(Glass[, -ncol(Glass)])
Yg <- ifelse(Glass$Type == '1', 1, -1)
run_projection_experiment(Xg, Yg, 'Glass', h_ref = 1.9)

# --- Experimento 4: Dataset BreastCancer ---
data(BreastCancer, package = 'mlbench')
# Remove linhas com valores ausentes (NA) e a coluna de ID
bc_data <- BreastCancer[complete.cases(BreastCancer), -1] 
# As features são as colunas de 1 a 9 na nova matriz de dados
X_bc <- data.matrix(bc_data[, 1:9])
# A classe 'benign' será +1 e 'malignant' será -1
Y_bc <- ifelse(bc_data$Class == 'benign', 1, -1)
# Um valor de referência razoável para este dataset
run_projection_experiment(X_bc, Y_bc, 'BreastCancer', h_ref = 25)

# --- Experimento 5: Dataset XOR Fictício ---
# Geração dos dados conforme especificado
N <- 100
n <- 2
C1_LABEL <- -1 # Definindo rótulos para o XOR
C2_LABEL <- 1

m1 <- c(2,2); m2 <- c(4,4) # Cantos para a Classe -1
m3 <- c(2,4); m4 <- c(4,2) # Cantos para a Classe +1
variancia <- 0.3

g1 <- matrix(rnorm(2 * N), ncol = n) * variancia + matrix(m1, nrow = N, ncol = n, byrow = TRUE)
g2 <- matrix(rnorm(2 * N), ncol = n) * variancia + matrix(m2, nrow = N, ncol = n, byrow = TRUE)
g3 <- matrix(rnorm(2 * N), ncol = n) * variancia + matrix(m3, nrow = N, ncol = n, byrow = TRUE)
g4 <- matrix(rnorm(2 * N), ncol = n) * variancia + matrix(m4, nrow = N, ncol = n, byrow = TRUE)

# Combina os clusters para formar as duas classes
xc1_xor <- rbind(g1, g2) # Classe -1
xc2_xor <- rbind(g3, g4) # Classe +1

# Monta a matriz de features e o vetor de rótulos final
X_xor <- rbind(xc1_xor, xc2_xor)
Y_xor <- c(rep(C1_LABEL, nrow(xc1_xor)), rep(C2_LABEL, nrow(xc2_xor)))
# Para o XOR, um h pequeno é esperado como ótimo
run_projection_experiment(X_xor, Y_xor, 'XOR Fictício', h_ref = 0.8)
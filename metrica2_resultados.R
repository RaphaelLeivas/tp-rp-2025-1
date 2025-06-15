# Limpeza inicial
rm(list = ls())
if (length(dev.list())) dev.off()

# Bibliotecas
library(kernlab)
library(mlbench)
library(cluster)

set.seed(203)

# Função de plotagem (sem alterações)
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

# ==========================================================================
# FUNÇÃO DE EXPERIMENTO PARA KERNEL RBF (IMPLEMENTAÇÃO COMPLETA)
# ==========================================================================
run_rbf_experiment <- function(X, Y, dataset_name, h_ref = NULL, h_min = 0.1, h_max = 200) {
  cat(sprintf("\n>>> Iniciando análise RBF para o dataset: %s\n", dataset_name))
  cat(sprintf("--- Usando faixa de h: [%.1f, %.1f]\n", h_min, h_max))
  
  C1_LABEL <- -1; C2_LABEL <- 1
  xc1 <- X[Y == C1_LABEL, ]; xc2 <- X[Y == C2_LABEL, ]
  nc1 <- nrow(xc1); nc2 <- nrow(xc2)
  if (nc1 < 2 || nc2 < 2) stop("Uma das classes não tem amostras suficientes.")
  X_ord <- rbind(xc1, xc2)

  dall <- as.matrix(dist(X_ord))
  
  h_list <- sort(unique(c(seq(h_min, h_max, length.out = 100), h_ref)))
  h_list <- h_list[h_list > 0]
  
  sil_arr <- numeric(length(h_list))
  proj_list <- vector("list", length(h_list))
  
  for (i in seq_along(h_list)) {
    h <- h_list[i]
    K <- exp(- (dall * dall) / (2 * h^2))
    
    k11 <- K[1:nc1, 1:nc1]; k12 <- K[1:nc1, (nc1+1):(nc1+nc2)]; k22 <- K[(nc1+1):(nc1+nc2), (nc1+1):(nc1+nc2)]
    p11 <- rowSums(k11)/nc1; p12 <- rowSums(k12)/nc2
    p21 <- colSums(k12)/nc1; p22 <- rowSums(k22)/nc2
    p1 <- cbind(p11, p12); p2 <- cbind(p21, p22)
    proj_list[[i]] <- list(p1=p1, p2=p2)
    
    pts <- rbind(p1, p2); lbls <- c(rep(C1_LABEL, nc1), rep(C2_LABEL, nc2))
    sil <- silhouette(as.integer(as.factor(lbls)), dist(pts))
    sil_arr[i] <- mean(sil[,'sil_width'], na.rm=T)
  }

  df_sil <- data.frame(h=h_list, sil=sil_arr); top10 <- head(df_sil[order(-df_sil$sil, na.last=NA), ], 10)
  print("Top 10 valores de h (RBF):"); print(top10)
  
  best_idx <- which.max(sil_arr); best_h <- h_list[best_idx]
  
  plot(h_list, sil_arr, type='l', lwd=2, col='darkgreen', xlab='Parâmetro h (RBF)', ylab='Silhueta Média', main=paste('RBF - Silhueta vs h -', dataset_name))
  points(h_list, sil_arr, pch=16)
  if(!is.na(best_h)) abline(v=best_h, col='red', lty=2, lwd=2)
  if(!is.null(h_ref)) abline(v=h_ref, col='purple', lty=3, lwd=2)
  legend('topright', bty="n", legend=c(if(!is.na(best_h)) sprintf('Melhor h = %.3f', best_h), if(!is.null(h_ref)) sprintf('h_ref = %.3f', h_ref)), col=c('red','purple'), lty=c(2,3), lwd=2)

  select_h <- sort(unique(c(h_list[1], h_list[20], if(!is.null(h_ref)) h_ref, best_h, h_list[length(h_list)-20], h_list[length(h_list)])))
  par(mfrow=c(2, ceiling(length(select_h)/2)), mar=c(4,4,3,1))
  for(h_val in select_h) {
    if(is.na(h_val)) next
    idx_to_plot <- which(h_list == h_val)[1]; pr <- proj_list[[idx_to_plot]]
    title <- sprintf('%s (RBF)\nh=%.3f%s', dataset_name, h_val, if(!is.null(h_ref)&&h_val==h_ref) ' (ref)' else if(h_val==best_h) ' (best)' else '')
    plot_similarity_space(pr$p1, pr$p2, title)
  }
  par(mfrow=c(1,1))
  cat(sprintf("<<< Análise RBF para %s concluída.\n", dataset_name))
}


# ==========================================================================
# FUNÇÃO DE EXPERIMENTO PARA KERNEL SIGMOIDAL (IMPLEMENTAÇÃO COMPLETA)
# ==========================================================================
run_sigmoid_experiment <- function(X, Y, dataset_name, h_ref = NULL, h_min = 0.1, h_max = 500) {
  cat(sprintf("\n>>> Iniciando análise SIGMOIDAL para o dataset: %s\n", dataset_name))
  cat(sprintf("--- Usando faixa de h: [%.1f, %.1f]\n", h_min, h_max))
  
  C1_LABEL <- -1; C2_LABEL <- 1
  xc1 <- X[Y == C1_LABEL, ]; xc2 <- X[Y == C2_LABEL, ]
  nc1 <- nrow(xc1); nc2 <- nrow(xc2)
  if (nc1 < 2 || nc2 < 2) stop("Uma das classes não tem amostras suficientes.")
  X_ord <- rbind(xc1, xc2)

  dot_product_matrix <- X_ord %*% t(X_ord)
  
  h_list <- sort(unique(c(seq(h_min, h_max, length.out = 100), h_ref)))
  h_list <- h_list[h_list > 0]
  
  sil_arr <- numeric(length(h_list))
  proj_list <- vector("list", length(h_list))
  
  for (i in seq_along(h_list)) {
    h <- h_list[i]
    K <- tanh(dot_product_matrix / (2 * h^2))
    
    k11 <- K[1:nc1, 1:nc1]; k12 <- K[1:nc1, (nc1+1):(nc1+nc2)]; k22 <- K[(nc1+1):(nc1+nc2), (nc1+1):(nc1+nc2)]
    p11 <- rowSums(k11)/nc1; p12 <- rowSums(k12)/nc2
    p21 <- colSums(k12)/nc1; p22 <- rowSums(k22)/nc2
    p1 <- cbind(p11, p12); p2 <- cbind(p21, p22)
    proj_list[[i]] <- list(p1 = p1, p2 = p2) 

    pts <- rbind(p1, p2); lbls <- c(rep(C1_LABEL, nc1), rep(C2_LABEL, nc2))
    sil <- silhouette(as.integer(as.factor(lbls)), dist(pts))
    sil_arr[i] <- mean(sil[, 'sil_width'], na.rm = TRUE)
  }

  df_sil <- data.frame(h = h_list, sil = sil_arr)
  top10 <- head(df_sil[order(-df_sil$sil, na.last = NA), ], 10)
  print("Top 10 valores de h (Sigmoidal):"); print(top10)

  best_idx <- which.max(sil_arr); best_h <- h_list[best_idx]
  
  plot(h_list, sil_arr, type = 'l', lwd = 2, col = 'blue',
       xlab = 'Parâmetro h (Sigmoidal)', ylab = 'Silhueta Média',
       main = paste('Sigmoidal - Silhueta vs h -', dataset_name))
  points(h_list, sil_arr, pch = 16)
  
  if (!is.na(best_h)) abline(v = best_h, col = 'red', lty = 2, lwd = 2)
  if (!is.null(h_ref)) abline(v = h_ref, col = 'purple', lty = 3, lwd = 2)
  
  legend('topright', bty="n",
         legend = c(if(!is.na(best_h)) sprintf('Melhor h = %.3f', best_h),
                    if(!is.null(h_ref)) sprintf('h_ref = %.3f', h_ref)),
         col = c('red', 'purple'), lty = c(2, 3), lwd = 2)

  select_h <- sort(unique(c(h_list[1], h_list[findInterval(0.25*max(h_list), h_list)], if (!is.null(h_ref)) h_ref, best_h, h_list[findInterval(0.75*max(h_list), h_list)], h_list[length(h_list)])))
  par(mfrow = c(2, ceiling(length(select_h)/2)), mar = c(4,4,3,1))
  for (h_val in select_h) {
    if(is.na(h_val)) next
    idx_to_plot <- which(h_list == h_val)[1] 
    pr <- proj_list[[idx_to_plot]]
    title <- sprintf('%s (Sigmoidal)\nh=%.3f%s', dataset_name, h_val,
                     if (!is.null(h_ref) && h_val == h_ref) ' (ref)' else if (h_val == best_h) ' (best)' else '')
    plot_similarity_space(pr$p1, pr$p2, title)
  }
  par(mfrow=c(1,1))
  cat(sprintf("<<< Análise SIGMOIDAL para %s concluída.\n", dataset_name))
}


# ==========================================================================
# BLOCOS DE EXECUÇÃO
# ==========================================================================
# --- MODIFICADO: Configurações centralizadas para todos os experimentos ---
datasets_info <- list(
  list(name="XOR", 
       h_ref_rbf=0.8, h_min_rbf=0.1, h_max_rbf=5,
       h_ref_sigmoid=2.5, h_min_sigmoid=0.1, h_max_sigmoid=5,
       loader=function() { N<-100;n<-2;C1_LABEL<--1;C2_LABEL<-1;m1<-c(2,2);m2<-c(4,4);m3<-c(2,4);m4<-c(4,2);v<-0.3;g1<-matrix(rnorm(2*N),ncol=n)*v+matrix(m1,nrow=N,ncol=n,byrow=T);g2<-matrix(rnorm(2*N),ncol=n)*v+matrix(m2,nrow=N,ncol=n,byrow=T);g3<-matrix(rnorm(2*N),ncol=n)*v+matrix(m3,nrow=N,ncol=n,byrow=T);g4<-matrix(rnorm(2*N),ncol=n)*v+matrix(m4,nrow=N,ncol=n,byrow=T);xc1<-rbind(g1,g2);xc2<-rbind(g3,g4);X<-rbind(xc1,xc2);Y<-c(rep(C1_LABEL,nrow(xc1)),rep(C2_LABEL,nrow(xc2))); list(X=X,Y=Y) }),
  list(name="BreastCancer", 
       h_ref_rbf=25, h_min_rbf=1, h_max_rbf=50,
       h_ref_sigmoid=8, h_min_sigmoid=0.1, h_max_sigmoid=20,
       loader=function() { data(BreastCancer, package='mlbench'); bc_data <- BreastCancer[complete.cases(BreastCancer),-1]; X <- data.matrix(bc_data[,1:9]); Y <- ifelse(bc_data$Class=='benign',1,-1); list(X=X,Y=Y) }),
  list(name="Zoo", 
       h_ref_rbf=1.7, h_min_rbf=0.1, h_max_rbf=5,
       h_ref_sigmoid=2.6, h_min_sigmoid=0.1, h_max_sigmoid=5,
       loader=function() { data(Zoo, package='mlbench'); X <- data.matrix(Zoo[,-c(1,17)]); Y <- ifelse(Zoo$type=='mammal',1,-1); list(X=X,Y=Y) }),
  list(name="Vehicle", 
       h_ref_rbf=50, h_min_rbf=1, h_max_rbf=200,
       h_ref_sigmoid=400, h_min_sigmoid=1, h_max_sigmoid=500,
       loader=function() { data(Vehicle, package='mlbench'); X <- data.matrix(Vehicle[,-ncol(Vehicle)]); Y <- ifelse(Vehicle$Class=='bus',1,-1); list(X=X,Y=Y) }),
  list(name="Glass", 
       h_ref_rbf=1.9, h_min_rbf=0.1, h_max_rbf=10,
       h_ref_sigmoid=60, h_min_sigmoid=1, h_max_sigmoid=100,
       loader=function() { data(Glass, package='mlbench'); X <- data.matrix(Glass[,-ncol(Glass)]); Y <- ifelse(Glass$Type=='1',1,-1); list(X=X,Y=Y) })
)

# --- Loop principal de execução ---
for (info in datasets_info) {
  data <- info$loader()
  
  # Executa a análise para o Kernel RBF
  run_rbf_experiment(data$X, data$Y, info$name, 
                     h_ref = info$h_ref_rbf,
                     h_min = info$h_min_rbf,
                     h_max = info$h_max_rbf)
  
  # Executa a análise para o Kernel Sigmoidal
  run_sigmoid_experiment(data$X, data$Y, info$name, 
                         h_ref = info$h_ref_sigmoid,
                         h_min = info$h_min_sigmoid,
                         h_max = info$h_max_sigmoid)
}
rm(list = ls())
if (length(dev.list())) {
  dev.off()
}

# --------------------------------------------------------------------------
# ANÁLISE DE HIPERPARÂMETROS SVM VIA PROJEÇÃO NO ESPAÇO DE SIMILARIDADES
# --------------------------------------------------------------------------
# Este script investiga o impacto do hiperparâmetro 'h' (relacionado ao 
# gamma) de um kernel RBF.
#
# Métrica de Qualidade Utilizada: Índice de Silhueta (Silhouette Score)
# Para cada valor de 'h', os dados são projetados em um espaço 2D onde os eixos
# representam a similaridade média com cada classe. A qualidade dessa projeção
# é medida pelo Índice de Silhueta médio. Um valor mais alto indica que as
# classes estão mais densas e bem separadas, sugerindo um bom valor para 'h'.

# O que o Índice de Silhueta mede: Um Índice de Silhueta mais alto (próximo de 1) 
# significa uma separação melhor e mais clara entre as classes. Um valor baixo
# (próximo de 0 ou negativo) significa que as classes estão misturadas.
# --------------------------------------------------------------------------


# Sessão de Carregamento de Bibliotecas ------------------------------------
library("kernlab")
library("mlbench")
library("cluster") 


# Sessão de Inicialização e Definições ------------------------------------
set.seed(203)
C1_LABEL = -1
C2_LABEL = 1


# Sessão de Carregamento e Pré-processamento dos Dados -------------------
# Pega os dados da package mlbench
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
  if (data2[i, label_column] == "benign") {
    Y[i] <- C1_LABEL
  } else {
    Y[i] <- C2_LABEL
  }
}

# Junta tudo na matriz dos dados
all_data <- cbind(X, Y)

# Embaralha a matriz dos dados de entrada para remover viés de coleta
all_data <- all_data[sample.int(nrow(all_data)), ]

# Separa os dados por classe novamente após o embaralhamento
xc1 <- all_data[which(all_data[,ncol(all_data)] == C1_LABEL), start_variables_column:end_variables_column]
xc2 <- all_data[which(all_data[,ncol(all_data)] == C2_LABEL), start_variables_column:end_variables_column]
nc1 <- nrow(xc1)
nc2 <- nrow(xc2)
nc_total <- nc1 + nc2

# Junta os dados de entrada na ordem correta (todos da classe 1, depois todos da classe 2)
X <- rbind(xc1, xc2)


# Sessão Principal: Projeção e Cálculo da Métrica -----------------------------
h_list <- seq(0.1, 40, 2)
sil_arr <- c() # Vetor para armazenar os resultados da silhueta

# Loop sobre diferentes valores do parâmetro 'h' do kernel
for (h in h_list) {
  # -- Cálculo da Matriz de Similaridade (Kernel RBF) --
  dall <- as.matrix(dist(X, diag = T, upper = T))
  kall <- exp(-(dall * dall) / (2 * (h)^2)) # Essa é a matriz de kernel
  
  # -- Particionamento da Matriz de Kernel --
  k11 <- kall[1:nc1, 1:nc1]
  k12 <- kall[1:nc1, (nc1+1):nc_total]
  k22 <- kall[(nc1+1):nc_total, (nc1+1):nc_total]
  
  # -- Projeção dos Dados no Espaço de Similaridade 2D --
  p11 <- rowSums(k11) / nc1
  p12 <- rowSums(k12) / nc2
  
  # p21 é a similaridade dos pontos da classe 2 com os da classe 1.
  # Usamos a transposta de k12 (colSums) para calcular isso eficientemente.
  p21 <- colSums(k12) / nc1 
  p22 <- rowSums(k22) / nc2
  
  p1 <- cbind(p11, p12) 
  p2 <- cbind(p21, p22)
  
  # Visualização da projeção para cada 'h'
  plot(
    NULL, main = paste("SVM: h = ", h), xlab = "Similaridade com Classe 1",
    ylab = "Similaridade com Classe 2", ylim = c(0, 1), xlim = c(0, 1)
  )
  points(p1[,1], p1[,2], col = "red", lwd = 2)
  points(p2[,1], p2[,2], col = "blue", lwd = 2)
  
  # -- Cálculo da Métrica de Qualidade (Índice de Silhueta) --
  # 1. Juntar todos os pontos projetados em uma única matriz
  proj_points <- rbind(p1, p2)
  
  # 2. Criar um vetor com os rótulos correspondentes aos pontos projetados
  proj_labels <- c(rep(C1_LABEL, nc1), rep(C2_LABEL, nc2))

  # 3. Calcular a matriz de distância para os pontos projetados
  dist_proj <- dist(proj_points)
  
  # 4. Calcular o objeto silhueta
  sil_obj <- silhouette(as.integer(as.factor(proj_labels)), dist_proj)
  
  # 5. Extrair e armazenar a largura média da silhueta
  if (!is.null(sil_obj)) {
    mean_sil <- mean(sil_obj[, "sil_width"])
  } else {
    mean_sil <- NA 
  }
  sil_arr <- c(sil_arr, mean_sil)
}


# Sessão de Visualização dos Resultados -------------------------------------
# Plota o Índice de Silhueta em função do parâmetro 'h'
plot(
  h_list,
  sil_arr,
  type = "b",
  lwd = 2,
  pch = 16,
  col = "blue",
  xlab = "Parâmetro do Kernel (h)",
  ylab = "Índice de Silhueta Médio",
  main = "Qualidade da Projeção vs. Parâmetro h",
  cex.main = 1.5,
  cex.axis = 1.5,
  cex.lab = 1.5,
  ylim = c(min(sil_arr, na.rm = TRUE) * 0.9, max(sil_arr, na.rm = TRUE) * 1.1)
)

# Adiciona linhas para destacar o melhor resultado encontrado
best_h_index <- which.max(sil_arr)
best_h <- h_list[best_h_index]
max_sil <- sil_arr[best_h_index]
abline(h = max_sil, col = "red", lty = 2)
abline(v = best_h, col = "red", lty = 2)
legend(
  "bottomright",
  legend = paste("Melhor h =", best_h, "\nSilhueta =", round(max_sil, 3)),
  cex = 1.2
)
# Limpeza inicial e carregamento de bibliotecas
rm(list = ls())
if (length(dev.list())) dev.off()

library(kernlab)
library(mlbench)

set.seed(203) # Garante a reprodutibilidade dos resultados

# ==========================================================================
# FUNÇÃO PARA CALCULAR A ACURÁCIA MÉDIA VIA VALIDAÇÃO CRUZADA (CV)
# ==========================================================================
# Esta função recebe os dados, o valor de h e realiza o treinamento e a validação.
calcular_acuracia_cv <- function(X, Y, h_param, n_folds = 5, C_param = 10) {
  
  # Junta os dados e embaralha para garantir que os folds sejam imparciais
  all_data <- cbind(as.matrix(X), Y)
  all_data <- all_data[sample(nrow(all_data)), ]
  
  fold_size <- floor(nrow(all_data) / n_folds)
  acc_array <- numeric(n_folds) # Vetor para guardar a acurácia de cada fold
  
  for (fold in 1:n_folds) {
    # Define os índices para treino e teste
    start_index <- (fold - 1) * fold_size + 1
    end_index <- fold * fold_size
    if (fold == n_folds) {
      end_index <- nrow(all_data) # Garante que o último fold pegue todos os dados restantes
    }
    
    test_indices <- start_index:end_index
    
    # Separa os dados
    data_for_test <- all_data[test_indices, ]
    data_for_train <- all_data[-test_indices, ]
    
    X_train <- data_for_train[, 1:(ncol(data_for_train) - 1)]
    Y_train <- data_for_train[, ncol(data_for_train)]
    
    X_test <- data_for_test[, 1:(ncol(data_for_test) - 1)]
    Y_test <- data_for_test[, ncol(data_for_test)]
    
    # Treina o modelo SVM com kernel RBF e o 'h' especificado
    svm_model <- ksvm(X_train, Y_train, 
                        type = 'C-bsvc', 
                        kernel = 'rbfdot',
                        kpar = list(sigma = 1 / (2 * h_param^2)), # Converte h para sigma
                        C = C_param)
    
    # Realiza a predição e calcula a acurácia do fold
    y_pred <- predict(svm_model, X_test)
    acc_array[fold] <- sum(y_pred == Y_test) / length(Y_test)
  }
  
  # Retorna a acurácia média de todos os folds
  return(mean(acc_array) * 100)
}


# ==========================================================================
# BLOCOS DE EXECUÇÃO PARA CADA DATASET
# ==========================================================================
C1_LABEL <- 1
C2_LABEL <- -1
resultados <- list()

# --- 1. Dataset XOR ---
cat("Processando Dataset: XOR...\n")
h_xor <- 5
N <- 100; n <- 2; variancia <- 0.3
m1 <- c(2,2); m2 <- c(4,4); m3 <- c(2,4); m4 <- c(4,2)
g1 <- matrix(rnorm(2*N), ncol=n)*variancia + matrix(m1, nrow=N, ncol=n, byrow=T)
g2 <- matrix(rnorm(2*N), ncol=n)*variancia + matrix(m2, nrow=N, ncol=n, byrow=T)
g3 <- matrix(rnorm(2*N), ncol=n)*variancia + matrix(m3, nrow=N, ncol=n, byrow=T)
g4 <- matrix(rnorm(2*N), ncol=n)*variancia + matrix(m4, nrow=N, ncol=n, byrow=T)
xc1_xor <- rbind(g1,g2); xc2_xor <- rbind(g3,g4)
X_xor <- rbind(xc1_xor, xc2_xor); Y_xor <- c(rep(C1_LABEL, nrow(xc1_xor)), rep(C2_LABEL, nrow(xc2_xor)))
resultados$XOR <- calcular_acuracia_cv(X_xor, Y_xor, h_param = h_xor)

# --- 2. Dataset BreastCancer ---
cat("Processando Dataset: BreastCancer...\n")
h_bc <- 13
data(BreastCancer, package = 'mlbench')
bc_data <- BreastCancer[complete.cases(BreastCancer), -1]
X_bc <- data.matrix(bc_data[, 1:9])
Y_bc <- ifelse(bc_data$Class == 'benign', C1_LABEL, C2_LABEL)
resultados$BreastCancer <- calcular_acuracia_cv(X_bc, Y_bc, h_param = h_bc)

# --- 3. Dataset Zoo ---
cat("Processando Dataset: Zoo...\n")
h_zoo <- 5
data(Zoo, package = 'mlbench')
X_zoo <- data.matrix(Zoo[, -c(1,17)])
Y_zoo <- ifelse(Zoo$type == 'mammal', C1_LABEL, C2_LABEL)
resultados$Zoo <- calcular_acuracia_cv(X_zoo, Y_zoo, h_param = h_zoo)

# --- 4. Dataset Vehicle ---
cat("Processando Dataset: Vehicle...\n")
h_vehicle <- 1
data(Vehicle, package = 'mlbench')
X_vehicle <- data.matrix(Vehicle[, -ncol(Vehicle)])
Y_vehicle <- ifelse(Vehicle$Class == 'bus', C1_LABEL, C2_LABEL)
resultados$Vehicle <- calcular_acuracia_cv(X_vehicle, Y_vehicle, h_param = h_vehicle)

# --- 5. Dataset Glass ---
cat("Processando Dataset: Glass...\n")
h_glass <- 1
data(Glass, package = 'mlbench')
X_glass <- data.matrix(Glass[, -ncol(Glass)])
Y_glass <- ifelse(Glass$Type == '1', C1_LABEL, C2_LABEL)
resultados$Glass <- calcular_acuracia_cv(X_glass, Y_glass, h_param = h_glass)

# ==========================================================================
# APRESENTAÇÃO DOS RESULTADOS FINAIS
# ==========================================================================
cat("\n--- Resultados Finais ---\n")
cat(sprintf("Dataset: XOR          | h = %-4d | Acurácia Média CV: %.2f%%\n", h_xor, resultados$XOR))
cat(sprintf("Dataset: BreastCancer   | h = %-4d | Acurácia Média CV: %.2f%%\n", h_bc, resultados$BreastCancer))
cat(sprintf("Dataset: Zoo            | h = %-4d | Acurácia Média CV: %.2f%%\n", h_zoo, resultados$Zoo))
cat(sprintf("Dataset: Vehicle        | h = %-4d | Acurácia Média CV: %.2f%%\n", h_vehicle, resultados$Vehicle))
cat(sprintf("Dataset: Glass          | h = %-4d | Acurácia Média CV: %.2f%%\n", h_glass, resultados$Glass))
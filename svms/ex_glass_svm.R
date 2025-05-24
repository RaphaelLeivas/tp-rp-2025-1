rm(list = ls())
if (length(dev.list())) {
  dev.off()
}

library("mlbench")
library("flextable")
library("kernlab")

set.seed(205)

# pega os dados da package mlbench
data("Glass")
data2 <- Glass

C1_LABEL <- 1
C2_LABEL <- -1

# length(which(data2$Type == 1)): 70
# length(which(data2$Type == 2)): 76
# 3: 17
# 4: 0
# 5: 13
# 6: 9
# 7: 29

# classe majotaria: 2 -> C1
# classe C2: todas as demais

# Realiza o tratamento dos dados para remoção de NA
data2 <- data2[complete.cases(data2),]

start_variables_column <- 1
end_variables_column <- 9
label_column <- 10

X <- data.matrix(data2[, start_variables_column:end_variables_column])

N <- nrow(X) # numero de amostras
n <- ncol(X) # numero de variaveis

Y <- matrix(NA, nrow = N, ncol = 1)

for (i in 1:N) {
  # ultima coluna é a coluna com a classe
  if (data2[i, label_column] == 2) {
    Y[i] <- C1_LABEL
  } else {
    Y[i] <- C2_LABEL
  }
}

# junta tudo na matriz dos dados
all_data <- cbind(X, Y)

# embaralha a matriz dos dados de entrada - remove bias de coleta
all_data <- all_data[sample.int(nrow(all_data)), ]

n_folds <- 10
fold_size <- floor(N / n_folds)

acc_array <- c()

acc_mean_array <- c()

h_list <- seq(0.1, 2, 0.1)
C <- 15

for (h in h_list) {
  for (fold in 1:n_folds) {
    num_of_corrects <- 0
    
    start_index <- fold_size * (fold - 1) + 1
    end_index <- start_index + fold_size
    
    data_for_test <- all_data[start_index:end_index, ]
    X_test <- data_for_test[, 1:n]
    Y_test <- data_for_test[, n+1]
    
    data_for_train <- all_data[-(start_index:end_index), ]
    X_train <- data_for_train[, 1:n]
    Y_train <- data_for_train[, n+1]
    
    # treina a rede
    svm_train <- ksvm(X_train, Y_train, type='C-bsvc', kernel='rbfdot',
                      kpar=list(sigma=h), C=C)
    
    # svm_train <- ksvm(X_train, Y_train, type='C-bsvc', kernel='besseldot',
    #                    C=C)
    
    yhat <- predict(svm_train, X_test, type="response")
    
    for (i in 1:length(Y_test)) {
      if (yhat[i] == Y_test[i]) {
        num_of_corrects <- num_of_corrects + 1
      }
    }
    
    acc_array <- c(acc_array, num_of_corrects / fold_size * 100)
    
    if (fold != 1) next
    
    # plota a projeção com o codigo do Braga
    dall <- as.matrix(dist(X_train, diag=TRUE, upper=TRUE))
    kall <- exp(-(dall*dall) / (2 * (h)^2))
    
    n_train <- nrow(X_train)
    k11 <- kall[(1:(n_train/2)), (1:(n_train/2))]
    k12 <- kall[(1:(n_train/2)), ((n_train / 2) + 1):n_train]
    k21 <- kall[((n_train / 2) + 1):n_train, (1:(n_train/2))]
    k22 <- kall[((n_train / 2) + 1):n_train, ((n_train / 2) + 1):n_train]
    
    p11 <- rowSums(k11)
    p12 <- rowSums(k12)
    p21 <- rowSums(k21)
    p22 <- rowSums(k22)
    
    p1 <- cbind(p11, p12) / (n_train / 2)
    p2 <- cbind(p21, p22) / (n_train / 2)
    
    pall <- cbind(p1, p2)
    
    plot(p1[,1], p1[,2], col = "red", xlim=c(0,1), ylim=c(0,1), lwd = 2)
    par(new=T)
    plot(p2[,1], p2[,2], col = "blue", xlim=c(0,1), ylim=c(0,1), lwd = 2)
  }
  
  # print(paste(mean(acc_array), " +/- ", sd(acc_array)))
  
  acc_mean_array <- c(acc_mean_array, mean(acc_array))
}


df <- data.frame(seq(1, 10, 1), round(acc_array, 2))
colnames(df) <- c("Fold", "Acurácia (%)")
ft <- flextable(df)
ft <- align(ft, align = "center", part = "all")

plot(h_list, acc_mean_array, xlab = "h", ylab="Acurácia Média", type="b", lty="dashed"
     , lwd = 2)


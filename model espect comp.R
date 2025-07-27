library(ir)
library(mixOmics)
library(ggplot2)

X <- readRDS("./terminamossiosi.rds")
x <- data

dats <- readRDS("./datos_analisis.rds")

respuesta <- dats[, c(1:2)]

sa <- reshape2::melt(scale(X))       # estandarizacion  de datos( Resta la media y divide por la desviación estándar para cada columna.)
ggplot(sa, aes(value)) + geom_density()   #melt modifica el fromato de los datos de ancho a largo·


idx <- sample(1:nrow(X), nrow(X) * 0.7)

trainx <- X[idx, ]            # el 70 %espectros mas representativos de todas las señales
y.train <- respuesta[idx, ]  # 70%  y mas concentraciones de los espectros mas repesentativso#

testx <- X[-idx, ]               # 30 % de los esepctro para la pureba de rmns
y.test <- respuesta[-idx, ]      # 30 % de las concentraciones  para la pureba de rmns

# componetes


cargas <- plotLoadings(mod, comp = 1, contrib = "max")    #cargas de la primera componete que solo muestra los valores maximos

ondas1 <- colnames(trainx)[abs(cargas$X) > 0.025]


cargas <- plotLoadings(mod, comp = 2, contrib = "max")    #cargas de la primera componete que solo muestra los valores maximos

ondas2 <- colnames(trainx)[abs(cargas$X) > 0.01]


cargas <- plotLoadings(mod, comp = 3, contrib = "max")    #cargas de la primera componete que solo muestra los valores maximos

ondas3 <- colnames(trainx)[abs(cargas$X) > 0.015]


cargas <- plotLoadings(mod, comp = 4, contrib = "max")    #cargas de la primera componete que solo muestra los valores maximos
ondas4 <- colnames(trainx)[abs(cargas$X) > 0.015]


unicosonda <- unique(c(ondas1, ondas2, ondas3, ondas4))

# los analisis mas repesntativos de cada uno de las componentes #

## girnd por  cuadirula









### moldeo pedicto mejor

###====spls=====

# modelo  de entrenamiento ###

ls.liver <- spls(
  X = trainx ,
  Y = y.train,
  ncomp = 3,
  mode = 'regression',
  scale = T
)

p <- perf(
  ls.liver,
  validation = "Mfold",
  folds = 5,
  progressBar = T,
  nrepeat = 3
)
plot(p)

####==tueno datos==#
tuenado <- tune.spls(
  trainx,
  y.train,
  ncomp = 3,
  mode = "regression",
  test.keepX =  c(100, 150, 200),
  validation = "Mfold",
  folds = 4,
  measure = "cor",
  scale = TRUE,
  progressBar = TRUE,
  tol = 1e-06,
  max.iter = 100,
  near.zero.var = FALSE,
  nrepeat = 2,
  multilevel = NULL
)

plot(tuenado)
tuenado$choice.keepX


###############



tuenado <- tune.spls(
  trainx,
  y.train,
  ncomp = 3,
  mode = "regression",
  test.keepX =  c(25, 85, 175),
  validation = "Mfold",
  folds = 4,
  measure = "cor",
  scale = TRUE,
  progressBar = TRUE,
  tol = 1e-06,
  max.iter = 100,
  near.zero.var = FALSE,
  nrepeat = 2,
  multilevel = NULL
)

plot(tuenado)
tuenado$choice.keepX



####============================== cambiao por mi =========================================================================================

library(mixOmics)

X <- readRDS("./terminamossiosi.rds")
reshape2::melt(scale(X))

Y <- readRDS("./datos_analisis.rds")
Y <- Y[, c(1:2)]

# Configuración básica
set.seed(123)
keepX_values <- seq(7, 50, by = 5)         # Valores de 10 a 100 en pasos de 10
folds <- 5                                    # Número de folds para validación cruzada




# Función optimizada para bajo consumo de recursos
eval_keepX_combination <- function(keepX_vec, X_, Y_)
{
  # Modelo ligero
  X <- X_
  Y <- Y_
  
  
  
  model <- spls(
    X,
    Y ,
    ncomp = 3,
    keepX = keepX_vec,
    mode = "regression",
    scale = TRUE,
    max.iter = 50,
    tol = 1e-04
  )                                      # Parámetros reducidos
  
  # Validación cruzada simplificada
  fold_indices <- split(sample(1:nrow(X)), rep(1:folds, length.out = nrow(X)))
  #
  # ######3 debugear
  # test_idx <-(40)
  # train_X <- X[-test_idx, , drop = FALSE]
  # train_Y <- Y[-test_idx,]
  # test_X <- X[test_idx, , drop = FALSE]
  # test_Y <- Y[test_idx,]
  #
  # fit <- spls(train_X, train_Y, ncomp = 3, keepX = keepX_vec,
  #             mode = "regression", scale = TRUE,
  #             max.iter = 30, tol = 1e-08)
  #
  # pred <- predict(fit, newdata = test_X)$predict[,,3]
  # # mean((test_Y - pred)^2)
  #
  #
  # ######========================================================================================
  
  
  
  
  mse <- (lapply(fold_indices, function(test_idx)
    
  {
    # test_idx<-34
    
    
    train_X <- X[-test_idx, , drop = FALSE]
    train_Y <- Y[-test_idx, ]
    test_X <- X[test_idx, , drop = FALSE]
    test_Y <- Y[test_idx, ]
    
    fit <- spls(
      train_X,
      train_Y,
      ncomp = 3,
      keepX = keepX_vec,
      mode = "regression",
      scale = TRUE,
      max.iter = 30,
      tol = 1e-04
    )
    
    
    pred1_ET <- predict(fit, newdata = test_X)$predict[, 1, 1]
    pred2_ET <- predict(fit, newdata = test_X)$predict[, 1, 2]
    pred3_ET <- predict(fit, newdata = test_X)$predict[, 1, 3]
    
    pred1_MT <- predict(fit, newdata = test_X)$predict[, 2, 1]
    pred2_MT <- predict(fit, newdata = test_X)$predict[, 2, 2]
    pred3_MT <- predict(fit, newdata = test_X)$predict[, 2, 3]
    
    
    MEANet1 <- mean((test_Y[, 1] - pred1_ET)^2)
    MEANet2 <- mean((test_Y [, 1] - pred2_ET)^2)
    MEANet3 <-  mean((test_Y [, 1] - pred3_ET)^2)
    
    MEANmt1 <- mean((test_Y[, 2] - pred1_MT)^2)
    MEANmt2 <- mean((test_Y[, 2] - pred2_MT)^2)
    MEANmt3 <-  mean((test_Y[, 2] - pred3_MT)^2)
    
    
    retunmET <- c(MEANet1, MEANet2, MEANet3)
    retunmMT <- c(MEANmt1, MEANmt2, MEANmt3)
    
    MS <- cbind(ETANOL = retunmET, METANOL = retunmMT)
  }))
  
  ######## la hace una lsita y deslista y hace una media de los 5 fols de cada componete #####
  
  etanol_comp1 <- mean(unlist(lapply(mse, function(x)
    x[1, 1])))
  etanol_comp2 <- mean(unlist(lapply(mse, function(x)
    x[2, 1])))
  etanol_comp3 <- mean(unlist(lapply(mse, function(x)
    x[3, 1])))
  
  metanol_comp1 <- mean(unlist(lapply(mse, function(x)
    x[1, 2])))
  metanol_comp2 <- mean(unlist(lapply(mse, function(x)
    x[2, 2])))
  metanol_comp3 <- mean(unlist(lapply(mse, function(x)
    x[3, 2])))
  
  mat_final <- cbind(
    c(etanol_comp1, etanol_comp2, etanol_comp3),
    c(metanol_comp1, metanol_comp2, metanol_comp3)
  )
  
  return(mat_final)
}
##etanol###

compET_1 <- mean(c(mse[[1]][1, 1], mse[[2]][1, 1], mse[[3]][1, 1], mse[[4]][1, 1], mse[[5]][1, 1]))


compET_2 <- mean(c(mse[[1]][2, 1], mse[[2]][2, 1], mse[[3]][2, 1], mse[[4]][2, 1], mse[[5]][2, 1]))

compET_3 <- mean(c(mse[[1]][3, 1], mse[[2]][3, 1], mse[[3]][3, 1], mse[[4]][3, 1], mse[[5]][3, 1]))


####metanol##

compMT_1 <- mean(c(mse[[1]][1, 2], mse[[2]][1, 2], mse[[3]][1, 2], mse[[4]][1, 2], mse[[5]][1, 2]))


compMT_2 <- mean(c(mse[[1]][2, 2], mse[[2]][2, 2], mse[[3]][2, 2], mse[[4]][2, 2], mse[[5]][2, 2]))

compMT_3 <- mean(c(mse[[1]][3, 2], mse[[2]][3, 2], mse[[3]][3, 2], mse[[4]][3, 2], mse[[5]][3, 2]))








# Grid search por niveles (para reducir memoria)
best_mse <- Inf
best_keepX <- c(10, 10, 10)  # Valores iniciales
keepX_values <- seq(7, 50, by = 5)

# Bucle por niveles (optimizado para RAM)
for (keep1 in keepX_values) {
  for (keep2 in keepX_values) {
    for (keep3 in keepX_values) {
      current_keepX <- c(keep1, keep2, keep3)
      current_mse <- eval_keepX_combination(current_keepX, X, Y)
      
      current_mse <- mean(current_mse[, 1])
      if (current_mse < best_mse) {
        best_mse <- current_mse
        best_keepX <- current_keepX
        cat(
          sprintf(
            "Nuevo óptimo: keepX = [%d, %d, %d], MSE = %.4f\n",
            keep1,
            keep2,
            keep3,
            best_mse
          )
        )
      }
    }
  }
}

# Resultados finales
cat("\nMEJOR COMBINACIÓN ENCONTRADA:\n")
cat("keepX para comp1:", best_keepX[1], "\n")
cat("keepX para comp2:", best_keepX[2], "\n")
cat("keepX para comp3:", best_keepX[3], "\n")
cat("MSE mínimo alcanzado:", best_mse, "\n")


vsfdgsfdg


#########================================con dos componetes ============================================================================================================#########


library(mixOmics)

X <- readRDS("./terminamossiosi.rds")
reshape2::melt(scale(X))

Y <- readRDS("./datos_analisis.rds")
Y <- Y[, c(1:2)]

# Configuración básica
set.seed(123)
folds <- 5                                    # Número de folds para validación cruzada




# Función optimizada para bajo consumo de recursos
eval_keepX_combination <- function(keepX_vec, X_, Y_)
{
  # Modelo ligero
  X <- X_
  Y <- Y_
  

  
  model <- spls(
    X,
    Y ,
    ncomp = 2,
    keepX = keepX_vec,
    mode = "regression",
    scale = TRUE,
    max.iter = 1000,
    tol = 1e-04
  )                                      # Parámetros reducidos
  
  # Validación cruzada simplificada
  fold_indices <- split(sample(1:nrow(X)), rep(1:folds, length.out = nrow(X)))
  
  mse <- (lapply(fold_indices, function(test_idx)
    
  {
    # test_idx<-34
    
    
    train_X <- X[-test_idx, , drop = FALSE]
    train_Y <- Y[-test_idx, ]
    test_X <- X[test_idx, , drop = FALSE]
    test_Y <- Y[test_idx, ]
    
    fit <- spls(
      train_X,
      train_Y,
      ncomp = 2,
      keepX = keepX_vec,
      mode = "regression",
      scale = TRUE,
      max.iter = 1000,
      tol = 1e-04
    )
    
    
    pred1_ET <- predict(fit, newdata = test_X)$predict[, 1, 1]
    pred2_ET <- predict(fit, newdata = test_X)$predict[, 1, 2]
    
    
    pred1_MT <- predict(fit, newdata = test_X)$predict[, 2, 1]
    pred2_MT <- predict(fit, newdata = test_X)$predict[, 2, 2]
    
    
    
    MEANet1 <- mean((test_Y[, 1] - pred1_ET)^2)
    MEANet2 <- mean((test_Y [, 1] - pred2_ET)^2)
    
    MEANmt1 <- mean((test_Y[, 2] - pred1_MT)^2)
    MEANmt2 <- mean((test_Y[, 2] - pred2_MT)^2)
    
    
    
    retunmET <- c(MEANet1, MEANet2)
    retunmMT <- c(MEANmt1, MEANmt2)
    
    MS <- cbind(ETANOL = retunmET, METANOL = retunmMT)
  }))
  
  ######## la hace una lsita y deslista y hace una media de los 5 fols de cada componete #####
  
  etanol_comp1 <- mean(unlist(lapply(mse, function(x)
    x[1, 1])))
  etanol_comp2 <- mean(unlist(lapply(mse, function(x)
    x[2, 1])))
  
  metanol_comp1 <- mean(unlist(lapply(mse, function(x)
    x[1, 2])))
  metanol_comp2 <- mean(unlist(lapply(mse, function(x)
    x[2, 2])))
  
  
  mat_final <- cbind(c(etanol_comp1, etanol_comp2),
                     c(metanol_comp1, metanol_comp2))
  
  return(mat_final)
}
##etanol###

#compET_1 <- mean(c(mse[[1]][1, 1],
#                   mse[[2]][1, 1],
#                  mse[[3]][1, 1],
#                 mse[[4]][1, 1],
#                mse[[5]][1, 1]))


#compET_2 <- mean(c(mse[[1]][2, 1],
#                   mse[[2]][2, 1],
#                   mse[[3]][2, 1],
#                   mse[[4]][2, 1],
#                  mse[[5]][2, 1]))


####metanol##

#compMT_1 <- mean(c(mse[[1]][1, 2],
#                   mse[[2]][1, 2],
#                   mse[[3]][1, 2],
#                   mse[[4]][1, 2],
#                   mse[[5]][1, 2]))


#compMT_2 <- mean(c(mse[[1]][2, 2],
#                   mse[[2]][2, 2],
#                   mse[[3]][2, 2],
#                   mse[[4]][2, 2],
#                   mse[[5]][2, 2]))









# Grid search por niveles (para reducir memoria)
best_mse <- Inf
best_keepX <- c(10, 10, 10)  # Valores iniciales
keepX_values <- seq(20, 80, by = 5)         # Valores de 10 a 100 en pasos de 10

LISTA <- list()
# Bucle por niveles (optimizado para RAM)
for (keep1 in keepX_values) {
  for (keep2 in keepX_values) {
    current_keepX <- c(keep1, keep2)
    current_mse_ <- eval_keepX_combination(current_keepX, X, Y)
    
    current_mse <- mean(current_mse_)
    if (current_mse < best_mse) {
      best_mse <- current_mse
      best_keepX <- current_keepX
      cat(sprintf(
        "Nuevo óptimo: keepX = [%d, %d], MSE = %.4f\n",
        keep1,
        keep2,
        best_mse
      ))
    }
  }
  LISTA[[keep1]] <- current_mse_
}


# Resultados finales

cat("\nMEJOR COMBINACIÓN ENCONTRADA:\n")
cat("keepX para comp1:", best_keepX[1], "\n")
cat("keepX para comp2:", best_keepX[2], "\n")
cat("MSE mínimo alcanzado:", best_mse, "\n")




####============================== =========================================================================================

x.melt <- reshape2::melt(scale(trainx))
ggplot(x.melt, aes(x = value)) + geom_density()



x.melt <- reshape2::melt(scale(y.train))

ggplot(x.melt, aes(x = value)) + geom_density()

x.melt <- reshape2::melt(y.train)
ggplot(x.melt, aes(x = value)) + geom_density()

y.train_scaled <- scale(y.train)

# Guarda los parámetros de escalado
mean_y <- attr(y.train_scaled, "scaled:center")
sd_y <- attr(y.train_scaled, "scaled:scale")

### MODELO TUNEADO ####
ls.liver <- spls(
  # modelo de regrecion
  X = trainx ,
  Y = y.train,
  ncomp = 3,
  mode = 'regression',
  scale = T,
  keepX = tuenado$choice.keepX
)

p <- perf(
  ls.liver,
  validation = "Mfold",
  folds = 5,
  progressBar = T,
  nrepeat = 3
)
plot(p)

### predicion ###

ped <- predict(ls.liver, newdata = testx)

ped$predict
rmse <- function(y_real, y_pred) {
  sqrt(mean((y_real - y_pred)^2))
}

(et1_CM <- rmse(y.test[, 1], ped$predict[, 1, 1] + mean_y[1]))
(et2_CM <- rmse(y.test[, 1], ped$predict[, 1, 2] + mean_y[1]))
(et3_CM <- rmse(y.test[, 1], ped$predict[, 1, 3] + mean_y[1]))

(et1_C <- rmse(y.test[, 1], ped$predict[, 1, 1]))
(et2_C <- rmse(y.test[, 1], ped$predict[, 1, 2]))
(et3_C <- rmse(y.test[, 1], ped$predict[, 1, 3]))

(mt1_CM <- rmse(y.test_corte[, 2], ped$predict[, 2, 1] + mean_y[2]))
(mt2_CM <- rmse(y.test_corte[, 2], ped$predict[, 2, 2] + mean_y[2]))
(mt3_CM <- rmse(y.test_corte[, 2], ped$predict[, 2, 3] + mean_y[2]))

(mt1_C <- rmse(y.test_corte[, 2], ped$predict[, 2, 1]))
(mt2_C <- rmse(y.test_corte[, 2], ped$predict[, 2, 2]))
(mt3_C <- rmse(y.test_corte[, 2], ped$predict[, 2, 3]))




#recojer el mejor numoero de componetes para el modelo en cuestion#

cargascomponete1 <- plotLoadings(moldeon2compun, comp = 1, contrib = "max")
ondascomponete1 <- colnames(trainx)[abs(cargascomponete1$X) > 0.025]

cargascomponete2 <- plotLoadings(moldeon2compun, comp = 2, contrib = "max")
ondascomponete2 <- colnames(trainx)[abs(cargascomponete2$X) > 0.025]

ONDASUNICAS <- unique(ondascomponete1, ondascomponete2)

ONDaspls <- pls(trainx[, ONDASUNICAS], y.train, ncomp = 2, scale = T)  # sin estandarizacion
p <- perf(ONDaspls,
          folds = 5,
          nrepeat = 3,
          progressBar = T)
plot(p)




mod <- spls(scale(trainx[, ondas]), y.train, ncomp = 15, scale = F)
p <- perf(mod,
          folds = 5,
          nrepeat = 3,
          progressBar = T)
plot(p)







ondas2 <- colnames(y.train)[abs(cargas$X) > 0.025]

tuenado <- tune.spls(
  X,
  Y,
  ncomp = 14,
  test.keepX = c(5, 15, 30),
  already.tested.X,
  validation = "Mfold",
  folds = 10,
  measure = "MSE",
  scale = TRUE,
  progressBar = TRUE,
  tol = 1e-06,
  max.iter = 100,
  near.zero.var = FALSE,
  nrepeat = 1,
  multilevel = NULL,
  light.output = TRUE,
  cpus
)


vars_finalkes <- tuneado$choice_keep_x

mod <- spls(scale(trainx[, ondas]), y.train, ncomp = 14, scale = F)


mod <- predict(mod, newdata = testx)


##### tuenado ####

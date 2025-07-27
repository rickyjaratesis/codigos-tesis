library(ir)
library(mixOmics)
library(ggplot2)

#X <- readRDS("./terminamossiosi.rds")

dats_corte <- readRDS("./datos_analisis_corte.rds")

respuesta_corte <- dats_corte[,c(1:2)]

espectr_corte<-fun_to_spectra(dats_corte)
X <- espectr_corte
sa <- reshape2::melt(scale(X))       # estandarizacion  de datos( Resta la media y divide por la desviación estándar para cada columna.)
ggplot(sa,aes(value))+geom_density()   #melt modifica el fromato de los datos de ancho a largo·


idx <- sample(1:nrow(X),nrow(X) *0.7)  

trainx <- X[idx,]           # el 70 %espectros mas representativos de todas las señales 
y.train_corte <- respuesta_corte[idx,]  # 70%  y mas concentraciones de los espectros mas repesentativso#

testx <- X[-idx,]               # 30 % de los esepctro para la pureba de rmns
y.test_corte <- respuesta_corte[-idx,]      # 30 % de las concentraciones  para la pureba de rmns


# Configuración básica
set.seed(123)
keepX_values <- seq(60, 110, by = 5)         # Valores de 10 a 100 en pasos de 10
folds <- 5                                    # Número de folds para validación cruzada

# Función optimizada para bajo consumo de recursos
eval_keepX_combination_3comps <- function(keepX_vec, X_, Y_)
{ # Modelo ligero
  X<-X_
  Y<-Y_
  
  
  model <- spls(X,Y , ncomp = 3, keepX = keepX_vec,
                mode = "regression", scale = TRUE,
                max.iter = 50, tol = 1e-04)                                      # Parámetros reducidos
  
  # Validación cruzada simplificada
  fold_indices <- split(sample(1:nrow(X)), rep(1:folds, length.out = nrow(X)))
 
  
  
  
  
  mse <- (lapply(fold_indices, function(test_idx) 
    
  {
    # test_idx<-34
    
    
    train_X <- X[-test_idx, , drop = FALSE]
    train_Y <- Y[-test_idx,]
    test_X <- X[test_idx, , drop = FALSE]
    test_Y <- Y[test_idx,]
    
    fit <- spls(train_X, train_Y, ncomp = 3, keepX = keepX_vec,
                mode = "regression", scale = TRUE,
                max.iter = 30, tol = 1e-04) 
    
    
    pred1_ET <- predict(fit, newdata = test_X)$predict[,1,1]
    pred2_ET <- predict(fit, newdata = test_X)$predict[,1,2]
    pred3_ET<- predict(fit, newdata = test_X)$predict[,1,3]
    
    pred1_MT<- predict(fit, newdata = test_X)$predict[,2,1]
    pred2_MT<- predict(fit, newdata = test_X)$predict[,2,2]
    pred3_MT<- predict(fit, newdata = test_X)$predict[,2,3]
    
    
    MEANet1<- mean((test_Y[,1] - pred1_ET)^2)
    MEANet2<- mean((test_Y [,1]- pred2_ET)^2)
    MEANet3<-  mean((test_Y [,1]- pred3_ET)^2)
    
    MEANmt1<- mean((test_Y[,2] - pred1_MT)^2)
    MEANmt2<- mean((test_Y[,2] - pred2_MT)^2)
    MEANmt3<-  mean((test_Y[,2]- pred3_MT)^2)
    
    
    retunmET<-c(MEANet1,MEANet2,MEANet3)
    retunmMT <-c(MEANmt1,MEANmt2,MEANmt3)
    
    MS<- cbind(ETANOL=retunmET,METANOL=retunmMT)
  }
  ))
  
  ######## la hace una lsita y deslista y hace una media de los 5 fols de cada componete #####
  
  etanol_comp1 <-mean(unlist( lapply(mse,function(x) x[1,1])))
  etanol_comp2 <-mean(unlist(lapply(mse,function(x) x[2,1])))
  etanol_comp3 <- mean(unlist(lapply(mse,function(x) x[3,1])))
  
  metanol_comp1 <- mean(unlist(lapply(mse,function(x) x[1,2])))
  metanol_comp2 <- mean(unlist(lapply(mse,function(x) x[2,2])))
  metanol_comp3 <- mean(unlist(lapply(mse,function(x) x[3,2])))
  
  mat_final <- cbind(c(etanol_comp1,etanol_comp2,etanol_comp3),c(metanol_comp1,metanol_comp2,metanol_comp3))
  
  return(mat_final)
}




# Grid search por niveles (para reducir memoria)
best_mse <- Inf
best_keepX <- c(10, 10, 10)  # Valores iniciales

# Bucle por niveles (optimizado para RAM)
for(keep1 in keepX_values) {
  for(keep2 in keepX_values) {
    for(keep3 in keepX_values) {
      current_keepX <- c(keep1, keep2, keep3)
      current_mse <- eval_keepX_combination_3comps(current_keepX, X, Y)
      
      current_mse <- mean(current_mse)
      if(current_mse < best_mse) {
        best_mse <- current_mse
        best_keepX <- current_keepX
        cat(sprintf("Nuevo óptimo: keepX = [%d, %d, %d], MSE = %.4f\n",
                    keep1, keep2, keep3, best_mse))
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



#########================================con dos componetes ============================================================================================================#########

# Función optimizada para bajo consumo de recursos
eval_keepX_combination2comps <- function(keepX_vec, X_, Y_)
{ # Modelo ligero
  X<-X_
  Y<-Y_
  
  # keepX_vec<-c(10,100)
  
  model <- spls(X,Y , ncomp = 2, keepX = keepX_vec,
                mode = "regression", scale = TRUE,
                max.iter = 1000, tol = 1e-04)                                      # Parámetros reducidos
  
  # Validación cruzada simplificada
  fold_indices <- split(sample(1:nrow(X)), rep(1:folds, length.out = nrow(X)))
  
  mse <- (lapply(fold_indices, function(test_idx) 
    
  {
    
    
    
    train_X <- X[-test_idx, , drop = FALSE]
    train_Y <- Y[-test_idx,]
    test_X <- X[test_idx, , drop = FALSE]
    test_Y <- Y[test_idx,]
    
    fit <- spls(train_X, train_Y, ncomp = 2, keepX = keepX_vec,
                mode = "regression", scale = TRUE,
                max.iter = 1000, tol = 1e-04) 
    
    
    pred1_ET <- predict(fit, newdata = test_X)$predict[,1,1]
    pred2_ET <- predict(fit, newdata = test_X)$predict[,1,2]
    
    
    pred1_MT<- predict(fit, newdata = test_X)$predict[,2,1]
    pred2_MT<- predict(fit, newdata = test_X)$predict[,2,2]
    
    
    
    MEANet1<- mean((test_Y[,1] - pred1_ET)^2)
    MEANet2<- mean((test_Y [,1]- pred2_ET)^2)
    
    MEANmt1<- mean((test_Y[,2] - pred1_MT)^2)
    MEANmt2<- mean((test_Y[,2] - pred2_MT)^2)
    
    
    
    retunmET<-c(MEANet1,MEANet2)
    retunmMT <-c(MEANmt1,MEANmt2)
    
    MS<- cbind(ETANOL=retunmET,METANOL=retunmMT)
  }
  ))
  
  ######## la hace una lsita y deslista y hace una media de los 5 fols de cada componete #####
  
  etanol_comp1 <-mean(unlist( lapply(mse,function(x) x[1,1])))
  etanol_comp2 <-mean(unlist(lapply(mse,function(x) x[2,1])))
  
  metanol_comp1 <- mean(unlist(lapply(mse,function(x) x[1,2])))
  metanol_comp2 <- mean(unlist(lapply(mse,function(x) x[2,2])))
  
  
  mat_final <- cbind(c(etanol_comp1,etanol_comp2),c(metanol_comp1,metanol_comp2))
  
  return(mat_final)
}



# Grid search por niveles (para reducir memoria)
best_mse <- Inf
best_keepX <- c(10, 10, 10)  # Valores iniciales

# Bucle por niveles (optimizado para RAM)
for(keep1 in keepX_values) {
  for(keep2 in keepX_values){
    current_keepX <- c(keep1, keep2)
    current_mse <- eval_keepX_combination2comps(current_keepX, X, Y)
    
    current_mse <- mean(current_mse)
    if(current_mse < best_mse) {
      best_mse <- current_mse
      best_keepX <- current_keepX
      cat(sprintf("Nuevo óptimo: keepX = [%d, %d], MSE = %.4f\n",
                  keep1, keep2, best_mse))
    }
  }
}


# Resultados finales
cat("\nMEJOR COMBINACIÓN ENCONTRADA:\n")
cat("keepX para comp1:", best_keepX[1], "\n")
cat("keepX para comp2:", best_keepX[2], "\n")
cat("MSE mínimo alcanzado:", best_mse, "\n")




##===================================================================================================================================


keepX_vec_final3 <- LOQUETEDE
# modelo  de entrenamiento ###

ls.liver <-spls(        
  X =trainx ,
  Y = y.train,
  ncomp = 3,
  mode = 'regression',scale = T, keepX = keepX_vec_final3
)


keepX_vec_final2 <- LOQUETEDE
# modelo  de entrenamiento ###

ls.liver <-spls(        
  X =trainx ,
  Y = y.train,
  ncomp = 2,
  mode = 'regression',scale = T, keepX = keepX_vec_final2
)



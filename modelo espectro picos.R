library(ir)

library(mixOmics)
library(ggplot2)

X <- readRDS("./terminamossiosi.rds")

dats_corte <- readRDS("./datos_analisis_corte .rds")
espectr_corte<-fun_to_spectra(dats_corte)
X<-espectr_corte


respuesta_corte <- dats_corte[,c(1:2)]

sa <- reshape2::melt(scale(X))       # estandarizacion  de datos( Resta la media y divide por la desviación estándar para cada columna.)
ggplot(sa,aes(value))+geom_density()   #melt modifica el fromato de los datos de ancho a largo·


idx <- sample(1:nrow(X),nrow(X) *0.7)  

trainx <- X[idx,]           # el 70 %espectros mas representativos de todas las señales 
y.train_corte <- respuesta_corte[idx,]  # 70%  y mas concentraciones de los espectros mas repesentativso#

testx <- X[-idx,]               # 30 % de los esepctro para la pureba de rmns
y.test_corte <- respuesta_corte[-idx,]      # 30 % de las concentraciones  para la pureba de rmns


### en estas lienas de arriba se puede ver que la preducion del modelo cuadno se lo compara con la el 30 % de test
## en el cual para el etanol se un alto margen de error en cabio del metnaol no existe cambio lo que indica que no esta apredniedo el modelos

plotLoadings(mod)

# componetes 
###================================================================================================================================================
#                   GRID SEARCH DE MODELO PICOS 
#


library(mixOmics)


X <- readRDS("./terminamossiosi.rds")
reshape2::melt(scale(X))

Y <- readRDS("./datos_analisis_corte .rds")
Y <- Y[,c(1:2)]




# Configuración básica
set.seed(123)
keepX_values <- seq(7, 50, by = 5)         # Valores de 10 a 100 en pasos de 10
folds <- 5                                    # Número de folds para validación cruzada




# Función optimizada para bajo consumo de recursos
eval_keepX_combination <- function(keepX_vec, X_, Y_)
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
##etanol###

compET_1 <- mean(c(mse[[1]][1, 1],
                   mse[[2]][1, 1],
                   mse[[3]][1, 1],
                   mse[[4]][1, 1],
                   mse[[5]][1, 1]))


compET_2 <- mean(c(mse[[1]][2, 1],
                   mse[[2]][2, 1],
                   mse[[3]][2, 1],
                   mse[[4]][2, 1],
                   mse[[5]][2, 1]))

compET_3 <- mean(c(mse[[1]][3, 1],
                   mse[[2]][3, 1],
                   mse[[3]][3, 1],
                   mse[[4]][3, 1],
                   mse[[5]][3, 1]))


####metanol##

compMT_1 <- mean(c(mse[[1]][1, 2],
                   mse[[2]][1, 2],
                   mse[[3]][1, 2],
                   mse[[4]][1, 2],
                   mse[[5]][1, 2]))


compMT_2 <- mean(c(mse[[1]][2, 2],
                   mse[[2]][2, 2],
                   mse[[3]][2, 2],
                   mse[[4]][2, 2],
                   mse[[5]][2, 2]))

compMT_3 <- mean(c(mse[[1]][3, 2],
                   mse[[2]][3, 2],
                   mse[[3]][3, 2],
                   mse[[4]][3, 2],
                   mse[[5]][3, 2]))








# Grid search por niveles (para reducir memoria)
best_mse <- Inf
best_keepX <- c(10, 10, 10)  # Valores iniciales

# Bucle por niveles (optimizado para RAM)
for(keep1 in keepX_values) {
  for(keep2 in keepX_values) {
    for(keep3 in keepX_values) {
      current_keepX <- c(keep1, keep2, keep3)
      current_mse <- eval_keepX_combination(current_keepX, X, Y)
      
      current_mse <- mean(current_mse[,1])
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


library(mixOmics)

X <- readRDS("./terminamossiosi.rds")
reshape2::melt(scale(X))

Y <- readRDS("./datos_analisis_corte .rds")
Y <- Y[,c(1:2)]

# Configuración básica
set.seed(123)
keepX_values <- seq(20,80, by = 5)         # Valores de 10 a 100 en pasos de 10
folds <- 5                                    # Número de folds para validación cruzada




# Función optimizada para bajo consumo de recursos
eval_keepX_combination <- function(keepX_vec, X_, Y_)
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
    # test_idx<-34
    
    
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
##etanol###

#compET_1 <- mean(c(mse[[1]][1, 1],
 #                  mse[[2]][1, 1],
  #                 mse[[3]][1, 1],
   #                mse[[4]][1, 1],
    #               mse[[5]][1, 1]))


#compET_2 <- mean(c(mse[[1]][2, 1],
 #                  mse[[2]][2, 1],
  #                 mse[[3]][2, 1],
   #                mse[[4]][2, 1],
    #               mse[[5]][2, 1]))


####metanol##
#
#compMT_1 <- mean(c(mse[[1]][1, 2],
3#                  mse[[2]][1, 2],
#                mse[[3]][1, 2],
 #                  mse[[4]][1, 2],
#                   mse[[5]][1, 2]))
#

#compMT_2 <- mean(c(mse[[1]][2, 2],
 #                  mse[[2]][2, 2],
  #                 mse[[3]][2, 2],
   #                mse[[4]][2, 2],
    #               mse[[5]][2, 2]))







# Grid search por niveles (para reducir memoria)
best_mse <- Inf
best_keepX <- c(10, 10, 10)  # Valores iniciales

# Bucle por niveles (optimizado para RAM)
for(keep1 in keepX_values) {
  for(keep2 in keepX_values){
    current_keepX <- c(keep1, keep2)
    current_mse <- eval_keepX_combination(current_keepX, X, Y)
    
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






~##===================================================================================================================================



# modelo  de entrenamiento ###

ls.liver <-spls(        
  X =trainx ,
  Y = y.train,
  ncomp = 4,
  mode = 'regression',scale = T
)

p<-perf(ls.liver,validation = "Mfold",folds = 5,progressBar =T,nrepeat = 3)
plot(p)

####==tueno datos==#
tuenado_corte <- tune.spls(trainx,y.train_corte, ncomp = 3, mode="regression",
                     test.keepX =  c(100,150,200),
                     validation = "Mfold", folds = 4, measure = "cor", scale = TRUE,
                     progressBar = TRUE, tol = 1e-06, max.iter = 100, near.zero.var = FALSE,
                     nrepeat = 2, multilevel = NULL)

plot(tuenado_corte)
tuenado_corte$choice.keepX


x.melt <- reshape2::melt(scale(trainx))
ggplot(x.melt,aes(x=value))+geom_density()



x.melt <- reshape2::melt(scale(y.train_corte))
ggplot(x.melt,aes(x=value))+geom_density()

x.melt <- reshape2::melt(y.train_corte)
ggplot(x.melt,aes(x=value))+geom_density()

y.train_scaled_corte <- scale(y.train_corte)

# Guarda los parámetros de escalado
mean_yC <- attr(y.train_scaled_corte, "scaled:center")
sd_yC<- attr(y.train_scaled_corte, "scaled:scale")

### MODELO TUNEADO ####
ls.liver_corte<-spls(        # modelo de regrecion 
  X =trainx ,
  Y = y.train_corte,
  ncomp = 3,
  mode = 'regression',scale = T,keepX = tuenado$choice.keepX
)

p<-perf(ls.liver_corte,validation = "Mfold",folds = 5,progressBar =T,nrepeat = 3)
plot(p)

### predicion ###

ped<-predict(ls.liver_corte,newdata=testx)

ped$predict
rmse<-function(y_real, y_pred){sqrt(mean((y_real-y_pred)^2))}

(et1_PM <- rmse(y.test_corte[,1],ped$predict[,1,1]+mean_yC[1]))
(et2_PM<- rmse(y.test_corte[,1],ped$predict[,1,2]+mean_yC[1]))
(et3_PM<- rmse(y.test_corte[,1],ped$predict[,1,3]+mean_yC[1]))

(et1_P <- rmse(y.test_corte[,1],ped$predict[,1,1]))
(et2_P <- rmse(y.test_corte[,1],ped$predict[,1,2]))
(et3_P<- rmse(y.test_corte[,1],ped$predict[,1,3]))



(mt1_PM <- rmse(y.test[,2],ped$predict[,2,1]+mean_yC[2]))
(mt2_PM<- rmse(y.test[,2],ped$predict[,2,2]+mean_yC[2]))
(mt3_PM <- rmse(y.test[,2],ped$predict[,2,3]+mean_yC[2]))

(mt1_P <- rmse(y.test[,2],ped$predict[,2,1]))
(mt2_P<- rmse(y.test[,2],ped$predict[,2,2]))
(mt3_P <- rmse(y.test[,2],ped$predict[,2,3]))



#recojer el mejor numoero de componetes para el modelo en cuestion#

cargascomponete1 <- plotLoadings(moldeon2compun,comp =1,contrib = "max") 
ondascomponete1 <- colnames(trainx)[abs(cargascomponete1$X)>0.025]

cargascomponete2 <- plotLoadings(moldeon2compun,comp =2,contrib = "max") 
ondascomponete2 <- colnames(trainx)[abs(cargascomponete2$X)>0.025]

ONDASUNICAS <- unique(ondascomponete1,ondascomponete2 )

ONDaspls <- pls(trainx[,ONDASUNICAS],y.train,ncomp = 2,scale = T)  # sin estandarizacion 
p <- perf(ONDaspls,folds = 5,nrepeat = 3,progressBar = T)
plot(p)




mod <- spls(scale(trainx[,ondas]),y.train,ncomp = 15,scale = F)
p <- perf(mod,folds = 5,nrepeat = 3,progressBar = T)
plot(p)




mod <- spls(scale(trainx[,ondas]),y.train,ncomp = 15,scale = F)
p <- perf(mod,folds = 5,nrepeat = 3,progressBar = T)
plot(p)





ondas2 <- colnames(y.train)[abs(cargas$X)>0.025]

tuenado <- tune.spls(X, Y, ncomp = 14,
                     test.keepX = c(5, 15, 30), already.tested.X,
                     validation = "Mfold", folds = 10, measure = "MSE", scale = TRUE,
                     progressBar = TRUE, tol = 1e-06, max.iter = 100, near.zero.var = FALSE,
                     nrepeat = 1, multilevel = NULL, light.output = TRUE, cpus)


vars_finalkes <- tuneado$choice_keep_x

mod <- spls(scale(trainx[,ondas]),y.train,ncomp = 14,scale = F)


mod <- predict(mod,newdata = testx)

####==modelo solo con los picos mas importantes==####





archivos <- list.files("./tragoanalisis",pattern = ".CSV")






###modelo tuneado con los picos cortados####


tuenado2 <- tune.spls(coratdotrain,y.train, ncomp = 3, mode="regression",
                      test.keepX =  c(100,150,200),
                      validation = "Mfold", folds = 4, measure = "cor", scale = TRUE,
                      progressBar = TRUE, tol = 1e-06, max.iter = 100, near.zero.var = FALSE,
                      nrepeat = 2, multilevel = NULL)

plot(tuenado2)
tuenado$choice.keepX

p<-perf(ls.liver2,validation = "Mfold",folds = 5,progressBar =T,nrepeat = 3)
plot(p)


x.melt <- reshape2::melt(scale(coratdotrain))
ggplot(x.melt,aes(x=value))+geom_density()


ls.liver2 <-spls(        # modelo de regrecion 
  X =coratdotrain ,
  Y = y.train,
  ncomp = 3,
  mode = 'regression',scale = T,keepX = tuenado$choice.keepX
)

p<-perf(ls.liver2,validation = "Mfold",folds = 5,progressBar =T,nrepeat = 3)
plot(p)


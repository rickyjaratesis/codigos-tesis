library(ir)
library(dplyr)





install_package <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      message(paste("Instalando paquete:", package))
      
      # Verificar si es un paquete de Bioconductor
      if (package %in% c("mixOmics")) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")  # Instala BiocManager si no está instalado
        }
        BiocManager::install(package)  # Instala el paquete de Bioconductor
      } else {
        install.packages(package, dependencies = TRUE)  # Instala paquetes de CRAN
      }
      
      library(package, character.only = TRUE)  # Cargar el paquete después de instalarlo
    } else {
      message(paste("El paquete", package, "ya está instalado y cargado."))
    }
  }
}


install_package(
  c(
    "ggplot2",
    "tidyr",
    "dplyr",
    "factoextra",
    "ir",
    "signal",
    "baseline",
    "tibble",
    "rsm",
    "plotly",
    "BiocManager",
    "mixOmics"
  )
)

##funcion ALS #####

to_spectra <- function(X) {
  spectra <- as.data.frame((ir_flatten(X)))
  wn <- spectra[, 1]
  spectra <- spectra[, -1]
  rownames(spectra) <- wn
  colnames(spectra) <- X$sample_id
  spectra.t <- t(spectra)
}
fun_to_spectra <- function(X) {
  spectra <- as.data.frame((ir_flatten(X)))
  wn <- spectra[, 1]
  spectra <- spectra[, -1]
  rownames(spectra) <- wn
  colnames(spectra) <- X$sample_id
  spectra.t <- t(spectra)
  return(spectra.t)
}
fun_to_ir <- function(s, X) {
  l <- nrow(s)
  wn <- colnames(s)
  for (i in 1:l) {
    s_i <- as.numeric(s[i, ])
    df_i <- data.frame(x = rev(wn), y = rev(s_i))
    df_i <- as.data.frame(lapply(df_i, as.numeric))
    X$spectra[[i]] <- df_i
  }
  return(X)
}


pre_procesado <- function(espectro, out, derivada) {
  ### correcion de linea base ###3
  Correcionbbase2 <- baseline(espectro, method = "als")
  Correcionbbase2 <- Correcionbbase2@corrected
  Correcionbbase2[Correcionbbase2 < 0] <- 0
  contaminado_ir <- fun_to_ir(Correcionbbase2, out)
  plot(contaminado_ir) + geom_path(aes(color = sample_id)) + theme(legend.position = "none") +
    labs(x = "Longitud de onda", y = "Absorbancia", title = "Correcion de Linea Base (ALS)")
  
  #### suavizado de espectro ###3
  
  suavizado_ir <- contaminado_ir %>% ir_smooth(method = "sg",
                                               p = 3,
                                               n = 91,
                                               m = 0)
  suavizado_ir.matriz <- fun_to_spectra(suavizado_ir)
  
  plot(suavizado_ir) + geom_path(aes(color = sample_id)) + theme(legend.position = "none") +
    labs(x = "Longitud de onda", y = "Absorbancia", title = "Suavizado ")
  
  
  ###==normalizacion ==#
  normalized_spectra <- t(apply(suavizado_ir.matriz, 1, function(row)
    row / sum(row)))
  datos_bc.filtered.normalized <- fun_to_ir(s = normalized_spectra, X = contaminado_ir)
  
  
  plot(datos_bc.filtered.normalized) + geom_path(aes(color = sample_id)) + theme(legend.position = "none") +
    labs(x = "Longitud de onda", y = "Absorbancia", title = "Normalizado ")
  
  ##====filtro S-G===###
  
  
  datos_bc.spectra.filtered <-  as.data.frame(t(apply(normalized_spectra, 1, function(x)
    sgolayfilt(x, p = 2, n = 7, m = derivada))))
  
  colnames(datos_bc.spectra.filtered) <- colnames(suavizado_ir.matriz)
  
  datos_ir <- fun_to_ir(datos_bc.spectra.filtered, datos_bc.filtered.normalized)
  
  
  
  
  plot(datos_ir) + geom_path(aes(color = sample_id)) + theme(legend.position = "none") +
    labs(x = "Longitud de onda",
         y = "Absorbancia",
         title = paste("1", "Derivada"))
  
  
  return(datos_bc.spectra.filtered)
  
}

cortar_espectro <- function(espectro_,respuesta){
  pcx <- pls(X = espectro_,Y = respuesta)
  
  loads <- plotLoadings(pcx,block = 1,comp = 1)
  contrib <- abs(loads$`1`$importance)
  contrib95 <- quantile(contrib,.97)
  
  columnas1 <- colnames(espectro_)[contrib>contrib95]
  loads <- plotLoadings(pcx,block = 1,comp = 2)
  contrib <- abs(loads$`1`$importance)
  contrib95 <- quantile(contrib,.97)
  
  columnas2 <- colnames(espectro_)[contrib>contrib95]
  columnas_finales <- unique(c(columnas1,columnas2))
  espectro_final <- espectro_[,columnas_finales]
  
  return(espectro_final)
}


# Función optimizada para bajo consumo de recursos
eval_keepX_combination <- function(keepX_vec, X_, Y_,folds,scalar)
{
  # Modelo ligero
  X <- X_
  Y <- Y_
                          # Parámetros reducidos
  
  # Validación cruzada simplificada
  fold_indices <- split(sample(1:nrow(X)), rep(1:folds, length.out = nrow(X)))
 
  
  
  
  mse <- (lapply(fold_indices, function(test_idx)
    
  {
    # test_idx<-34
    
    
    train_X <- X[-test_idx, , drop = FALSE]
    train_Y <- Y[-test_idx, ]
    test_X <- X[test_idx, , drop = FALSE]
    test_Y <- Y[test_idx, ]
    
    if(scalar ==T){
      
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
      
    }else{
      
      
      
      fit <- spls(
        scale(train_X),
        train_Y,
        ncomp = 2,
        keepX = keepX_vec,
        mode = "regression",
        scale = F,
        max.iter = 1000,
        tol = 1e-04
      )
      
    }
 
    
    pred1_ET <- predict(fit, newdata = test_X)$predict[, 1, 1]
    pred2_ET <- predict(fit, newdata = test_X)$predict[, 1, 2]

    pred1_MT <- predict(fit, newdata = test_X)$predict[, 2, 1]
    pred2_MT <- predict(fit, newdata = test_X)$predict[, 2, 2]

    
    MEANet1 <- (mean((test_Y[, 1] - pred1_ET)^2))
    MEANet2 <- (mean((test_Y [, 1] - pred2_ET)^2))

    MEANmt1 <- (mean((test_Y[, 2] - pred1_MT)^2))
    MEANmt2 <- (mean((test_Y[, 2] - pred2_MT)^2))

    
    retunmET <- c(MEANet1, MEANet2)
    retunmMT <- c(MEANmt1, MEANmt2)
    
    MS <- cbind(ETANOL = retunmET, METANOL = retunmMT)
  }))
  

  etanol_comp1 <- mean(unlist(lapply(mse, function(x)
    x[1, 1])))
  etanol_comp2 <- mean(unlist(lapply(mse, function(x)
    x[2, 1])))

  metanol_comp1 <- mean(unlist(lapply(mse, function(x)
    x[1, 2])))
  metanol_comp2 <- mean(unlist(lapply(mse, function(x)
    x[2, 2])))

  
  mat_final <- cbind(
    c(etanol_comp1, etanol_comp2 ),
    c(metanol_comp1, metanol_comp2 )
  )
  
  return(mat_final)
}

fun_evaluar_modelo <- function(espectro_,keepX_values,respuesta_,folds = 5,scalar){
  # Grid search por niveles (para reducir memoria)
  best_mse <- Inf
  best_keepX <- c(10, 10, 10)  # Valores iniciales
  # Bucle por niveles (optimizado para RAM)
  for (keep1 in keepX_values) {
    for (keep2 in keepX_values) {
        current_keepX <- c(keep1, keep2)
        current_mse_<- eval_keepX_combination(current_keepX, espectro_, respuesta_,folds,scalar)
        
        current_mse <- mean(current_mse_[,2])
        if (current_mse < best_mse) {
          best_mse <- current_mse
          best_keepX <- current_keepX
          cat(
            sprintf(
              "Nuevo óptimo: keepX = [%d, %d], MSE = %.4f\n",
              keep1,
              keep2,
              best_mse
            )
          )
      }
    }
  }
  
  # Resultados finales
  cat("\nMEJOR COMBINACIÓN ENCONTRADA:\n")
  cat("keepX para comp1:", best_keepX[1], "\n")
  cat("keepX para comp2:", best_keepX[2], "\n")
  cat("MSE mínimo alcanzado:", best_mse, "\n")
  
  
  return(list(best_keepX,best_mse =best_mse,BEST = current_mse_))
  
}

fun_trainable <- function(X){
  
  set.seed(123)
  
  idx <- sample(1:nrow(X),nrow(X)*0.7)
  
  X.train <- X[idx,]
  X.test <- X[-idx,]
  
  return(list(train = X.train,test = X.test))
}


### PREPARAMOS DATOS


archivos <- list.files("./spects_para_analisis/", full.names = T)
aux <- ir::ir_import_csv(archivos)
aux$sample_id  <- sub("\\.[cC][sS][vV]$", "", aux$sample_id)

metadata <- read.csv("./datos de concentraciones y diluciones.csv", sep = ";")

colnames(metadata)[3] <- "sample_id"

metadata$sample_id <- as.character(metadata$sample_id)
aux$sample_id <- as.character(aux$sample_id)

out <- dplyr::inner_join(metadata, aux, by = "sample_id")




outsinespectro <- out[-22, ]

out1 <- ir_as_ir(outsinespectro)

espectro <- as.matrix(fun_to_spectra(out1))




nums_ondastrain <- as.numeric(colnames(espectro))
nums_ondastrain2 <- nums_ondastrain[nums_ondastrain > 800 &
                                      nums_ondastrain < 1250]
nums_ondastrain_2 <- as.character(nums_ondastrain2)
espectro_cortado <- espectro[, nums_ondastrain_2]




completo1 <- pre_procesado(espectro,out = out1,derivada = 1)
completo2 <- pre_procesado(espectro,out = out1,derivada = 2)


espectro_cortado1 <- pre_procesado(espectro_cortado,out = out1,derivada = 1)
espectro_cortado2 <- pre_procesado(espectro_cortado,out = out1,derivada = 2)

respuesta <- as.matrix(metadata[-22,1:2])

completo1 <- cortar_espectro(completo1,respuesta)
completo2 <- cortar_espectro(completo2,respuesta)
espectro_cortado1 <- cortar_espectro(espectro_cortado1,respuesta)
espectro_cortado2 <- cortar_espectro(espectro_cortado2,respuesta)


respuestas_modelos <- fun_trainable(respuesta)

respuesta.train <- respuestas_modelos$train
respuesta.test <- respuestas_modelos$test

completo1_modelos <- fun_trainable(completo1)
completo2_modelos <- fun_trainable(completo2)
espectro_cortado1_modelos <- fun_trainable(espectro_cortado1)
espectro_cortado2_modelos <- fun_trainable(espectro_cortado2)


#### evaluacion modelo completo 1 derivada[5, 10, 35]

length(colnames(completo1_modelos$train))
keepX_values_ <- seq(115,150,5)
res_primera_derivada <- fun_evaluar_modelo(espectro_ = completo1_modelos$train,
                                           keepX_values = keepX_values_,
                                           respuesta_ = respuesta.train,folds = 5,scalar = T)


#### evaluacion modelo completo 2 derivada [5, 30, 5], MSE = 3.3754
keepX_values_ <- seq(170,190,2)

res_segunda_derivada <- fun_evaluar_modelo(espectro_ = completo2_modelos$train,
                                          keepX_values = keepX_values_,
                                          respuesta_ = respuesta.train,folds = 5,scalar = T)



#### evaluacion modelo CORTADO 1 derivada[5, 20, 30]


keepX_values_ <- seq(1,20,2)


res_primera_derivada.cortado <- fun_evaluar_modelo(espectro_ = espectro_cortado1_modelos$train,
                                          keepX_values = keepX_values_,
                                          respuesta_ = respuesta.train,folds = 5,scalar = T)


#### evaluacion modelo CORTADO 2 derivada[5, 20, 30]
keepX_values_ <- seq(5,20,2)

res_segunda_derivada.cortado <- fun_evaluar_modelo(espectro_ = espectro_cortado2_modelos$train,
                                           keepX_values = keepX_values_,
                                           respuesta_ = respuesta.train,folds = 5,scalar = T)





tabla.normal <- data.frame(primera_Derivada = res_primera_derivada$best_mse,
                      primera_Derivada.cortado = res_primera_derivada.cortado$best_mse,
                      segunda_Derivada = res_segunda_derivada$best_mse,
                      psegunda_Derivada.cortado = res_primera_derivada.cortado$best_mse)



optimo_keep <- res_segunda_derivada[[1]]


fit <- spls(
  completo2_modelos$train,
  respuesta.train,
  ncomp = 2,
  keepX = optimo_keep,
  mode = "regression",
  scale = TRUE,
  max.iter = 1000,
  tol = 1e-04
)

pred <- predict(fit,newdata = (completo1_modelos$test))



pred1_ET <-pred$predict[, 1, 1]
pred2_ET <-pred$predict[, 1, 2]

pred1_MT <-pred$predict[, 2, 1]
pred2_MT <- pred$predict[, 2, 2]

# Suponiendo que ya tienes:
# - respuesta.train: matriz o data.frame con columnas Etanol y Metanol (training)
# - respuesta.test:  matriz o data.frame con las mismas dos columnas (test)
# - pred1_ET, pred2_ET, pred1_MT, pred2_MT: vectores de predicción en escala Z

# 1) Calcula media y desviación estándar de Y en el entrenamiento
mu    <- colMeans(respuesta.train)            
sigma <- apply(respuesta.train, 2, sd)        

# 2) “Desescala” las predicciones
pred1_ET_orig <- pred1_ET * sigma[1] + mu[1]
pred2_ET_orig <- pred2_ET * sigma[1] + mu[1]

pred1_MT_orig <- pred1_MT * sigma[2] + mu[2]
pred2_MT_orig <- pred2_MT * sigma[2] + mu[2]

# 3) Funciones para MSE y RMSE
mse  <- function(actual, pred) mean((actual - pred)^2)
rmse <- function(actual, pred) sqrt(mse(actual, pred))

# 4) Calcula MSE y RMSE en escala original
# — Etanol
mse_et1  <- mse(respuesta.test[,1], pred1_ET_orig)
rmse_et1 <- rmse(respuesta.test[,1], pred1_ET_orig)

mse_et2  <- mse(respuesta.test[,1], pred2_ET_orig)
rmse_et2 <- rmse(respuesta.test[,1], pred2_ET_orig)

# — Metanol
mse_mt1  <- mse(respuesta.test[,2], pred1_MT_orig)
rmse_mt1 <- rmse(respuesta.test[,2], pred1_MT_orig)

mse_mt2  <- mse(respuesta.test[,2], pred2_MT_orig)
rmse_mt2 <- rmse(respuesta.test[,2], pred2_MT_orig)

# 5) Organiza y muestra resultados
rownames <- c("Modelo 1", "Modelo 2")

results_mse <- cbind(
  ETANOL  = c(mse_et1,  mse_et2),
  METANOL = c(mse_mt1,  mse_mt2)
)
rownames(results_mse) <- rownames

results_rmse <- cbind(
  ETANOL  = c(rmse_et1, rmse_et2),
  METANOL = c(rmse_mt1, rmse_mt2)
)
rownames(results_rmse) <- rownames

# Imprime tablas
print("MSE en escala original:")
print(results_mse)

print("RMSE en escala original:")
print(results_rmse)

# ——————————————————————————————————————————————
# Opcional: calcular errores en la escala Z directamente
# (escalamos Y_test con mu y sigma de entrenamiento)
y1_test_z <- scale(respuesta.test[,1], center = mu[1], scale = F)
y2_test_z <- scale(respuesta.test[,2], center = mu[2], scale =F)

mse_z_et1  <- mse(y1_test_z, pred1_ET)
rmse_z_et1 <- rmse(y1_test_z, pred1_ET)

mse_z_et2  <- mse(y1_test_z, pred2_ET)
rmse_z_et2 <- rmse(y1_test_z, pred2_ET)

mse_z_mt1  <- mse(y2_test_z, pred1_MT)
rmse_z_mt1 <- rmse(y2_test_z, pred1_MT)

mse_z_mt2  <- mse(y2_test_z, pred2_MT)
rmse_z_mt2 <- rmse(y2_test_z, pred2_MT)

results_mse_z <- cbind(
  ETANOL  = c(mse_z_et1,  mse_z_et2),
  METANOL = c(mse_z_mt1,  mse_z_mt2)
)
rownames(results_mse_z) <- rownames

results_rmse_z <- cbind(
  ETANOL  = c(rmse_z_et1, rmse_z_et2),
  METANOL = c(rmse_z_mt1, rmse_z_mt2)
)
rownames(results_rmse_z) <- rownames

print("MSE en escala Z:")
print(results_mse_z)

print("RMSE en escala Z:")
print(results_rmse_z)


### yNUGUILLA

ARCHS <- list.files("./trago_para_analsiis/",full.names = T)

spectro <- ir_import_csv(ARCHS)
yunq <- fun_to_spectra(spectro)
yunqw <- pre_procesado(yunq,out = spectro,derivada = 2)
yunqw <- yunqw[,colnames(completo1)]
pred <- predict(fit,newdata = yunqw)

pred$predict[,2,2]* sigma[2] + mu[2]
pred$predict[,2,1]* sigma[2] + mu[2]



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

contaminado_mx <- as.matrix(fun_to_spectra(out1))

### correcion de linea base ###3
Correcionbbase2 <- baseline(contaminado_mx, method = "als")
Correcionbbase2 <- Correcionbbase2@corrected
Correcionbbase2[Correcionbbase2 < 0] <- 0
contaminado_ir <- fun_to_ir(Correcionbbase2, out1)
plot(contaminado_ir) + geom_path(aes(color = sample_id)) + theme(legend.position = "none")+labs(x = "Longitud de onda",y = "Absorbancia", title ="Correcion de Linea Base (ALS)")

#### suavizado de espectro ###3

suavizado_ir <- contaminado_ir %>% ir_smooth(method = "sg",
                                             p = 3,
                                             n = 91,
                                             m = 0)
suavizado_ir.matriz <- fun_to_spectra(suavizado_ir)

plot(suavizado_ir) + geom_path(aes(color = sample_id)) + theme(legend.position = "none")+labs(x = "Longitud de onda",y = "Absorbancia", title ="Suavizado ")


###==normalizacion ==#
normalized_spectra <- t(apply(suavizado_ir.matriz, 1, function(row)
  row / sum(row)))
datos_bc.filtered.normalized <- fun_to_ir(s = normalized_spectra, X = contaminado_ir)


plot(datos_bc.filtered.normalized) + geom_path(aes(color = sample_id)) + theme(legend.position = "none")+labs(x = "Longitud de onda",y = "Absorbancia", title ="Normalizado ")

##====filtro S-G===###
datos_bc.spectra.filtered <-  as.data.frame(t(apply(normalized_spectra, 1, function(x)
  sgolayfilt(x, p = 2, n = 7, m = 2))))
colnames(datos_bc.spectra.filtered) <- colnames(suavizado_ir.matriz)

datos_ir <- fun_to_ir(datos_bc.spectra.filtered, datos_bc.filtered.normalized)


plot(datos_ir) + geom_path(aes(color = sample_id)) + theme(legend.position = "none")+labs(x = "Longitud de onda",y = "Absorbancia", title ="Segunda Derivada ")

aux <- plot(datos_ir)
ggsave(plot = aux,filename= "./prueba.jpeg")

saveRDS(datos_ir,"./datos_analisis.rds")
##=== filtro S-G primera derivada===###

datos_bc.spectra.filtered_1derivada <-  as.data.frame(t(apply(normalized_spectra, 1, function(x)
  sgolayfilt(x, p = 2, n = 7, m = 1))))
colnames(datos_bc.spectra.filtered_1derivada) <- colnames(suavizado_ir.matriz)

datos_ir_1derivada <- fun_to_ir(datos_bc.spectra.filtered_1derivada, datos_bc.filtered.normalized)

plot(datos_ir_1derivada) + geom_path(aes(color = sample_id)) + theme(legend.position = "none")+labs(x = "Longitud de onda",y = "Absorbancia", title ="Primera  Derivada ")

aux <- plot(datos_ir)
ggsave(plot = aux,filename= "./prueba.jpeg")

saveRDS(datos_ir_1derivada,"./datos_analisis_1derivada.rds")




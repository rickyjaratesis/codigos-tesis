library(ir)
library(dplyr)



archivos <- list.files("./spects/",pattern = ".CSV$", full.names = T)



milista <-  lapply(archivos, function(x)
  read.csv(x, sep = ";", header = F))
milista <- lapply(milista, function(df) {
  df[,2] <- as.numeric(df[,2])
  df
})
archivos_names <- list.files("./spects/")

if(!dir.exists("./spects_para_analisis")){
  dir.create("./spects_para_analisis")
}



for (l in seq_along(milista)) {
  write.table(
    milista[[l]],
    file = file.path("./spects_para_analisis/", paste0(archivos_names[l], ".csv")),
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}



archivos <- list.files("./spects_para_analisis/",full.names = T)

print(archivos)          # Muestra las rutas
file.exists(archivos)    # Verifica si existen los archivos




aux <-ir::ir_import_csv(archivos)
aux$sample_id  <- sub("\\.[cC][sS][vV]$", "", aux$sample_id)

metadata <- read.csv("./datos de concentraciones y diluciones.csv",sep = ";")

colnames(metadata)[3] <- "sample_id"

metadata$sample_id <- as.character(metadata$sample_id)
aux$sample_id <- as.character(aux$sample_id)

out <- dplyr::inner_join(metadata, aux, by = "sample_id")


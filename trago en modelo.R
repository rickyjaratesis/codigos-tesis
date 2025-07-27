library(ir)
library(dplyr)
library(baseline)


archivos <- list.files("./tragoanalisis",pattern = ".CSV",full.names = T)

milista <-  lapply(archivos, function(x)
  
  read.csv(x, sep = ";", header = F))

milista <- lapply(milista, function(df) {
  df[,2] <- as.numeric(df[,2])
  df
})


if(!dir.exists("./trago_para_analsiis")){
  dir.create("./trago_para_analsiis")
}


archivos_names<-c("yung1","yung2","yung3","yung4")
  
for (l in seq_along(milista)) {
  write.table(
    milista[[l]],
    file = file.path("./trago_para_analsiis/", paste0(archivos_names[l], ".csv")),
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}



archivos <- list.files("./trago_para_analsiis/",full.names = T)

print(archivos)          # Muestra las rutas
file.exists(archivos)    # Verifica si existen los archivos




aux <-ir::ir_import_csv(archivos)
aux$sample_id  <- sub("\\.[cC][sS][vV]$", "", aux$sample_id)

#######funciones ###
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

#####baseline ####
contamido_trago <- ir_as_ir(aux)


contamido_trago<- as.matrix(fun_to_spectra(contamido_trago))

Correcionbbase2 <- baseline(contamido_trago, method = "als")
Correcionbbase2 <- Correcionbbase2@corrected
Correcionbbase2[Correcionbbase2 < 0] <- 0
contaminado_ir <- fun_to_ir(Correcionbbase2, aux)
plot(contaminado_ir) + geom_path(aes(color = sample_id)) + theme(legend.position = "none")


##suavizado===##

suavizado_ir <- contaminado_ir %>% ir_smooth(method = "sg",
                                             p = 3,
                                             n = 91,
                                             m = 0)
suavizado_ir.matriz <- fun_to_spectra(suavizado_ir)

normalized_spectra <- t(apply(suavizado_ir.matriz, 1, function(row)
  row / sum(row)))
datos_bc.filtered.normalized <- fun_to_ir(s = normalized_spectra, X = contaminado_ir)


plot(datos_bc.filtered.normalized) + geom_path(aes(color = sample_id)) + theme(legend.position = "none")+labs(x = "Longitud de onda",y = "Absorbancia", title ="Suavizado del Alcohol Artesanal")




###==normalizacion ==#
normalized_spectra <- t(apply(suavizado_ir.matriz, 1, function(row) row / sum(row)))
datos_bc.filtered.normalized <- fun_to_ir(s = normalized_spectra, X = contaminado_ir)


plot(datos_bc.filtered.normalized) + geom_path(aes(color = sample_id))+labs(x = "Longitud de onda",y = "Absorbancia", title ="Normalizacion del Alcohol Artesanal")

##====filtro S-G===###

library(signal)

datos_bc.spectra.filtered <-  as.data.frame(t(apply(normalized_spectra, 1, function(x)
  sgolayfilt(x, p = 2, n = 7, m = 2))))
colnames(datos_bc.spectra.filtered) <- colnames(suavizado_ir.matriz)
datos_ir <- fun_to_ir(datos_bc.spectra.filtered, datos_bc.filtered.normalized)

plot(datos_ir)
aux <- plot(datos_ir)+ geom_path(aes(color = sample_id))+labs(x = "Longitud de onda",y = "Absorbancia", title ="Segunda Derivada")


datos_ir <- fun_to_ir(datos_bc.spectra.filtered, datos_bc.filtered.normalized)

plot(datos_ir)
aux <- plot(datos_ir)
ggsave(plot = aux,filename= "./prueba2yung.jpeg")

saveRDS(datos_ir,"./datos_analisis2yung.rds")

##### predecir  el modelo ####

espcetroyungulla <- fun_to_spectra(datos_ir)

espcetroyungulla[espcetroyungulla < 0] <- 0


#### predcion espectro completo ###


ped<-predict(ls.liver,newdata=(espcetroyungulla))

ped$predict
(et <-ped$predict[,1,1])
(et <-ped$predict[,1,2])
(et <-ped$predict[,1,3])


(et <-ped$predict[,1,1] + mean_y[1])
(et <-ped$predict[,1,2] + mean_y[1])
(et <-ped$predict[,1,3] + mean_y[1])


(metanol <- ped$predict[,2,1])
(metanol <- ped$predict[,2,2])
(metanol <- ped$predict[,2,3])



(metanol <- ped$predict[,2,1]+ mean_y[2])
(metanol <- ped$predict[,2,2]+ mean_y[2])
(metanol <- ped$predict[,2,3]+ mean_y[2])


ped$predict

##### solo los picos####
ped<-predict(ls.liver_corte,newdata=(espcetroyungulla))
ped$predict

(et <-ped$predict[,1,1])
(et <-ped$predict[,1,2])
(et <-ped$predict[,1,3])



(et <-ped$predict[,1,1] + mean_yC[1])
(et <-ped$predict[,1,2] + mean_yC[1])
(et <-ped$predict[,1,3] + mean_yC[1])

(metanol <- ped$predict[,2,1])
(metanol <- ped$predict[,2,2])
(metanol <- ped$predict[,2,3])



(metanol <- ped$predict[,2,1]+ mean_yC[2])
(metanol <- ped$predict[,2,2]+ mean_yC[2])
(metanol <- ped$predict[,2,3]+ mean_yC[2])



##### evaluacion de los rmes entre los moodelos para evaluar el mejor #####
###====    espectro metanol  completo      ==###

et_sepectro_completo<- data.frame(etanol_todo_sinmedias=c(comp1=et1_C,
                                                  comp2=et2_C,
                                                   comp3=et3_C))

metanaol_sepectro_completo<- data.frame(metanol_todo_sin_medias=c(comp1=mt1_C,
                                                      comp2=mt2_C,
                                                      comp3=mt3_C ))

###====    espectro completo+media     ==###

et_spect_comp_media<- data.frame(etanol_todo_con_medias=c(comp1=et1_CM,
                                                    comp2=et2_CM,
                                                    comp3=et3_CM))

met_spect_comp_media<- data.frame(metanol_todo_media=c(comp1=mt1_CM,
                                                       comp2=mt2_CM,
                                                       comp3=mt3_CM ))


######====    espectro (picos)    ==###

etnaol_sepectro_picos<- data.frame(etanol_picos_sin_medias_=c(comp1=et1_P,
                                                  comp2=et2_P,
                                                  comp3=et3_P))

metanaol_sepectro_picos<- data.frame(metanol_picos_sin_medias=c(comp1=mt1_P,
                                                    comp2=mt2_P,
                                                    comp3=mt3_P ))

######====    espectro (picos) + media     ==###

etnaol_sepectro_picos_media<- data.frame(etanol_picos_con_media=c(comp1=et1_PM,
                                                 comp2=et2_PM,
                                                 comp3=et3_PM))

metanaol_sepectro_picos_media<- data.frame(metanol_picos_con_media=c(comp1=mt1_PM,
                                                    comp2=mt2_PM,
                                                    comp3=mt3_PM ))



resultadosmetanol <-bind_cols(metanaol_sepectro_completo, met_spect_comp_media, metanaol_sepectro_picos, metanaol_sepectro_picos_media, )
colnames(resultadosmetanol)<-c("Completo","Completo_media","pico", "pico_media")

resultadosEtanol <-bind_cols(etnaol_sepectro_completo, met_spect_comp_media,etnaol_sepectro_picos,etnaol_sepectro_picos_media )
colnames(resultadosEtanol)<-c("Completo","Completo_media","pico", "pico_media")

resultadosmetanol$Alcohol<-"METANOL"
resultadosEtanol$Alcohol<-"ETANOL"

df_final<-rbind(resultadosmetanol,resultadosEtanol)


write.csv(df_final, "data_final.CSV")

library(tidyr)
library(ggplot2)

# Convertir a formato largo (tidy)
df_long <- df_final %>%
  pivot_longer(cols = c(Completo, Completo_media, pico, pico_media),
               names_to = "Metodo",
               values_to = "Valor")

# Ver estructura
head(df_long)

ggplot(df_long, aes(x = Metodo, y = Valor, fill = Alcohol)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.7) +
  scale_fill_manual(values = c("METANOL" = "#1f77b4", "ETANOL" = "#ff7f0e")) +
  labs(title = "Comparación RMSE de Métodos entre Metanol y Etanol",
       x = "Método de Análisis",
       y = "Valor de Medición") +
  theme_bw()



todos_Dfs <- reshape2::melt(df_final)

head(todos_Dfs)

p <- ggplot(data=todos_Dfs, aes(x= Metodo, y=rmse, fill=alcohol)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()
# Use custom colors
p + scale_fill_manual(values=c('#999999','#E69F00'))
# Use brewer color palettes
p + scale_fill_brewer(palette="Blues")
    

#### === GRAFICO  ====####
library(ggplot2)
library(dplyr)

# Calcular mínimos por Alcohol y Método
minimos <- df_final %>%
  group_by(Alcohol) %>%
  summarise(
    min_completo = min(Completo, na.rm = TRUE),
    min_pico = min(pico, na.rm = TRUE)
  ) %>% 
  reshape2::melt(id.vars = "Alcohol", variable.name = "Metodo", value.name = "minimo") %>%
  mutate(Metodo = gsub("min_", "", Metodo))

# Ordenar componentes por valor (opcional)
df_final$df_final <- factor(df_final$df_final, levels = df_final$df_final[order(-df_final$Completo)])

# Gráfico
ggplot(df_final, aes(x = X)) +
  geom_col(aes(y = Completo, fill = Alcohol), position = position_dodge(0.8), width = 0.7) +
  geom_col(aes(y = pico, fill = Alcohol), position = position_dodge(0.8), width = 0.7, alpha = 0.7) +
  
  # Líneas horizontales por Alcohol y Método
  geom_hline(data = minimos, 
             aes(yintercept = minimo, color = Alcohol, linetype = Metodo),
             linewidth = 0.8, show.legend = TRUE) +
  
  # Esquema de colores
  scale_fill_manual(values = c('ETANOL' = '#E69F00', 'METANOL' = '#56B4E9')) +
  scale_color_manual(values = c('ETANOL' = '#D55E00', 'METANOL' = '#0072B2')) +
  scale_linetype_manual(values = c('completo' = "solid", 'pico' = "dashed")) +
  
  # Temas y etiquetas
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top") +
  labs(x = "Componente", y = "Valor",
       title = "Comparación por tipo de Alcohol",
       subtitle = "Líneas sólidas: mínimo en método Completo\nLíneas discontinuas: mínimo en método Pico") +
  
  # Opcional: etiquetar mínimos en eje Y
  geom_text(data = minimos, 
            aes(x = 0.5, y = minimo + 0.5, 
                label = paste("Mín", Alcohol, "=", round(minimo, 2)),
                color = Alcohol),
            hjust = 0, vjust = 0, size = 3, show.legend = FALSE)



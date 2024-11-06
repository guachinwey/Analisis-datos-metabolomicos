# Descarga de datos

# Indicamos los enlaces de descarga del dataset
link_DataValues <- "https://raw.githubusercontent.com/nutrimetabolomics/metaboData/main/Datasets/2018-MetabotypingPaper/DataValues_S013.csv"
link_DataInfo <- "https://raw.githubusercontent.com/nutrimetabolomics/metaboData/main/Datasets/2018-MetabotypingPaper/DataInfo_S013.csv"

# Ponemos nombre a los archivos
DataValues <- "DataValues_S013.csv"
DataInfo <- "DataInfo_S013.csv"

# Procedemos a la descarga
download.file(link_DataValues, destfile = DataValues)
download.file(link_DataInfo, destfile = DataInfo)

# Cargamos los archivos
data_values <- read.csv(DataValues, row.names = 1)
data_info <- read.csv(DataInfo, row.names = 1)

# Creación del objeto Summarized Experiment

# Transponemos data_values
data_values <- t(data_values)

# Separamos los metadatos y los datos de los ensayos
metadata <- data_values[1:9,]
metadata <- t(metadata)
metadata <- as.data.frame(metadata)
assay_data <- data_values[-c(1:9), ]

# Asignamos los sujetos como identificativo de fila en metadata y de las columnas de assay_data
rownames(metadata) <- metadata$SUBJECTS
colnames(assay_data) <- metadata$SUBJECTS

# Convertimos a numéricas las columnas de metadatos correspondientes
metadata$SUBJECTS <- as.numeric(metadata$SUBJECTS)
metadata$AGE <- as.numeric(metadata$AGE)
metadata$MEDDM_T0 <- as.numeric(metadata$MEDDM_T0)
metadata$MEDCOL_T0 <- as.numeric(metadata$MEDCOL_T0)
metadata$MEDINF_T0 <- as.numeric(metadata$MEDINF_T0)
metadata$MEDHTA_T0 <- as.numeric(metadata$MEDHTA_T0)

# Forzamos los metadatos a ser reconocidos como tal
library(S4Vectors)
colData <- DataFrame(Subjects = metadata$SUBJECTS,
                     Surgery = metadata$SURGERY,
                     Age = metadata$AGE, 
                     Gender = metadata$GENDER,
                     Group = metadata$Group,
                     MEDMM_T0 = metadata$MEDDM_T0, 
                     MEDCOL_T0 = metadata$MEDCOL_T0, 
                     MEDINF_T0 = metadata$MEDINF_T0, 
                     MEDHTA_T0 = metadata$MEDHTA_T0)

# Convertimos assay_data a numérico y manejamos los NA
assay_data[is.na(assay_data)] <- 0
assay_data <- as.matrix(assay_data)
assay_data <- matrix(as.numeric(assay_data), nrow = 686, ncol = 39)

# Asignamos los identificativos de las filas y columnas
rownames(assay_data) <- rownames(data_values)[-c(1:9)]
colnames(assay_data) <- metadata$SUBJECTS

# Comprobamos que las dimensiones de assay_data y colData son compatibles
print(dim(assay_data))
print(dim(colData))
print(all(rownames(colData) == colnames(assay_data)))

# Cargamos la librería SummarizedExperiment
library(SummarizedExperiment)

# Creamos el objeto SummarizedExperiment
summarized_experiment <- SummarizedExperiment(assays = list(counts = assay_data), colData = colData)

# Guardamos el objeto SummarizedExperiment
save(summarized_experiment, file = "summarized_experiment.RDA")

# Análisis del dataset

# Revisión del objeto SummarizedExperiment
summarized_experiment

# Revisión de los metadatos en colData
colData(summarized_experiment)

# Revisión de las filas de assay
head(rownames(assay(summarized_experiment)))

# Resumen inicial del dataset
summary(assay(summarized_experiment))

# Boxplot del assay
boxplot(assay(summarized_experiment))

# Eliminamos valores negativos
assay(summarized_experiment)[assay(summarized_experiment) < 0] <- NA

# Normalización de los datos
library(DESeq2)
sum(is.na(assay(summarized_experiment)))
assay(summarized_experiment)[is.na(assay(summarized_experiment))] <- 0
assay(summarized_experiment) <- round(assay(summarized_experiment))
DESDataset <- DESeqDataSet(summarized_experiment, design = ~ Group)

# Aplicamos la normalización
DES_vst <- rlog(DESDataset)
assay(summarized_experiment) <- assay(DES_vst)

# Boxplot tras normalización
boxplot(assay(summarized_experiment))

# Estudio de componentes principales
data_for_pcs <- assay(summarized_experiment)
pcs <- prcomp(data_for_pcs)
barplot(pcs$sdev)

# Representación de los dos primeros componentes principales
plot(pcs$rotation[, 1], pcs$rotation[, 2], main = "Representación de los dos primeros componentes principales")
text(pcs$rotation[, 1], pcs$rotation[, 2], metadata$Group, cex = 0.8, pps = 3)

# Diagrama de dispersión de los componentes principales
plot(pcs$x[,1], pcs$x[,2])

# Estudio de la distribución de las características
hist(colData(summarized_experiment)$Age, main = "Distribución de la edad", xlab = "Edad", ylab = "Frecuencia de muestras")

# Distribución por grupo
library(ggplot2)
ggplot(as.data.frame(colData(summarized_experiment)), aes(x = Group)) +
  geom_bar() + 
  xlab("Grupos") +
  ylab("Número de muestras") +
  ggtitle("Número de muestras por grupo")

# Distribución por tipo de cirugía
ggplot(as.data.frame(colData(summarized_experiment)), aes(x = Surgery)) +
  geom_bar() + 
  xlab("Tipo de cirugía") +
  ylab("Número de muestras") +
  ggtitle("Número de muestras por tipo de cirugía")

# Generación de un heatmap
library(pheatmap)
pheatmap(assay(summarized_experiment), cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = FALSE, show_colnames = FALSE, main = "Heatmap de todas las muestras")

# Heatmap de la correlación entre los metabolitos
cor_matrix <- cor(assay(summarized_experiment))
pheatmap(cor_matrix, main = "Correlación entre metabolitos", display_numbers = TRUE)

# Filtrado de metabolitos
valor_na <- 0.5 * ncol(summarized_experiment)
na_correctos <- rowSums(is.na(assay(summarized_experiment))) <= valor_na
varianza <- apply(assay(summarized_experiment), 1, var, na.rm = TRUE)
var_correctos <- varianza > 0.1
metabolitos_filtrados <- na_correctos & var_correctos

# Filtramos los metabolitos
summarized_experiment_filtered <- summarized_experiment[metabolitos_filtrados, ]

# Revisamos el número de metabolitos filtrados
dim(summarized_experiment_filtered)

# Heatmap con los metabolitos filtrados
pheatmap(assay(summarized_experiment_filtered), cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = FALSE, show_colnames = FALSE, main = "Heatmap con metabolitos filtrados")

# Heatmap de la correlación entre metabolitos filtrados
cor_matrix_filtered <- cor(assay(summarized_experiment_filtered))
pheatmap(cor_matrix_filtered, main = "Correlación entre metabolitos filtrados", display_numbers = TRUE)
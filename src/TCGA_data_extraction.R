if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instalar TCGAbiolinks
BiocManager::install("TCGAbiolinks")

# -----------------------------------------------------------------------------
# Script para descargar datos de TCGA-LUAD (Cáncer de Pulmón - Adenocarcinoma)
# -----------------------------------------------------------------------------

# 1. Cargar la librería
library(TCGAbiolinks)
library(SummarizedExperiment) # Para manejar el objeto de datos

# 2. Definir el proyecto y la ruta de salida


project_name <- "TCGA-LUAD"
output_dir <- "TCGA_data" 

# 3. Construir la consulta para los datos de expresión génica
# Queremos la matriz de conteos generada por el pipeline HTSeq
query_expression <- GDCquery(
  project = project_name,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# 4. Descargar los datos
GDCdownload(
  query = query_expression,
  method = "api",
  directory = output_dir,
  files.per.chunk = 10 
)

# 5. Preparar los datos en un objeto R
# Esto lee todos los ficheros descargados y los une en una única matriz
luad_data <- GDCprepare(
  query = query_expression,
  directory = output_dir
)

# 6. Extraer la matriz de conteos y los datos clínicos
# El objeto 'luad_data' es un SummarizedExperiment que contiene todo
counts_matrix <- assay(luad_data)
clinical_data <- colData(luad_data)

# 7. Limpieza de la matriz de conteos
rownames(counts_matrix) <- gsub("\\..*$", "", rownames(counts_matrix))

# 8. Inspección y Limpieza de los Datos Clínicos
# Convertimos a un data.frame estándar para facilitar la manipulación
clinical_df <- as.data.frame(clinical_data)

# Identificamos las columnas que son listas
is_list_col <- sapply(clinical_df, is.list)
list_cols <- names(is_list_col[is_list_col])

print("Columnas con tipo de dato 'list' que serán excluidas:")
print(list_cols)

# Creamos un nuevo data.frame solo con las columnas que NO son listas
clinical_df_clean <- clinical_df[, !is_list_col]


key_clinical_vars <- c(
  # Identificadores
  "barcode", "patient", "definition",
  # Variables del tumor
  "ajcc_pathologic_stage", "primary_diagnosis", "tissue_or_organ_of_origin",
  "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m",
  # Datos demográficos
  "age_at_index", "gender", "race",
  # Datos de supervivencia
  "vital_status", "days_to_death", "days_to_last_follow_up",
  # Opcional: Exposición
  "tobacco_smoking_status", "pack_years_smoked"
)
# Filtramos para quedarnos solo con las columnas clave que existen en nuestro dataframe
# Usamos `any_of` para evitar errores si alguna columna no existe
library(dplyr)
clinical_df_final <- select(clinical_df, all_of(key_clinical_vars))

print("Dimensiones de la tabla clínica final:")
print(dim(clinical_df_final))
print("Primeras filas de la tabla clínica final:")
print(head(clinical_df_final))


# 9. Guardar los datos limpios en formato CSV
write.csv(t(counts_matrix), file = "TCGA-LUAD_star_counts.csv")
write.csv(clinical_df_final, file = "TCGA-LUAD_clinical_data_clean.csv", row.names = FALSE)


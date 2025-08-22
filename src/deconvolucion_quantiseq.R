

if (!require("remotes", quietly = TRUE)) install.packages("remotes")
if (!require("immunedeconv", quietly = TRUE)) remotes::install_github("grst/immunedeconv")
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!require("tibble", quietly = TRUE)) install.packages("tibble")
if (!require("tidyr", quietly = TRUE)) install.packages("tidyr")
# Librerías para anotación y normalización
if (!require("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
if (!require("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi")
if (!require("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")


library(immunedeconv)
library(dplyr)
library(tibble)
library(tidyr)
library(biomaRt)

#Carga de datos

print("--- Cargando datos de conteos crudos ---")
bulk_path <- "data/processed/TCGA-LUAD_bulk_for_R.tsv"
bulk_df <- read.delim(bulk_path, sep = "\t") %>% tibble::column_to_rownames("sample_id")
bulk_counts_matrix <- t(as.matrix(bulk_df)) # Genes x Muestras



print("--- Iniciando traducción de Ensembl IDs a Símbolos de Gen ---")

# 1. TRADUCCIÓN DE IDS
# Usamos el paquete `org.Hs.eg.db` para el mapeo
# Primero, limpiamos la versión del Ensembl ID (si la tuviera)
ensembl_ids <- gsub("\\..*$", "", rownames(bulk_counts_matrix))

# Hacemos el mapeo
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first") 

# Creamos un dataframe con el mapeo y eliminamos los IDs sin símbolo
mapping_df <- data.frame(ENSEMBL = ensembl_ids, SYMBOL = gene_symbols, stringsAsFactors = FALSE)
mapping_df <- na.omit(mapping_df)

# Filtramos la matriz de conteos original para quedarnos solo con los genes que tienen símbolo
bulk_counts_matrix_mapped <- bulk_counts_matrix[mapping_df$ENSEMBL, ]
# Y le asignamos los nuevos nombres de fila (símbolos de gen)
rownames(bulk_counts_matrix_mapped) <- mapping_df$SYMBOL

# Manejamos los símbolos de gen duplicados (nos quedamos con el que tiene mayor expresión promedio)
gene_means <- rowMeans(bulk_counts_matrix_mapped)
bulk_counts_matrix_final <- bulk_counts_matrix_mapped[
  !duplicated(rownames(bulk_counts_matrix_mapped)) | !rev(duplicated(rev(rownames(bulk_counts_matrix_mapped)))),
]
print("Traducción de IDs completada.")


print("--- Iniciando la deconvolución con quanTIseq estándar ---")
deconv_results <- deconvolute(
  gene_expression = bulk_counts_matrix_final,
  method = "quantiseq",
  tumor = TRUE,
  arrays = FALSE,
  scale_mrna = TRUE
)

deconv_results_wide <- deconv_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample")

# Guardamos el dataframe de resultados en un fichero CSV para poder cargarlo
# de vuelta en nuestro notebook de Python para el análisis final.

print("--- Guardando los resultados ---")

# Ruta de salida
results_path <- "data/processed/deconv_results_quantisec_new.csv"

# Usamos `write.csv` para guardar los resultados
write.csv(deconv_results, file = results_path, row.names = FALSE)

print(paste("Resultados de la deconvolución guardados en:", results_path))

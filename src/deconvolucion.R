# --- APARTADO 1: INSTALACIÓN Y CARGA DE LIBRERÍAS ---
# Instalamos y cargamos todas las librerías necesarias.
# `remotes` se usa para instalar paquetes desde GitHub.
# `immunedeconv` es nuestra herramienta principal.
# `dplyr` y `tibble` son para la manipulación de datos.

if (!require("remotes", quietly = TRUE)) install.packages("remotes")
if (!require("immunedeconv", quietly = TRUE)) remotes::install_github("grst/immunedeconv")
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!require("tibble", quietly = TRUE)) install.packages("tibble")

library(immunedeconv)
library(dplyr)
library(tibble)

# --- APARTADO 2: CARGA Y FORMATEO DE DATOS ---
#
# Cargamos los ficheros .tsv que preparamos en Python y los formateamos
# para que sean compatibles con immunedeconv.

print("--- Cargando y formateando los datos de entrada ---")

# Rutas a los ficheros
bulk_path <- "data/processed/TCGA-LUAD_bulk_for_R.tsv"
signature_path <- "data/processed/optimized_signature_for_R.tsv"

# Cargar los datos
bulk_df <- read.delim(bulk_path, sep = "\t", stringsAsFactors = FALSE)
signature_df <- read.delim(signature_path, sep = "\t", stringsAsFactors = FALSE)

# --- Formateo de la Matriz de Bulk ---
# immunedeconv espera:
# - Genes como filas (con nombres de gen como nombres de fila).
# - Muestras como columnas.
#
# Primero, convertimos la primera columna ('sample_id') en los nombres de las filas
bulk_df <- bulk_df %>% tibble::column_to_rownames("sample_id")

# Luego, transponemos el dataframe
bulk_matrix <- t(as.matrix(bulk_df))

print(paste("Matriz de bulk formateada con", nrow(bulk_matrix), "genes y", ncol(bulk_matrix), "muestras."))


# --- Formateo de la Matriz de Firmas ---
# immunedeconv la espera como una matriz simple con genes como nombres de fila.
signature_matrix <- signature_df %>% tibble::column_to_rownames("gene") %>% as.matrix()

print(paste("Matriz de firmas formateada con", nrow(signature_matrix), "genes y", ncol(signature_matrix), "tipos celulares."))

# --- APARTADO 3: EJECUCIÓN DE LA DECONVOLUCIÓN ---
#
# Ejecutamos la deconvolución usando el método EPIC con nuestra firma personalizada.
# EPIC está diseñado para datos de RNA-seq de tumores y es una elección robusta.

print("--- Iniciando la deconvolución con el método EPIC ---")

# 1. Definir los argumentos para la función
# `gene_expression_matrix` es nuestra matriz de bulk
# `signature_matrix` es nuestra matriz de firmas
# `signature_genes` es simplemente la lista de todos los genes en nuestra firma
signature_gene_list <- rownames(signature_matrix)

# 2. Llamar a la función epic_custom
deconv_results <- deconvolute_epic_custom(
  gene_expression_matrix = bulk_matrix,
  signature_matrix = signature_matrix,
  signature_genes = signature_gene_list
)

# El resultado de esta función es una matriz (muestras x tipos celulares)
deconv_results_wide <- deconv_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample")

print("--- Deconvolución completada ---")

# --- APARTADO 4: GUARDADO DE RESULTADOS ---
#
# Guardamos el dataframe de resultados en un fichero CSV para poder cargarlo
# de vuelta en nuestro notebook de Python para el análisis final.

print("--- Guardando los resultados ---")

# Ruta de salida
results_path <- "data/processed/deconv_results_epic.csv"

# Usamos `write.csv` para guardar los resultados
write.csv(deconv_results_wide, file = results_path, row.names = FALSE)

print(paste("Resultados de la deconvolución guardados en:", results_path))

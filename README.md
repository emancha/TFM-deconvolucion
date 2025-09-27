# Deconvolución del Microambiente Tumoral y Análisis Pronóstico en Cáncer de Pulmón mediante scRNA-seq y Machine Learning

## Resumen del Proyecto

Este repositorio contiene el código y los análisis realizados para el Trabajo de Fin de Máster `Deconvolución del Microambiente Tumoral y Análisis Pronóstico en Cáncer de Pulmón mediante scRNA-seq y Machine Learning`. El objetivo principal es desarrollar y validar un pipeline de bioinformática para estimar la composición de tipos celulares en muestras de adenocarcinoma de pulmón (LUAD) a partir de datos de bulk RNA-seq, utilizando un atlas de single-cell RNA-seq como referencia.

El proyecto abarca desde el pre-procesamiento de datos de scRNA-seq, la construcción y comparación de matrices de firmas genéticas mediante técnicas de Machine Learning, hasta la aplicación de algoritmos de deconvolución y el análisis clínico de los resultados en la cohorte de TCGA, incluyendo análisis de supervivencia.

---

## Estructura del Repositorio

El proyecto está organizado en la siguiente estructura de carpetas:

-   **/data/**: Contiene los datos. **Nota:** Debido a su tamaño, los ficheros de datos principales no están incluidos en el repositorio. Se deben descargar de las fuentes originales (ver "Instalación").
    -   `/raw/`: Datos originales descargados.
    -   `/processed/`: Datos intermedios y finales generados por los notebooks.
-   **/notebooks/**: Jupyter Notebooks que contienen el flujo de trabajo completo del análisis, organizados de forma secuencial.
    -   `1_Exploracion_y_Filtrado_de_Datos.ipynb`: Carga, pre-procesa y prepara el dataset de referencia de scRNA-seq.
    -   `2_Construccion_y_Validacion_de_Firmas.ipynb`: Construye, compara y valida las matrices de firmas genéticas candidatas.
    -   `3_Exploracion_TCGA.ipynb`: Carga, limpia y caracteriza la cohorte de TCGA-LUAD.
    -   `4_Deconvolucion_y_Analisis_Clinico.ipynb`: Realiza la deconvolución comparativa y los análisis clínicos y de supervivencia finales.
-   **/src/**: Módulos de Python y R con funciones de utilidad reutilizables para el pre-procesamiento, normalización, modelado y visualización.
    -   `ensembl_api.py`: Funciones para interactuar con la API de Ensembl y obtener nombres de genes.
    -   `preprocessing.py`: Funciones para preprocesar los datos.
    -   `counts_to_tpm.py`: Script de Python 
    -   `models.py`: Definición de los modelos de Machine Learning.
    -   `plotting.py`: Funciones para generar las matrices de confusión
    -   `deconvolucion.R`: Script de R utilizado para ejecutar los algoritmos de deconvolución EPIC.
    -   `deconvolucion_quantiseq.R`: Script de R utilizado para ejecutar los algoritmos de deconvolución quanTIseq.
    -   `TCGA_data_extraction.R`: Script de R utilizado para extraer los datos de TCGA.
    -   `plot_kaplan_meier.py`: Script de Python utilizado para generar las curvas de supervivencia de Kaplan-Meier.

-   **/outputs/**: Ficheros generados por los notebooks de análisis.
    -   `/figures/`: Figuras y gráficos (matrices de confusión, curvas de supervivencia, etc.).
    -   `/models/`: Modelos de Machine Learning entrenados y guardados.
-   `environment.yml`: Fichero de entorno de Conda para garantizar la reproducibilidad completa del entorno de Python y R.

---

## Instalación y Reproducibilidad

Para replicar este análisis, sigue los siguientes pasos.

### 1. Requisitos Previos

-   [Conda](https://docs.conda.io/en/latest/miniconda.html) para la gestión de entornos.
-   [Git](https://git-scm.com/) para clonar el repositorio.

### 2. Clonar el Repositorio

```bash
git clone https://github.com/emancha/TFM-deconvolucion
cd TFM_repo
```

### 3. Crear el Entorno de Conda

El fichero `environment.yml` contiene todas las dependencias de Python y R necesarias.

```bash
conda env create -f environment.yml
conda activate entorno_TFM
```

### 4. Descargar los Datos

Los datasets principales no están incluidos en este repositorio debido a su tamaño.

-   **Dataset de scRNA-seq:** Descargar el fichero `.h5ad` de la colección de CZI CELLxGENE: [Enlace a la colección](https://cellxgene.cziscience.com/collections/3f7c572c-cd73-4b51-a313-207c7f20f188) y colocarlo en la carpeta `data/raw/`.
-   **Dataset de TCGA:** Los datos de TCGA se descargan automáticamente mediante el script `deconvolution.R` (ver paso 5).

### 5. Ejecutar el Pipeline

El análisis completo se puede reproducir ejecutando los notebooks en su orden numérico:

1.  **Ejecutar los notebooks de Python:** Abre Jupyter Lab/Notebook en el entorno `entorno_TFM` y ejecuta los notebooks de la carpeta `/notebooks/` en orden, desde el 1 al 4.
    -   *Nota:* El Notebook 4 tiene un punto de parada donde se debe ejecutar el script de R.
2.  **Ejecutar el script de R:** En el punto indicado del Notebook 4, abre `deconvolution.R` con RStudio y ejecútalo. Esto generará los ficheros de resultados de la deconvolución en `data/processed/`.
3.  **Continuar con el notebook de Python:** Una vez finalizado el script de R, puedes continuar con la ejecución del resto de celdas del Notebook 4 para generar los análisis clínicos finales.

---
## Palabras Clave

Adenocarcinoma de Pulmón; Microambiente Tumoral; Transcriptoma; Aprendizaje Automático; Deconvolución Celular; Análisis de Supervivencia; Single-Cell Analysis; TCGA.

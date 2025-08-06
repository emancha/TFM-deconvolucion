# src/preprocessing2.py

import scanpy as sc
import numpy as np

# ==============================================================================
# Función de Diagnóstico / Validación
# ==============================================================================

def validate_adata_state(adata):
    """
    Diagnostica el estado de pre-procesamiento de un objeto AnnData en memoria.
    Retorna un diccionario con el informe de estado.
    """
    print("--- Diagnosticando estado del objeto AnnData ---")
    report = {}

    # Comprobación 1: Métricas de QC
    report['qc_metrics_calculated'] = all(c in adata.obs for c in ['n_genes_by_counts', 'total_counts'])
    print(f"Métricas de QC calculadas: {report['qc_metrics_calculated']}")

    # Comprobación 2: Normalización y Logaritmo
    # Comprobamos si los datos son de tipo float y no enteros.
    is_raw = np.issubdtype(adata.X.dtype, np.integer)
    report['is_normalized_and_logged'] = not is_raw
    print(f"Datos normalizados y logaritmizados: {report['is_normalized_and_logged']}")

    # Comprobación 3: Genes Altamente Variables (HVG)
    report['hvg_calculated'] = 'highly_variable' in adata.var.columns
    print(f"Genes Altamente Variables calculados: {report['hvg_calculated']}")
    
    # Comprobación 4: Reducción de Dimensionalidad
    report['pca_calculated'] = 'X_pca' in adata.obsm
    report['neighbors_calculated'] = 'neighbors' in adata.uns
    report['umap_calculated'] = 'X_umap' in adata.obsm
    print(f"PCA calculado: {report['pca_calculated']}")
    print(f"Grafo de vecinos calculado: {report['neighbors_calculated']}")
    print(f"UMAP calculado: {report['umap_calculated']}")
    
    return report


# ==============================================================================
# Funciones de Pre-procesamiento Atómicas
# ==============================================================================

def calculate_and_filter_qc(adata_in, min_genes, min_cells, mt_qc_threshold):
    """
    Calcula métricas de QC y filtra células y genes de baja calidad.
    Retorna un nuevo objeto AnnData filtrado.
    """
    print("--- Ejecutando Paso 1: Control de Calidad y Filtrado ---")
    
    adata = adata_in.copy()
    
    # Identifica genes mitocondriales
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    
    # Calcula métricas de QC
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Guarda el número de células/genes antes de filtrar
    n_obs_before, n_vars_before = adata.shape
    
    # Filtrado
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[adata.obs.pct_counts_mt < mt_qc_threshold, :].copy()
    
    print(f"Filtrado QC: {n_obs_before - adata.n_obs} células y {n_vars_before - adata.n_vars} genes eliminados.")
    return adata


def normalize_and_log_transform(adata_in, target_sum):
    """
    Normaliza por tamaño de librería y aplica una transformación logarítmica.
    Retorna un nuevo objeto AnnData procesado.
    """
    print("--- Ejecutando Paso 2: Normalización y Transformación Logarítmica ---")
    
    adata = adata_in.copy()
    
    # Guardamos los conteos crudos en una capa por si se necesitan más tarde
    adata.layers['counts'] = adata.X.copy()
    
    sc.pp.normalize_total(adata, target_sum)
    sc.pp.log1p(adata)
    
    # Guardamos los datos logaritmizados antes de cualquier escalado posterior
    adata.layers['log_normalized'] = adata.X.copy()
    
    return adata


def find_highly_variable_genes(adata_in, n_top_genes, flavor):
    """
    Identifica Genes Altamente Variables (HVG). No filtra el objeto, solo añade la anotación.
    Retorna un nuevo objeto AnnData con la columna .var['highly_variable'].
    """
    print(f"--- Ejecutando Paso 3: Encontrando {n_top_genes} Genes Altamente Variables ---")
    
    adata = adata_in.copy()
    
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        subset=False,  # No se filtra nada, sino que se marca
        flavor=flavor
    )
    
    return adata


def annotate_for_visualization(adata_in, scale_max_value, pca_svd_solver, n_pcs):
    """
    Calcula PCA, vecinos y UMAP para visualización sin alterar la matriz .X principal.

    Esta función realiza los cálculos que requieren escalado en una copia temporal
    y luego transfiere los resultados (embeddings y grafos) de vuelta al objeto
    original.

    Retorna:
    -------
    AnnData
        El objeto AnnData original, anotado con .obsm['X_pca'], .obsm['X_umap'], etc.
    """
    print("--- Ejecutando Anotación para Visualización (PCA, Vecinos, UMAP) ---")
    
    # Hacemos una copia para no modificar el objeto de entrada directamente
    adata = adata_in.copy()
    
    # Verificamos si los HVG han sido marcados, lo cual es un prerrequisito
    if 'highly_variable' not in adata.var.columns:
        raise ValueError("Error: Se deben calcular los Genes Altamente Variables antes de este paso.")

    # 1. Creamos una copia temporal que filtra por HVG
    print(f"Creando copia temporal con {adata.var.highly_variable.sum()} genes para PCA...")
    adata_pca = adata[:, adata.var.highly_variable].copy()
    
    # 2. Escalamos y calculamos PCA en la copia temporal
    print("Escalando datos y ejecutando PCA...")
    sc.pp.scale(adata_pca, max_value=scale_max_value)
    sc.tl.pca(adata_pca, svd_solver=pca_svd_solver, n_comps=n_pcs)
    
    # 3. Calculamos vecinos y UMAP en la copia temporal
    print("Calculando grafo de vecinos y UMAP...")
    sc.pp.neighbors(adata_pca, n_pcs=n_pcs)
    sc.tl.umap(adata_pca)
    
    # 4. Transferimos los resultados clave de vuelta al objeto principal
    print("Transfiriendo resultados al objeto principal...")
    adata.obsm['X_pca'] = adata_pca.obsm['X_pca']
    adata.obsm['X_umap'] = adata_pca.obsm['X_umap']
    adata.uns['neighbors'] = adata_pca.uns['neighbors']
    
    # No necesitamos devolver la copia temporal, se eliminará al salir de la función
    print("Anotación para visualización completada.")
    
    return adata
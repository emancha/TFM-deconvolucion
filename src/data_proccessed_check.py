import scanpy as sc
import pandas as pd
import numpy as np

RUTA = "C:/Users/emanc/OneDrive - Universidad de Extremadura (1)/MUBinf/09.TFM/Datos/"
ARCHIVO = "lung_cancer_filtrado.h5ad"

adata= sc.read_h5ad(f"{RUTA}{ARCHIVO}", backed="r")

print("\nRecuento de tipos celulares en el dataset de Cáncer de Pulmón:")
cell_counts = adata.obs['cell_type'].value_counts()
print(cell_counts)

#Verificación de transformación o normalización

#print(type(adata.X), adata.X.shape)
#print(adata.X[:20, :20])

#RESULTADO: Datos transformados (log1p) y normalizados por célula
#CONFIRMACIÖN:
    #subset = adata[:50000, :50000].X

    #print(type(subset))  # Debería ser scipy.sparse
    #print(subset.max())  # Verifica si los valores están <~15


#RESULTADO: EL VALOR MAYOR ES DE 8.84 (50000 VALORES)
'''
#genes altamente variables ya seleccionados
print("Genes altamente variables")
print("highly_variable" in adata.var)
#RESULTADO: NO TIENE GENES ALTAMENTE VARIABLES SELECCIONADOS

#grafo de vecinos computado
print("grafo de vecinos computado:")
print("neighbors" in adata.uns)
#RESULTADO: NO TIENE GRAFOS DE VECINOS COMPUTADOS

#reducción de dimensionalidad (PCA, UMAP, t-SNE)
print("PCA:")
print("X_pca" in adata.obsm)
print("UMAP:")
print("X_umap" in adata.obsm)
print("tSNE:")
print("X_tsne" in adata.obsm)
#RESULTADO: SOLO TIENE UMAP CALCULADO
#GRAFICO
    # sc.pl.umap(adata, color='leiden' if 'leiden' in adata.obs else adata.obs.columns[0])

#CLUSTERING
print("leiden:")
print([k for k in adata.obs.columns if 'leiden' in k or 'louvain' in k])
#RESULTADO: LEIDEN [0.2-1.4]

#ANNOTATIONS
print("genes anotados:")
print([col for col in adata.obs.columns if 'cell_type' in col or 'annotation' in col])
#RESULTADOS: ['author_first_cell_type', 'author_cell_type', 'cell_type_ontology_term_id', 'cell_type']

#GENES MARCADORES IDENTIFICADOS
print("Genes marcadores identificados:")
print("rank_genes_groups" in adata.uns)
#RESULTADO: NO

#VERSIÓN CRUDA
print("Versión cruda de los datos:")
print(adata.raw is not None)
'''

print(adata.uns)

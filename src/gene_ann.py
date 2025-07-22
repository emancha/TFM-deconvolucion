import scanpy as sc
import pandas as pd
import numpy as np
from ensembl_api import fetch_gene_symbols_ensembl

RUTA = "C:/Users/emanc/OneDrive - Universidad de Extremadura (1)/MUBinf/09.TFM/Datos/"
ARCHIVO = "lung_cancer_filtrado.h5ad"

lung_cancer_adata= sc.read_h5ad(f"{RUTA}{ARCHIVO}")

ensembl_ids = lung_cancer_adata.var.index.tolist()
#Llamada a la función por lotes
gene_symbol_dict = fetch_gene_symbols_ensembl(ensembl_ids)
# Añadir como nueva columna
lung_cancer_adata.var['gene_name'] = lung_cancer_adata.var.index.map(gene_symbol_dict)
# Si el gen no tiene nombre, se usa el 'feature_name'
lung_cancer_adata.var['gene_name'] = lung_cancer_adata.var['gene_name'].fillna(lung_cancer_adata.var['feature_name'])


print(lung_cancer_adata.var[['feature_name', 'gene_name']].head(10))
print(f"Genes sin nombre (NaN): {lung_cancer_adata.var['gene_name'].isna().sum()}")
#Guardar los datos nuevos
lung_cancer_adata.write(f"{RUTA}lung_cancer_gene_name.h5ad")
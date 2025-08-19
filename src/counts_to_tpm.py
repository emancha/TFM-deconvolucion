import numpy as np
from scipy.sparse import csr_matrix, diags


def counts_to_tpm(counts_matrix, lengths_series):
    """
    Convierte una matriz de conteos (densa o dispersa) a TPM.

    Parámetros:
    ----------
    counts_matrix : np.ndarray or scipy.sparse.spmatrix
        Matriz de conteos (células/muestras x genes).
    lengths_series : pd.Series
        Serie con la longitud de cada gen, con el índice de genes alineado.

    Retorna:
    -------
    np.ndarray or scipy.sparse.spmatrix
        Matriz de datos en TPM, manteniendo el tipo de entrada si es posible.
    """
    
    # 1. Asegurarse de que la matriz de entrada es scipy.sparse para operaciones eficientes
    if not isinstance(counts_matrix, csr_matrix):
        counts_matrix = csr_matrix(counts_matrix)
        
    # 2. Normalizar por longitud de gen en kilobases (RPK)
    lengths_kb = lengths_series.values / 1000
    inv_lengths_kb = 1 / lengths_kb
    rpk_matrix = counts_matrix.dot(diags(inv_lengths_kb))

    # 3. Calcular el "per million" scaling factor
    sum_rpk_per_cell = np.asarray(rpk_matrix.sum(axis=1))
    per_million_scalers = sum_rpk_per_cell / 1e6
    per_million_scalers[per_million_scalers == 0] = 1 # Evitar división por cero
    
    # 4. Dividir RPK por el factor de escala para obtener TPM
    inv_scalers = 1 / per_million_scalers
    tpm_matrix = diags(inv_scalers.flatten()).dot(rpk_matrix)
    
    return tpm_matrix.tocsr()
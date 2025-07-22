import requests
import pandas as pd
import time

def fetch_gene_symbols_ensembl(ensembl_ids):
    server = "https://rest.ensembl.org"
    endpoint = "/lookup/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    batch_size = 1000
    result = {}

    for i in range(0, len(ensembl_ids), batch_size):
        batch_ids = ensembl_ids[i:i+batch_size]
        data = {"ids": batch_ids}

        response = requests.post(server + endpoint, headers=headers, json=data)

        if not response.ok:
            print(f"Error en la petición: {response.status_code}")
            continue

        decoded = response.json()

        for ensembl_id, info in decoded.items():
            if info is not None:
                result[ensembl_id] = info.get("display_name", None)
            else:
                result[ensembl_id] = None

        time.sleep(0.1)

    return result
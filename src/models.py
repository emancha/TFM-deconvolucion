import os
import joblib
import pandas as pd
import numpy as np

# Importaciones de Modelos y Métricas
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from xgboost import XGBClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score

# --- Parámetros Comunes ---
RANDOM_STATE = 42


def train_and_evaluate_model(
    model_type, 
    X_train, y_train, 
    X_test, y_test, 
    output_dir, 
    model_name,
    **kwargs
):
    """
    Función centralizada para entrenar, evaluar y guardar un modelo de clasificación.

    Parámetros:
    ----------
    model_type : str
        El tipo de modelo a entrenar ('rf', 'xgb', 'mlp').
    X_train, y_train, X_test, y_test : array-like
        Los datos de entrenamiento y prueba.
    output_dir : str
        El directorio base para guardar modelos y figuras.
    model_name : str
        Un nombre único para esta ejecución del modelo (ej. 'rf_marker_genes').
    **kwargs : dict
        Argumentos adicionales para pasar al constructor del modelo.

    Retorna:
    -------
    tuple
        Contiene el modelo entrenado, el reporte de clasificación (str) y la matriz de confusión.
    """
    
    # --- Selección y Entrenamiento del Modelo ---
    model, le = None, None # le para el LabelEncoder de XGBoost
    
    if model_type == 'rf':
        print(f"--- Entrenando RandomForestClassifier: {model_name} ---")
        model = RandomForestClassifier(random_state=RANDOM_STATE, n_jobs=-1, **kwargs)
        model.fit(X_train, y_train)
        
    elif model_type == 'xgb':
        print(f"--- Entrenando XGBClassifier: {model_name} ---")
        le = LabelEncoder()
        y_train_encoded = le.fit_transform(y_train)
        
        model = XGBClassifier(
            random_state=RANDOM_STATE, 
            n_jobs=-1, 
            use_label_encoder=False, # Recomendado para evitar warnings
            eval_metric='mlogloss', 
            **kwargs
        )
        model.fit(X_train, y_train_encoded)

    elif model_type == 'mlp':
        print(f"--- Entrenando MLPClassifier: {model_name} ---")
        # Asegurarse de que los datos son densos para MLP
        if hasattr(X_train, 'toarray'):
            X_train = X_train.toarray()
            
        default_params = {
            'hidden_layer_sizes': (100, 50),
            'max_iter': 300,
            'early_stopping': True,
            'alpha': 0.001
        }
        default_params.update(kwargs)
        
        model = MLPClassifier(random_state=RANDOM_STATE, **default_params)
        model.fit(X_train, y_train)
        
    else:
        raise ValueError("Tipo de modelo no soportado. Elige entre 'rf', 'xgb', 'mlp'.")

    # --- Evaluación del Modelo ---
    print("Evaluando el modelo...")
    
    y_pred = None
    if model_type == 'xgb':
        y_pred_encoded = model.predict(X_test)
        y_pred = le.inverse_transform(y_pred_encoded)
    elif model_type == 'mlp':
        if hasattr(X_test, 'toarray'):
            X_test = X_test.toarray()
        y_pred = model.predict(X_test)
    else: # rf
        y_pred = model.predict(X_test)
        
    accuracy = accuracy_score(y_test, y_pred)
    report = classification_report(y_test, y_pred)
    cm = confusion_matrix(y_test, y_pred, labels=model.classes_ if model_type != 'xgb' else le.classes_)
    
    print(f"Precisión (Accuracy): {accuracy * 100:.2f}%")

    # --- Guardado del Modelo ---
    model_dir = os.path.join(output_dir, 'models')
    os.makedirs(model_dir, exist_ok=True)
    model_path = os.path.join(model_dir, f"{model_name}.joblib")
    
    if le:
        joblib.dump({'model': model, 'label_encoder': le}, model_path)
    else:
        joblib.dump(model, model_path)
    print(f"Modelo guardado en: {model_path}")
    
    return model, report, cm
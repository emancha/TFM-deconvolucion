from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
import joblib

def train_and_evaluate_rf(X_train, y_train, X_test, y_test, model_name="random_forest"):
    """Entrena y evalúa un modelo Random Forest."""
    
    print(f"--- Entrenando el modelo: {model_name} ---")
    model = RandomForestClassifier(random_state=42)
    model.fit(X_train, y_train)
    
    print("Evaluando el modelo...")
    y_pred = model.predict(X_test)
    
    report = classification_report(y_test, y_pred)
    cm = confusion_matrix(y_test, y_pred)
    
    print("Guardando el modelo...")
    joblib.dump(model, f"models/{model_name}.joblib")
    
    return model, report, cm

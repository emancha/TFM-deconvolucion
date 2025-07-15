import scanpy as sc
import joblib
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report
import os


final_adata= sc.read_h5ad(r"C:\Users\emanc\OneDrive - Universidad de Extremadura (1)\MUBinf\09.TFM\Datos\lung_cancer_filtrado.h5ad", backed="r")

X = final_adata.X
y = final_adata.obs['cell_type']

# 2. Dividir los datos
X_train, X_test, y_train, y_test = train_test_split(
    X, 
    y, 
    test_size=0.2,       # 20% para test
    random_state=42,     # Para reproducibilidad
    stratify=y           # ¡CRÍTICO! Mantiene la proporción de tipos celulares
)

#Carga del modelo
filename = r'C:\Users\emanc\OneDrive - Universidad de Extremadura (1)\MUBinf\09.TFM\Datos\random_forest_lung_cancer_v1.joblib'

nombre = os.path.splitext(os.path.basename(filename))[0]

modelo = joblib.load(filename)

print("Modelo cargado exitosamente.")

##Overall Accuracy

#Predicción del modelo de los tipos celulares para el set de test
y_pred = modelo.predict(X_test)

#Comparar las predicciones (y_pred) con las etiquetas reales (y_test)
accuracy = accuracy_score(y_test, y_pred)

print(f"Precisión del modelo (Accuracy): {accuracy * 100:.2f}%")

##Confusion Matrix

# Genera la matriz de confusión
cm = confusion_matrix(y_test, y_pred)

# Para visualizarla de forma más bonita
plt.figure(figsize=(10, 8))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
            xticklabels=modelo.classes_, yticklabels=modelo.classes_)
plt.ylabel('Etiqueta Real (True Label)')
plt.xlabel('Etiqueta Predicha (Predicted Label)')
plt.title(f'Matriz de Confusión_{nombre}')
print("Matriz de confuión generada")

plt.savefig(f'C:/Users/emanc/OneDrive/Documentos/MUBinf/Asignaturas/09. TFM/TFM_repo/data/processed/confusion_matrix_{nombre}.pdf')
plt.show()
print("Matriz de confusión guardada ")

##Classification Report

# Genera el reporte
report = classification_report(y_test, y_pred)
print("Reporte de Clasificación:\n")
print(report)

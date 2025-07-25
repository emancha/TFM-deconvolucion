import seaborn as sns
import matplotlib.pyplot as plt

def plot_confusion_matrix(cm, classes, title='Matriz de Confusión', save_path=None):
    """Plotea una matriz de confusión y la guarda."""
    plt.figure(figsize=(12, 10))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                xticklabels=classes, yticklabels=classes)
    plt.ylabel('Etiqueta Real')
    plt.xlabel('Etiqueta Predicha')
    plt.title(title)
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
        
    plt.show()


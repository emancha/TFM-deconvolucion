from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import numpy as np

def plot_kaplan_meier(df, cell_type, ax):
    """Función para generar una curva de Kaplan-Meier para un tipo celular."""
    
    median_val = df[cell_type].median()
    
    # Problema al analizar muestras de quantiseq (algunas tienen mediana 0)
    if median_val == 0.0:
        # Si la mediana es 0, los grupos son "cero" vs. "mayor que cero"
        df['group'] = np.where(df[cell_type] > 0, 'Presente', 'Ausente')
        
        presente_group = df[df['group'] == 'Presente']
        ausente_group = df[df['group'] == 'Ausente']
        
        # Verificamos si ambos grupos tienen muestras
        if presente_group.empty or ausente_group.empty:
            ax.set_title(f'{cell_type}\n(No hay suficiente variación para el análisis)')
            ax.text(0.5, 0.5, 'Datos insuficientes', ha='center', va='center', fontsize=12)
            ax.axis('off')
            return
        
        kmf_presente = KaplanMeierFitter()
        kmf_ausente = KaplanMeierFitter()
        
        kmf_presente.fit(presente_group['survival_time'], presente_group['event_status'], label=f'{cell_type} Presente (n={len(presente_group)})')
        kmf_ausente.fit(ausente_group['survival_time'], ausente_group['event_status'], label=f'{cell_type} Ausente (n={len(ausente_group)})')
        
        results = logrank_test(
            presente_group['survival_time'], ausente_group['survival_time'],
            event_observed_A=presente_group['event_status'], event_observed_B=ausente_group['event_status']
        )
        
        kmf_presente.plot_survival_function(ax=ax)
        kmf_ausente.plot_survival_function(ax=ax)
        
    else:
        # Lógica para cuando la mediana es > 0
        df['group'] = np.where(df[cell_type] >= median_val, 'Alta', 'Baja')
    
        alta_group = df[df['group'] == 'Alta']
        baja_group = df[df['group'] == 'Baja']
        
        kmf_alta = KaplanMeierFitter()
        kmf_baja = KaplanMeierFitter()
        
        kmf_alta.fit(alta_group['survival_time'], alta_group['event_status'], label=f'Alta {cell_type} (n={len(alta_group)})')
        kmf_baja.fit(baja_group['survival_time'], baja_group['event_status'], label=f'Baja {cell_type} (n={len(baja_group)})')
        
        results = logrank_test(
            alta_group['survival_time'], baja_group['survival_time'],
            event_observed_A=alta_group['event_status'], event_observed_B=baja_group['event_status']
        )
        
        kmf_alta.plot_survival_function(ax=ax)
        kmf_baja.plot_survival_function(ax=ax)
    
    ax.set_title(f'Supervivencia según la fracción de {cell_type}\nLog-Rank p-value: {results.p_value:.3f}')
    ax.set_xlabel('Tiempo (días)')
    ax.set_ylabel('Probabilidad de Supervivencia')
    ax.grid(True)

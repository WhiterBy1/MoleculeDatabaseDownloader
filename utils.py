# utils.py
"""
Utilidades generales del sistema
INPUT: various
OUTPUT: helper functions
"""

import os
import sys
import re
from pathlib import Path

def validate_molecule_name(name):
    """
    Valida si un nombre de molécula es válido
    
    INPUT:
    - name (str): Nombre a validar
    
    OUTPUT:
    - bool: True si es válido
    """
    if not name or not name.strip():
        return False
    
    # Muy corto o muy largo
    if len(name.strip()) < 2 or len(name.strip()) > 100:
        return False
    
    # Solo números (probablemente no es un nombre válido)
    if name.strip().isdigit():
        return False
    
    return True

def format_file_size(size_bytes):
    """
    Formatea un tamaño de archivo en bytes a formato legible
    
    INPUT:
    - size_bytes (int): Tamaño en bytes
    
    OUTPUT:
    - str: Tamaño formateado
    """
    if size_bytes == 0:
        return "0 B"
    
    size_names = ["B", "KB", "MB", "GB", "TB"]
    import math
    i = int(math.floor(math.log(size_bytes, 1024)))
    p = math.pow(1024, i)
    s = round(size_bytes / p, 2)
    return f"{s} {size_names[i]}"

def open_folder_in_explorer(folder_path):
    """
    Abre una carpeta en el explorador del sistema
    
    INPUT:
    - folder_path (str): Ruta de la carpeta
    
    OUTPUT:
    - bool: True si se abrió correctamente
    """
    try:
        folder_path = Path(folder_path)
        if not folder_path.exists():
            return False
        
        if sys.platform == "win32":
            os.startfile(folder_path)
        elif sys.platform == "darwin":  # macOS
            os.system(f"open '{folder_path}'")
        else:  # Linux
            os.system(f"xdg-open '{folder_path}'")
        
        return True
    except Exception:
        return False

def clean_molecule_list(molecule_list):
    """
    Limpia una lista de nombres de moléculas
    
    INPUT:
    - molecule_list (list): Lista de nombres sin limpiar
    
    OUTPUT:
    - list: Lista limpia y validada
    """
    cleaned = []
    
    for molecule in molecule_list:
        if isinstance(molecule, str):
            clean_name = molecule.strip()
            if validate_molecule_name(clean_name):
                cleaned.append(clean_name)
    
    # Remover duplicados manteniendo orden
    return list(dict.fromkeys(cleaned))

def estimate_download_time(num_molecules, avg_time_per_molecule=3):
    """
    Estima el tiempo de descarga para un lote de moléculas
    
    INPUT:
    - num_molecules (int): Número de moléculas
    - avg_time_per_molecule (int): Tiempo promedio por molécula en segundos
    
    OUTPUT:
    - str: Tiempo estimado formateado
    """
    total_seconds = num_molecules * avg_time_per_molecule
    
    if total_seconds < 60:
        return f"{total_seconds} segundos"
    elif total_seconds < 3600:
        minutes = total_seconds // 60
        return f"{minutes} minutos"
    else:
        hours = total_seconds // 3600
        minutes = (total_seconds % 3600) // 60
        return f"{hours}h {minutes}m"

def validate_excel_file(file_path):
    """
    Valida si un archivo Excel es válido y accesible
    
    INPUT:
    - file_path (str): Ruta al archivo Excel
    
    OUTPUT:
    - tuple: (is_valid: bool, error_message: str)
    """
    try:
        file_path = Path(file_path)
        
        if not file_path.exists():
            return False, "El archivo no existe"
        
        if file_path.suffix.lower() not in ['.xlsx', '.xls']:
            return False, "El archivo no es un Excel válido"
        
        if file_path.stat().st_size == 0:
            return False, "El archivo está vacío"
        
        # Intentar abrirlo brevemente
        import pandas as pd
        try:
            pd.read_excel(file_path, nrows=1)
        except Exception as e:
            return False, f"No se puede leer el archivo: {str(e)}"
        
        return True, ""
        
    except Exception as e:
        return False, f"Error validando archivo: {str(e)}"

def create_backup_name(original_name):
    """
    Crea un nombre de respaldo único
    
    INPUT:
    - original_name (str): Nombre original
    
    OUTPUT:
    - str: Nombre único con timestamp
    """
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"{original_name}_backup_{timestamp}"

def parse_molecule_identifier(identifier):
    """
    Intenta determinar qué tipo de identificador es
    
    INPUT:
    - identifier (str): Identificador de molécula
    
    OUTPUT:
    - dict: Información del identificador
    """
    identifier = identifier.strip()
    
    result = {
        "original": identifier,
        "type": "name",  # name, cid, smiles, inchi
        "cleaned": identifier
    }
    
    # Es un CID de PubChem?
    if identifier.isdigit() and len(identifier) <= 10:
        result["type"] = "cid"
        result["cleaned"] = identifier
    
    # Es un SMILES?
    elif any(char in identifier for char in ['(', ')', '=', '#', '@', '[', ']']):
        result["type"] = "smiles"
    
    # Es un InChI?
    elif identifier.startswith("InChI="):
        result["type"] = "inchi"
    
    # Limpiar nombre si es tipo "name"
    elif result["type"] == "name":
        # Remover caracteres problemáticos pero mantener estructura química
        cleaned = re.sub(r'[^\w\s\-\(\),\.]', '', identifier)
        result["cleaned"] = cleaned.strip()
    
    return result

class ProgressTracker:
    """Clase simple para trackear progreso"""
    
    def __init__(self, total_items):
        self.total_items = total_items
        self.completed_items = 0
        self.failed_items = 0
    
    def mark_completed(self):
        self.completed_items += 1
    
    def mark_failed(self):
        self.failed_items += 1
    
    def get_progress_percent(self):
        if self.total_items == 0:
            return 0
        return (self.completed_items + self.failed_items) / self.total_items * 100
    
    def get_summary(self):
        return {
            "total": self.total_items,
            "completed": self.completed_items,
            "failed": self.failed_items,
            "pending": self.total_items - self.completed_items - self.failed_items,
            "progress_percent": self.get_progress_percent()
        }

def log_system_info():
    """
    Registra información del sistema para debugging
    
    OUTPUT:
    - dict: Información del sistema
    """
    import platform
    
    info = {
        "platform": platform.system(),
        "python_version": platform.python_version(),
        "working_directory": os.getcwd()
    }
    
    # Verificar dependencias críticas
    try:
        import rdkit
        info["rdkit_version"] = rdkit.__version__
    except ImportError:
        info["rdkit_version"] = "Not installed"
    
    try:
        import pandas
        info["pandas_version"] = pandas.__version__
    except ImportError:
        info["pandas_version"] = "Not installed"
    
    return info
# config.py
"""
Configuraciones centralizadas del sistema
INPUT: None
OUTPUT: Constantes y configuraciones
"""

from pathlib import Path

# Directorios
BASE_DIR = Path(__file__).parent
DEFAULT_OUTPUT_DIR = "./molecules"
PROGRESS_DIR = "./progress"
TEMP_DIR = "./temp"

# APIs disponibles
SOURCES = {
    "pubchem": {
        "base_url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug",
        "search_endpoint": "/compound/name/{query}/cids/JSON",
        "download_3d": "/compound/cid/{id}/record/SDF/?record_type=3d",
        "download_2d": "/compound/cid/{id}/SDF"
    },
    "chemspider": {
        "base_url": "http://www.chemspider.com/API",
        # Agregar endpoints según necesidad
    }
}

# Configuración de descarga
TIMEOUT = 30
RETRY_ATTEMPTS = 3
DELAY_BETWEEN_REQUESTS = 1

# Formatos soportados
EXCEL_FORMATS = ['.xlsx', '.xls']
TEXT_FORMATS = ['.txt', '.csv']

# Categorización
CATEGORIZE_OPTIONS = ["formula", "weight", "atoms", "none"]
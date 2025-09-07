# sources.py
"""
Manejo de fuentes de datos para búsqueda de moléculas
INPUT: query (str), source_name (str)
OUTPUT: list of molecule_ids
"""

import requests
import time
from config import SOURCES, TIMEOUT, RETRY_ATTEMPTS, DELAY_BETWEEN_REQUESTS

class SourceManager:
    """Gestor unificado para todas las fuentes de datos"""
    
    def search_molecule(self, query, source_name="pubchem"):
        """
        Busca una molécula en la fuente especificada
        
        INPUT:
        - query (str): Nombre o ID de la molécula
        - source_name (str): Nombre de la fuente ("pubchem", "chemspider")
        
        OUTPUT:
        - list: Lista de IDs encontrados o lista vacía si no hay resultados
        """
        if source_name not in SOURCES:
            raise ValueError(f"Fuente no soportada: {source_name}")
        
        if source_name == "pubchem":
            return self._search_pubchem(query)
        elif source_name == "chemspider":
            return self._search_chemspider(query)
        
        return []
    
    def _search_pubchem(self, query):
        """Búsqueda específica en PubChem"""
        source_config = SOURCES["pubchem"]
        search_url = source_config["base_url"] + source_config["search_endpoint"].format(query=query)
        
        for attempt in range(RETRY_ATTEMPTS):
            try:
                response = requests.get(search_url, timeout=TIMEOUT)
                
                if response.status_code == 200:
                    data = response.json()
                    if "IdentifierList" in data and "CID" in data["IdentifierList"]:
                        return data["IdentifierList"]["CID"]
                elif response.status_code == 404:
                    return []  # No encontrado
                
                time.sleep(DELAY_BETWEEN_REQUESTS)
                
            except Exception as e:
                if attempt == RETRY_ATTEMPTS - 1:
                    print(f"Error en búsqueda PubChem: {e}")
                time.sleep(DELAY_BETWEEN_REQUESTS * (attempt + 1))
        
        return []
    
    def _search_chemspider(self, query):
        """Búsqueda específica en ChemSpider (implementar según API)"""
        # Placeholder para ChemSpider
        print(f"ChemSpider search not implemented yet for: {query}")
        return []
    
    def get_download_url(self, molecule_id, source_name="pubchem", format_3d=True):
        """
        Obtiene la URL de descarga para una molécula
        
        INPUT:
        - molecule_id (str/int): ID de la molécula
        - source_name (str): Fuente de datos
        - format_3d (bool): True para 3D, False para 2D
        
        OUTPUT:
        - str: URL de descarga o None si no está disponible
        """
        if source_name not in SOURCES:
            return None
        
        source_config = SOURCES[source_name]
        
        if format_3d and "download_3d" in source_config:
            return source_config["base_url"] + source_config["download_3d"].format(id=molecule_id)
        elif "download_2d" in source_config:
            return source_config["base_url"] + source_config["download_2d"].format(id=molecule_id)
        
        return None
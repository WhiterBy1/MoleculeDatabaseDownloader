# downloader.py
"""
Lógica principal de descarga de moléculas
INPUT: molecule_id, source_name
OUTPUT: sdf_content (str) or None
"""

import requests
import time
from rdkit import Chem
from rdkit.Chem import AllChem
from config import TIMEOUT, RETRY_ATTEMPTS, DELAY_BETWEEN_REQUESTS
from storage import StorageManager

class MoleculeDownloader:
    """Descargador principal de moléculas"""
    
    def __init__(self, source_manager):
        """
        INPUT:
        - source_manager: Instancia de SourceManager
        """
        self.source_manager = source_manager
        self.storage_manager = StorageManager()
        
        # Configuración
        self.output_dir = "./molecules"
        self.categorize_by = "formula"
        self.progress_callback = None
        self.log_callback = None
    
    def configure(self, output_dir=None, categorize_by=None, progress_callback=None, log_callback=None):
        """
        Configura el descargador
        
        INPUT:
        - output_dir (str): Directorio de salida
        - categorize_by (str): Criterio de categorización
        - progress_callback (function): Callback para progreso
        - log_callback (function): Callback para logging
        """
        if output_dir:
            self.output_dir = output_dir
            self.storage_manager.set_output_dir(output_dir)
        
        if categorize_by:
            self.categorize_by = categorize_by
            self.storage_manager.set_categorize_by(categorize_by)
        
        self.progress_callback = progress_callback
        self.log_callback = log_callback
    
    def download_molecule(self, query, source_name="pubchem"):
        """
        Descarga una molécula completa: buscar + descargar + guardar
        
        INPUT:
        - query (str): Nombre o ID de la molécula
        - source_name (str): Fuente de datos
        
        OUTPUT:
        - bool: True si se descargó exitosamente
        """
        try:
            # 1. Buscar la molécula
            self._log(f"🔍 Buscando '{query}' en {source_name}...")
            molecule_ids = self.source_manager.search_molecule(query, source_name)
            
            if not molecule_ids:
                self._log(f"❌ No se encontraron resultados para '{query}'")
                return False
            
            # Usar el primer resultado
            molecule_id = molecule_ids[0]
            self._log(f"✅ Encontrado ID: {molecule_id}")
            
            # 2. Descargar estructura
            sdf_content = self.download_structure(molecule_id, source_name)
            
            if not sdf_content:
                self._log(f"❌ No se pudo descargar estructura para ID {molecule_id}")
                return False
            
            # 3. Procesar y guardar
            success = self.storage_manager.save_molecule(sdf_content, query, molecule_id)
            
            if success:
                self._log(f"✅ Molécula '{query}' guardada exitosamente")
                return True
            else:
                self._log(f"❌ Error al guardar molécula '{query}'")
                return False
                
        except Exception as e:
            self._log(f"❌ Error procesando '{query}': {e}")
            return False
    
    def download_structure(self, molecule_id, source_name="pubchem"):
        """
        Descarga la estructura SDF de una molécula
        
        INPUT:
        - molecule_id (str/int): ID de la molécula
        - source_name (str): Fuente de datos
        
        OUTPUT:
        - str: Contenido SDF o None si hay error
        """
        self._log(f"⬇️ Descargando estructura para ID {molecule_id}...")
        
        # Intentar 3D primero
        sdf_content = self._download_sdf(molecule_id, source_name, format_3d=True)
        
        if sdf_content:
            self._log("✅ Estructura 3D descargada")
            return self._process_structure(sdf_content)
        
        # Si no hay 3D, intentar 2D y convertir
        self._log("⚠️ Estructura 3D no disponible, intentando 2D...")
        sdf_content = self._download_sdf(molecule_id, source_name, format_3d=False)
        
        if sdf_content:
            self._log("✅ Estructura 2D descargada, convirtiendo a 3D...")
            return self._convert_2d_to_3d(sdf_content)
        
        self._log("❌ No se pudo descargar ninguna estructura")
        return None
    
    def _download_sdf(self, molecule_id, source_name, format_3d=True):
        """Descarga el archivo SDF desde la fuente"""
        download_url = self.source_manager.get_download_url(molecule_id, source_name, format_3d)
        
        if not download_url:
            return None
        
        for attempt in range(RETRY_ATTEMPTS):
            try:
                response = requests.get(download_url, timeout=TIMEOUT)
                
                if response.status_code == 200:
                    return response.text
                elif response.status_code == 404:
                    return None  # No disponible
                
                time.sleep(DELAY_BETWEEN_REQUESTS * (attempt + 1))
                
            except Exception as e:
                if attempt == RETRY_ATTEMPTS - 1:
                    self._log(f"❌ Error en descarga: {e}")
                time.sleep(DELAY_BETWEEN_REQUESTS * (attempt + 1))
        
        return None
    
    def _process_structure(self, sdf_content):
        """
        Procesa la estructura descargada (añadir hidrógenos si es necesario)
        
        INPUT:
        - sdf_content (str): Contenido SDF original
        
        OUTPUT:
        - str: Contenido SDF procesado
        """
        try:
            mol = Chem.MolFromMolBlock(sdf_content)
            if mol is None:
                self._log("⚠️ No se pudo leer la estructura")
                return sdf_content
            
            # Verificar si tiene hidrógenos explícitos
            h_count_explicit = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
            h_count_implicit = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())
            
            if h_count_explicit == 0 and h_count_implicit > 0:
                self._log("🔧 Añadiendo hidrógenos explícitos...")
                mol = Chem.AddHs(mol)
                
                # Verificar si tiene coordenadas 3D
                try:
                    conf = mol.GetConformer()
                    if not conf.Is3D():
                        self._log("🔧 Optimizando geometría 3D...")
                        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                        AllChem.UFFOptimizeMolecule(mol)
                except:
                    pass
                
                return Chem.MolToMolBlock(mol)
            
            return sdf_content
            
        except Exception as e:
            self._log(f"⚠️ Error procesando estructura: {e}")
            return sdf_content
    
    def _convert_2d_to_3d(self, sdf_content):
        """
        Convierte estructura 2D a 3D usando RDKit
        
        INPUT:
        - sdf_content (str): Contenido SDF 2D
        
        OUTPUT:
        - str: Contenido SDF 3D o None si falla
        """
        try:
            mol = Chem.MolFromMolBlock(sdf_content)
            if mol is None:
                self._log("❌ No se pudo leer estructura 2D")
                return None
            
            # Añadir hidrógenos
            mol = Chem.AddHs(mol)
            
            # Generar conformación 3D
            params = AllChem.ETKDG()
            params.randomSeed = 42  # Para reproducibilidad
            embed_result = AllChem.EmbedMolecule(mol, params)
            
            if embed_result == -1:
                self._log("⚠️ Intentando con parámetros alternativos...")
                params.useRandomCoords = True
                embed_result = AllChem.EmbedMolecule(mol, params)
            
            if embed_result != -1:
                # Optimizar geometría
                try:
                    AllChem.UFFOptimizeMolecule(mol, maxIters=500)
                    self._log("✅ Estructura 3D generada y optimizada")
                except:
                    self._log("⚠️ Optimización parcialmente exitosa")
                
                return Chem.MolToMolBlock(mol)
            else:
                self._log("❌ No se pudo generar conformación 3D")
                return None
                
        except Exception as e:
            self._log(f"❌ Error convirtiendo 2D a 3D: {e}")
            return None
    
    def _log(self, message):
        """Helper para logging"""
        print(message)
        if self.log_callback:
            self.log_callback(message)
# downloader.py
"""
Lógica principal de descarga de moléculas - VERSIÓN ROBUSTA
INPUT: molecule_id, source_name
OUTPUT: sdf_content (str) or None
"""

import requests
import time
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, rdMolDescriptors
from rdkit.Chem.rdDistGeom import EmbedMolecule
import numpy as np
from config import TIMEOUT, RETRY_ATTEMPTS, DELAY_BETWEEN_REQUESTS
from storage import StorageManager
#velocity optimizing
import concurrent.futures
import threading
from functools import lru_cache
import pickle
import os
import hashlib
import gzip
import json
from pathlib import Path

class MoleculeDownloader:
    """Descargador principal de moléculas"""
    
    def __init__(self, source_manager):
        """
        INPUT:
        - source_manager: Instancia de SourceManager
        """
        self.source_manager = source_manager
        self.storage_manager = StorageManager()
        
        #parallel config
        self.cache_dir = Path("./cache")
        self.cache_dir.mkdir(exist_ok=True)
        self.max_workers = 6  # Configurable
        self.use_parallel = True  # Flag para activar/desactivar
        
        # Session reutilizable para velocidad
        self.session = self._setup_optimized_session()
        
        # Configuración
        self.output_dir = "./molecules"
        self.categorize_by = "formula"
        self.progress_callback = None
        self.log_callback = None
    
    def _setup_optimized_session(self):
        """NUEVA: Sesión HTTP optimizada"""
        session = requests.Session()
        session.headers.update({
            'User-Agent': 'Mozilla/5.0 (compatible; MoleculeDownloader/2.0)',
            'Accept': 'text/plain,application/sdf',
            'Connection': 'keep-alive'
        })
        
        # Configuración de reintentos más agresiva
        from requests.adapters import HTTPAdapter
        from urllib3.util.retry import Retry
        
        retry_strategy = Retry(
            total=2,
            status_forcelist=[429, 500, 502, 503, 504],
            backoff_factor=0.1
        )
        
        adapter = HTTPAdapter(
            max_retries=retry_strategy,
            pool_connections=10,
            pool_maxsize=10
        )
        
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        return session
    
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
    
    
    
    # NUEVO: Método para lotes paralelos
    def download_molecules_batch(self, queries, source_name="pubchem"):
        """
        NUEVO: Descarga por lotes con paralelización automática
        Se activa automáticamente para >3 moléculas
        """
        if len(queries) <= 3 or not self.use_parallel:
            # Para pocas moléculas, usar método secuencial
            return self._download_sequential(queries, source_name)
        else:
            # Para muchas moléculas, usar paralelo
            return self._download_parallel(queries, source_name)
    
    def _download_sequential(self, queries, source_name):
        """Descarga secuencial (método original)"""
        results = {}
        for query in queries:
            results[query] = self.download_molecule(query, source_name)
        return results
    
    def _download_parallel(self, queries, source_name):
        """NUEVA: Descarga paralela optimizada"""
        self._log(f"🚀 Iniciando descarga paralela: {len(queries)} moléculas con {self.max_workers} hilos")
        
        results = {}
        start_time = time.time()
        
        # Fase 1: Búsqueda paralela de IDs
        molecule_ids = self._parallel_search(queries, source_name)
        
        # Fase 2: Descarga paralela
        valid_queries = [(query, ids[0]) for query, ids in molecule_ids.items() if ids]
        failed_searches = [query for query, ids in molecule_ids.items() if not ids]
        
        # Marcar búsquedas fallidas
        for query in failed_searches:
            results[query] = False
        
        # Descargar las válidas
        if valid_queries:
            download_results = self._parallel_download_structures(valid_queries, source_name)
            results.update(download_results)
        
        # Estadísticas
        end_time = time.time()
        successful = sum(1 for success in results.values() if success)
        self._log(f"⚡ Completado en {end_time - start_time:.2f}s - {successful}/{len(queries)} exitosas")
        
        return results
    
    def _parallel_search(self, queries, source_name):
        """Búsqueda de IDs en paralelo con cache"""
        results = {}
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_query = {
                executor.submit(self._search_with_cache, query, source_name): query 
                for query in queries
            }
            
            for future in concurrent.futures.as_completed(future_to_query):
                query = future_to_query[future]
                try:
                    molecule_ids = future.result()
                    results[query] = molecule_ids
                except Exception as e:
                    self._log(f"❌ Error buscando {query}: {e}")
                    results[query] = []
        
        return results
    
    @lru_cache(maxsize=500)
    def _search_with_cache(self, query, source_name):
        """Búsqueda con cache para evitar repeticiones"""
        cache_file = self.cache_dir / f"search_{hashlib.md5((query + source_name).encode()).hexdigest()}.pkl"
        
        # Verificar cache
        if cache_file.exists():
            try:
                with open(cache_file, 'rb') as f:
                    cached_result = pickle.load(f)
                    if time.time() - cached_result['timestamp'] < 86400:  # 24 horas
                        return cached_result['ids']
            except:
                pass
        
        # Búsqueda real
        molecule_ids = self.source_manager.search_molecule(query, source_name)
        
        # Guardar en cache
        try:
            with open(cache_file, 'wb') as f:
                pickle.dump({
                    'ids': molecule_ids,
                    'timestamp': time.time()
                }, f)
        except:
            pass
        
        return molecule_ids
    
    def _parallel_download_structures(self, query_id_pairs, source_name):
        """Descarga de estructuras en paralelo"""
        results = {}
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_query = {
                executor.submit(self._download_single_structure, query, mol_id, source_name): query 
                for query, mol_id in query_id_pairs
            }
            
            for future in concurrent.futures.as_completed(future_to_query):
                query = future_to_query[future]
                try:
                    success = future.result()
                    results[query] = success
                except Exception as e:
                    self._log(f"❌ Error descargando {query}: {e}")
                    results[query] = False
        
        return results
    
    
    def _download_single_structure(self, query, molecule_id, source_name):
        """Descarga individual optimizada"""
        try:
            # Verificar cache
            cache_file = self.cache_dir / f"structure_{molecule_id}_{source_name}.sdf"
            if cache_file.exists():
                try:
                    with open(cache_file, 'r') as f:
                        sdf_content = f.read()
                        if self._validate_3d_structure_fast(sdf_content):
                            # Guardar desde cache
                            success = self.storage_manager.save_molecule(sdf_content, query, molecule_id)
                            if success:
                                self._log(f"💾 {query} (desde cache)")
                            return success
                except:
                    pass
            
            # Descargar nueva
            sdf_content = self._download_structure_optimized(molecule_id, source_name)
            
            if sdf_content:
                # Guardar en cache
                try:
                    with open(cache_file, 'w') as f:
                        f.write(sdf_content)
                except:
                    pass
                
                # Guardar molécula
                success = self.storage_manager.save_molecule(sdf_content, query, molecule_id)
                if success:
                    self._log(f"✅ {query}")
                return success
            
            return False
            
        except Exception as e:
            self._log(f"❌ Error con {query}: {e}")
            return False
    
    def _download_structure_optimized(self, molecule_id, source_name):
        """Descarga optimizada usando session reutilizable"""
        # Usar 2D para velocidad y confiabilidad
        download_url = self.source_manager.get_download_url(molecule_id, source_name, format_3d=False)
        
        if not download_url:
            return None
        
        try:
            response = self.session.get(download_url, timeout=10)
            if response.status_code == 200:
                return self._convert_2d_to_3d_fast(response.text)
        except:
            pass
        
        return None
    
    def _convert_2d_to_3d_fast(self, sdf_content):
        """Conversión 2D→3D optimizada para velocidad"""
        try:
            mol = Chem.MolFromMolBlock(sdf_content)
            if mol is None:
                return None
            
            mol = Chem.RemoveHs(mol)
            mol = Chem.AddHs(mol)
            
            # Parámetros optimizados para velocidad
            params = AllChem.ETKDG()
            params.randomSeed = 42
            params.maxAttempts = 25  # Reducido para velocidad
            params.useExpTorsionAnglePrefs = False  # Más rápido
            
            if AllChem.EmbedMolecule(mol, params) != -1:
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
                result_sdf = Chem.MolToMolBlock(mol)
                
                if self._validate_3d_structure_fast(result_sdf):
                    return result_sdf
            
            return None
            
        except Exception:
            return None
    
    def _validate_3d_structure_fast(self, sdf_content):
        """Validación rápida"""
        try:
            mol = Chem.MolFromMolBlock(sdf_content)
            if mol is None or mol.GetNumConformers() == 0:
                return False
            
            conf = mol.GetConformer()
            zero_coords = 0
            
            # Solo verificar primeros 5 átomos para velocidad
            check_atoms = min(5, mol.GetNumAtoms())
            for i in range(check_atoms):
                pos = conf.GetAtomPosition(i)
                if abs(pos.x) < 0.001 and abs(pos.y) < 0.001 and abs(pos.z) < 0.001:
                    zero_coords += 1
            
            return zero_coords <= 1
            
        except Exception:
            return False
    
    # Configuración de paralelización
    def set_parallel_config(self, max_workers=6, use_parallel=True):
        """Configura parámetros de paralelización"""
        self.max_workers = max_workers
        self.use_parallel = use_parallel
        self._log(f"🔧 Configuración paralela: {max_workers} hilos, activo: {use_parallel}")
    
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
        
        # Estrategia: Siempre descargar 2D primero para tener estructura limpia
        self._log("📥 Descargando estructura 2D para base limpia...")
        sdf_2d = self._download_sdf(molecule_id, source_name, format_3d=False)
        
        if sdf_2d:
            self._log("✅ Estructura 2D descargada, generando 3D optimizada...")
            result_3d = self._convert_2d_to_3d_robust(sdf_2d)
            if result_3d:
                return result_3d
        
        # Fallback: intentar 3D de PubChem pero procesarla
        self._log("⚠️ Intentando 3D de PubChem como fallback...")
        sdf_3d = self._download_sdf(molecule_id, source_name, format_3d=True)
        
        if sdf_3d:
            result_cleaned = self._clean_broken_3d_structure(sdf_3d)
            if result_cleaned and self._validate_3d_structure(result_cleaned):
                return result_cleaned
        
        self._log("❌ No se pudo obtener estructura válida")
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
    
    def _clean_broken_3d_structure(self, sdf_content):
        """
        Limpia estructuras 3D con hidrógenos en (0,0,0)
        
        INPUT:
        - sdf_content (str): Contenido SDF con hidrógenos problemáticos
        
        OUTPUT:
        - str: Contenido SDF corregido o None
        """
        try:
            self._log("🧹 Limpiando estructura 3D problemática...")
            
            mol = Chem.MolFromMolBlock(sdf_content)
            if mol is None:
                return None
            
            # Estrategia: Mantener átomos pesados, regenerar hidrógenos
            mol_heavy = Chem.RemoveHs(mol)  # Quitar todos los hidrógenos
            
            if mol.GetNumConformers() > 0:
                # Extraer coordenadas de átomos pesados
                old_conf = mol.GetConformer()
                heavy_coords = []
                
                heavy_idx = 0
                for atom_idx in range(mol.GetNumAtoms()):
                    atom = mol.GetAtomWithIdx(atom_idx)
                    if atom.GetAtomicNum() != 1:  # No es hidrógeno
                        pos = old_conf.GetAtomPosition(atom_idx)
                        heavy_coords.append((pos.x, pos.y, pos.z))
                        heavy_idx += 1
                
                # Crear nueva molécula con hidrógenos
                mol_with_h = Chem.AddHs(mol_heavy)
                
                # Aplicar coordenadas a átomos pesados
                new_conf = Chem.Conformer(mol_with_h.GetNumAtoms())
                heavy_idx = 0
                
                for atom_idx in range(mol_with_h.GetNumAtoms()):
                    atom = mol_with_h.GetAtomWithIdx(atom_idx)
                    if atom.GetAtomicNum() != 1:  # Átomo pesado
                        if heavy_idx < len(heavy_coords):
                            x, y, z = heavy_coords[heavy_idx]
                            new_conf.SetAtomPosition(atom_idx, [x, y, z])
                            heavy_idx += 1
                    else:
                        # Hidrógeno: posición temporal (se optimizará)
                        new_conf.SetAtomPosition(atom_idx, [0, 0, 0])
                
                mol_with_h.AddConformer(new_conf, assignId=True)
                
                # Optimizar solo hidrógenos manteniendo átomos pesados fijos
                self._optimize_hydrogen_positions(mol_with_h)
                
                result_sdf = Chem.MolToMolBlock(mol_with_h)
                
                if self._validate_3d_structure(result_sdf):
                    self._log("✅ Estructura 3D limpiada exitosamente")
                    return result_sdf
            
            return None
            
        except Exception as e:
            self._log(f"❌ Error limpiando estructura 3D: {e}")
            return None
    
    def _convert_2d_to_3d_robust(self, sdf_content):
        """
        Conversión robusta de 2D a 3D con múltiples estrategias
        
        INPUT:
        - sdf_content (str): Contenido SDF 2D
        
        OUTPUT:
        - str: Contenido SDF 3D o None
        """
        try:
            mol = Chem.MolFromMolBlock(sdf_content)
            if mol is None:
                self._log("❌ No se pudo leer estructura 2D")
                return None
            
            # Preparar molécula limpia
            mol = Chem.RemoveHs(mol)  # Quitar hidrógenos si los hay
            mol = Chem.AddHs(mol)     # Añadir hidrógenos limpios
            
            # Estrategia 1: ETKDG mejorado con múltiples intentos
            self._log("🔄 Generando 3D con ETKDG mejorado...")
            for attempt in range(5):
                try:
                    mol_copy = Chem.Mol(mol)  # Copia para cada intento
                    
                    params = AllChem.ETKDG()
                    params.randomSeed = 42 + attempt * 123
                    params.maxAttempts = 100
                    params.numThreads = 1
                    params.useExpTorsionAnglePrefs = True
                    params.useBasicKnowledge = True
                    params.enforceChirality = True
                    
                    embed_result = AllChem.EmbedMolecule(mol_copy, params)
                    
                    if embed_result != -1:
                        # Optimización en etapas
                        self._multi_stage_optimization(mol_copy)
                        
                        result_sdf = Chem.MolToMolBlock(mol_copy)
                        if self._validate_3d_structure(result_sdf):
                            self._log(f"✅ ETKDG exitoso en intento {attempt + 1}")
                            return result_sdf
                
                except Exception as e:
                    self._log(f"⚠️ ETKDG intento {attempt + 1} falló: {e}")
                    continue
            
            # Estrategia 2: Distance Geometry con conformaciones múltiples
            self._log("🔄 Intentando Distance Geometry...")
            try:
                mol_copy = Chem.Mol(mol)
                
                # Generar múltiples conformaciones y elegir la mejor
                conf_ids = rdDistGeom.EmbedMultipleConfs(
                    mol_copy, 
                    numConfs=10, 
                    randomSeed=42,
                    clearConfs=True,
                    useExpTorsionAnglePrefs=True,
                    useBasicKnowledge=True
                )
                
                if conf_ids:
                    best_conf_id = None
                    best_energy = float('inf')
                    
                    for conf_id in conf_ids:
                        try:
                            # Optimizar conformación
                            ff = AllChem.UFFGetMoleculeForceField(mol_copy, confId=conf_id)
                            if ff:
                                ff.Minimize(maxIts=500)
                                energy = ff.CalcEnergy()
                                
                                if energy < best_energy:
                                    best_energy = energy
                                    best_conf_id = conf_id
                        except:
                            continue
                    
                    if best_conf_id is not None:
                        result_sdf = Chem.MolToMolBlock(mol_copy, confId=best_conf_id)
                        if self._validate_3d_structure(result_sdf):
                            self._log("✅ Distance Geometry exitoso")
                            return result_sdf
            
            except Exception as e:
                self._log(f"⚠️ Distance Geometry falló: {e}")
            
            # Estrategia 3: Generación básica forzada
            self._log("🔄 Generación básica como último recurso...")
            try:
                mol_copy = Chem.Mol(mol)
                
                params = AllChem.ETKDG()
                params.useRandomCoords = True
                params.randomSeed = -1
                params.maxAttempts = 200
                params.enforceChirality = False
                params.useExpTorsionAnglePrefs = False
                
                if AllChem.EmbedMolecule(mol_copy, params) != -1:
                    # Optimización básica
                    try:
                        AllChem.UFFOptimizeMolecule(mol_copy, maxIters=1000)
                    except:
                        pass
                    
                    result_sdf = Chem.MolToMolBlock(mol_copy)
                    self._log("⚠️ Generada conformación básica")
                    return result_sdf
            
            except Exception as e:
                self._log(f"❌ Generación básica falló: {e}")
            
            return None
            
        except Exception as e:
            self._log(f"❌ Error en conversión 2D→3D: {e}")
            return None
    
    def _multi_stage_optimization(self, mol):
        """
        Optimización en múltiples etapas para mejor geometría
        
        INPUT:
        - mol: Molécula RDKit con conformación
        """
        try:
            # Etapa 1: Optimización suave
            ff = AllChem.UFFGetMoleculeForceField(mol)
            if ff:
                ff.Minimize(maxIts=100, forceTol=1e-3)
            
            # Etapa 2: Optimización normal
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            
            # Etapa 3: Optimización fina
            ff = AllChem.UFFGetMoleculeForceField(mol)
            if ff:
                ff.Minimize(maxIts=200, forceTol=1e-4)
        
        except Exception as e:
            self._log(f"⚠️ Optimización parcial: {e}")
    
    def _optimize_hydrogen_positions(self, mol):
        """
        Optimiza solo las posiciones de hidrógenos manteniendo átomos pesados fijos
        
        INPUT:
        - mol: Molécula RDKit con conformación
        """
        try:
            ff = AllChem.UFFGetMoleculeForceField(mol)
            if not ff:
                return
            
            # Fijar todos los átomos pesados
            for atom_idx in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() != 1:  # No es hidrógeno
                    ff.AddFixedPoint(atom_idx)
            
            # Optimizar solo hidrógenos
            ff.Minimize(maxIts=500)
            
        except Exception as e:
            self._log(f"⚠️ Optimización de hidrógenos parcial: {e}")
    
    def _validate_3d_structure(self, sdf_content):
        """
        Valida si la estructura 3D tiene coordenadas válidas
        
        INPUT:
        - sdf_content (str): Contenido SDF
        
        OUTPUT:
        - bool: True si tiene coordenadas 3D válidas
        """
        try:
            mol = Chem.MolFromMolBlock(sdf_content)
            if mol is None or mol.GetNumConformers() == 0:
                return False
            
            conf = mol.GetConformer()
            total_atoms = mol.GetNumAtoms()
            
            # Contar átomos en (0,0,0)
            zero_coords = 0
            all_coords = []
            
            for i in range(total_atoms):
                pos = conf.GetAtomPosition(i)
                all_coords.append([pos.x, pos.y, pos.z])
                
                if abs(pos.x) < 0.001 and abs(pos.y) < 0.001 and abs(pos.z) < 0.001:
                    zero_coords += 1
            
            # Test 1: No más de 1 átomo en (0,0,0)
            if zero_coords > 1:
                self._log(f"❌ Validación fallida: {zero_coords} átomos en (0,0,0)")
                return False
            
            # Test 2: Verificar dispersión espacial
            coords_array = np.array(all_coords)
            
            # Calcular desviación estándar para cada dimensión
            std_x = np.std(coords_array[:, 0])
            std_y = np.std(coords_array[:, 1]) 
            std_z = np.std(coords_array[:, 2])
            
            # Debe haber dispersión en las 3 dimensiones
            if std_x < 0.1 or std_y < 0.1 or std_z < 0.1:
                self._log(f"❌ Validación fallida: Baja dispersión 3D (σx={std_x:.3f}, σy={std_y:.3f}, σz={std_z:.3f})")
                return False
            
            # Test 3: Verificar que es realmente 3D
            if hasattr(conf, 'Is3D') and not conf.Is3D():
                self._log("❌ Validación fallida: No es estructura 3D")
                return False
            
            self._log(f"✅ Estructura válida: {total_atoms} átomos, dispersión 3D adecuada")
            return True
            
        except Exception as e:
            self._log(f"⚠️ Error en validación: {e}")
            return False
    
    def _log(self, message):
        """Helper para logging"""
        print(message)
        if self.log_callback:
            self.log_callback(message)
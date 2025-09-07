# storage.py
"""
Gestión de almacenamiento y organización de archivos
INPUT: sdf_content, molecule_properties
OUTPUT: saved_file_path or success status
"""

import re
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors

class StorageManager:
    """Gestor de almacenamiento y organización de moléculas"""
    
    def __init__(self, output_dir="./molecules", categorize_by="formula"):
        """
        INPUT:
        - output_dir (str): Directorio base de salida
        - categorize_by (str): Criterio de categorización
        """
        self.output_dir = Path(output_dir)
        self.categorize_by = categorize_by
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def set_output_dir(self, output_dir):
        """
        Cambia el directorio de salida
        
        INPUT:
        - output_dir (str): Nuevo directorio
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def set_categorize_by(self, categorize_by):
        """
        Cambia el criterio de categorización
        
        INPUT:
        - categorize_by (str): Nuevo criterio ("formula", "weight", "atoms", "none")
        """
        self.categorize_by = categorize_by
    
    def save_molecule(self, sdf_content, query_name, molecule_id):
        """
        Guarda una molécula con organización automática
        
        INPUT:
        - sdf_content (str): Contenido del archivo SDF
        - query_name (str): Nombre original de búsqueda
        - molecule_id (str): ID de la molécula
        
        OUTPUT:
        - bool: True si se guardó exitosamente
        """
        try:
            # 1. Extraer propiedades
            properties = self.extract_properties(sdf_content)
            
            # 2. Generar nombre de archivo
            filename = self.generate_filename(query_name, molecule_id, properties)
            
            # 3. Determinar categoría y crear directorio
            category = self.get_category(properties)
            category_dir = self.output_dir / category
            category_dir.mkdir(exist_ok=True)
            
            # 4. Guardar archivo
            file_path = category_dir / filename
            
            # Evitar sobreescribir archivos existentes
            if file_path.exists():
                counter = 1
                base_name = file_path.stem
                extension = file_path.suffix
                while file_path.exists():
                    new_name = f"{base_name}_{counter}{extension}"
                    file_path = category_dir / new_name
                    counter += 1
            
            # Guardar el archivo
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(sdf_content)
            
            print(f"✅ Archivo guardado: {file_path}")
            return True
            
        except Exception as e:
            print(f"❌ Error guardando archivo: {e}")
            return False
    
    def extract_properties(self, sdf_content):
        """
        Extrae propiedades químicas de la molécula
        
        INPUT:
        - sdf_content (str): Contenido SDF
        
        OUTPUT:
        - dict: Propiedades extraídas
        """
        properties = {
            "name": "unknown",
            "formula": "unknown", 
            "weight": 0,
            "atoms": 0
        }
        
        try:
            mol = Chem.MolFromMolBlock(sdf_content)
            if mol:
                # Fórmula molecular
                properties["formula"] = rdMolDescriptors.CalcMolFormula(mol)
                
                # Peso molecular
                properties["weight"] = round(Descriptors.MolWt(mol), 2)
                
                # Número de átomos
                properties["atoms"] = mol.GetNumAtoms()
                
                # Intentar extraer nombre desde propiedades del SDF
                if mol.HasProp("PUBCHEM_IUPAC_NAME"):
                    properties["name"] = mol.GetProp("PUBCHEM_IUPAC_NAME")
                elif mol.HasProp("PUBCHEM_COMPOUND_CID"):
                    properties["name"] = f"CID_{mol.GetProp('PUBCHEM_COMPOUND_CID')}"
                
        except Exception as e:
            print(f"⚠️ Error extrayendo propiedades: {e}")
        
        return properties
    
    def generate_filename(self, query_name, molecule_id, properties):
        """
        Genera un nombre de archivo apropiado
        
        INPUT:
        - query_name (str): Nombre de búsqueda original
        - molecule_id (str): ID de la molécula  
        - properties (dict): Propiedades de la molécula
        
        OUTPUT:
        - str: Nombre de archivo sanitizado
        """
        # Usar el nombre IUPAC si está disponible, sino el query original
        if properties["name"] != "unknown" and len(properties["name"]) > 3:
            base_name = properties["name"]
        else:
            base_name = query_name
        
        # Añadir fórmula si está disponible
        if properties["formula"] != "unknown":
            base_name += f"_{properties['formula']}"
        
        # Añadir ID como respaldo
        base_name += f"_ID{molecule_id}"
        
        # Sanitizar nombre
        filename = self.sanitize_filename(base_name) + ".sdf"
        
        return filename
    
    def sanitize_filename(self, name):
        """
        Sanitiza un nombre para uso como archivo
        
        INPUT:
        - name (str): Nombre original
        
        OUTPUT:
        - str: Nombre sanitizado
        """
        # Reemplazar caracteres problemáticos
        sanitized = re.sub(r'[\\/*?:"<>|,\[\]()]', "_", name)
        
        # Reemplazar espacios múltiples
        sanitized = re.sub(r'\s+', '_', sanitized)
        
        # Remover caracteres especiales adicionales
        sanitized = re.sub(r'[^\w\-_.]', '', sanitized)
        
        # Limitar longitud
        if len(sanitized) > 60:
            sanitized = sanitized[:57] + "..."
        
        # Asegurar que no esté vacío
        if not sanitized or sanitized == "...":
            sanitized = "molecule"
        
        return sanitized
    
    def get_category(self, properties):
        """
        Determina la categoría para organización
        
        INPUT:
        - properties (dict): Propiedades de la molécula
        
        OUTPUT:
        - str: Nombre de la categoría
        """
        if self.categorize_by == "none":
            return "all_molecules"
        
        elif self.categorize_by == "formula":
            if properties["formula"] != "unknown":
                # Categorizar por primer elemento
                match = re.match(r'^([A-Z][a-z]?)', properties["formula"])
                if match:
                    element = match.group(1)
                    return f"element_{element}"
                else:
                    return "formula_other"
            return "formula_unknown"
        
        elif self.categorize_by == "weight":
            weight = properties["weight"]
            if weight == 0:
                return "weight_unknown"
            elif weight < 100:
                return "weight_0-100"
            elif weight < 200:
                return "weight_100-200"
            elif weight < 300:
                return "weight_200-300"
            elif weight < 500:
                return "weight_300-500"
            else:
                return "weight_500+"
        
        elif self.categorize_by == "atoms":
            atoms = properties["atoms"]
            if atoms == 0:
                return "atoms_unknown"
            elif atoms < 10:
                return "atoms_0-10"
            elif atoms < 20:
                return "atoms_10-20"
            elif atoms < 50:
                return "atoms_20-50"
            else:
                return "atoms_50+"
        
        else:
            return "misc"
    
    def get_statistics(self):
        """
        Obtiene estadísticas del almacenamiento actual
        
        OUTPUT:
        - dict: Estadísticas de archivos guardados
        """
        stats = {
            "total_files": 0,
            "categories": {},
            "total_size_mb": 0
        }
        
        try:
            for sdf_file in self.output_dir.rglob("*.sdf"):
                stats["total_files"] += 1
                
                # Categoría
                category = sdf_file.parent.name
                if category not in stats["categories"]:
                    stats["categories"][category] = 0
                stats["categories"][category] += 1
                
                # Tamaño
                stats["total_size_mb"] += sdf_file.stat().st_size / (1024 * 1024)
            
            stats["total_size_mb"] = round(stats["total_size_mb"], 2)
            
        except Exception as e:
            print(f"Error calculando estadísticas: {e}")
        
        return stats
    
    def cleanup_empty_directories(self):
        """
        Limpia directorios vacíos en el almacenamiento
        
        OUTPUT:
        - int: Número de directorios eliminados
        """
        removed_count = 0
        
        try:
            for category_dir in self.output_dir.iterdir():
                if category_dir.is_dir():
                    # Verificar si está vacío
                    if not any(category_dir.iterdir()):
                        category_dir.rmdir()
                        removed_count += 1
                        print(f"🗑️ Directorio vacío eliminado: {category_dir.name}")
        
        except Exception as e:
            print(f"Error limpiando directorios: {e}")
        
        return removed_count
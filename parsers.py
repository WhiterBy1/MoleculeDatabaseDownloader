# parsers.py
"""
Parsers para diferentes tipos de archivos
INPUT: file_path (str), options (dict)
OUTPUT: list of queries (list)
"""

import pandas as pd
from pathlib import Path
from config import EXCEL_FORMATS, TEXT_FORMATS

class FileParser:
    """Parser unificado para diferentes tipos de archivos"""
    
    def parse_file(self, file_path, **options):
        """
        Parser principal que detecta el tipo de archivo y lo procesa
        
        INPUT:
        - file_path (str): Ruta al archivo
        - options (dict): Opciones específicas del parser
        
        OUTPUT:
        - list: Lista de consultas de moléculas extraídas
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"Archivo no encontrado: {file_path}")
        
        extension = file_path.suffix.lower()
        
        if extension in EXCEL_FORMATS:
            return self.parse_excel(file_path, **options)
        elif extension in TEXT_FORMATS:
            return self.parse_text(file_path, **options)
        else:
            raise ValueError(f"Formato de archivo no soportado: {extension}")
    
    def parse_excel(self, file_path, column_name=None, sheet_name=0):
        """
        Parser para archivos Excel
        
        INPUT:
        - file_path (Path): Ruta al archivo Excel
        - column_name (str): Nombre de la columna con las moléculas
        - sheet_name (str/int): Nombre o índice de la hoja
        
        OUTPUT:
        - list: Lista de nombres de moléculas
        """
        try:
            # Leer el archivo Excel
            df = pd.read_excel(file_path, sheet_name=sheet_name)
            
            # Si no se especifica columna, usar la primera
            if column_name is None:
                column_name = df.columns[0]
            
            if column_name not in df.columns:
                raise ValueError(f"Columna '{column_name}' no encontrada en el archivo")
            
            # Extraer valores y limpiar
            molecules = df[column_name].dropna().astype(str).str.strip()
            molecules = molecules[molecules != ''].tolist()
            
            return molecules
            
        except Exception as e:
            raise Exception(f"Error al leer archivo Excel: {e}")
    
    def parse_text(self, file_path, delimiter=None):
        """
        Parser para archivos de texto y CSV
        
        INPUT:
        - file_path (Path): Ruta al archivo
        - delimiter (str): Delimitador para CSV (None para auto-detectar)
        
        OUTPUT:
        - list: Lista de nombres de moléculas
        """
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read().strip()
            
            if delimiter:
                # Es un CSV con delimitador específico
                lines = content.split('\n')
                molecules = []
                for line in lines:
                    if line.strip():
                        # Tomar el primer elemento de cada línea
                        molecules.append(line.split(delimiter)[0].strip())
            else:
                # Es un archivo de texto simple (una molécula por línea)
                molecules = [line.strip() for line in content.split('\n') if line.strip()]
            
            return molecules
            
        except Exception as e:
            raise Exception(f"Error al leer archivo de texto: {e}")
    
    def get_excel_columns(self, file_path, sheet_name=0):
        """
        Obtiene las columnas disponibles en un archivo Excel
        
        INPUT:
        - file_path (str): Ruta al archivo Excel
        - sheet_name (str/int): Nombre o índice de la hoja
        
        OUTPUT:
        - list: Lista de nombres de columnas
        """
        try:
            df = pd.read_excel(file_path, sheet_name=sheet_name, nrows=0)
            return df.columns.tolist()
        except Exception as e:
            raise Exception(f"Error al leer columnas del Excel: {e}")
    
    def get_excel_sheets(self, file_path):
        """
        Obtiene las hojas disponibles en un archivo Excel
        
        INPUT:
        - file_path (str): Ruta al archivo Excel
        
        OUTPUT:
        - list: Lista de nombres de hojas
        """
        try:
            with pd.ExcelFile(file_path) as xls:
                return xls.sheet_names
        except Exception as e:
            raise Exception(f"Error al leer hojas del Excel: {e}")
# progress.py
"""
Manejo de progreso y recuperación de sesiones
INPUT: session_data (dict)
OUTPUT: save/load state operations
"""

import json
import os
from datetime import datetime
from pathlib import Path
from config import PROGRESS_DIR

class ProgressManager:
    """Gestor de progreso y recuperación de sesiones"""
    
    def __init__(self):
        self.progress_dir = Path(PROGRESS_DIR)
        self.progress_dir.mkdir(exist_ok=True)
    
    def create_session(self, session_name, total_items, metadata=None):
        """
        Crea una nueva sesión de descarga
        
        INPUT:
        - session_name (str): Nombre único de la sesión
        - total_items (int): Total de elementos a procesar
        - metadata (dict): Información adicional (archivo origen, configuración, etc.)
        
        OUTPUT:
        - str: ID de la sesión creada
        """
        session_id = self._generate_session_id(session_name)
        
        session_data = {
            "session_id": session_id,
            "session_name": session_name,
            "created_at": datetime.now().isoformat(),
            "total_items": total_items,
            "completed_items": 0,
            "failed_items": 0,
            "completed_queries": [],
            "failed_queries": [],
            "pending_queries": [],
            "metadata": metadata or {},
            "status": "active"
        }
        
        self._save_session(session_data)
        return session_id
    
    def load_session(self, session_id):
        """
        Carga una sesión existente
        
        INPUT:
        - session_id (str): ID de la sesión
        
        OUTPUT:
        - dict: Datos de la sesión o None si no existe
        """
        session_file = self.progress_dir / f"{session_id}.json"
        
        if not session_file.exists():
            return None
        
        try:
            with open(session_file, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception as e:
            print(f"Error al cargar sesión: {e}")
            return None
    
    def update_progress(self, session_id, completed_query=None, failed_query=None):
        """
        Actualiza el progreso de una sesión
        
        INPUT:
        - session_id (str): ID de la sesión
        - completed_query (str): Query completada exitosamente
        - failed_query (str): Query que falló
        
        OUTPUT:
        - bool: True si se actualizó correctamente
        """
        session_data = self.load_session(session_id)
        if not session_data:
            return False
        
        if completed_query:
            if completed_query not in session_data["completed_queries"]:
                session_data["completed_queries"].append(completed_query)
                session_data["completed_items"] = len(session_data["completed_queries"])
        
        if failed_query:
            if failed_query not in session_data["failed_queries"]:
                session_data["failed_queries"].append(failed_query)
                session_data["failed_items"] = len(session_data["failed_queries"])
        
        # Actualizar timestamp
        session_data["updated_at"] = datetime.now().isoformat()
        
        self._save_session(session_data)
        return True
    
    def get_pending_queries(self, session_id, all_queries):
        """
        Obtiene las queries pendientes de una sesión
        
        INPUT:
        - session_id (str): ID de la sesión
        - all_queries (list): Lista completa de queries originales
        
        OUTPUT:
        - list: Queries que aún no se han procesado
        """
        session_data = self.load_session(session_id)
        if not session_data:
            return all_queries
        
        processed = set(session_data["completed_queries"] + session_data["failed_queries"])
        pending = [q for q in all_queries if q not in processed]
        
        return pending
    
    def complete_session(self, session_id):
        """
        Marca una sesión como completada
        
        INPUT:
        - session_id (str): ID de la sesión
        
        OUTPUT:
        - bool: True si se completó correctamente
        """
        session_data = self.load_session(session_id)
        if not session_data:
            return False
        
        session_data["status"] = "completed"
        session_data["completed_at"] = datetime.now().isoformat()
        
        self._save_session(session_data)
        return True
    
    def list_sessions(self, status_filter=None):
        """
        Lista todas las sesiones disponibles
        
        INPUT:
        - status_filter (str): Filtro por estado ("active", "completed", "failed")
        
        OUTPUT:
        - list: Lista de sesiones con información básica
        """
        sessions = []
        
        for session_file in self.progress_dir.glob("*.json"):
            try:
                with open(session_file, 'r', encoding='utf-8') as f:
                    session_data = json.load(f)
                
                if status_filter is None or session_data.get("status") == status_filter:
                    sessions.append({
                        "session_id": session_data["session_id"],
                        "session_name": session_data["session_name"],
                        "created_at": session_data["created_at"],
                        "status": session_data.get("status", "unknown"),
                        "progress": f"{session_data['completed_items']}/{session_data['total_items']}"
                    })
            except Exception:
                continue
        
        return sorted(sessions, key=lambda x: x["created_at"], reverse=True)
    
    def delete_session(self, session_id):
        """
        Elimina una sesión
        
        INPUT:
        - session_id (str): ID de la sesión
        
        OUTPUT:
        - bool: True si se eliminó correctamente
        """
        session_file = self.progress_dir / f"{session_id}.json"
        
        try:
            if session_file.exists():
                session_file.unlink()
                return True
        except Exception as e:
            print(f"Error al eliminar sesión: {e}")
        
        return False
    
    def _generate_session_id(self, session_name):
        """Genera un ID único para la sesión"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        clean_name = "".join(c for c in session_name if c.isalnum() or c in "._-")[:20]
        return f"{clean_name}_{timestamp}"
    
    def _save_session(self, session_data):
        """Guarda los datos de la sesión"""
        session_file = self.progress_dir / f"{session_data['session_id']}.json"
        
        with open(session_file, 'w', encoding='utf-8') as f:
            json.dump(session_data, f, indent=2, ensure_ascii=False)
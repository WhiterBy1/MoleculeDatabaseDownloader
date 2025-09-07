# main.py
"""
Punto de entrada principal del sistema
Orquesta todos los módulos y maneja la GUI
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import threading
from pathlib import Path

# Importar módulos locales
from sources import SourceManager
from parsers import FileParser
from downloader import MoleculeDownloader
from progress import ProgressManager
from storage import StorageManager
from config import DEFAULT_OUTPUT_DIR, CATEGORIZE_OPTIONS

class MoleculeDownloaderApp:
    """Aplicación principal con GUI mejorada"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("Descargador de Moléculas 3D - Versión Modular")
        self.root.geometry("900x800")
        
        # Inicializar gestores
        self.source_manager = SourceManager()
        self.file_parser = FileParser()
        self.downloader = MoleculeDownloader(self.source_manager)
        self.progress_manager = ProgressManager()
        self.storage_manager = StorageManager()
        
        # Variables de configuración
        self.output_dir = tk.StringVar(value=DEFAULT_OUTPUT_DIR)
        self.categorize_by = tk.StringVar(value="formula")
        self.source = tk.StringVar(value="pubchem")
        
        # Variables para Excel
        self.excel_file = tk.StringVar()
        self.selected_sheet = tk.StringVar()
        self.selected_column = tk.StringVar()
        
        # Estado de la aplicación
        self.current_session_id = None
        self.is_processing = False
        
        self.create_gui()
    
    def create_gui(self):
        """Crea la interfaz gráfica modular"""
        
        # Notebook principal
        notebook = ttk.Notebook(self.root)
        notebook.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Tab 1: Configuración
        self.create_config_tab(notebook)
        
        # Tab 2: Descarga Individual
        self.create_single_tab(notebook)
        
        # Tab 3: Descarga por Lotes (Texto)
        self.create_batch_tab(notebook)
        
        # Tab 4: Descarga desde Excel (NUEVA)
        self.create_excel_tab(notebook)
        
        # Tab 5: Gestión de Sesiones (NUEVA)
        self.create_sessions_tab(notebook)
        
        # Tab 6: Monitor
        self.create_monitor_tab(notebook)
    
    def create_config_tab(self, notebook):
        """Tab de configuración general"""
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="⚙️ Configuración")
        
        # Directorio de salida
        ttk.Label(frame, text="Directorio de salida:").grid(row=0, column=0, sticky='w', padx=5, pady=5)
        ttk.Entry(frame, textvariable=self.output_dir, width=50).grid(row=0, column=1, padx=5, pady=5)
        ttk.Button(frame, text="📁", command=self.select_output_dir).grid(row=0, column=2, padx=5, pady=5)
        
        # Categorización
        ttk.Label(frame, text="Categorizar por:").grid(row=1, column=0, sticky='w', padx=5, pady=5)
        ttk.Combobox(frame, textvariable=self.categorize_by, values=CATEGORIZE_OPTIONS, 
                    state="readonly").grid(row=1, column=1, sticky='w', padx=5, pady=5)
        
        # Fuente
        ttk.Label(frame, text="Fuente de datos:").grid(row=2, column=0, sticky='w', padx=5, pady=5)
        ttk.Combobox(frame, textvariable=self.source, values=["pubchem", "chemspider"], 
                    state="readonly").grid(row=2, column=1, sticky='w', padx=5, pady=5)
    
    def create_single_tab(self, notebook):
        """Tab para descarga individual"""
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="🔬 Individual")
        
        ttk.Label(frame, text="Nombre de la molécula:").grid(row=0, column=0, sticky='w', padx=5, pady=5)
        self.single_query = tk.StringVar()
        ttk.Entry(frame, textvariable=self.single_query, width=50).grid(row=0, column=1, padx=5, pady=5)
        ttk.Button(frame, text="⬇️ Descargar", command=self.download_single).grid(row=0, column=2, padx=5, pady=5)
    
    def create_batch_tab(self, notebook):
        """Tab para descarga por lotes"""
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="📋 Por Lotes")
        
        ttk.Label(frame, text="Lista de moléculas (una por línea):").pack(anchor='w', padx=5, pady=5)
        
        self.batch_text = tk.Text(frame, height=15, width=70)
        self.batch_text.pack(fill='both', expand=True, padx=5, pady=5)
        
        ttk.Button(frame, text="▶️ Procesar Lista", command=self.download_from_text).pack(pady=10)
    
    def create_excel_tab(self, notebook):
        """Tab para descarga desde Excel (NUEVA FUNCIONALIDAD)"""
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="📊 Desde Excel")
        
        # Selección de archivo
        ttk.Label(frame, text="Archivo Excel:").grid(row=0, column=0, sticky='w', padx=5, pady=5)
        ttk.Entry(frame, textvariable=self.excel_file, width=50).grid(row=0, column=1, padx=5, pady=5)
        ttk.Button(frame, text="📁", command=self.select_excel_file).grid(row=0, column=2, padx=5, pady=5)
        
        # Selección de hoja
        ttk.Label(frame, text="Hoja:").grid(row=1, column=0, sticky='w', padx=5, pady=5)
        self.sheet_combo = ttk.Combobox(frame, textvariable=self.selected_sheet, state="readonly")
        self.sheet_combo.grid(row=1, column=1, sticky='w', padx=5, pady=5)
        self.sheet_combo.bind('<<ComboboxSelected>>', self.on_sheet_selected)
        
        # Selección de columna
        ttk.Label(frame, text="Columna con moléculas:").grid(row=2, column=0, sticky='w', padx=5, pady=5)
        self.column_combo = ttk.Combobox(frame, textvariable=self.selected_column, state="readonly")
        self.column_combo.grid(row=2, column=1, sticky='w', padx=5, pady=5)
        
        # Vista previa
        preview_frame = ttk.LabelFrame(frame, text="Vista previa")
        preview_frame.grid(row=3, column=0, columnspan=3, sticky='ew', padx=5, pady=10)
        
        self.preview_text = tk.Text(preview_frame, height=8, width=70)
        self.preview_text.pack(padx=5, pady=5)
        
        # Botones
        button_frame = ttk.Frame(frame)
        button_frame.grid(row=4, column=0, columnspan=3, pady=10)
        
        ttk.Button(button_frame, text="🔄 Cargar Archivo", 
                  command=self.load_excel_preview).pack(side='left', padx=5)
        ttk.Button(button_frame, text="▶️ Iniciar Descarga", 
                  command=self.start_excel_download).pack(side='left', padx=5)
    
    def create_sessions_tab(self, notebook):
        """Tab para gestión de sesiones (NUEVA FUNCIONALIDAD)"""
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="📋 Sesiones")
        
        # Lista de sesiones
        ttk.Label(frame, text="Sesiones guardadas:").pack(anchor='w', padx=5, pady=5)
        
        # Treeview para mostrar sesiones
        columns = ('name', 'status', 'progress', 'created')
        self.sessions_tree = ttk.Treeview(frame, columns=columns, show='headings', height=10)
        
        self.sessions_tree.heading('name', text='Nombre')
        self.sessions_tree.heading('status', text='Estado')
        self.sessions_tree.heading('progress', text='Progreso')
        self.sessions_tree.heading('created', text='Creado')
        
        self.sessions_tree.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Botones de gestión
        button_frame = ttk.Frame(frame)
        button_frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Button(button_frame, text="🔄 Actualizar", 
                  command=self.refresh_sessions).pack(side='left', padx=5)
        ttk.Button(button_frame, text="▶️ Continuar Seleccionada", 
                  command=self.resume_selected_session).pack(side='left', padx=5)
        ttk.Button(button_frame, text="🗑️ Eliminar Seleccionada", 
                  command=self.delete_selected_session).pack(side='right', padx=5)
        
        # Cargar sesiones al inicio
        self.refresh_sessions()
    
    def create_monitor_tab(self, notebook):
        """Tab de monitoreo"""
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="📊 Monitor")
        
        # Barra de progreso
        ttk.Label(frame, text="Progreso actual:").pack(anchor='w', padx=5, pady=5)
        self.progress = ttk.Progressbar(frame, mode='determinate')
        self.progress.pack(fill='x', padx=5, pady=5)
        
        self.progress_label = ttk.Label(frame, text="Listo")
        self.progress_label.pack(anchor='w', padx=5, pady=5)
        
        # Log
        ttk.Label(frame, text="Log de actividad:").pack(anchor='w', padx=5, pady=(20,5))
        self.log_text = tk.Text(frame, height=20, width=80)
        self.log_text.pack(fill='both', expand=True, padx=5, pady=5)
    
    # =========================================================================
    # MÉTODOS PARA EXCEL
    # =========================================================================
    
    def select_excel_file(self):
        """Selecciona archivo Excel"""
        file_path = filedialog.askopenfilename(
            title="Seleccionar archivo Excel",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )
        if file_path:
            self.excel_file.set(file_path)
            self.load_excel_structure()
    
    def load_excel_structure(self):
        """Carga la estructura del Excel (hojas y columnas)"""
        try:
            file_path = self.excel_file.get()
            if not file_path:
                return
            
            # Cargar hojas
            sheets = self.file_parser.get_excel_sheets(file_path)
            self.sheet_combo['values'] = sheets
            if sheets:
                self.selected_sheet.set(sheets[0])
                self.on_sheet_selected()
                
        except Exception as e:
            messagebox.showerror("Error", f"Error al cargar archivo Excel: {e}")
    
    def on_sheet_selected(self, event=None):
        """Se ejecuta cuando se selecciona una hoja"""
        try:
            file_path = self.excel_file.get()
            sheet_name = self.selected_sheet.get()
            
            if file_path and sheet_name:
                columns = self.file_parser.get_excel_columns(file_path, sheet_name)
                self.column_combo['values'] = columns
                if columns:
                    self.selected_column.set(columns[0])
                    
        except Exception as e:
            messagebox.showerror("Error", f"Error al cargar columnas: {e}")
    
    def load_excel_preview(self):
        """Carga vista previa de los datos del Excel"""
        try:
            file_path = self.excel_file.get()
            sheet_name = self.selected_sheet.get()
            column_name = self.selected_column.get()
            
            if not all([file_path, sheet_name, column_name]):
                messagebox.showwarning("Advertencia", "Seleccione archivo, hoja y columna")
                return
            
            # Obtener datos de la columna
            molecules = self.file_parser.parse_excel(
                file_path, 
                column_name=column_name, 
                sheet_name=sheet_name
            )
            
            # Mostrar vista previa (primeros 20)
            preview_text = f"Total de moléculas encontradas: {len(molecules)}\n\n"
            preview_text += "Vista previa (primeras 20):\n"
            preview_text += "\n".join(molecules[:20])
            
            if len(molecules) > 20:
                preview_text += f"\n... y {len(molecules) - 20} más"
            
            self.preview_text.delete(1.0, tk.END)
            self.preview_text.insert(1.0, preview_text)
            
        except Exception as e:
            messagebox.showerror("Error", f"Error al cargar vista previa: {e}")
    
    def start_excel_download(self):
        """Inicia descarga desde Excel"""
        if self.is_processing:
            messagebox.showwarning("Advertencia", "Ya hay un proceso en ejecución")
            return
        
        try:
            file_path = self.excel_file.get()
            sheet_name = self.selected_sheet.get()
            column_name = self.selected_column.get()
            
            if not all([file_path, sheet_name, column_name]):
                messagebox.showwarning("Advertencia", "Complete todos los campos")
                return
            
            # Extraer moléculas
            molecules = self.file_parser.parse_excel(
                file_path, 
                column_name=column_name, 
                sheet_name=sheet_name
            )
            
            if not molecules:
                messagebox.showwarning("Advertencia", "No se encontraron moléculas válidas")
                return
            
            # Crear sesión de progreso
            session_name = f"Excel_{Path(file_path).stem}_{column_name}"
            self.current_session_id = self.progress_manager.create_session(
                session_name=session_name,
                total_items=len(molecules),
                metadata={
                    "source_file": file_path,
                    "sheet_name": sheet_name,
                    "column_name": column_name,
                    "source_type": "excel"
                }
            )
            
            # Iniciar descarga en thread separado
            self.is_processing = True
            thread = threading.Thread(
                target=self._download_batch_thread, 
                args=(molecules, self.current_session_id)
            )
            thread.daemon = True
            thread.start()
            
            messagebox.showinfo("Inicio", f"Descarga iniciada para {len(molecules)} moléculas")
            
        except Exception as e:
            messagebox.showerror("Error", f"Error al iniciar descarga: {e}")
    
    # =========================================================================
    # MÉTODOS PARA SESIONES
    # =========================================================================
    
    def refresh_sessions(self):
        """Actualiza la lista de sesiones"""
        try:
            # Limpiar árbol
            for item in self.sessions_tree.get_children():
                self.sessions_tree.delete(item)
            
            # Cargar sesiones
            sessions = self.progress_manager.list_sessions()
            
            for session in sessions:
                self.sessions_tree.insert('', 'end', values=(
                    session['session_name'],
                    session['status'],
                    session['progress'],
                    session['created_at'][:19]  # Solo fecha y hora
                ))
                
        except Exception as e:
            messagebox.showerror("Error", f"Error al cargar sesiones: {e}")
    
    def resume_selected_session(self):
        """Continúa una sesión seleccionada"""
        if self.is_processing:
            messagebox.showwarning("Advertencia", "Ya hay un proceso en ejecución")
            return
        
        selected = self.sessions_tree.selection()
        if not selected:
            messagebox.showwarning("Advertencia", "Seleccione una sesión")
            return
        
        try:
            # Obtener datos de la sesión
            item = self.sessions_tree.item(selected[0])
            session_name = item['values'][0]
            
            # Buscar sesión por nombre
            sessions = self.progress_manager.list_sessions()
            target_session = None
            for session in sessions:
                if session['session_name'] == session_name:
                    target_session = session
                    break
            
            if not target_session:
                messagebox.showerror("Error", "No se pudo encontrar la sesión")
                return
            
            session_id = target_session['session_id']
            session_data = self.progress_manager.load_session(session_id)
            
            if session_data['status'] == 'completed':
                messagebox.showinfo("Info", "Esta sesión ya está completada")
                return
            
            # Reconstruir lista de moléculas según el tipo de fuente
            molecules = self._reconstruct_molecule_list(session_data)
            
            if not molecules:
                messagebox.showerror("Error", "No se pudo reconstruir la lista de moléculas")
                return
            
            # Obtener moléculas pendientes
            pending_molecules = self.progress_manager.get_pending_queries(session_id, molecules)
            
            if not pending_molecules:
                messagebox.showinfo("Info", "No hay moléculas pendientes en esta sesión")
                return
            
            # Continuar descarga
            self.current_session_id = session_id
            self.is_processing = True
            
            thread = threading.Thread(
                target=self._download_batch_thread,
                args=(pending_molecules, session_id)
            )
            thread.daemon = True
            thread.start()
            
            messagebox.showinfo("Reanudado", f"Continuando con {len(pending_molecules)} moléculas pendientes")
            
        except Exception as e:
            messagebox.showerror("Error", f"Error al reanudar sesión: {e}")
    
    def delete_selected_session(self):
        """Elimina una sesión seleccionada"""
        selected = self.sessions_tree.selection()
        if not selected:
            messagebox.showwarning("Advertencia", "Seleccione una sesión")
            return
        
        if messagebox.askyesno("Confirmar", "¿Está seguro de eliminar esta sesión?"):
            try:
                item = self.sessions_tree.item(selected[0])
                session_name = item['values'][0]
                
                # Buscar y eliminar sesión
                sessions = self.progress_manager.list_sessions()
                for session in sessions:
                    if session['session_name'] == session_name:
                        self.progress_manager.delete_session(session['session_id'])
                        break
                
                self.refresh_sessions()
                messagebox.showinfo("Eliminado", "Sesión eliminada correctamente")
                
            except Exception as e:
                messagebox.showerror("Error", f"Error al eliminar sesión: {e}")
    
    def _reconstruct_molecule_list(self, session_data):
        """Reconstruye la lista de moléculas desde los metadatos de la sesión"""
        try:
            metadata = session_data.get('metadata', {})
            source_type = metadata.get('source_type', 'unknown')
            
            if source_type == 'excel':
                # Reconstruir desde Excel
                return self.file_parser.parse_excel(
                    metadata['source_file'],
                    column_name=metadata['column_name'],
                    sheet_name=metadata['sheet_name']
                )
            elif source_type == 'text':
                # Reconstruir desde archivo de texto
                return self.file_parser.parse_text(metadata['source_file'])
            else:
                # Para listas manuales, usar las queries registradas
                completed = session_data.get('completed_queries', [])
                failed = session_data.get('failed_queries', [])
                return completed + failed
                
        except Exception as e:
            print(f"Error reconstruyendo lista: {e}")
            return []
    
    # =========================================================================
    # MÉTODOS PARA DESCARGAS
    # =========================================================================
    
    def download_single(self):
        """Descarga una molécula individual"""
        query = self.single_query.get().strip()
        if not query:
            messagebox.showerror("Error", "Ingrese el nombre de una molécula")
            return
        
        if self.is_processing:
            messagebox.showwarning("Advertencia", "Ya hay un proceso en ejecución")
            return
        
        self.is_processing = True
        thread = threading.Thread(target=self._download_single_thread, args=(query,))
        thread.daemon = True
        thread.start()
    
    def download_from_text(self):
        """Descarga desde el área de texto"""
        text_content = self.batch_text.get(1.0, tk.END).strip()
        if not text_content:
            messagebox.showerror("Error", "Ingrese algunas moléculas")
            return
        
        queries = [line.strip() for line in text_content.split('\n') if line.strip()]
        if not queries:
            messagebox.showerror("Error", "No se encontraron consultas válidas")
            return
        
        if self.is_processing:
            messagebox.showwarning("Advertencia", "Ya hay un proceso en ejecución")
            return
        
        # Crear sesión
        session_name = f"Manual_List_{len(queries)}_items"
        session_id = self.progress_manager.create_session(
            session_name=session_name,
            total_items=len(queries),
            metadata={"source_type": "manual"}
        )
        
        self.current_session_id = session_id
        self.is_processing = True
        
        thread = threading.Thread(target=self._download_batch_thread, args=(queries, session_id))
        thread.daemon = True
        thread.start()
    
    def _download_single_thread(self, query):
        """Thread para descarga individual"""
        try:
            self.log(f"Iniciando descarga de: {query}")
            
            # Configurar descargador
            self.downloader.configure(
                output_dir=self.output_dir.get(),
                categorize_by=self.categorize_by.get()
            )
            
            success = self.downloader.download_molecule(query, self.source.get())
            
            if success:
                self.log("✓ Descarga completada")
                messagebox.showinfo("Éxito", f"Molécula '{query}' descargada")
            else:
                self.log("✗ Error en la descarga")
                messagebox.showerror("Error", f"No se pudo descargar '{query}'")
                
        except Exception as e:
            self.log(f"Error: {e}")
            messagebox.showerror("Error", str(e))
        finally:
            self.is_processing = False
    
    def _download_batch_thread(self, queries, session_id):
        """Thread para descarga por lotes con manejo de sesiones"""
        try:
            self.log(f"Iniciando descarga por lotes: {len(queries)} moléculas")
            
            # Configurar descargador
            self.downloader.configure(
                output_dir=self.output_dir.get(),
                categorize_by=self.categorize_by.get(),
                progress_callback=self.update_progress,
                log_callback=self.log
            )
            
            completed = 0
            failed = 0
            
            for i, query in enumerate(queries):
                try:
                    self.log(f"Procesando: {query} ({i+1}/{len(queries)})")
                    self.update_progress(i, len(queries))
                    
                    success = self.downloader.download_molecule(query, self.source.get())
                    
                    if success:
                        completed += 1
                        self.progress_manager.update_progress(session_id, completed_query=query)
                        self.log(f"✓ {query} completada")
                    else:
                        failed += 1
                        self.progress_manager.update_progress(session_id, failed_query=query)
                        self.log(f"✗ {query} falló")
                        
                except Exception as e:
                    failed += 1
                    self.progress_manager.update_progress(session_id, failed_query=query)
                    self.log(f"✗ Error con {query}: {e}")
            
            # Completar sesión
            self.progress_manager.complete_session(session_id)
            self.update_progress(len(queries), len(queries))
            
            message = f"Proceso completado:\n✓ Exitosas: {completed}\n✗ Fallidas: {failed}\n📊 Total: {len(queries)}"
            self.log(message)
            messagebox.showinfo("Completado", message)
            
        except Exception as e:
            self.log(f"Error crítico: {e}")
            messagebox.showerror("Error", str(e))
        finally:
            self.is_processing = False
            self.current_session_id = None
    
    def select_output_dir(self):
        """Selecciona directorio de salida"""
        directory = filedialog.askdirectory(initialdir=self.output_dir.get())
        if directory:
            self.output_dir.set(directory)
    
    def update_progress(self, current, total):
        """Actualiza barra de progreso"""
        if total > 0:
            progress_percent = (current / total) * 100
            self.progress['value'] = progress_percent
            self.progress_label.config(text=f"Procesando: {current}/{total} ({progress_percent:.1f}%)")
        else:
            self.progress['value'] = 0
            self.progress_label.config(text="Listo")
    
    def log(self, message):
        """Añade mensaje al log"""
        import time
        timestamp = time.strftime('%H:%M:%S')
        self.log_text.insert(tk.END, f"{timestamp} - {message}\n")
        self.log_text.see(tk.END)
        self.root.update_idletasks()

def main():
    """Función principal"""
    root = tk.Tk()
    app = MoleculeDownloaderApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()

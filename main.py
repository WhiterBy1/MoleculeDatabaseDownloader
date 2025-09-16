# main.py
"""
Punto de entrada principal del sistema
Orquesta todos los módulos y maneja la GUI
"""

import customtkinter as ctk
from PIL import Image
import tkinter as tk
from tkinter import filedialog, messagebox
import threading
from pathlib import Path
import time
import pandas as pd
import os
import subprocess

# Importar módulos locales
from sources import SourceManager
from parsers import FileParser
from downloader import MoleculeDownloader
from progress import ProgressManager
from storage import StorageManager
from config import DEFAULT_OUTPUT_DIR, CATEGORIZE_OPTIONS

# Configurar tema moderno
ctk.set_appearance_mode("dark")  # "system", "dark", "light"
ctk.set_default_color_theme("blue")  # "blue", "green", "dark-blue"

class MoleculeApp:
    """Aplicación moderna con interfaz elegante y funcionalidades completas"""
    
    def __init__(self):
        self.root = ctk.CTk()
        self.root.title("🧪 Molecular Structure Downloader Pro")
        self.root.geometry("1400x900")
        self.root.minsize(1200, 800)
        
        # Configurar icono de la ventana
        try:
            self.root.iconbitmap("assets/molecule_icon.ico")
        except:
            pass
        
        # Inicializar gestores (mantener los existentes)
        self.source_manager = SourceManager()
        self.file_parser = FileParser()
        self.downloader = MoleculeDownloader(self.source_manager)
        self.progress_manager = ProgressManager()
        self.storage_manager = StorageManager()
        
        # Variables de configuración
        self.output_dir = ctk.StringVar(value=DEFAULT_OUTPUT_DIR)
        self.categorize_by = ctk.StringVar(value="formula")
        self.source = ctk.StringVar(value="pubchem")
        
        # Variables para Excel
        self.excel_file = ctk.StringVar()
        self.selected_sheet = ctk.StringVar()
        self.selected_column = ctk.StringVar()
        self.excel_sheets = []
        self.excel_columns = []
        self.excel_preview_data = []
        
        # Estado
        self.current_session_id = None
        self.is_processing = False
        
        # Widgets de referencia
        self.nav_buttons = []
        self.current_view = None
        self.progress_var = ctk.DoubleVar()
        self.status_var = ctk.StringVar(value="Ready")
        
        # Cargar iconos y configurar GUI
        self.load_icons()
        self.setup_gui()
        
        # Configurar callbacks del downloader
        self.configure_downloader()
    
    def load_icons(self):
        """Carga iconos para la interfaz"""
        self.icons = {}
        
        # Iconos por defecto (emoji/texto) - funciona sin archivos externos
        self.icons = {
            'download': "⬇️",
            'settings': "⚙️", 
            'molecules': "🔬",
            'excel': "📊",
            'sessions': "📋",
            'monitor': "📈",
            'folder': "📁",
            'play': "▶️",
            'refresh': "🔄",
            'delete': "🗑️",
            'pause': "⏸️",
            'stop': "⏹️",
            'batch': "📋",
            'upload': "📤",
            'save': "💾"
        }
    
    def setup_gui(self):
        """Configura la interfaz principal"""
        # Frame principal con padding
        main_frame = ctk.CTkFrame(self.root, fg_color="transparent")
        main_frame.pack(fill="both", expand=True, padx=15, pady=15)
        
        # Header elegante
        self.create_header(main_frame)
        
        # Navegación lateral y contenido principal
        content_frame = ctk.CTkFrame(main_frame, fg_color="transparent")
        content_frame.pack(fill="both", expand=True, pady=(15, 0))
        
        # Sidebar de navegación
        self.create_sidebar(content_frame)
        
        # Área de contenido principal
        self.create_main_content(content_frame)
        
        # Footer con estadísticas y progreso
        self.create_footer(main_frame)
    
    def create_header(self, parent):
        """Crea el header elegante"""
        header_frame = ctk.CTkFrame(parent, height=80)
        header_frame.pack(fill="x", pady=(0, 10))
        header_frame.pack_propagate(False)
        
        # Logo y título
        title_frame = ctk.CTkFrame(header_frame, fg_color="transparent")
        title_frame.pack(side="left", fill="y", padx=20, pady=10)
        
        title_label = ctk.CTkLabel(
            title_frame,
            text="🧪 Molecular Downloader Pro",
            font=ctk.CTkFont(size=26, weight="bold")
        )
        title_label.pack(anchor="w")
        
        subtitle_label = ctk.CTkLabel(
            title_frame,
            text="Advanced 3D molecular structure acquisition system",
            font=ctk.CTkFont(size=11),
            text_color=("gray60", "gray40")
        )
        subtitle_label.pack(anchor="w")
        
        # Botones de acción rápida
        actions_frame = ctk.CTkFrame(header_frame, fg_color="transparent")
        actions_frame.pack(side="right", fill="y", padx=20, pady=15)
        
        quick_download_btn = ctk.CTkButton(
            actions_frame,
            text=f"{self.icons['download']} Quick Download",
            command=self.show_quick_download,
            height=35,
            width=140
        )
        quick_download_btn.pack(side="right", padx=(10, 0))
        
        open_folder_btn = ctk.CTkButton(
            actions_frame,
            text=f"{self.icons['folder']} Open Results",
            command=self.open_results_folder,
            height=35,
            width=120,
            fg_color=("gray75", "gray25"),
            hover_color=("gray70", "gray30")
        )
        open_folder_btn.pack(side="right")
    
    def create_sidebar(self, parent):
        """Crea la barra lateral de navegación"""
        sidebar_frame = ctk.CTkFrame(parent, width=280)
        sidebar_frame.pack(side="left", fill="y", padx=(0, 15))
        sidebar_frame.pack_propagate(False)
        
        # Título de navegación
        nav_title = ctk.CTkLabel(
            sidebar_frame,
            text="Navigation",
            font=ctk.CTkFont(size=16, weight="bold")
        )
        nav_title.pack(pady=(20, 15), padx=20)
        
        # Botones de navegación
        self.nav_buttons = []
        nav_items = [
            ("settings", "Configuration", self.show_settings),
            ("molecules", "Single Download", self.show_single_download),
            ("batch", "Batch from Text", self.show_batch_download),
            ("excel", "Batch from Excel", self.show_excel_download),
            ("sessions", "Session Manager", self.show_sessions),
            ("monitor", "Progress Monitor", self.show_monitor)
        ]
        
        for icon_key, text, command in nav_items:
            btn = ctk.CTkButton(
                sidebar_frame,
                text=f"{self.icons[icon_key]} {text}",
                command=command,
                anchor="w",
                height=45,
                width=240,
                fg_color="transparent",
                text_color=("gray40", "gray60"),
                hover_color=("gray90", "gray20")
            )
            btn.pack(pady=3, padx=20, fill="x")
            self.nav_buttons.append(btn)
        
        # Separador
        separator = ctk.CTkFrame(sidebar_frame, height=2)
        separator.pack(fill="x", padx=20, pady=20)
        
        # Estado del sistema
        self.create_system_status(sidebar_frame)
        
        # Activar primer botón
        self.activate_nav_button(0)
    
    def create_system_status(self, parent):
        """Crea el panel de estado del sistema"""
        status_frame = ctk.CTkFrame(parent)
        status_frame.pack(fill="x", padx=20, pady=(0, 20))
        
        status_title = ctk.CTkLabel(
            status_frame,
            text="System Status",
            font=ctk.CTkFont(size=14, weight="bold")
        )
        status_title.pack(pady=(15, 10))
        
        # Status text
        self.status_label = ctk.CTkLabel(
            status_frame,
            textvariable=self.status_var,
            font=ctk.CTkFont(size=11)
        )
        self.status_label.pack(pady=(0, 10))
        
        # Progress bar
        self.progress_bar = ctk.CTkProgressBar(
            status_frame,
            variable=self.progress_var,
            width=200
        )
        self.progress_bar.pack(pady=(0, 15), padx=15)
        
        # Progress percentage
        self.progress_label = ctk.CTkLabel(
            status_frame,
            text="0%",
            font=ctk.CTkFont(size=10)
        )
        self.progress_label.pack(pady=(0, 15))
    
    def create_main_content(self, parent):
        """Crea el área de contenido principal"""
        self.content_frame = ctk.CTkFrame(parent)
        self.content_frame.pack(side="right", fill="both", expand=True)
        
        # Por defecto mostrar configuración
        self.current_view = None
        self.show_settings()
    
    def create_footer(self, parent):
        """Crea el footer con estadísticas"""
        footer_frame = ctk.CTkFrame(parent, height=40)
        footer_frame.pack(fill="x", pady=(10, 0))
        footer_frame.pack_propagate(False)
        
        # Statistics
        self.stats_label = ctk.CTkLabel(
            footer_frame,
            text="Ready • 0 molecules downloaded • 0 sessions active",
            font=ctk.CTkFont(size=11)
        )
        self.stats_label.pack(side="left", padx=20, pady=12)
        
        # Action buttons
        button_frame = ctk.CTkFrame(footer_frame, fg_color="transparent")
        button_frame.pack(side="right", padx=20, pady=8)
        
        self.emergency_stop_btn = ctk.CTkButton(
            button_frame,
            text=f"{self.icons['stop']} Emergency Stop",
            command=self.emergency_stop,
            height=25,
            width=120,
            fg_color="red",
            hover_color="darkred"
        )
        self.emergency_stop_btn.pack(side="right", padx=(10, 0))
        
        # Update stats initially
        self.update_footer_stats()
    
    def activate_nav_button(self, index):
        """Activa un botón de navegación"""
        for i, btn in enumerate(self.nav_buttons):
            if i == index:
                btn.configure(fg_color=("gray80", "gray20"))
            else:
                btn.configure(fg_color="transparent")
    
    def clear_content(self):
        """Limpia el contenido actual"""
        for widget in self.content_frame.winfo_children():
            widget.destroy()
    
    # =========================================================================
    # VISTAS DE CONTENIDO
    # =========================================================================
    
    def show_settings(self):
        """Muestra la vista de configuración"""
        if self.current_view == "settings":
            return
        
        self.activate_nav_button(0)
        self.current_view = "settings"
        self.clear_content()
        
        # Scrollable frame
        scrollable_frame = ctk.CTkScrollableFrame(self.content_frame)
        scrollable_frame.pack(fill="both", expand=True, padx=25, pady=25)
        
        # Título de sección
        title = ctk.CTkLabel(
            scrollable_frame,
            text=f"{self.icons['settings']} Configuration",
            font=ctk.CTkFont(size=22, weight="bold")
        )
        title.pack(pady=(0, 20), anchor="w")
        
        # Tarjeta de directorio de salida
        output_card = ctk.CTkFrame(scrollable_frame)
        output_card.pack(fill="x", pady=(0, 15))
        
        ctk.CTkLabel(
            output_card,
            text="📁 Output Directory",
            font=ctk.CTkFont(size=16, weight="bold")
        ).pack(anchor="w", padx=20, pady=(20, 10))
        
        dir_frame = ctk.CTkFrame(output_card, fg_color="transparent")
        dir_frame.pack(fill="x", padx=20, pady=(0, 20))
        
        self.output_entry = ctk.CTkEntry(
            dir_frame,
            textvariable=self.output_dir,
            height=35,
            font=ctk.CTkFont(size=12)
        )
        self.output_entry.pack(side="left", fill="x", expand=True, padx=(0, 10))
        
        browse_btn = ctk.CTkButton(
            dir_frame,
            text="Browse",
            command=self.select_output_dir,
            width=80,
            height=35
        )
        browse_btn.pack(side="right")
        
        # Tarjeta de categorización
        cat_card = ctk.CTkFrame(scrollable_frame)
        cat_card.pack(fill="x", pady=(0, 15))
        
        ctk.CTkLabel(
            cat_card,
            text="🗂️ Organization Method",
            font=ctk.CTkFont(size=16, weight="bold")
        ).pack(anchor="w", padx=20, pady=(20, 10))
        
        categorize_combo = ctk.CTkComboBox(
            cat_card,
            variable=self.categorize_by,
            values=CATEGORIZE_OPTIONS,
            height=35,
            width=300
        )
        categorize_combo.pack(anchor="w", padx=20, pady=(0, 20))
        
        # Tarjeta de fuente de datos
        source_card = ctk.CTkFrame(scrollable_frame)
        source_card.pack(fill="x", pady=(0, 15))
        
        ctk.CTkLabel(
            source_card,
            text="🌐 Data Source",
            font=ctk.CTkFont(size=16, weight="bold")
        ).pack(anchor="w", padx=20, pady=(20, 10))
        
        source_combo = ctk.CTkComboBox(
            source_card,
            variable=self.source,
            values=["pubchem", "chemspider"],
            height=35,
            width=300
        )
        source_combo.pack(anchor="w", padx=20, pady=(0, 20))
        
        # Botón de guardar configuración
        save_btn = ctk.CTkButton(
            scrollable_frame,
            text=f"{self.icons['save']} Save Configuration",
            command=self.save_configuration,
            height=40,
            width=200
        )
        save_btn.pack(pady=20)
    
    def show_single_download(self):
        """Muestra la vista de descarga individual"""
        if self.current_view == "single":
            return
        
        self.activate_nav_button(1)
        self.current_view = "single"
        self.clear_content()
        
        # Scrollable frame
        scrollable_frame = ctk.CTkScrollableFrame(self.content_frame)
        scrollable_frame.pack(fill="both", expand=True, padx=25, pady=25)
        
        # Título
        title = ctk.CTkLabel(
            scrollable_frame,
            text=f"{self.icons['molecules']} Single Molecule Download",
            font=ctk.CTkFont(size=22, weight="bold")
        )
        title.pack(pady=(0, 20), anchor="w")
        
        # Tarjeta principal
        main_card = ctk.CTkFrame(scrollable_frame)
        main_card.pack(fill="x", pady=(0, 20))
        
        # Campo de entrada
        ctk.CTkLabel(
            main_card,
            text="Enter molecule name or identifier:",
            font=ctk.CTkFont(size=14, weight="bold")
        ).pack(anchor="w", padx=30, pady=(30, 10))
        
        input_frame = ctk.CTkFrame(main_card, fg_color="transparent")
        input_frame.pack(fill="x", padx=30, pady=(0, 20))
        
        self.single_entry = ctk.CTkEntry(
            input_frame,
            placeholder_text="e.g., aspirin, caffeine, C8H10N4O2",
            height=40,
            font=ctk.CTkFont(size=14)
        )
        self.single_entry.pack(side="left", fill="x", expand=True, padx=(0, 15))
        
        download_btn = ctk.CTkButton(
            input_frame,
            text=f"{self.icons['download']} Download",
            command=self.download_single,
            height=40,
            width=120,
            font=ctk.CTkFont(size=14, weight="bold")
        )
        download_btn.pack(side="right")
        
        # Info adicional
        info_text = """
Tips for better results:
• Use common names: aspirin, caffeine, glucose
• Use IUPAC names for precise identification
• PubChem CIDs work great: 2244 (aspirin)
• SMILES strings are supported
• InChI identifiers are accepted
        """
        
        info_label = ctk.CTkLabel(
            main_card,
            text=info_text.strip(),
            font=ctk.CTkFont(size=11),
            justify="left"
        )
        info_label.pack(anchor="w", padx=30, pady=(0, 30))
        
        # Bind Enter key
        self.single_entry.bind("<Return>", lambda e: self.download_single())
    
    def show_batch_download(self):
        """Muestra la vista de descarga por lotes desde texto"""
        if self.current_view == "batch":
            return
        
        self.activate_nav_button(2)
        self.current_view = "batch"
        self.clear_content()
        
        # Scrollable frame
        scrollable_frame = ctk.CTkScrollableFrame(self.content_frame)
        scrollable_frame.pack(fill="both", expand=True, padx=25, pady=25)
        
        # Título
        title = ctk.CTkLabel(
            scrollable_frame,
            text=f"{self.icons['batch']} Batch Download from Text",
            font=ctk.CTkFont(size=22, weight="bold")
        )
        title.pack(pady=(0, 20), anchor="w")
        
        # Tarjeta principal
        main_card = ctk.CTkFrame(scrollable_frame)
        main_card.pack(fill="both", expand=True)
        
        # Área de texto
        ctk.CTkLabel(
            main_card,
            text="Enter molecule list (one per line):",
            font=ctk.CTkFont(size=14, weight="bold")
        ).pack(anchor="w", padx=30, pady=(30, 10))
        
        self.batch_textbox = ctk.CTkTextbox(
            main_card,
            height=300,
            font=ctk.CTkFont(size=12)
        )
        self.batch_textbox.pack(fill="both", expand=True, padx=30, pady=(0, 20))
        
        # Botones de acción
        button_frame = ctk.CTkFrame(main_card, fg_color="transparent")
        button_frame.pack(fill="x", padx=30, pady=(0, 30))
        
        load_file_btn = ctk.CTkButton(
            button_frame,
            text=f"{self.icons['upload']} Load from File",
            command=self.load_text_file,
            height=40
        )
        load_file_btn.pack(side="left", padx=(0, 10))
        
        clear_btn = ctk.CTkButton(
            button_frame,
            text="Clear List",
            command=self.clear_batch_text,
            height=40,
            fg_color="gray",
            hover_color="darkgray"
        )
        clear_btn.pack(side="left", padx=(0, 10))
        
        process_btn = ctk.CTkButton(
            button_frame,
            text=f"{self.icons['play']} Start Batch Download",
            command=self.start_batch_download,
            height=40,
            width=200
        )
        process_btn.pack(side="right")
        
        # Sample text
        sample_text = """aspirin
caffeine
morphine
penicillin
glucose
insulin"""
        self.batch_textbox.insert("0.0", sample_text)
    
    def show_excel_download(self):
        """Muestra la vista de descarga desde Excel"""
        if self.current_view == "excel":
            return
        
        self.activate_nav_button(3)
        self.current_view = "excel"
        self.clear_content()
        
        # Scrollable frame
        scrollable_frame = ctk.CTkScrollableFrame(self.content_frame)
        scrollable_frame.pack(fill="both", expand=True, padx=25, pady=25)
        
        # Título
        title = ctk.CTkLabel(
            scrollable_frame,
            text=f"{self.icons['excel']} Batch Download from Excel",
            font=ctk.CTkFont(size=22, weight="bold")
        )
        title.pack(pady=(0, 20), anchor="w")
        
        # Selección de archivo
        file_card = ctk.CTkFrame(scrollable_frame)
        file_card.pack(fill="x", pady=(0, 15))
        
        ctk.CTkLabel(
            file_card,
            text="📄 Excel File Selection",
            font=ctk.CTkFont(size=16, weight="bold")
        ).pack(anchor="w", padx=20, pady=(20, 10))
        
        file_frame = ctk.CTkFrame(file_card, fg_color="transparent")
        file_frame.pack(fill="x", padx=20, pady=(0, 20))
        
        self.excel_entry = ctk.CTkEntry(
            file_frame,
            textvariable=self.excel_file,
            placeholder_text="Select Excel file...",
            height=35
        )
        self.excel_entry.pack(side="left", fill="x", expand=True, padx=(0, 10))
        
        browse_excel_btn = ctk.CTkButton(
            file_frame,
            text="Browse",
            command=self.select_excel_file,
            width=80,
            height=35
        )
        browse_excel_btn.pack(side="right")
        
        # Configuración de Excel
        config_card = ctk.CTkFrame(scrollable_frame)
        config_card.pack(fill="x", pady=(0, 15))
        
        ctk.CTkLabel(
            config_card,
            text="⚙️ Excel Configuration",
            font=ctk.CTkFont(size=16, weight="bold")
        ).pack(anchor="w", padx=20, pady=(20, 10))
        
        config_frame = ctk.CTkFrame(config_card, fg_color="transparent")
        config_frame.pack(fill="x", padx=20, pady=(0, 20))
        
        # Sheet selection
        sheet_frame = ctk.CTkFrame(config_frame, fg_color="transparent")
        sheet_frame.pack(fill="x", pady=(0, 10))
        
        ctk.CTkLabel(sheet_frame, text="Sheet:", width=80).pack(side="left")
        self.sheet_combo = ctk.CTkComboBox(
            sheet_frame,
            variable=self.selected_sheet,
            values=["Select file first..."],
            width=200,
            command=self.on_sheet_selected
        )
        self.sheet_combo.pack(side="left", padx=(10, 0))
        
        # Column selection
        column_frame = ctk.CTkFrame(config_frame, fg_color="transparent")
        column_frame.pack(fill="x", pady=(0, 10))
        
        ctk.CTkLabel(column_frame, text="Column:", width=80).pack(side="left")
        self.column_combo = ctk.CTkComboBox(
            column_frame,
            variable=self.selected_column,
            values=["Select sheet first..."],
            width=200,
            command=self.on_column_selected
        )
        self.column_combo.pack(side="left", padx=(10, 0))
        
        # Preview
        preview_card = ctk.CTkFrame(scrollable_frame)
        preview_card.pack(fill="both", expand=True, pady=(0, 15))
        
        ctk.CTkLabel(
            preview_card,
            text="👁️ Data Preview",
            font=ctk.CTkFont(size=16, weight="bold")
        ).pack(anchor="w", padx=20, pady=(20, 10))
        
        self.excel_preview = ctk.CTkTextbox(
            preview_card,
            height=200,
            font=ctk.CTkFont(size=11)
        )
        self.excel_preview.pack(fill="both", expand=True, padx=20, pady=(0, 20))
        
        # Botones de acción
        action_frame = ctk.CTkFrame(scrollable_frame, fg_color="transparent")
        action_frame.pack(fill="x", pady=(0, 20))
        
        load_preview_btn = ctk.CTkButton(
            action_frame,
            text=f"{self.icons['refresh']} Load Preview",
            command=self.load_excel_preview,
            height=40
        )
        load_preview_btn.pack(side="left", padx=(0, 10))
        
        start_excel_btn = ctk.CTkButton(
            action_frame,
            text=f"{self.icons['play']} Start Excel Download",
            command=self.start_excel_download,
            height=40,
            width=200
        )
        start_excel_btn.pack(side="right")
    
    def show_sessions(self):
        """Muestra el gestor de sesiones"""
        if self.current_view == "sessions":
            return
        
        self.activate_nav_button(4)
        self.current_view = "sessions"
        self.clear_content()
        
        # Scrollable frame
        scrollable_frame = ctk.CTkScrollableFrame(self.content_frame)
        scrollable_frame.pack(fill="both", expand=True, padx=25, pady=25)
        
        # Título y botón de refresh
        header_frame = ctk.CTkFrame(scrollable_frame, fg_color="transparent")
        header_frame.pack(fill="x", pady=(0, 20))
        
        title = ctk.CTkLabel(
            header_frame,
            text=f"{self.icons['sessions']} Session Manager",
            font=ctk.CTkFont(size=22, weight="bold")
        )
        title.pack(side="left")
        
        refresh_btn = ctk.CTkButton(
            header_frame,
            text=f"{self.icons['refresh']} Refresh",
            command=self.refresh_sessions,
            height=35,
            width=100
        )
        refresh_btn.pack(side="right")
        
        # Lista de sesiones
        sessions_card = ctk.CTkFrame(scrollable_frame)
        sessions_card.pack(fill="both", expand=True)
        
        ctk.CTkLabel(
            sessions_card,
            text="📋 Available Sessions",
            font=ctk.CTkFont(size=16, weight="bold")
        ).pack(anchor="w", padx=20, pady=(20, 10))
        
        # Headers
        headers_frame = ctk.CTkFrame(sessions_card)
        headers_frame.pack(fill="x", padx=20, pady=(0, 10))
        
        ctk.CTkLabel(headers_frame, text="Name", width=200, font=ctk.CTkFont(weight="bold")).pack(side="left", padx=5)
        ctk.CTkLabel(headers_frame, text="Status", width=100, font=ctk.CTkFont(weight="bold")).pack(side="left", padx=5)
        ctk.CTkLabel(headers_frame, text="Progress", width=100, font=ctk.CTkFont(weight="bold")).pack(side="left", padx=5)
        ctk.CTkLabel(headers_frame, text="Created", width=150, font=ctk.CTkFont(weight="bold")).pack(side="left", padx=5)
        ctk.CTkLabel(headers_frame, text="Actions", width=200, font=ctk.CTkFont(weight="bold")).pack(side="left", padx=5)
        
        # Sessions container
        self.sessions_container = ctk.CTkScrollableFrame(sessions_card, height=300)
        self.sessions_container.pack(fill="both", expand=True, padx=20, pady=(0, 20))
        
        # Load sessions
        self.refresh_sessions()
    
    def show_monitor(self):
        """Muestra el monitor de progreso"""
        if self.current_view == "monitor":
            return
        
        self.activate_nav_button(5)
        self.current_view = "monitor"
        self.clear_content()
        
        # Scrollable frame
        scrollable_frame = ctk.CTkScrollableFrame(self.content_frame)
        scrollable_frame.pack(fill="both", expand=True, padx=25, pady=25)
        
        # Título
        title = ctk.CTkLabel(
            scrollable_frame,
            text=f"{self.icons['monitor']} Progress Monitor",
            font=ctk.CTkFont(size=22, weight="bold")
        )
        title.pack(pady=(0, 20), anchor="w")
        
        # Current progress card
        progress_card = ctk.CTkFrame(scrollable_frame)
        progress_card.pack(fill="x", pady=(0, 15))
        
        ctk.CTkLabel(
            progress_card,
            text="📊 Current Operation",
            font=ctk.CTkFont(size=16, weight="bold")
        ).pack(anchor="w", padx=20, pady=(20, 10))
        
        # Progress details
        progress_details_frame = ctk.CTkFrame(progress_card, fg_color="transparent")
        progress_details_frame.pack(fill="x", padx=20, pady=(0, 20))
        
        self.current_operation_label = ctk.CTkLabel(
            progress_details_frame,
            text="No active operation",
            font=ctk.CTkFont(size=14)
        )
        self.current_operation_label.pack(anchor="w", pady=(0, 10))
        
        self.main_progress_bar = ctk.CTkProgressBar(
            progress_details_frame,
            width=400,
            height=20
        )
        self.main_progress_bar.pack(anchor="w", pady=(0, 5))
        
        self.main_progress_label = ctk.CTkLabel(
            progress_details_frame,
            text="0 / 0 (0%)",
            font=ctk.CTkFont(size=12)
        )
        self.main_progress_label.pack(anchor="w")
        
        # Statistics card
        stats_card = ctk.CTkFrame(scrollable_frame)
        stats_card.pack(fill="x", pady=(0, 15))
        
        ctk.CTkLabel(
            stats_card,
            text="📈 Session Statistics",
            font=ctk.CTkFont(size=16, weight="bold")
        ).pack(anchor="w", padx=20, pady=(20, 10))
        
        stats_grid = ctk.CTkFrame(stats_card, fg_color="transparent")
        stats_grid.pack(fill="x", padx=20, pady=(0, 20))
        
        # Stats items
        stats_items = [
            ("Total Processed:", "total_processed"),
            ("Successful:", "successful"),
            ("Failed:", "failed"),
            ("Success Rate:", "success_rate")
        ]
        
        self.stats_labels = {}
        for i, (label_text, key) in enumerate(stats_items):
            row = i // 2
            col = i % 2
            
            item_frame = ctk.CTkFrame(stats_grid, fg_color="transparent")
            item_frame.grid(row=row, column=col, sticky="w", padx=(0, 40), pady=2)
            
            ctk.CTkLabel(item_frame, text=label_text, width=100).pack(side="left")
            self.stats_labels[key] = ctk.CTkLabel(item_frame, text="0", font=ctk.CTkFont(weight="bold"))
            self.stats_labels[key].pack(side="left")
        
        # Log card
        log_card = ctk.CTkFrame(scrollable_frame)
        log_card.pack(fill="both", expand=True)
        
        log_header = ctk.CTkFrame(log_card, fg_color="transparent")
        log_header.pack(fill="x", padx=20, pady=(20, 10))
        
        ctk.CTkLabel(
            log_header,
            text="📝 Activity Log",
            font=ctk.CTkFont(size=16, weight="bold")
        ).pack(side="left")
        
        clear_log_btn = ctk.CTkButton(
            log_header,
            text="Clear Log",
            command=self.clear_log,
            height=30,
            width=80
        )
        clear_log_btn.pack(side="right")
        
        self.log_textbox = ctk.CTkTextbox(
            log_card,
            height=250,
            font=ctk.CTkFont(size=11)
        )
        self.log_textbox.pack(fill="both", expand=True, padx=20, pady=(0, 20))
    
    # =========================================================================
    # FUNCIONALIDADES DE CONFIGURACIÓN
    # =========================================================================
    
    def configure_downloader(self):
        """Configura el descargador con callbacks"""
        self.downloader.configure(
            output_dir=self.output_dir.get(),
            categorize_by=self.categorize_by.get(),
            progress_callback=self.update_progress,
            log_callback=self.log_message
        )
    
    def save_configuration(self):
        """Guarda la configuración actual"""
        try:
            self.configure_downloader()
            self.log_message("✅ Configuration saved successfully")
            messagebox.showinfo("Success", "Configuration saved successfully!")
        except Exception as e:
            self.log_message(f"❌ Error saving configuration: {e}")
            messagebox.showerror("Error", f"Error saving configuration: {e}")
    
    def select_output_dir(self):
        """Selecciona directorio de salida"""
        directory = filedialog.askdirectory(initialdir=self.output_dir.get())
        if directory:
            self.output_dir.set(directory)
            self.log_message(f"📁 Output directory changed to: {directory}")
    
    # =========================================================================
    # FUNCIONALIDADES DE DESCARGA
    # =========================================================================
    
    def show_quick_download(self):
        """Muestra diálogo de descarga rápida"""
        self.show_single_download()
        if hasattr(self, 'single_entry'):
            self.single_entry.focus()
    
    def download_single(self):
        """Descarga una molécula individual"""
        if self.is_processing:
            messagebox.showwarning("Warning", "Another download is in progress")
            return
        
        query = self.single_entry.get().strip()
        if not query:
            messagebox.showerror("Error", "Please enter a molecule name")
            return
        
        # Iniciar descarga en thread separado
        self.is_processing = True
        self.update_status("Starting download...")
        
        thread = threading.Thread(
            target=self._download_single_thread,
            args=(query,),
            daemon=True
        )
        thread.start()
    
    def _download_single_thread(self, query):
        """Thread para descarga individual"""
        try:
            self.log_message(f"🔍 Starting download for: {query}")
            
            # Configurar descargador
            self.configure_downloader()
            
            success = self.downloader.download_molecule(query, self.source.get())
            
            if success:
                self.log_message(f"✅ Successfully downloaded: {query}")
                self.root.after(0, lambda: messagebox.showinfo("Success", f"Molecule '{query}' downloaded successfully!"))
            else:
                self.log_message(f"❌ Failed to download: {query}")
                self.root.after(0, lambda: messagebox.showerror("Error", f"Could not download '{query}'"))
                
        except Exception as e:
            self.log_message(f"❌ Error downloading {query}: {e}")
            self.root.after(0, lambda: messagebox.showerror("Error", f"Error: {e}"))
        finally:
            self.is_processing = False
            self.root.after(0, lambda: self.update_status("Ready"))
            self.root.after(0, lambda: self.update_progress(0, 1))
    
    def load_text_file(self):
        """Carga archivo de texto con moléculas"""
        file_path = filedialog.askopenfilename(
            title="Select text file",
            filetypes=[("Text files", "*.txt"), ("CSV files", "*.csv"), ("All files", "*.*")]
        )
        
        if file_path:
            try:
                molecules = self.file_parser.parse_text(file_path)
                self.batch_textbox.delete("0.0", "end")
                self.batch_textbox.insert("0.0", "\n".join(molecules))
                self.log_message(f"📄 Loaded {len(molecules)} molecules from {Path(file_path).name}")
            except Exception as e:
                messagebox.showerror("Error", f"Error loading file: {e}")
    
    def clear_batch_text(self):
        """Limpia el área de texto por lotes"""
        self.batch_textbox.delete("0.0", "end")
    
    def start_batch_download(self):
        """Inicia descarga por lotes"""
        if self.is_processing:
            messagebox.showwarning("Warning", "Another download is in progress")
            return
        
        text_content = self.batch_textbox.get("0.0", "end").strip()
        if not text_content:
            messagebox.showerror("Error", "Please enter some molecules")
            return
        
        queries = [line.strip() for line in text_content.split('\n') if line.strip()]
        if not queries:
            messagebox.showerror("Error", "No valid molecules found")
            return
        
        # Crear sesión
        session_name = f"Batch_Text_{len(queries)}_molecules"
        session_id = self.progress_manager.create_session(
            session_name=session_name,
            total_items=len(queries),
            metadata={"source_type": "manual", "queries": queries}
        )
        
        self.current_session_id = session_id
        self.is_processing = True
        
        # Iniciar descarga en thread separado
        thread = threading.Thread(
            target=self._download_batch_thread,
            args=(queries, session_id),
            daemon=True
        )
        thread.start()
        
        messagebox.showinfo("Started", f"Batch download started for {len(queries)} molecules")
    
    # =========================================================================
    # FUNCIONALIDADES DE EXCEL
    # =========================================================================
    
    def select_excel_file(self):
        """Selecciona archivo Excel"""
        file_path = filedialog.askopenfilename(
            title="Select Excel file",
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
            self.excel_sheets = self.file_parser.get_excel_sheets(file_path)
            self.sheet_combo.configure(values=self.excel_sheets)
            
            if self.excel_sheets:
                self.selected_sheet.set(self.excel_sheets[0])
                self.on_sheet_selected()
                
            self.log_message(f"📊 Loaded Excel file: {Path(file_path).name}")
                
        except Exception as e:
            messagebox.showerror("Error", f"Error loading Excel file: {e}")
    
    def on_sheet_selected(self, *args):
        """Se ejecuta cuando se selecciona una hoja"""
        try:
            file_path = self.excel_file.get()
            sheet_name = self.selected_sheet.get()
            
            if file_path and sheet_name:
                self.excel_columns = self.file_parser.get_excel_columns(file_path, sheet_name)
                self.column_combo.configure(values=self.excel_columns)
                
                if self.excel_columns:
                    self.selected_column.set(self.excel_columns[0])
                    
        except Exception as e:
            messagebox.showerror("Error", f"Error loading columns: {e}")
    
    def on_column_selected(self, *args):
        """Se ejecuta cuando se selecciona una columna"""
        self.load_excel_preview()
    
    def load_excel_preview(self):
        """Carga vista previa de los datos del Excel"""
        try:
            file_path = self.excel_file.get()
            sheet_name = self.selected_sheet.get()
            column_name = self.selected_column.get()
            
            if not all([file_path, sheet_name, column_name]):
                return
            
            # Obtener datos de la columna
            molecules = self.file_parser.parse_excel(
                file_path, 
                column_name=column_name, 
                sheet_name=sheet_name
            )
            
            self.excel_preview_data = molecules
            
            # Mostrar vista previa (primeros 20)
            preview_text = f"Total molecules found: {len(molecules)}\n\n"
            preview_text += "Preview (first 20):\n"
            preview_text += "\n".join(f"{i+1}. {mol}" for i, mol in enumerate(molecules[:20]))
            
            if len(molecules) > 20:
                preview_text += f"\n... and {len(molecules) - 20} more"
            
            self.excel_preview.delete("0.0", "end")
            self.excel_preview.insert("0.0", preview_text)
            
        except Exception as e:
            messagebox.showerror("Error", f"Error loading preview: {e}")
    
    def start_excel_download(self):
        """Inicia descarga desde Excel"""
        if self.is_processing:
            messagebox.showwarning("Warning", "Another download is in progress")
            return
        
        if not self.excel_preview_data:
            messagebox.showwarning("Warning", "Please load and preview Excel data first")
            return
        
        molecules = self.excel_preview_data
        
        # Crear sesión de progreso
        session_name = f"Excel_{Path(self.excel_file.get()).stem}_{self.selected_column.get()}"
        session_id = self.progress_manager.create_session(
            session_name=session_name,
            total_items=len(molecules),
            metadata={
                "source_file": self.excel_file.get(),
                "sheet_name": self.selected_sheet.get(),
                "column_name": self.selected_column.get(),
                "source_type": "excel"
            }
        )
        
        self.current_session_id = session_id
        self.is_processing = True
        
        # Iniciar descarga en thread separado
        thread = threading.Thread(
            target=self._download_batch_thread,
            args=(molecules, session_id),
            daemon=True
        )
        thread.start()
        
        messagebox.showinfo("Started", f"Excel download started for {len(molecules)} molecules")
    
    # =========================================================================
    # FUNCIONALIDADES DE SESIONES
    # =========================================================================
    
    def refresh_sessions(self):
        """Actualiza la lista de sesiones"""
        try:
            # Limpiar container
            for widget in self.sessions_container.winfo_children():
                widget.destroy()
            
            # Cargar sesiones
            sessions = self.progress_manager.list_sessions()
            
            if not sessions:
                no_sessions_label = ctk.CTkLabel(
                    self.sessions_container,
                    text="No sessions found",
                    font=ctk.CTkFont(size=14),
                    text_color="gray"
                )
                no_sessions_label.pack(pady=50)
                return
            
            for session in sessions:
                self.create_session_item(session)
                
            self.log_message(f"🔄 Refreshed sessions: {len(sessions)} found")
                
        except Exception as e:
            pass #debo explicar el error
    
    def create_session_item(self, session):
        """Crea un item de sesión en la lista"""
        item_frame = ctk.CTkFrame(self.sessions_container)
        item_frame.pack(fill="x", pady=2)
        
        # Session info
        info_frame = ctk.CTkFrame(item_frame, fg_color="transparent")
        info_frame.pack(fill="x", padx=10, pady=5)
        
        ctk.CTkLabel(info_frame, text=session['session_name'], width=200).pack(side="left", padx=5)
        ctk.CTkLabel(info_frame, text=session['status'], width=100).pack(side="left", padx=5)
        ctk.CTkLabel(info_frame, text=session['progress'], width=100).pack(side="left", padx=5)
        ctk.CTkLabel(info_frame, text=session['created_at'][:19], width=150).pack(side="left", padx=5)
        
        # Action buttons
        actions_frame = ctk.CTkFrame(info_frame, fg_color="transparent")
        actions_frame.pack(side="left", padx=5)
        
        if session['status'] != 'completed':
            resume_btn = ctk.CTkButton(
                actions_frame,
                text="Resume",
                command=lambda s=session: self.resume_session(s),
                height=25,
                width=60
            )
            resume_btn.pack(side="left", padx=2)
        
        delete_btn = ctk.CTkButton(
            actions_frame,
            text="Delete",
            command=lambda s=session: self.delete_session(s),
            height=25,
            width=60,
            fg_color="red",
            hover_color="darkred"
        )
        delete_btn.pack(side="left", padx=2)
    
    def resume_session(self, session):
        """Continúa una sesión"""
        if self.is_processing:
            messagebox.showwarning("Warning", "Another download is in progress")
            return
        
        try:
            session_id = session['session_id']
            session_data = self.progress_manager.load_session(session_id)
            
            if not session_data:
                messagebox.showerror("Error", "Could not load session data")
                return
            
            if session_data['status'] == 'completed':
                messagebox.showinfo("Info", "This session is already completed")
                return
            
            # Reconstruir lista de moléculas
            molecules = self._reconstruct_molecule_list(session_data)
            
            if not molecules:
                messagebox.showerror("Error", "Could not reconstruct molecule list")
                return
            
            # Obtener moléculas pendientes
            pending_molecules = self.progress_manager.get_pending_queries(session_id, molecules)
            
            if not pending_molecules:
                messagebox.showinfo("Info", "No pending molecules in this session")
                return
            
            # Continuar descarga
            self.current_session_id = session_id
            self.is_processing = True
            
            thread = threading.Thread(
                target=self._download_batch_thread,
                args=(pending_molecules, session_id),
                daemon=True
            )
            thread.start()
            
            messagebox.showinfo("Resumed", f"Continuing with {len(pending_molecules)} pending molecules")
            
        except Exception as e:
            messagebox.showerror("Error", f"Error resuming session: {e}")
    
    def delete_session(self, session):
        """Elimina una sesión"""
        if messagebox.askyesno("Confirm", f"Delete session '{session['session_name']}'?"):
            try:
                self.progress_manager.delete_session(session['session_id'])
                self.refresh_sessions()
                self.log_message(f"🗑️ Deleted session: {session['session_name']}")
                messagebox.showinfo("Deleted", "Session deleted successfully")
            except Exception as e:
                messagebox.showerror("Error", f"Error deleting session: {e}")
    
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
            elif source_type == 'manual':
                # Para listas manuales, usar las queries guardadas
                return metadata.get('queries', [])
            else:
                # Usar las queries completadas y fallidas
                completed = session_data.get('completed_queries', [])
                failed = session_data.get('failed_queries', [])
                return completed + failed
                
        except Exception as e:
            self.log_message(f"❌ Error reconstructing molecule list: {e}")
            return []
    
    # =========================================================================
    # FUNCIONALIDADES DE DESCARGA EN LOTES
    # =========================================================================
    
    def _download_batch_thread(self, queries, session_id):
        """Thread para descarga por lotes con manejo de sesiones"""
        try:
            self.root.after(0, lambda: self.update_status(f"Processing {len(queries)} molecules"))
            self.log_message(f"🚀 Starting batch download: {len(queries)} molecules")
            
            # Configurar descargador
            self.configure_downloader()
            
            completed = 0
            failed = 0
            
            for i, query in enumerate(queries):
                try:
                    # Update progress on main thread
                    self.root.after(0, lambda i=i, total=len(queries): self.update_progress(i, total))
                    self.root.after(0, lambda q=query, i=i, total=len(queries): 
                                  self.update_status(f"Processing: {q} ({i+1}/{total})"))
                    
                    self.log_message(f"🔍 Processing: {query} ({i+1}/{len(queries)})")
                    
                    success = self.downloader.download_molecule(query, self.source.get())
                    
                    if success:
                        completed += 1
                        self.progress_manager.update_progress(session_id, completed_query=query)
                        self.log_message(f"✅ {query} completed")
                    else:
                        failed += 1
                        self.progress_manager.update_progress(session_id, failed_query=query)
                        self.log_message(f"❌ {query} failed")
                        
                    # Update stats on main thread
                    self.root.after(0, lambda: self.update_monitor_stats(completed, failed, len(queries)))
                        
                except Exception as e:
                    failed += 1
                    self.progress_manager.update_progress(session_id, failed_query=query)
                    self.log_message(f"❌ Error with {query}: {e}")
            
            # Completar sesión
            self.progress_manager.complete_session(session_id)
            
            # Final updates on main thread
            self.root.after(0, lambda: self.update_progress(len(queries), len(queries)))
            self.root.after(0, lambda: self.update_status("Completed"))
            
            message = f"Process completed:\n✅ Successful: {completed}\n❌ Failed: {failed}\n📊 Total: {len(queries)}"
            self.log_message(message)
            
            self.root.after(0, lambda: messagebox.showinfo("Completed", message))
            
        except Exception as e:
            self.log_message(f"❌ Critical error: {e}")
            self.root.after(0, lambda: messagebox.showerror("Error", f"Critical error: {e}"))
        finally:
            self.is_processing = False
            self.current_session_id = None
            self.root.after(0, lambda: self.update_status("Ready"))
    
    # =========================================================================
    # FUNCIONALIDADES DE MONITOR Y PROGRESO
    # =========================================================================
    
    def update_progress(self, current, total):
        """Actualiza barra de progreso"""
        if total > 0:
            progress_value = current / total
            self.progress_var.set(progress_value)
            percent = progress_value * 100
            self.progress_label.configure(text=f"{percent:.1f}%")
            
            # Update main progress in monitor
            if hasattr(self, 'main_progress_bar'):
                self.main_progress_bar.set(progress_value)
                self.main_progress_label.configure(text=f"{current} / {total} ({percent:.1f}%)")
        else:
            self.progress_var.set(0)
            self.progress_label.configure(text="0%")
    
    def update_status(self, status):
        """Actualiza el estado del sistema"""
        self.status_var.set(status)
        if hasattr(self, 'current_operation_label'):
            self.current_operation_label.configure(text=status)
    
    def update_monitor_stats(self, completed, failed, total):
        """Actualiza estadísticas del monitor"""
        if hasattr(self, 'stats_labels'):
            self.stats_labels['total_processed'].configure(text=str(completed + failed))
            self.stats_labels['successful'].configure(text=str(completed))
            self.stats_labels['failed'].configure(text=str(failed))
            
            if total > 0:
                success_rate = (completed / (completed + failed)) * 100 if (completed + failed) > 0 else 0
                self.stats_labels['success_rate'].configure(text=f"{success_rate:.1f}%")
    
    def log_message(self, message):
        """Añade mensaje al log"""
        timestamp = time.strftime('%H:%M:%S')
        log_entry = f"[{timestamp}] {message}\n"
        
        if hasattr(self, 'log_textbox'):
            self.log_textbox.insert("end", log_entry)
            self.log_textbox.see("end")
        
        print(log_entry.strip())  # También imprimir en consola
    
    def clear_log(self):
        """Limpia el log de actividad"""
        if hasattr(self, 'log_textbox'):
            self.log_textbox.delete("0.0", "end")
    
    def update_footer_stats(self):
        """Actualiza estadísticas del footer"""
        try:
            # Get storage stats
            stats = self.storage_manager.get_statistics()
            sessions = self.progress_manager.list_sessions()
            active_sessions = len([s for s in sessions if s['status'] == 'active'])
            
            status_text = f"Ready • {stats['total_files']} molecules downloaded • {active_sessions} sessions active"
            self.stats_label.configure(text=status_text)
            
            # Schedule next update
            self.root.after(5000, self.update_footer_stats)  # Update every 5 seconds
        except:
            pass
    
    # =========================================================================
    # FUNCIONALIDADES ADICIONALES
    # =========================================================================
    
    def open_results_folder(self):
        """Abre la carpeta de resultados"""
        try:
            folder_path = self.output_dir.get()
            if os.path.exists(folder_path):
                if os.name == 'nt':  # Windows
                    os.startfile(folder_path)
                elif os.name == 'posix':  # macOS y Linux
                    if os.uname().sysname == 'Darwin':  # macOS
                        subprocess.call(['open', folder_path])
                    else:  # Linux
                        subprocess.call(['xdg-open', folder_path])
                        
                self.log_message(f"📁 Opened results folder: {folder_path}")
            else:
                messagebox.showwarning("Warning", "Results folder does not exist yet")
        except Exception as e:
            messagebox.showerror("Error", f"Cannot open folder: {e}")
    
    def emergency_stop(self):
        """Parada de emergencia"""
        if self.is_processing:
            if messagebox.askyesno("Emergency Stop", "Are you sure you want to stop the current operation?"):
                self.is_processing = False
                self.update_status("Stopped by user")
                self.log_message("🛑 Emergency stop activated by user")
                messagebox.showinfo("Stopped", "Operation stopped")
        else:
            messagebox.showinfo("Info", "No operation is currently running")
    
    def run(self):
        """Ejecuta la aplicación"""
        # Start footer stats updates
        self.update_footer_stats()
        
        # Show welcome message
        self.log_message("🧪 Molecular Downloader Pro initialized")
        self.log_message("✨ Ready for molecular structure downloads")
        
        self.root.mainloop()

def main():
    """Función principal"""
    try:
        app = MoleculeApp()
        app.run()
    except Exception as e:
        messagebox.showerror("Critical Error", f"Failed to start application: {e}")
        print(f"Critical error: {e}")

if __name__ == "__main__":
    main()
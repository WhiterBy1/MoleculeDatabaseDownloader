# 📋 Guía Paso a Paso: Cómo Usar el Descargador de Moléculas

## 🚀 PARTE 1: Instalación desde GitHub

### Paso 1: Descargar el código
1. Ve al repositorio de GitHub donde está alojado el proyecto
2. Haz clic en el botón verde **"Code"** 
3. Selecciona **"Download ZIP"**
4. Extrae el archivo ZIP en tu computadora (por ejemplo, en `Escritorio/molecule_downloader/`)



### Paso 2: Navegar a la carpeta del proyecto en tu IDE (Visual studio Code)
desde el terminal ya estando en la carpeta madre se puede abrir escribiendo:
```bash
code Escritorio/molecule_downloader
```
*(Ajusta la ruta según donde extrajiste el ZIP)*

### Paso 3: Instalar dependencias

```bash
pip install -r requirements.txt
```
⏳ *Espera a que se instalen todas las librerías (puede tomar 2-5 minutos)*

### Paso 4: Ejecutar la aplicación
```bash
python main.py
```
o darle a play al archivo con nombre main desde VScode
🎉 *¡La aplicación debería abrirse!*

---

## 📊 PARTE 2: Usar la Descarga desde Excel

### Paso 1: Preparar tu archivo Excel
- Tu Excel debe tener una columna con nombres de moléculas o los ID

### Paso 2: Cargar tu archivo Excel
1. Ve a la pestaña **"📊 Desde Excel"**
2. Haz clic en **📁** junto a "Archivo Excel:"
3. Navega y selecciona tu archivo `.xlsx` o `.xls`
4. Haz clic en **"Abrir"**

### Paso 3: Seleccionar hoja y columna (si no es la por defecto que se autoselecciona)
1. En el menú desplegable **"Hoja:"**, selecciona la hoja que contiene tus moléculas
2. En **"Columna con moléculas:"**, selecciona la columna que tiene los nombres
3. Haz clic en **"🔄 Cargar Archivo"**

### Paso 4: Verificar vista previa
- En el área "Vista previa" aparecerán las primeras 20 moléculas
- Verifica que sean correctas
- Si ves errores, cambia la columna seleccionada

### Paso 5: Iniciar descarga
1. Haz clic en **"▶️ Iniciar Descarga"**
2. Aparecerá un mensaje confirmando el inicio
3. Ve a la pestaña **"📊 Monitor"** para ver el progreso

### Paso 6: Monitorear progreso
- La barra de progreso muestra el avance
- El log muestra cada molécula procesada:
  - ✅ = Descarga exitosa
  - ❌ = Error en la descarga

---

## 🔄 PARTE 3: Recuperar Sesión Interrumpida

### Si la descarga se interrumpe:
1. Vuelve a ejecutar `python main.py`
2. Ve a la pestaña **"📋 Sesiones"**
3. Haz clic en **"🔄 Actualizar"**
4. Selecciona tu sesión interrumpida
5. Haz clic en **"▶️ Continuar Seleccionada"**
6. La descarga continuará desde donde se quedó

---

## 📁 PARTE 4: Encontrar tus moléculas descargadas

### Las moléculas se organizan automáticamente:
```
molecules/
├── 
```

### Para abrir la carpeta:
1. Ve al menú de configuración
2. Haz clic en **"Abrir Carpeta"** (si está disponible)
3. O navega manualmente al directorio que configuraste

---

## ⚠️ Solución de Problemas Comunes

### Error: "ModuleNotFoundError"
```bash
pip install --upgrade pip
pip install -r requirements.txt
```

### Error: "No se puede leer el archivo Excel"
- Verifica que el archivo no esté abierto en Excel
- Asegúrate de que sea `.xlsx` o `.xls`

### Error: "No internet connection"
- Verifica tu conexión a internet
- Las moléculas se descargan de PubChem (requiere internet)

### La aplicación no se abre
- Verifica que instalaste Python correctamente
- Ejecuta: `python --version` (debe mostrar Python 3.8+)

---

## 🎯 Ejemplo Completo


**Pasos:**
1. ✅ Ejecutar `python main.py`
2. ✅ Ir a **"📊 Desde Excel"**
3. ✅ Cargar archivo Excel
4. ✅ Seleccionar columna **"compound_name"**
5. ✅ Clic en **"🔄 Cargar Archivo"**
6. ✅ Verificar vista previa
7. ✅ Clic en **"▶️ Iniciar Descarga"**
8. ✅ Ir a **"📊 Monitor"** para ver progreso
9. ✅ ¡Listo! Moléculas en la carpeta configurada

**Resultado esperado:**
- Archivos SDF descargados
  Organizados automáticamente por categoría
- Sesión guardada para recuperación

¡Tu sistema está listo para descargar cientos de moléculas automáticamente! 🧬

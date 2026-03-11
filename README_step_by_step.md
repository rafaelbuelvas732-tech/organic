# 🧪 Extractor de Alcoholes desde PubChem (Versión Ultra Robusta)

Este script te permite construir una base de datos de **alcoholes** (compuestos con grupo –OH) obteniendo información desde la API de PubChem y enriqueciéndola con **descriptores moleculares calculados con RDKit**. Está diseñado para ser **robusto, reanudable y eficiente**, pudiendo extraer **miles de compuestos** en un tiempo razonable (entre 30 y 60 minutos, dependiendo de tu conexión y configuración).

## 🎯 ¿Para qué sirve?

- Obtener una colección grande y diversa de alcoholes, clasificados por tipo (primarios, secundarios, terciarios, dioles, polioles, aromáticos, etc.).
- Descargar propiedades fisicoquímicas y topológicas desde PubChem (peso molecular, XLogP, TPSA, etc.).
- Calcular descriptores adicionales con RDKit (índices de conectividad, forma, recuento de anillos, etc.) para análisis más profundos.
- Almacenar todo en una **base de datos SQLite** local, con **checkpoints** para poder reanudar si se interrumpe.
- Generar al final un archivo CSV listo para usar en Python, Excel o cualquier herramienta de análisis.

## 🧠 ¿Qué conceptos importantes debo conocer?

### 🔹 SMARTS vs SMILES
- **SMILES** es una forma de escribir una molécula como texto (ej. `CCO` para etanol).
- **SMARTS** es un lenguaje más potente para **buscar patrones estructurales** dentro de moléculas. Por ejemplo, `[CH2][OH]` busca un carbono con dos hidrógenos unido a un OH (alcohol primario).  
  El script usa SMARTS para encontrar compuestos que coincidan con cada tipo de alcohol.

### 🔹 Cuotas (quotas)
- El script tiene un objetivo de cuántos compuestos queremos de cada categoría (ej. 30.000 primarios).  
- Esto nos permite controlar el tamaño final de la base y equilibrar la representación de cada tipo.

### 🔹 Checkpoint y reanudación
- El script guarda periódicamente un archivo `checkpoint.json` con los IDs de compuestos ya procesados.  
- Si la ejecución se interrumpe (por error, por apagado, o porque pulsas Ctrl+C), al volver a ejecutar continuará desde donde lo dejó, sin repetir trabajo.

### 🔹 SQLite
- Los datos se guardan en un archivo local (`alcoholes_robusto.db`).  
- Esto permite almacenar millones de registros de forma eficiente y consultarlos después sin necesidad de volver a descargar.

### 🔹 Procesamiento paralelo
- Se usan **hilos** (ThreadPoolExecutor) para descargar lotes de propiedades desde PubChem más rápido.
- Se usan **procesos** (ProcessPoolExecutor) para calcular los descriptores RDKit, aprovechando múltiples núcleos de la CPU.

### 🔹 Descriptores RDKit
- Son valores numéricos que describen aspectos de la molécula:  
  - **Chi0v, Chi1v, Chi2v**: índices de conectividad de Kier-Hall (con valencia), relacionados con la ramificación y el tamaño.  
  - **Kappa1, Kappa2, Kappa3**: índices de forma.  
  - **LabuteASA**: área superficial accesible aproximada.  
  - **QED**: "quantitative estimate of drug-likeness" (estimación de similitud a fármacos).  
  - **FpDensityMorgan**: densidad de fingerprints circulares, útil para búsqueda de similitud.  
- También se calculan propiedades equivalentes a las de PubChem (MolLogP, TPSA, etc.) para comparar o reemplazar nulos.

## 📦 Requisitos

- Python 3.8 o superior.
- Librerías necesarias (instalar con pip):

```bash
pip install requests pandas tqdm rdkit-pypi openpyxl
```

(Opcional: `openpyxl` solo si quieres guardar también en Excel).

## 🚀 ¿Cómo usar el script?

### 1. Configuración inicial

Abre el script y, si lo deseas, ajusta estos parámetros al principio (están en las primeras líneas después de los imports):

- **Cuotas** (`quota` en `ALCOHOL_CATEGORIES`): define cuántos compuestos quieres de cada categoría. Para una base rápida (~10.000 compuestos totales), pon valores como:
  ```python
  "primary":    {"smarts": "[CH2][OH]",                           "quota": 2000},
  "secondary":  {"smarts": "[CH1]([OH])([#6])[#6]",               "quota": 1500},
  ...
  ```
  Si quieres todos los que puedas, deja las cuotas altas (30.000, 25.000, etc.), pero la extracción durará más.

- **`calls_per_second`**: controla la velocidad de peticiones a PubChem para no saturar el servidor. El valor por defecto `0.33` significa 1 petición cada 3 segundos. Puedes aumentarlo a `0.5` si tu conexión es buena, pero respeta los límites de PubChem.

- **`BATCH_SIZE`**: número de compuestos por lote al descargar propiedades. `50` está bien.

- **`FETCH_WORKERS`** y **`RDKIT_WORKERS`**: número de hilos/procesos paralelos. Ajusta según los núcleos de tu CPU (ej. 4 para un i5).

### 2. Ejecución

Simplemente ejecuta el script:

```bash
python extractor_ultra_robusto.py
```

Verás en pantalla el progreso: cuántos compuestos se van encontrando, lotes procesados, etc. También se crea un archivo `extractor.log` con el registro detallado.

### 3. Si se interrumpe

Puedes parar con `Ctrl+C`. El script guardará un checkpoint y al volver a ejecutarlo continuará desde donde lo dejó. No perderás los datos ya descargados.

### 4. Resultados

Al finalizar, tendrás:

- Una base de datos SQLite: `alcoholes_robusto.db` (puedes consultarla con cualquier cliente SQLite o desde Python).
- Un archivo CSV: `alcoholes_robusto_<fecha>_<hora>.csv` con todos los datos.
- (Opcional) un archivo Excel si tienes `openpyxl`.

Además, en la consola se mostrará un resumen con:
- Número de compuestos por categoría y porcentaje de cuota alcanzado.
- Clasificación según RDKit (primarios, secundarios, etc.).
- Estadísticas básicas de peso molecular.

## 📁 Estructura de los datos

Las columnas del CSV final se dividen en:

### Identificadores (de PubChem)
- `cid`, `molecularformula`, `canonicalsmiles`, `isomericsmiles`, `inchi`, `inchikey`, `iupacname`
- `category` (categoría asignada durante la búsqueda), `extraction_date`

### Propiedades fisicoquímicas (PubChem)
- `molecularweight`, `xlogp`, `exactmass`, `monoisotopicmass`, `tpsa`, `hbonddonorcount`, `hbondacceptorcount`, `rotatablebondcount`, `heavyatomcount`, `volume3d`

### Propiedades topológicas (PubChem)
- `complexity`, `charge`, `isotopeatomcount`, `atomstereocount`, `definedatomstereocount`, `undefinedatomstereocount`, `bondstereocount`, `definedbondstereocount`, `undefinedbondstereocount`, `covalentunitcount`

### Descriptores RDKit (calculados)
- `RDKit_AlcoholClass` (clasificación según grupos OH)
- `RDKit_NumOH` (número de grupos OH)
- `RDKit_MolLogP`, `RDKit_MolWt`, `RDKit_NumHeavyAtoms`, `RDKit_NumRotatableBonds`, `RDKit_NumHDonors`, `RDKit_NumHAcceptors`, `RDKit_TPSA`, `RDKit_NumStereocenters`, `RDKit_MolMR`, `RDKit_NumAtoms`
- `RDKit_NumRings`, `RDKit_NumAromaticRings`, `RDKit_NumAliphaticRings`, `RDKit_NumSaturatedRings`, `RDKit_NumHeterocycles`
- `RDKit_FractionCsp3`, `RDKit_NumValenceElectrons`, `RDKit_NumRadicalElectrons`
- Índices topológicos: `RDKit_Chi0v`, `RDKit_Chi1v`, `RDKit_Chi2v`, `RDKit_Kappa1`, `RDKit_Kappa2`, `RDKit_Kappa3`, `RDKit_LabuteASA`, `RDKit_BalabanJ`
- Desglose de anillos: `RDKit_NumSaturatedCarbocycles`, `RDKit_NumSaturatedHeterocycles`, `RDKit_NumAromaticCarbocycles`, `RDKit_NumAromaticHeterocycles`, `RDKit_NumAliphaticCarbocycles`, `RDKit_NumAliphaticHeterocycles`
- `RDKit_QED`, `RDKit_FpDensityMorgan1`, `RDKit_FpDensityMorgan2`, `RDKit_FpDensityMorgan3`

## 🔬 Posibles análisis con estos datos

- **Estudio de diversidad**: ¿qué tipos de alcoholes están más representados?
- **Relación estructura-propiedad**: ¿cómo influye el número de anillos en la lipofilicidad (XLogP)?
- **Comparación de métodos**: ¿coinciden los valores de TPSA de PubChem con los de RDKit?
- **Clasificación**: ¿podemos predecir la clase de alcohol a partir de los descriptores?
- **Búsqueda de fármacos potenciales**: usar QED y otros índices para filtrar compuestos con propiedades favorables.

## 🧩 Personalización y ampliación

### Añadir más categorías
Puedes agregar nuevas categorías al diccionario `ALCOHOL_CATEGORIES` siguiendo el mismo formato. Por ejemplo, para alcoholes con grupos funcionales adicionales:

```python
"amino_alcohol": {"smarts": "[NX3][CH2][OH]", "quota": 5000}
```

### Cambiar a SMILES en lugar de InChI para RDKit
Si prefieres usar SMILES (puede ser más rápido y fiable), modifica la función `_rdkit_worker_inchi` para que reciba SMILES y cambia su nombre a `_rdkit_worker_smiles`. Luego, en `enhance_with_rdkit_from_inchi`, usa la columna `canonicalsmiles` en lugar de `inchi`.

### Dividir en grupos temáticos
Al final del script, puedes añadir una llamada a un función que genere los cuatro archivos (identificadores, fisicoquímicas, topológicas, anillos) como hicimos en versiones anteriores. Así tendrás subconjuntos listos para análisis específicos.

## 🛠️ Solución de problemas

- **Error de conexión**: si la API de PubChem no responde, el script reintentará hasta 5 veces con esperas crecientes. Si persiste, revisa tu conexión a internet.
- **RDKit no instalado**: si no tienes RDKit, el script lo detectará y continuará sin los descriptores calculados (solo con datos de PubChem).
- **Columna InChI no encontrada**: el script usa la columna `inchi` (en minúsculas) para RDKit. Si tu DataFrame la tiene con otro nombre, ajústalo en `posibles_inchi`.
- **Memoria insuficiente**: si el número de compuestos es muy grande, la fusión final puede consumir mucha RAM. Puedes reducir las cuotas o procesar por partes.

## 📚 Referencias

- [PubChem REST API](https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest)
- [RDKit Documentación](https://www.rdkit.org/docs/)
- [SMARTS tutorial](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)

## 👨‍🔬 Autor

Este script fue desarrollado como parte de un aprendizaje en quimioinformática, con el objetivo de construir una base de datos de alcoholes para explorar relaciones estructura–propiedad.

---

**¡Ahora a explorar la química de los alcoholes!** 🧪🔬

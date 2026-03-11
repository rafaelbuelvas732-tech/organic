# Base de Datos de Alcoholes – Proyecto de Quimioinformática

Este repositorio contiene una base de datos curada y enriquecida de **alcoholes** (compuestos con grupo –OH unido a carbono saturado), construida mediante consultas a la API de PubChem y cálculos de descriptores con RDKit. El objetivo es proporcionar un recurso estructurado y de alta calidad para estudios de relación estructura–propiedad, modelado QSAR, clasificación, aprendizaje automático y exploración química.

## 📌 Objetivos

- Obtener una colección representativa y amplia de alcoholes (decenas de miles de compuestos) distribuidos en categorías estructurales (primarios, secundarios, terciarios, dioles, polioles, aromáticos, alílicos, grasos, etc.).
- Extraer propiedades fisicoquímicas y topológicas directamente desde PubChem.
- Calcular descriptores moleculares avanzados con RDKit (índices de conectividad, forma, recuento de anillos, densidad de fingerprints, QED, etc.).
- Organizar la información en subconjuntos temáticos para facilitar análisis específicos.
- Implementar un pipeline robusto y escalable que permita reanudar extracciones interrumpidas (checkpoint, base de datos SQLite, reintentos con backoff).
- Sentar las bases para ampliar el proyecto a otros grupos funcionales (aldehídos, cetonas, ácidos carboxílicos, aminas, etc.).

## 🧪 Metodología

El proceso de construcción de la base de datos consta de varias etapas:

1. **Definición de categorías mediante SMARTS**:  
   Se utilizan patrones SMARTS (más potentes y precisos que SMILES) para definir cada tipo estructural de alcohol. Por ejemplo:
   - Primarios: `[CH2][OH]`
   - Secundarios: `[CH1]([OH])([#6])[#6]`
   - Terciarios: `[C;X4;H0]([OH])([#6])([#6])[#6]`
   - Aromáticos (fenoles): `c[OH]`
   - Dioles, polioles, alcoholes grasos, etc.

2. **Búsqueda de compuestos en PubChem**:  
   Para cada categoría se consulta el endpoint `/compound/fastsubstructure/smarts` de PubChem, obteniendo listas de CIDs (Compound IDs). Se respetan cuotas máximas por categoría (configurables).

3. **Extracción de propiedades básicas**:  
   Se recuperan 26 propiedades fisicoquímicas y topológicas mediante el endpoint `/compound/cid/property`. Estas incluyen:
   - Identificadores: `CID`, `MolecularFormula`, `CanonicalSMILES`, `IsomericSMILES`, `InChI`, `InChIKey`, `IUPACName`
   - Fisicoquímicas: `MolecularWeight`, `XLogP`, `ExactMass`, `MonoisotopicMass`, `TPSA`, `HBondDonorCount`, `HBondAcceptorCount`, `RotatableBondCount`, `HeavyAtomCount`, `Volume3D`
   - Topológicas: `Complexity`, `Charge`, `IsotopeAtomCount`, `AtomStereoCount`, `DefinedAtomStereoCount`, `UndefinedAtomStereoCount`, `BondStereoCount`, `DefinedBondStereoCount`, `UndefinedBondStereoCount`, `CovalentUnitCount`

4. **Almacenamiento progresivo en SQLite**:  
   Los datos se guardan en una base de datos SQLite local (`alcoholes_robusto.db`) a medida que se obtienen, con un sistema de checkpoint que permite reanudar la extracción en caso de interrupción.

5. **Cálculo de descriptores RDKit**:  
   A partir del `InChI` (o alternativamente del `CanonicalSMILES`) se calculan más de 30 descriptores moleculares usando RDKit, entre ellos:
   - **Clasificación estructural**: `RDKit_AlcoholClass` (primario, secundario, terciario, diol, poliol, aromático, etc.) y `RDKit_NumOH` (número de grupos hidroxilo).
   - **Fisicoquímicos**: `RDKit_MolLogP`, `RDKit_MolWt`, `RDKit_NumHeavyAtoms`, `RDKit_NumRotatableBonds`, `RDKit_NumHDonors`, `RDKit_NumHAcceptors`, `RDKit_TPSA`, `RDKit_MolMR`.
   - **Topológicos**: índices de conectividad (`RDKit_Chi0v`, `RDKit_Chi1v`, `RDKit_Chi2v`), índices de forma (`RDKit_Kappa1`, `RDKit_Kappa2`, `RDKit_Kappa3`), `RDKit_LabuteASA`, `RDKit_BalabanJ`.
   - **Recuento de anillos**: `RDKit_NumRings`, `RDKit_NumAromaticRings`, `RDKit_NumAliphaticRings`, `RDKit_NumSaturatedRings`, `RDKit_NumHeterocycles`, y desglose por tipo de anillo (carbociclos/heterociclos saturados/aromáticos/alifáticos).
   - **Otros**: `RDKit_FractionCsp3`, `RDKit_NumValenceElectrons`, `RDKit_NumRadicalElectrons`, `RDKit_QED` (drug-likeness), `RDKit_FpDensityMorgan1/2/3` (densidad de fingerprints).

6. **Limpieza y filtrado**:  
   - Se eliminan compuestos sin grupos hidroxilo según PubChem (`HBondDonorCount = 0`) o según RDKit (`RDKit_NumOH = 0`), garantizando que solo se conserven alcoholes verdaderos.
   - Se descartan columnas que resulten completamente nulas (sin información para ningún compuesto).

7. **División en subconjuntos temáticos**:  
   Para facilitar el análisis, la base final se divide en cuatro archivos CSV, cada uno con un grupo específico de propiedades (siempre incluyendo `CID` para poder recombinarlos):
   - **Identificadores**: `alcoholes_identificadores.csv`
   - **Fisicoquímicas**: `alcoholes_fisicoquimicas.csv`
   - **Topológicas**: `alcoholes_topologicas.csv`
   - **Anillos y aromaticidad**: `alcoholes_anillos.csv`

## 📁 Estructura de los datos

Los archivos generados se encuentran en la carpeta `bases_divididas/`. Cada archivo incluye la columna `CID` como clave primaria.

### 1. `alcoholes_identificadores.csv`

| Columna | Descripción |
|---------|-------------|
| CID | PubChem Compound ID |
| MolecularFormula | Fórmula molecular |
| SMILES | Representación SMILES canónica |
| ConnectivitySMILES | SMILES de conectividad (sin estereoquímica) |
| InChI | International Chemical Identifier |
| InChIKey | Hash del InChI |
| IUPACName | Nombre sistemático IUPAC |
| Category | Categoría asignada en la búsqueda original |
| Category_Description | Descripción textual de la categoría |
| Extraction_Date | Fecha de extracción de PubChem |
| RDKit_AlcoholClass | Clasificación estructural según RDKit |

### 2. `alcoholes_fisicoquimicas.csv`

| Columna | Descripción |
|---------|-------------|
| CID | ... |
| MolecularWeight | Peso molecular (PubChem) |
| XLogP | Coeficiente de partición octanol/agua (PubChem) |
| ExactMass | Masa exacta (PubChem) |
| MonoisotopicMass | Masa del isótopo más ligero (PubChem) |
| TPSA | Área polar superficial topológica (PubChem) |
| HBondDonorCount | Número de donores de H (PubChem) |
| HBondAcceptorCount | Número de aceptores de H (PubChem) |
| RotatableBondCount | Enlaces simples rotables (PubChem) |
| HeavyAtomCount | Átomos pesados (PubChem) |
| Volume3D | Volumen molecular 3D (PubChem) |
| RDKit_MolLogP | LogP calculado por RDKit |
| RDKit_MolWt | Peso molecular calculado por RDKit |
| RDKit_NumHeavyAtoms | Átomos pesados (RDKit) |
| RDKit_NumRotatableBonds | Enlaces rotables (RDKit) |
| RDKit_NumHDonors | Donores de H (RDKit) |
| RDKit_NumHAcceptors | Aceptores de H (RDKit) |
| RDKit_TPSA | TPSA calculada por RDKit |
| RDKit_NumOH | Número de grupos hidroxilo (RDKit) |
| RDKit_MolMR | Refractividad molar (RDKit) |
| RDKit_QED | Drug-likeness (Quantitative Estimate of Drug-likeness) |
| RDKit_FpDensityMorgan1/2/3 | Densidad de fingerprints Morgan (radio 1,2,3) |

### 3. `alcoholes_topologicas.csv`

| Columna | Descripción |
|---------|-------------|
| CID | ... |
| Complexity | Complejidad molecular (PubChem) |
| Charge | Carga formal neta |
| IsotopeAtomCount | Número de átomos isotópicos no estándar |
| AtomStereoCount | Centros estereogénicos totales |
| DefinedAtomStereoCount | Centros con estereoquímica definida |
| UndefinedAtomStereoCount | Centros con estereoquímica no definida |
| BondStereoCount | Enlaces con isomería geométrica |
| DefinedBondStereoCount | Enlaces con estereoquímica definida |
| UndefinedBondStereoCount | Enlaces con estereoquímica no definida |
| CovalentUnitCount | Unidades covalentes (sales/complejos) |
| RDKit_Chi0v, RDKit_Chi1v, RDKit_Chi2v | Índices de conectividad de Kier–Hall (con valencia) |
| RDKit_Kappa1, RDKit_Kappa2, RDKit_Kappa3 | Índices de forma de Kier |
| RDKit_LabuteASA | Área superficial accesible aproximada (Labute) |
| RDKit_BalabanJ | Índice de Balaban |
| RDKit_NumAtoms | Número total de átomos |
| RDKit_NumValenceElectrons | Electrones de valencia |
| RDKit_NumRadicalElectrons | Electrones desapareados (radicales) |
| RDKit_FractionCsp3 | Fracción de carbonos con hibridación sp3 |

### 4. `alcoholes_anillos.csv`

| Columna | Descripción |
|---------|-------------|
| CID | ... |
| RDKit_NumRings | Número total de anillos |
| RDKit_NumAromaticRings | Anillos aromáticos |
| RDKit_NumAliphaticRings | Anillos alifáticos |
| RDKit_NumSaturatedRings | Anillos saturados |
| RDKit_NumHeterocycles | Anillos con heteroátomos |
| RDKit_NumSaturatedCarbocycles | Carbociclos saturados |
| RDKit_NumSaturatedHeterocycles | Heterociclos saturados |
| RDKit_NumAromaticCarbocycles | Carbociclos aromáticos |
| RDKit_NumAromaticHeterocycles | Heterociclos aromáticos |
| RDKit_NumAliphaticCarbocycles | Carbociclos alifáticos |
| RDKit_NumAliphaticHeterocycles | Heterociclos alifáticos |

## 🔧 Requisitos e instalación

Para ejecutar los scripts y trabajar con los datos necesitas:

- Python ≥ 3.8
- Librerías: `pandas`, `numpy`, `requests`, `rdkit`, `matplotlib`, `seaborn`, `scikit-learn`, `tqdm`, `openpyxl` (opcional para Excel)

Instalación rápida:

```bash
pip install pandas numpy requests matplotlib seaborn scikit-learn tqdm openpyxl
pip install rdkit-pypi

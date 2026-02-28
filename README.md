# Base de Datos de Alcoholes – Proyecto de Quimioinformática

Este repositorio contiene una base de datos curada de **alcoholes** (compuestos con grupo –OH unido a carbono saturado), construida a partir de consultas a la API de PubChem y enriquecida con descriptores calculados con RDKit. El objetivo es proporcionar un recurso estructurado para estudios de relación estructura–propiedad, modelado QSAR, clasificación y exploración química.

## 📌 Objetivos

- Obtener una colección representativa de alcoholes (~1000 compuestos) distribuidos en categorías estructurales (primarios, secundarios, terciarios, dioles, polioles, grasos, aromáticos).
- Extraer propiedades fisicoquímicas y topológicas desde PubChem.
- Calcular descriptores moleculares adicionales con RDKit (índices de conectividad, forma, recuento de anillos, etc.).
- Organizar la información en subconjuntos temáticos para facilitar análisis específicos.
- Sentar las bases para ampliar el proyecto a otros grupos funcionales (aldehídos, cetonas, ácidos, etc.).

## 🧪 Metodología

1. **Búsqueda de compuestos**: se utilizaron patrones SMILES (y listas predefinidas para terciarios) para recuperar CIDs de PubChem mediante el endpoint `/compound/fastsubstructure/smiles`.
2. **Extracción de propiedades básicas**: se consultaron 26 propiedades fisicoquímicas y topológicas desde el endpoint `/compound/cid/property`.
3. **Cálculo de descriptores RDKit**: a partir del SMILES canónico se generaron:
   - Propiedades equivalentes a PubChem (para cubrir posibles nulos): MolLogP, MolWt, NumHeavyAtoms, etc.
   - Descriptores únicos: índices de conectividad (Chi0, Chi1, Chi2), índices de forma (Kappa1-3), LabuteASA, BalabanJ, FractionCsp3, etc.
   - Recuento detallado de anillos (aromáticos, alifáticos, heterociclos).
4. **Clasificación con RDKit**: mediante patrones SMARTS se asignó una categoría adicional (`RDKit_Category`) para validar/contrastar con la categoría original.
5. **Limpieza**:
   - Eliminación de columnas totalmente nulas.
   - Filtrado de compuestos no alcoholes (`HBondDonorCount = 0`).
6. **División en subconjuntos**: se crearon cuatro archivos CSV temáticos para trabajar de forma modular.

## 📁 Estructura de los datos

Los archivos generados se encuentran en la carpeta `bases_divididas/` (o el nombre que hayas usado). Cada archivo incluye la columna `CID` para poder fusionarlos.

### 1. `alcoholes_identificadores.csv`
| Columna | Descripción |
|---------|-------------|
| CID | PubChem Compound ID (clave primaria) |
| MolecularFormula | Fórmula molecular |
| SMILES | Representación SMILES canónica |
| ConnectivitySMILES | SMILES de conectividad (sin estereoquímica) |
| InChI | International Chemical Identifier |
| InChIKey | Hash del InChI |
| IUPACName | Nombre sistemático IUPAC |
| Category | Categoría asignada en la búsqueda original |
| Category_Description | Descripción textual de la categoría |
| Extraction_Date | Fecha de extracción de PubChem |
| RDKit_Category | Categoría calculada con RDKit |

### 2. `alcoholes_fisicoquimicas.csv`
| Columna | Descripción |
|---------|-------------|
| CID | … |
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
| RDKit_MolLogP | LogP calculado por RDKit (Crippen) |
| RDKit_MolWt | Peso molecular calculado por RDKit |
| RDKit_NumHeavyAtoms | Átomos pesados (RDKit) |
| RDKit_NumRotatableBonds | Enlaces rotables (RDKit) |
| RDKit_NumHDonors | Donores de H (RDKit) |
| RDKit_NumHAcceptors | Aceptores de H (RDKit) |
| RDKit_TPSA | TPSA calculada por RDKit |

### 3. `alcoholes_topologicas.csv`
| Columna | Descripción |
|---------|-------------|
| CID | … |
| Complexity | Complejidad molecular (PubChem) |
| Charge | Carga formal neta |
| IsotopeAtomCount | Número de átomos isotópicos no estándar |
| AtomStereoCount | Centros estereogénicos totales |
| DefinedAtomStereoCount | Centros con estereoquímica definida |
| UndefinedAtomStereoCount | Centros con estereoquímica no definida |
| BondStereoCount | Enlaces con isomería geométrica |
| DefinedBondStereoCount | Enlaces con estereoquímica definida |
| UndefinedBondStereoCount | Enlaces con estereoquímica no definida |
| CovalentUnitCount | Unidades covalentes (para sales/complejos) |
| RDKit_Chi0, RDKit_Chi1, RDKit_Chi2 | Índices de conectividad de Kier–Hall |
| RDKit_Kappa1, RDKit_Kappa2, RDKit_Kappa3 | Índices de forma de Kier |
| RDKit_LabuteASA | Área superficial accesible aproximada (Labute) |
| RDKit_BalabanJ | Índice de Balaban |
| RDKit_NumAtoms | Número total de átomos |
| RDKit_NumValenceElectrons | Electrones de valencia |
| RDKit_NumRadicalElectrons | Electrones desapareados (radicales) |

### 4. `alcoholes_anillos.csv`
| Columna | Descripción |
|---------|-------------|
| CID | … |
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
- Librerías: `pandas`, `numpy`, `requests`, `rdkit`, `matplotlib`, `seaborn`, `scikit-learn`

Instalación rápida:

```bash
pip install pandas numpy requests matplotlib seaborn scikit-learn
pip install rdkit-pypi

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extractor de alcoholes desde PubChem con enriquecimiento RDKit.
Versión ULTRA ROBUSTA: reintentos con backoff, rate limiting, checkpoint,
logging profesional y corrección de descriptores (FractionCSP3).
"""

import requests
import time
import pandas as pd
import urllib.parse
import sqlite3
import logging
import json
import signal
import sys
from typing import List, Dict, Set, Optional, Any
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from datetime import datetime
from tqdm import tqdm
import os

# ── CONFIGURACIÓN DE LOGGING ─────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("extractor.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# ── CATEGORÍAS DE ALCOHOLES ────────────────────────────────────────
ALCOHOL_CATEGORIES = {
    "primary":       {"smarts": "[CH2][OH]",                           "quota": 30000},
    "secondary":     {"smarts": "[CH1]([OH])([#6])[#6]",               "quota": 25000},
    "tertiary":      {"smarts": "[C;X4;H0]([OH])([#6])([#6])[#6]",     "quota": 20000},
    "diol":          {"smarts": "[OH][CX4][CX4][OH]",                  "quota": 15000},
    "polyol":        {"smarts": "[OH][CX4][CX4][CX4][OH]",             "quota": 10000},
    "aromatic":      {"smarts": "c[OH]",                               "quota": 15000},
    "allylic":       {"smarts": "[OH][CX4][CX3]=[CX3]",                "quota": 10000},
    "fatty":         {"smarts": "CCCCCC[CH2][OH]",                     "quota": 10000},
    "cycloalcohol":  {"smarts": "[OH][C;R]",                           "quota": 10000},
    "unsaturated":   {"smarts": "[OH][C]=[C]",                         "quota": 5000},
    "general":       {"smarts": "[#6][OH]",                            "quota": 5000},
}

BASE_PROPERTIES = [
    "MolecularWeight", "MolecularFormula", "CanonicalSMILES", "IsomericSMILES",
    "InChI", "InChIKey", "IUPACName", "XLogP", "ExactMass", "MonoisotopicMass",
    "TPSA", "Complexity", "Charge", "HBondDonorCount", "HBondAcceptorCount",
    "RotatableBondCount", "HeavyAtomCount", "IsotopeAtomCount", "AtomStereoCount",
    "DefinedAtomStereoCount", "UndefinedAtomStereoCount", "BondStereoCount",
    "DefinedBondStereoCount", "UndefinedBondStereoCount", "CovalentUnitCount", "Volume3D",
]

NUMERIC_PROPS = {"MolecularWeight", "XLogP", "ExactMass", "MonoisotopicMass", "TPSA", "Volume3D"}
INT_PROPS = {"Charge", "HBondDonorCount", "HBondAcceptorCount", "RotatableBondCount",
             "HeavyAtomCount", "IsotopeAtomCount", "AtomStereoCount", "DefinedAtomStereoCount",
             "UndefinedAtomStereoCount", "BondStereoCount", "DefinedBondStereoCount",
             "UndefinedBondStereoCount", "CovalentUnitCount", "Complexity"}

BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
REQUEST_TIMEOUT = 90
MAX_RETRIES = 5
INITIAL_BACKOFF = 1

# ── WORKER RDKit (CORREGIDO) ────────────────────────────────────────
def _rdkit_worker_inchi(inchi: str) -> Dict[str, Any]:
    """
    Calcula descriptores RDKit a partir de InChI.
    Versión corregida: usa Chi0v, Chi1v, Chi2v (con valencia) y los demás descriptores según tu lista.
    """
    rdkit_columns = [
        "RDKit_AlcoholClass", "RDKit_NumOH",
        "RDKit_MolLogP", "RDKit_MolWt", "RDKit_NumHeavyAtoms",
        "RDKit_NumRotatableBonds", "RDKit_NumHDonors", "RDKit_NumHAcceptors",
        "RDKit_TPSA", "RDKit_NumStereocenters", "RDKit_MolMR",
        "RDKit_NumAtoms", "RDKit_NumRings", "RDKit_NumAromaticRings",
        "RDKit_NumAliphaticRings", "RDKit_NumSaturatedRings", "RDKit_NumHeterocycles",
        "RDKit_FractionCsp3", "RDKit_NumValenceElectrons", "RDKit_NumRadicalElectrons",
        "RDKit_Chi0v", "RDKit_Chi1v", "RDKit_Chi2v",
        "RDKit_Kappa1", "RDKit_Kappa2", "RDKit_Kappa3",
        "RDKit_LabuteASA", "RDKit_BalabanJ",
        "RDKit_NumSaturatedCarbocycles", "RDKit_NumSaturatedHeterocycles",
        "RDKit_NumAromaticCarbocycles", "RDKit_NumAromaticHeterocycles",
        "RDKit_NumAliphaticCarbocycles", "RDKit_NumAliphaticHeterocycles",
        "RDKit_QED",
        "RDKit_FpDensityMorgan1", "RDKit_FpDensityMorgan2", "RDKit_FpDensityMorgan3"
    ]
    
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski, MolSurf, QED
        from rdkit.Chem import rdMolDescriptors
        
        if not inchi or pd.isna(inchi) or not inchi.strip():
            return {col: None for col in rdkit_columns}
        
        mol = Chem.MolFromInchi(inchi)
        if mol is None:
            logger.warning(f"InChI inválido para RDKit: {inchi[:50]}...")
            return {col: None for col in rdkit_columns}
        
        # Contar grupos OH
        oh_pattern = Chem.MolFromSmarts("[OH]")
        n_oh = len(mol.GetSubstructMatches(oh_pattern)) if oh_pattern else 0
        
        # Clasificación informativa
        if mol.HasSubstructMatch(Chem.MolFromSmarts("c[OH]")):
            cls = "aromatic"
        elif mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2][OH]")) and n_oh == 1:
            cls = "primary"
        elif mol.HasSubstructMatch(Chem.MolFromSmarts("[CH1]([OH])")) and n_oh == 1:
            cls = "secondary"
        elif mol.HasSubstructMatch(Chem.MolFromSmarts("[C;H0]([OH])([#6])([#6])[#6]")) and n_oh == 1:
            cls = "tertiary"
        elif n_oh > 2:
            cls = "polyol"
        elif n_oh == 2:
            cls = "diol"
        elif n_oh == 0:
            cls = "no_alcohol"
        else:
            cls = "other"
        
        def safe_desc(func, *args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                logger.debug(f"Error en descriptor {func.__name__}: {e}")
                return None
        
        return {
            "RDKit_AlcoholClass": cls,
            "RDKit_NumOH": n_oh,
            "RDKit_MolLogP": safe_desc(Descriptors.MolLogP, mol),
            "RDKit_MolWt": safe_desc(Descriptors.MolWt, mol),
            "RDKit_NumHeavyAtoms": mol.GetNumHeavyAtoms(),
            "RDKit_NumRotatableBonds": safe_desc(Descriptors.NumRotatableBonds, mol),
            "RDKit_NumHDonors": safe_desc(Lipinski.NumHDonors, mol),
            "RDKit_NumHAcceptors": safe_desc(Lipinski.NumHAcceptors, mol),
            "RDKit_TPSA": safe_desc(MolSurf.TPSA, mol),
            "RDKit_NumStereocenters": safe_desc(rdMolDescriptors.CalcNumAtomStereoCenters, mol),
            "RDKit_MolMR": safe_desc(Descriptors.MolMR, mol),
            "RDKit_NumAtoms": mol.GetNumAtoms(),
            "RDKit_NumRings": safe_desc(rdMolDescriptors.CalcNumRings, mol),
            "RDKit_NumAromaticRings": safe_desc(rdMolDescriptors.CalcNumAromaticRings, mol),
            "RDKit_NumAliphaticRings": safe_desc(rdMolDescriptors.CalcNumAliphaticRings, mol),
            "RDKit_NumSaturatedRings": safe_desc(rdMolDescriptors.CalcNumSaturatedRings, mol),
            "RDKit_NumHeterocycles": safe_desc(rdMolDescriptors.CalcNumHeterocycles, mol),
            "RDKit_FractionCsp3": safe_desc(Descriptors.FractionCSP3, mol),
            "RDKit_NumValenceElectrons": safe_desc(Descriptors.NumValenceElectrons, mol),
            "RDKit_NumRadicalElectrons": safe_desc(Descriptors.NumRadicalElectrons, mol),
            # CORREGIDO: usamos las versiones con valencia (Chi0v, Chi1v, Chi2v)
            "RDKit_Chi0v": safe_desc(Descriptors.Chi0v, mol),
            "RDKit_Chi1v": safe_desc(Descriptors.Chi1v, mol),
            "RDKit_Chi2v": safe_desc(Descriptors.Chi2v, mol),
            "RDKit_Kappa1": safe_desc(Descriptors.Kappa1, mol),
            "RDKit_Kappa2": safe_desc(Descriptors.Kappa2, mol),
            "RDKit_Kappa3": safe_desc(Descriptors.Kappa3, mol),
            "RDKit_LabuteASA": safe_desc(Descriptors.LabuteASA, mol),
            "RDKit_BalabanJ": safe_desc(Descriptors.BalabanJ, mol),
            "RDKit_NumSaturatedCarbocycles": safe_desc(rdMolDescriptors.CalcNumSaturatedCarbocycles, mol),
            "RDKit_NumSaturatedHeterocycles": safe_desc(rdMolDescriptors.CalcNumSaturatedHeterocycles, mol),
            "RDKit_NumAromaticCarbocycles": safe_desc(rdMolDescriptors.CalcNumAromaticCarbocycles, mol),
            "RDKit_NumAromaticHeterocycles": safe_desc(rdMolDescriptors.CalcNumAromaticHeterocycles, mol),
            "RDKit_NumAliphaticCarbocycles": safe_desc(rdMolDescriptors.CalcNumAliphaticCarbocycles, mol),
            "RDKit_NumAliphaticHeterocycles": safe_desc(rdMolDescriptors.CalcNumAliphaticHeterocycles, mol),
            "RDKit_QED": safe_desc(QED.qed, mol),
            "RDKit_FpDensityMorgan1": safe_desc(Descriptors.FpDensityMorgan1, mol),
            "RDKit_FpDensityMorgan2": safe_desc(Descriptors.FpDensityMorgan2, mol),
            "RDKit_FpDensityMorgan3": safe_desc(Descriptors.FpDensityMorgan3, mol),
        }
    except Exception as e:
        logger.error(f"Error inesperado en worker RDKit: {e}")
        return {col: None for col in rdkit_columns}
# ── CLASE EXTRACTORA ─────────────────────────────────────────────────
class PubChemExtractorRobusto:
    """
    Clase que maneja la extracción desde PubChem, el almacenamiento en SQLite
    y el enriquecimiento con RDKit.
    """
    def __init__(self, db_path="alcoholes_robusto.db", calls_per_second=0.33,
                 batch_size=50, max_workers=4, checkpoint_file="checkpoint.json"):
        self.db_path = db_path
        self.batch_size = batch_size
        self.max_workers = max_workers
        self.calls_per_second = calls_per_second
        self.checkpoint_file = checkpoint_file
        self.processed = set()
        self.shutdown_flag = False
        
        signal.signal(signal.SIGINT, self._signal_handler)
        
        self._init_db()
        self._load_checkpoint()
    
    def _signal_handler(self, sig, frame):
        logger.info("⚠️  Interrupción recibida. Guardando checkpoint y saliendo...")
        self.shutdown_flag = True
        self._save_checkpoint()
        sys.exit(0)
    
    def _init_db(self):
        conn = sqlite3.connect(self.db_path)
        cols = ["cid INTEGER PRIMARY KEY", "category TEXT", "extraction_date TEXT"]
        for p in BASE_PROPERTIES:
            if p in NUMERIC_PROPS:
                cols.append(f"{p} REAL")
            elif p in INT_PROPS:
                cols.append(f"{p} INTEGER")
            else:
                cols.append(f"{p} TEXT")
        conn.execute(f"CREATE TABLE IF NOT EXISTS compounds ({', '.join(cols)})")
        conn.commit()
        conn.close()
        self.processed = self._load_processed_cids()
        logger.info(f"Base de datos inicializada. {len(self.processed)} compuestos ya procesados.")
    
    def _load_processed_cids(self) -> Set[int]:
        conn = sqlite3.connect(self.db_path)
        cids = {row[0] for row in conn.execute("SELECT cid FROM compounds").fetchall()}
        conn.close()
        return cids
    
    def _load_checkpoint(self):
        if os.path.exists(self.checkpoint_file):
            with open(self.checkpoint_file, 'r') as f:
                checkpoint = json.load(f)
                self.processed.update(checkpoint.get('processed_cids', []))
                logger.info(f"Checkpoint cargado: {len(checkpoint.get('processed_cids', []))} CIDs adicionales.")
    
    def _save_checkpoint(self):
        checkpoint = {
            'processed_cids': list(self.processed),
            'timestamp': datetime.now().isoformat()
        }
        with open(self.checkpoint_file, 'w') as f:
            json.dump(checkpoint, f)
        logger.info("Checkpoint guardado.")
    
    def _request_with_backoff(self, url: str, method='GET', **kwargs) -> Optional[Dict]:
        for attempt in range(MAX_RETRIES):
            if self.shutdown_flag:
                return None
            try:
                time.sleep(1.0 / self.calls_per_second)
                resp = requests.request(method, url, timeout=REQUEST_TIMEOUT, **kwargs)
                if resp.status_code == 200:
                    return resp.json()
                elif resp.status_code in (429, 503):
                    wait = INITIAL_BACKOFF * (2 ** attempt) + 0.5
                    logger.warning(f"Límite de tasa alcanzado (código {resp.status_code}). Esperando {wait:.2f}s...")
                    time.sleep(wait)
                else:
                    logger.error(f"Error HTTP {resp.status_code} en {url}")
                    return None
            except requests.exceptions.Timeout:
                logger.warning(f"Timeout en intento {attempt+1}/{MAX_RETRIES}")
                if attempt == MAX_RETRIES - 1:
                    return None
                time.sleep(INITIAL_BACKOFF * (2 ** attempt))
            except requests.exceptions.RequestException as e:
                logger.error(f"Error de red: {e}")
                if attempt == MAX_RETRIES - 1:
                    return None
                time.sleep(INITIAL_BACKOFF * (2 ** attempt))
        return None
    
    def fetch_cids_for_category(self, smarts: str, quota: int) -> List[int]:
        all_cids = []
        start = 0
        max_per_request = 5000
        encoded = urllib.parse.quote(smarts, safe='')
        logger.info(f"Buscando CIDs con SMARTS: {smarts}")
        while len(all_cids) < quota and not self.shutdown_flag:
            url = f"{BASE_URL}/compound/fastsubstructure/smarts/{encoded}/cids/JSON?MaxRecords={max_per_request}&Start={start}"
            data = self._request_with_backoff(url)
            if not data or 'IdentifierList' not in data:
                break
            cids = data['IdentifierList'].get('CID', [])
            if not cids:
                break
            all_cids.extend(cids)
            logger.info(f"  → Página: {len(all_cids):,} CIDs obtenidos")
            start += len(cids)
            if len(cids) < max_per_request:
                break
        return all_cids[:quota]
    
    def fetch_properties_batch(self, cid_list: List[int], category: str) -> List[Dict]:
        if not cid_list:
            return []
        cids_str = ','.join(str(c) for c in cid_list)
        props_str = ','.join(BASE_PROPERTIES)
        url = f"{BASE_URL}/compound/cid/{cids_str}/property/{props_str}/JSON"
        data = self._request_with_backoff(url, timeout=REQUEST_TIMEOUT*2)
        if not data or "PropertyTable" not in data:
            return []
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        rows = []
        for p in data["PropertyTable"].get("Properties", []):
            inchi = p.get("InChI", "")
            if not inchi or pd.isna(inchi) or not inchi.strip():
                logger.debug(f"CID {p.get('CID')} omitido por falta de InChI")
                continue
            row = dict(p)
            row["cid"] = row.pop("CID", None)
            row["category"] = category
            row["extraction_date"] = now
            rows.append(row)
        return rows
    
    def insert_batch_db(self, rows: List[Dict]) -> int:
        if not rows:
            return 0
        conn = sqlite3.connect(self.db_path)
        existing_cols = {row[1] for row in conn.execute("PRAGMA table_info(compounds)").fetchall()}
        clean = []
        for r in rows:
            if "cid" not in r:
                continue
            clean.append({k: v for k, v in r.items() if k in existing_cols})
        if not clean:
            conn.close()
            return 0
        cols = list(clean[0].keys())
        sql = f"INSERT OR IGNORE INTO compounds ({','.join(cols)}) VALUES ({','.join(['?']*len(cols))})"
        data = [tuple(r.get(c) for c in cols) for r in clean]
        cursor = conn.executemany(sql, data)
        inserted = cursor.rowcount
        conn.commit()
        conn.close()
        return inserted
    
    def process_category(self, cat: str, info: Dict):
        logger.info(f"{'─'*60}\nProcesando {cat} – {info['smarts']} (cuota {info['quota']:,})")
        conn = sqlite3.connect(self.db_path)
        cur = conn.execute("SELECT COUNT(*) FROM compounds WHERE category=?", (cat,))
        already = cur.fetchone()[0]
        conn.close()
        remaining = info["quota"] - already
        if remaining <= 0:
            logger.info(f"Categoría '{cat}' completa ({already:,}/{info['quota']:,})")
            return 0
        logger.info(f"Faltan {remaining:,} compuestos (ya tenemos {already:,})")
        all_cids = self.fetch_cids_for_category(info["smarts"], info["quota"])
        new_cids = [c for c in all_cids if c not in self.processed][:remaining]
        logger.info(f"Encontrados {len(all_cids):,} CIDs, {len(new_cids):,} nuevos")
        if not new_cids:
            return 0
        batches = [new_cids[i:i+self.batch_size] for i in range(0, len(new_cids), self.batch_size)]
        saved = 0
        failed = 0
        with ThreadPoolExecutor(max_workers=self.max_workers) as ex:
            futures = {ex.submit(self.fetch_properties_batch, b, cat): b for b in batches}
            with tqdm(total=len(batches), desc=f"   {cat}", unit="lote") as pbar:
                for fut in as_completed(futures):
                    if self.shutdown_flag:
                        ex.shutdown(wait=False)
                        return saved
                    try:
                        rows = fut.result()
                        if rows:
                            inserted = self.insert_batch_db(rows)
                            if inserted:
                                for r in rows:
                                    if r.get("cid"):
                                        self.processed.add(r["cid"])
                                saved += inserted
                        else:
                            failed += 1
                    except Exception as e:
                        logger.error(f"Error procesando lote: {e}")
                        failed += 1
                    pbar.update(1)
                    pbar.set_postfix(guardados=f"{saved:,}", fallos=failed)
        logger.info(f"Categoría {cat}: {saved} guardados, {failed} lotes fallidos")
        self._save_checkpoint()
        return saved
    
    def build(self) -> pd.DataFrame:
        logger.info("="*70)
        total_quota = sum(c["quota"] for c in ALCOHOL_CATEGORIES.values())
        logger.info(f"🎯 META TOTAL: {total_quota:,} alcoholes")
        logger.info("="*70)
        for i, (cat, info) in enumerate(ALCOHOL_CATEGORIES.items()):
            if self.shutdown_flag:
                break
            self.process_category(cat, info)
            if i < len(ALCOHOL_CATEGORIES)-1 and not self.shutdown_flag:
                logger.info("⏸️ Pausa de 30 segundos entre categorías...")
                time.sleep(30)
        logger.info(f"✅ COMPLETADO: {len(self.processed):,} compuestos únicos")
        return self.load_df()
    
    def load_df(self) -> pd.DataFrame:
        conn = sqlite3.connect(self.db_path)
        df = pd.read_sql_query("SELECT * FROM compounds", conn)
        conn.close()
        return df
    
    def save(self, df: pd.DataFrame, prefix: str = "alcoholes_robusto"):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        csv = f"{prefix}_{len(df)}_{timestamp}.csv"
        df.to_csv(csv, index=False, encoding="utf-8")
        logger.info(f"💾 CSV guardado: {csv}")
        try:
            import openpyxl
            excel = f"{prefix}_{len(df)}_{timestamp}.xlsx"
            df.to_excel(excel, index=False)
            logger.info(f"💾 Excel guardado: {excel}")
        except ImportError:
            logger.warning("openpyxl no instalado. No se generó archivo Excel.")
    
    def summary(self, df: pd.DataFrame):
        logger.info(f"\n📊 Resumen final: {len(df):,} compuestos")
        if "category" in df.columns:
            logger.info("\nPor categoría (cuota):")
            for cat, n in df["category"].value_counts().items():
                quota = ALCOHOL_CATEGORIES.get(cat, {}).get("quota", 0)
                pct = n / quota * 100 if quota else 0
                bar = "█" * int(pct / 5)
                logger.info(f"  {cat:<14} {n:>7,}/{quota:>6,} ({pct:5.1f}%) {bar}")
        rdkit_cols = [c for c in df.columns if c.startswith("RDKit_")]
        if rdkit_cols:
            logger.info(f"\n🧪 Columnas RDKit: {len(rdkit_cols)}")
            if "RDKit_AlcoholClass" in df.columns:
                logger.info("\nClasificación RDKit (según grupos OH):")
                counts = df["RDKit_AlcoholClass"].value_counts(dropna=False)
                for cls, cnt in counts.items():
                    logger.info(f"  {cls}: {cnt}")
        if "MolecularWeight" in df.columns:
            mw = pd.to_numeric(df["MolecularWeight"], errors="coerce")
            logger.info(f"\n⚖️  Peso molecular: media={mw.mean():.1f}  min={mw.min():.1f}  max={mw.max():.1f}")

# ── FUNCIÓN DE ENRIQUECIMIENTO RDKit ────────────────────────────────
def enhance_with_rdkit_from_inchi(df: pd.DataFrame, n_workers: int = 4, chunksize: int = 100,
                                  filter_non_alcohols: bool = True) -> pd.DataFrame:
    logger.info("="*50)
    logger.info("VERIFICACIÓN DE RDKit")
    logger.info("="*50)
    try:
        import rdkit
        logger.info(f"✅ RDKit versión: {rdkit.__version__}")
    except ImportError:
        logger.error("❌ RDKit no instalado. Se omite enriquecimiento.")
        return df
    
    posibles_inchi = ["InChI", "inchi"]
    col_inchi = None
    for col in posibles_inchi:
        if col in df.columns:
            col_inchi = col
            break
    
    if not col_inchi:
        logger.error("❌ No se encontró columna InChI en el DataFrame.")
        return df
    
    logger.info(f"✅ Usando columna InChI: '{col_inchi}'")
    inchi_list = df[col_inchi].fillna("").tolist()
    total_con_datos = sum(1 for i in inchi_list if i.strip())
    logger.info(f"Filas con InChI no vacío: {total_con_datos}/{len(df)}")
    
    logger.info(f"🧪 Calculando descriptores RDKit con {n_workers} procesos...")
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        results = list(tqdm(ex.map(_rdkit_worker_inchi, inchi_list, chunksize=chunksize),
                            total=len(df), desc="RDKit"))
    
    no_vacios = sum(1 for r in results if any(v is not None for v in r.values()))
    logger.info(f"✅ Filas con al menos un descriptor calculado: {no_vacios}/{len(results)}")
    
    rdkit_df = pd.DataFrame(results)
    if rdkit_df.empty:
        logger.warning("⚠️ DataFrame RDKit vacío. No se añaden columnas.")
        return df
    
    for col in rdkit_df.columns:
        if col in df.columns:
            logger.debug(f"Columna {col} ya existe en df, se sobrescribirá.")
            df = df.drop(columns=[col])
    
    df_concat = pd.concat([df.reset_index(drop=True), rdkit_df], axis=1)
    logger.info("✅ Enriquecimiento RDKit completado.")
    
    if filter_non_alcohols and "RDKit_NumOH" in df_concat.columns:
        non_oh = (df_concat["RDKit_NumOH"] == 0).sum()
        if non_oh > 0:
            logger.info(f"🔍 Eliminando {non_oh} compuestos sin grupos OH (no alcoholes según RDKit)...")
            df_concat = df_concat[df_concat["RDKit_NumOH"] > 0].copy()
            logger.info(f"   Quedan {len(df_concat)} compuestos.")
    
    return df_concat

# ── EJECUCIÓN PRINCIPAL ────────────────────────────────────────────
if __name__ == "__main__":
    DB_PATH = "alcoholes_robusto.db"
    BATCH_SIZE = 50
    FETCH_WORKERS = 4
    RDKIT_WORKERS = 4
    CALLS_PER_SECOND = 0.33
    
    extractor = PubChemExtractorRobusto(
        db_path=DB_PATH,
        batch_size=BATCH_SIZE,
        max_workers=FETCH_WORKERS,
        calls_per_second=CALLS_PER_SECOND,
        checkpoint_file="checkpoint.json"
    )
    
    try:
        df = extractor.build()
    except KeyboardInterrupt:
        logger.info("⏸️ Interrupción por usuario. Cargando datos existentes...")
        df = extractor.load_df()
    except Exception as e:
        logger.exception(f"❌ Error inesperado durante la extracción: {e}")
        df = extractor.load_df()
    
    if not df.empty:
        if 'HBondDonorCount' in df.columns:
            non_alcohols_pubchem = (df['HBondDonorCount'] == 0).sum()
            if non_alcohols_pubchem > 0:
                logger.info(f"🔍 Eliminando {non_alcohols_pubchem} compuestos con HBondDonorCount = 0 (PubChem)...")
                df = df[df['HBondDonorCount'] > 0].copy()
                logger.info(f"   Quedan {len(df)} compuestos.")
        else:
            logger.warning("No se encuentra columna HBondDonorCount. No se puede filtrar por PubChem.")
        
        logger.info(f"📄 DataFrame base: {df.shape[0]} filas, {df.shape[1]} columnas")
        
        if not any(c.startswith("RDKit_") for c in df.columns):
            logger.info("🧪 Aplicando enriquecimiento RDKit usando InChI...")
            df = enhance_with_rdkit_from_inchi(df, n_workers=RDKIT_WORKERS, chunksize=100, filter_non_alcohols=True)
        else:
            logger.info("✅ Ya existen columnas RDKit en el DataFrame.")
        
        extractor.save(df)
        extractor.summary(df)
    else:
        logger.warning("No hay datos para guardar.")

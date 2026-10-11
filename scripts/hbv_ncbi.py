"""Helper akses NCBI E-utilities + logika kurasi rekaman HBV.

Catatan desain (hasil probe empiris, lihat docs/LAPORAN_PENGUMPULAN_DATASET.md):
- Qualifier `gene` TIDAK andal pada rekaman HBV (gene="P" hanya ada di ~39%
  sampel). Qualifier `product` jauh lebih konsisten ("polymerase" 99%,
  "large S protein" 100%). Ekstraksi gen karena itu berbasis `product`
  dengan `gene` sebagai cadangan.
- Genotipe HBV tidak ada di field terstruktur. Pada 100% sampel probe,
  genotipe tertulis di qualifier `note` pada feature `source`.
- Banyak CDS HBV memakai compound location join(x..end, 1..y) karena genom
  HBV sirkular. Ekstraksi WAJIB memakai feature.extract(), bukan slicing.
"""
from __future__ import annotations

import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field

from defusedxml import ElementTree as DET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UA = {"User-Agent": "hbv-dataset-collector/1.0 (tugas kuliah JST)"}

# NCBI: maksimum 3 request/detik tanpa API key.
_MIN_INTERVAL = 0.4
_last_call = [0.0]


def _throttle() -> None:
    delta = time.monotonic() - _last_call[0]
    if delta < _MIN_INTERVAL:
        time.sleep(_MIN_INTERVAL - delta)
    _last_call[0] = time.monotonic()


def _post(endpoint: str, params: dict, retries: int = 4) -> str:
    data = urllib.parse.urlencode(params).encode()
    last_err = None
    for attempt in range(retries):
        _throttle()
        try:
            req = urllib.request.Request(EUTILS + endpoint, data=data, headers=UA)
            with urllib.request.urlopen(req, timeout=300) as resp:
                return resp.read().decode("utf-8", "replace")
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last_err = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"gagal akses {endpoint}: {last_err}")


def esearch_count(term: str) -> int:
    root = DET.fromstring(_post("esearch.fcgi", {
        "db": "nuccore", "term": term, "retmax": 0, "retmode": "xml"}))
    node = root.find("Count")
    return int(node.text) if node is not None and node.text else 0


def esearch_ids(term: str, retmax: int) -> list[str]:
    """Ambil UID. Diurutkan by accession agar hasil reprodusibel."""
    root = DET.fromstring(_post("esearch.fcgi", {
        "db": "nuccore", "term": term, "retmax": retmax,
        "retmode": "xml", "sort": "accession"}))
    return [e.text for e in root.findall(".//IdList/Id") if e.text]


def efetch_genbank(uids: list[str]) -> str:
    return _post("efetch.fcgi", {
        "db": "nuccore", "id": ",".join(uids),
        "rettype": "gbwithparts", "retmode": "text"})


# --------------------------------------------------------------------------
# Parsing genotipe
# --------------------------------------------------------------------------
# "subgenotype C2" memuat substring "genotype C2" -> ikut tertangkap sebagai C.
# Itu memang diinginkan: subgenotipe mengandung huruf genotipe induknya.
_GENOTYPE_RE = re.compile(r"genotype\s*[:=]?\s*([A-Ja-j])(?![A-Za-z])")
_RECOMB_RE = re.compile(r"recombinant|genotype\s*[:=]?\s*[A-Ja-j]\s*/\s*[A-Ja-j]",
                        re.IGNORECASE)

VALID_GENOTYPES = tuple("ABCDEFGH")


def parse_genotype(notes: list[str], description: str) -> tuple[str | None, str]:
    """Kembalikan (genotipe, alasan). genotipe None bila tidak dapat dipastikan."""
    blob = " ; ".join(notes + [description])
    if _RECOMB_RE.search(blob):
        return None, "rekombinan/antar-genotipe"
    letters = {m.group(1).upper() for m in _GENOTYPE_RE.finditer(blob)}
    if not letters:
        return None, "genotipe tidak tercantum"
    if len(letters) > 1:
        return None, f"genotipe ambigu ({'/'.join(sorted(letters))})"
    g = letters.pop()
    if g not in VALID_GENOTYPES:
        return None, f"genotipe di luar cakupan A-H ({g})"
    return g, ""


def parse_subgenotype(notes: list[str]) -> str:
    m = re.search(r"sub(?:genotype|type)\s*[:=]?\s*([A-Ja-j]\d+)", " ; ".join(notes))
    return m.group(1).upper() if m else ""


# --------------------------------------------------------------------------
# Identifikasi & ekstraksi CDS gen S dan gen P
# --------------------------------------------------------------------------
# Gen S HBV = satu ORF dengan tiga start codon in-frame (preS1/preS2/S) yang
# menghasilkan protein L, M, S. Region gen S yang utuh = CDS "large S protein".
S_PRODUCTS = ("large s protein", "large surface protein", "pres1/pres2/s",
              "large envelope protein", "l protein")
P_PRODUCTS = ("polymerase", "dna polymerase", "p protein",
              "reverse transcriptase", "polymerase/reverse transcriptase")

# Rentang panjang yang wajar (bp). Gen S utuh ~1203, gen P ~2532.
S_LEN_RANGE = (1100, 1310)
P_LEN_RANGE = (2400, 2700)


def _norm(values: list[str]) -> list[str]:
    return [v.strip().lower() for v in values]


def find_cds(record, kind: str):
    """Cari CDS gen S atau gen P. Bila ada beberapa kandidat, ambil terpanjang
    (untuk S: membedakan 'large S' dari 'middle S'/'S' yang lebih pendek)."""
    products = S_PRODUCTS if kind == "S" else P_PRODUCTS
    gene_names = ("s", "pres/s", "pres1/pres2/s") if kind == "S" else ("p", "pol")
    cands = []
    for f in record.features:
        if f.type != "CDS":
            continue
        prods = _norm(f.qualifiers.get("product", []))
        genes = _norm(f.qualifiers.get("gene", []))
        hit = any(p in products for p in prods) or any(g in gene_names for g in genes)
        if hit:
            cands.append(f)
    if not cands:
        return None
    return max(cands, key=lambda f: len(f))


def is_partial(feature) -> bool:
    return "<" in str(feature.location) or ">" in str(feature.location)


def span_positions(feature) -> set[int]:
    pos = set()
    for part in feature.location.parts:
        pos.update(range(int(part.start), int(part.end)))
    return pos


@dataclass
class Curated:
    """Hasil pemeriksaan satu rekaman genom."""
    accession: str = ""
    version: str = ""
    organism: str = ""
    genotype: str = ""
    subgenotype: str = ""
    genome_len: int = 0
    topology: str = ""
    has_S: bool = False
    has_P: bool = False
    S_len: int = 0
    P_len: int = 0
    S_loc: str = ""
    P_loc: str = ""
    S_seq: str = ""
    P_seq: str = ""
    S_partial: bool = False
    P_partial: bool = False
    S_nonacgt: int = 0
    P_nonacgt: int = 0
    overlap_frac: float = 0.0
    country: str = ""
    collection_date: str = ""
    url: str = ""
    eligible: bool = False
    reasons: list[str] = field(default_factory=list)


_NON_ACGT = re.compile(r"[^ACGT]")


def curate(record) -> Curated:
    c = Curated()
    c.version = record.id
    c.accession = record.id.split(".")[0]
    c.genome_len = len(record.seq)
    c.topology = record.annotations.get("topology", "")
    c.organism = record.annotations.get("organism", "")
    c.url = f"https://www.ncbi.nlm.nih.gov/nuccore/{c.version}"

    src_notes: list[str] = []
    for f in record.features:
        if f.type == "source":
            src_notes += f.qualifiers.get("note", [])
            c.country = (f.qualifiers.get("geo_loc_name")
                         or f.qualifiers.get("country") or [""])[0]
            c.collection_date = f.qualifiers.get("collection_date", [""])[0]

    genotype, why = parse_genotype(src_notes, record.description)
    c.subgenotype = parse_subgenotype(src_notes)
    if genotype:
        c.genotype = genotype
    else:
        c.reasons.append(why)

    s_cds = find_cds(record, "S")
    p_cds = find_cds(record, "P")
    c.has_S, c.has_P = s_cds is not None, p_cds is not None
    if not c.has_S:
        c.reasons.append("anotasi gen S tidak ada")
    if not c.has_P:
        c.reasons.append("anotasi gen P tidak ada")

    for cds, tag in ((s_cds, "S"), (p_cds, "P")):
        if cds is None:
            continue
        partial = is_partial(cds)
        seq = str(cds.extract(record.seq)).upper()
        bad = len(_NON_ACGT.findall(seq))
        lo, hi = S_LEN_RANGE if tag == "S" else P_LEN_RANGE
        setattr(c, f"{tag}_seq", seq)
        setattr(c, f"{tag}_len", len(seq))
        setattr(c, f"{tag}_loc", str(cds.location))
        setattr(c, f"{tag}_partial", partial)
        setattr(c, f"{tag}_nonacgt", bad)
        if partial:
            c.reasons.append(f"CDS gen {tag} parsial")
        if not (lo <= len(seq) <= hi):
            c.reasons.append(f"panjang gen {tag} di luar rentang ({len(seq)} bp)")
        if bad:
            c.reasons.append(f"gen {tag} memuat {bad} basa non-ACGT")

    if s_cds is not None and p_cds is not None:
        sset, pset = span_positions(s_cds), span_positions(p_cds)
        if sset:
            c.overlap_frac = len(sset & pset) / len(sset)

    c.eligible = not c.reasons
    return c

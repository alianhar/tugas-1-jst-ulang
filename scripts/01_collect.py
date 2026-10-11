"""Tahap 1 — pencarian & pengunduhan rekaman genom HBV dari NCBI Nucleotide.

Strategi pencarian (per genotipe, bukan satu query global) dipilih karena
genotipe HBV tidak tersimpan di field terstruktur NCBI; satu query global
"Hepatitis B virus complete genome" menghasilkan 13.097 rekaman yang sangat
didominasi genotipe C dan D, sehingga kelas minoritas (G, H) hampir tidak
akan terambil.

Untuk menjaga reprodusibilitas sekaligus menghindari bias era submission
(accession lama vs baru), skrip mengambil daftar UID yang besar lalu
mengambil subsampel deterministik dengan seed tetap.

Output: data/raw_gb/<ACCESSION>.gb  (berkas GenBank asli, tidak diubah)
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from hbv_ncbi import VALID_GENOTYPES, efetch_genbank, esearch_count, esearch_ids

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw_gb"
RAW.mkdir(parents=True, exist_ok=True)

# Dua strategi pencarian dipakai dan digabung (union), diterapkan SERAGAM ke
# delapan genotipe. Kriteria inklusi tidak boleh berbeda antar kelas, karena
# itu akan menciptakan artefak teknis yang bisa "dicurangi" model klasifikasi.
#
# Uji empiris (lihat laporan): term "judul" memberi rekaman dengan mutu
# anotasi JAUH lebih baik (proporsi rekaman punya anotasi gen S dan P lebih
# tinggi), sedangkan term "luas" memberi jangkauan lebih besar dan menolong
# kelas minoritas. Keduanya digabung agar dapat kedua manfaat itu.
TERMS = {
    # Kata kunci awal sesuai arahan tugas.
    "judul": ('"Hepatitis B virus"[Organism] AND "complete genome"[Title] '
              'AND 3000:3400[SLEN]'),
    # Lebih luas: genom HBV sepanjang ~3,0-3,4 kb tanpa syarat frasa judul.
    "luas": ('"Hepatitis B virus"[Organism] AND 3000:3400[SLEN] '
             'AND biomol_genomic[PROP] NOT patent[Properties]'),
}
ID_POOL_PER_GENOTYPE = 900   # hanya UID (murah), untuk dasar subsampel
CANDIDATES_PER_GENOTYPE = 150  # per term, per genotipe
SEED = 20251011
FETCH_BATCH = 100


def main() -> None:
    rng = random.Random(SEED)
    manifest: dict[str, dict] = {}
    all_uids: list[str] = []
    seen: set[str] = set()

    print("== Tahap 1a: inventarisasi ketersediaan per genotipe ==")
    print(f"  {'gt':>2} | {'judul':>18} | {'luas':>18} | kandidat kumulatif")
    for g in VALID_GENOTYPES:
        manifest[g] = {}
        line = f"  {g:>2} |"
        for label, base in TERMS.items():
            term = f'{base} AND "genotype {g}"[All Fields]'
            total = esearch_count(term)
            ids = esearch_ids(term, retmax=ID_POOL_PER_GENOTYPE)
            take = min(CANDIDATES_PER_GENOTYPE, len(ids))
            picked = sorted(rng.sample(ids, take), key=int)
            manifest[g][label] = {"total_di_ncbi": total,
                                  "uid_terambil": len(ids),
                                  "uid_disubsampel": take, "term": term}
            for u in picked:
                if u not in seen:
                    seen.add(u)
                    all_uids.append(u)
            line += f" {total:6d} ada/{take:3d} ambil |"
        print(f"{line} {len(all_uids)}")

    # Lanjutkan unduhan: lewati UID yang berkasnya sudah ada dari sesi lalu.
    log_path = ROOT / "data" / "uid_terunduh.json"
    done = set(json.loads(log_path.read_text())) if log_path.exists() else set()
    todo = [u for u in all_uids if u not in done]
    print(f"\ntotal UID unik: {len(all_uids)} | sudah ada: "
          f"{len(all_uids) - len(todo)} | perlu diunduh: {len(todo)}")

    print("\n== Tahap 1b: pengunduhan GenBank ==")
    written, batches = 0, 0
    for i in range(0, len(todo), FETCH_BATCH):
        batch = todo[i:i + FETCH_BATCH]
        raw = efetch_genbank(batch)
        batches += 1
        for block in _split_genbank(raw):
            acc = _accession_of(block)
            if not acc:
                continue
            (RAW / f"{acc}.gb").write_text(block, encoding="utf-8")
            written += 1
        done.update(batch)
        print(f"  batch {batches}: {len(batch)} UID diminta, "
              f"kumulatif {written} berkas .gb baru")

    log_path.write_text(json.dumps(sorted(done)), encoding="utf-8")
    (ROOT / "data" / "manifest_pencarian.json").write_text(
        json.dumps({"seed": SEED, "terms": TERMS,
                    "candidates_per_genotype_per_term": CANDIDATES_PER_GENOTYPE,
                    "per_genotipe": manifest,
                    "uid_kandidat_total": len(all_uids),
                    "berkas_gb_baru": written,
                    "berkas_gb_total": len(list(RAW.glob('*.gb')))},
                   indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nselesai: {len(list(RAW.glob('*.gb')))} berkas GenBank di {RAW}")


def _split_genbank(raw: str) -> list[str]:
    """Pecah stream GenBank jadi blok per rekaman, byte asli dipertahankan."""
    blocks, cur = [], []
    for line in raw.splitlines(keepends=True):
        cur.append(line)
        if line.rstrip("\r\n") == "//":
            blocks.append("".join(cur))
            cur = []
    return blocks


def _accession_of(block: str) -> str | None:
    for line in block.splitlines():
        if line.startswith("VERSION"):
            return line.split()[1].split(".")[0]
        if line.startswith("ACCESSION"):
            parts = line.split()
            if len(parts) > 1:
                return parts[1]
    return None


if __name__ == "__main__":
    main()

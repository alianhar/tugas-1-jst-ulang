"""Tahap 2 — pemeriksaan, penyaringan, dan penyusunan dataset.

Alur:
  baca semua .gb -> kurasi per rekaman -> deduplikasi -> pembatasan kuota per
  genotipe -> tulis metadata (CSV+XLSX) + tiga berkas FASTA dataset.

Kriteria KELAYAKAN sebuah rekaman (semua harus terpenuhi):
  1. Genotipe dapat dipastikan tunggal dan termasuk A-H
  2. Bukan rekombinan antar-genotipe (label kelasnya tidak tunggal)
  3. Anotasi CDS gen S dan gen P keduanya ada
  4. Kedua CDS lengkap (bukan parsial / tidak ada tanda < atau >)
  5. Panjang gen S 1100-1310 bp dan gen P 2400-2700 bp
  6. Sekuens kedua gen hanya memuat basa A/C/G/T (tanpa N atau IUPAC ambigu)
  7. Pasangan sekuens (S,P) belum pernah muncul di rekaman lain (unik)
"""
from __future__ import annotations

import hashlib
import json
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
from hbv_ncbi import VALID_GENOTYPES, curate

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw_gb"
PROC = ROOT / "data" / "processed"
META = ROOT / "metadata"
PROC.mkdir(parents=True, exist_ok=True)
META.mkdir(parents=True, exist_ok=True)

TARGET_PER_GENOTYPE = 45   # batas atas agar kelas tidak terlalu tidak seimbang
MIN_TOTAL = 300            # ketentuan tugas
SEED = 20251011
WRAP = 60


def main() -> None:
    files = sorted(RAW.glob("*.gb"))
    print(f"== Tahap 2a: kurasi {len(files)} rekaman ==")

    curated, parse_fail = [], 0
    for path in files:
        try:
            rec = next(SeqIO.parse(path, "genbank"))
        except Exception:
            parse_fail += 1
            continue
        curated.append(curate(rec))
    print(f"  terparse: {len(curated)} | gagal parse: {parse_fail}")

    # --- deduplikasi ------------------------------------------------------
    # Keunikan dituntut pada gen S saja, gen P saja, DAN pasangannya.
    # Alasan: deduplikasi hanya pada pasangan (S,P) menyisakan rekaman dengan
    # gen S identik tetapi gen P berbeda. Rekaman seperti itu menjadi duplikat
    # persis di Dataset 1 (gen S saja) dan berisiko membocorkan data antara
    # himpunan latih dan uji. Karena ketiga dataset harus memakai himpunan
    # rekaman yang sama agar perbandingan skenario tidak terkonfounding,
    # keunikan ditegakkan untuk ketiganya sekaligus.
    print("\n== Tahap 2b: deduplikasi (gen S, gen P, dan pasangannya) ==")
    seen: dict[str, dict[str, str]] = {"S": {}, "P": {}, "SP": {}}
    counts = Counter()
    for c in sorted(curated, key=lambda x: x.accession):
        if not c.eligible:
            continue
        keys = {"S": c.S_seq, "P": c.P_seq, "SP": f"{c.S_seq}|{c.P_seq}"}
        for tag, payload in keys.items():
            h = hashlib.sha256(payload.encode()).hexdigest()
            if h in seen[tag]:
                label = {"S": "gen S", "P": "gen P", "SP": "pasangan gen S+P"}[tag]
                c.eligible = False
                c.reasons.append(f"{label} duplikat dengan {seen[tag][h]}")
                counts[tag] += 1
                break
        else:
            for tag, payload in keys.items():
                seen[tag][hashlib.sha256(payload.encode()).hexdigest()] = c.version
    print(f"  dibuang karena gen S duplikat          : {counts['S']}")
    print(f"  dibuang karena gen P duplikat          : {counts['P']}")
    print(f"  dibuang karena pasangan S+P duplikat   : {counts['SP']}")

    # --- statistik penolakan ---------------------------------------------
    reason_counter = Counter()
    for c in curated:
        for r in c.reasons:
            reason_counter[r.split(" (")[0].split(" dengan ")[0]] += 1
    print("\n== Tahap 2c: sebab penolakan (satu rekaman bisa >1 sebab) ==")
    for r, n in reason_counter.most_common():
        print(f"  {n:4d}  {r}")

    eligible = [c for c in curated if c.eligible]
    by_g = defaultdict(list)
    for c in eligible:
        by_g[c.genotype].append(c)
    print(f"\n  rekaman layak sebelum kuota: {len(eligible)}")
    for g in VALID_GENOTYPES:
        print(f"    genotipe {g}: {len(by_g[g])}")

    # --- pembatasan kuota per genotipe ------------------------------------
    rng = random.Random(SEED)
    selected = []
    for g in VALID_GENOTYPES:
        pool = sorted(by_g[g], key=lambda x: x.accession)
        if len(pool) > TARGET_PER_GENOTYPE:
            pool = sorted(rng.sample(pool, TARGET_PER_GENOTYPE),
                          key=lambda x: x.accession)
        selected.extend(pool)
    chosen = {c.version for c in selected}
    print(f"\n== Tahap 2d: terpilih {len(selected)} rekaman "
          f"(kuota {TARGET_PER_GENOTYPE}/genotipe) ==")
    dist = Counter(c.genotype for c in selected)
    for g in VALID_GENOTYPES:
        print(f"    genotipe {g}: {dist[g]}")
    if len(selected) < MIN_TOTAL:
        print(f"  !! PERINGATAN: total {len(selected)} < ketentuan {MIN_TOTAL}")

    # --- metadata ---------------------------------------------------------
    rows = []
    for c in sorted(curated, key=lambda x: (x.genotype or "Z", x.accession)):
        rows.append({
            "accession": c.accession,
            "accession_version": c.version,
            "nama_organisme": c.organism,
            "genotipe": c.genotype or "tidak dapat dipastikan",
            "subgenotipe": c.subgenotype,
            "panjang_genom_bp": c.genome_len,
            "topologi": c.topology,
            "anotasi_gen_S": "Ada" if c.has_S else "Tidak ada",
            "panjang_gen_S_bp": c.S_len or "",
            "lokasi_gen_S": c.S_loc,
            "anotasi_gen_P": "Ada" if c.has_P else "Tidak ada",
            "panjang_gen_P_bp": c.P_len or "",
            "lokasi_gen_P": c.P_loc,
            "fraksi_S_overlap_P": round(c.overlap_frac, 4) if c.has_S and c.has_P else "",
            "negara": c.country,
            "tanggal_koleksi": c.collection_date,
            "tautan_ncbi": c.url,
            "status_kelayakan": "Layak" if c.eligible else "Tidak layak",
            "alasan_tidak_layak": "; ".join(c.reasons),
            "masuk_dataset_final": "Ya" if c.version in chosen else "Tidak",
        })
    df = pd.DataFrame(rows)
    df.to_csv(META / "metadata_rekaman.csv", index=False, encoding="utf-8-sig")
    with pd.ExcelWriter(META / "metadata_rekaman.xlsx", engine="openpyxl") as xw:
        df.to_excel(xw, sheet_name="semua_rekaman", index=False)
        df[df.masuk_dataset_final == "Ya"].to_excel(
            xw, sheet_name="dataset_final", index=False)
    print(f"\n  metadata ditulis: {len(df)} baris "
          f"-> metadata/metadata_rekaman.{{csv,xlsx}}")

    # --- berkas FASTA -----------------------------------------------------
    def fasta_block(c, gene: str) -> str:
        seq = c.S_seq if gene == "S" else c.P_seq
        head = f">{c.version}|Genotype={c.genotype}|Gene={gene}"
        body = "\n".join(seq[i:i + WRAP] for i in range(0, len(seq), WRAP))
        return f"{head}\n{body}\n"

    d1 = "".join(fasta_block(c, "S") for c in selected)
    d2 = "".join(fasta_block(c, "P") for c in selected)
    d3 = "".join(fasta_block(c, "S") + fasta_block(c, "P") for c in selected)
    (PROC / "dataset1_genS.fasta").write_text(d1, encoding="utf-8")
    (PROC / "dataset2_genP.fasta").write_text(d2, encoding="utf-8")
    (PROC / "dataset3_genS_genP.fasta").write_text(d3, encoding="utf-8")
    print(f"  dataset1_genS.fasta        : {len(selected)} sekuens")
    print(f"  dataset2_genP.fasta        : {len(selected)} sekuens")
    print(f"  dataset3_genS_genP.fasta   : {len(selected) * 2} sekuens "
          f"({len(selected)} pasang)")

    # --- ringkasan untuk laporan -----------------------------------------
    s_lens = [c.S_len for c in selected]
    p_lens = [c.P_len for c in selected]
    overlaps = [c.overlap_frac for c in selected]
    summary = {
        "rekaman_diunduh": len(files),
        "rekaman_layak": len(eligible),
        "rekaman_dataset_final": len(selected),
        "distribusi_genotipe": {g: dist[g] for g in VALID_GENOTYPES},
        "sebab_penolakan": dict(reason_counter.most_common()),
        "gen_S_panjang": {"min": min(s_lens), "max": max(s_lens),
                          "modus": Counter(s_lens).most_common(1)[0][0]},
        "gen_P_panjang": {"min": min(p_lens), "max": max(p_lens),
                          "modus": Counter(p_lens).most_common(1)[0][0]},
        "overlap_S_dalam_P": {
            "min": round(min(overlaps), 4), "max": round(max(overlaps), 4),
            "rerata": round(sum(overlaps) / len(overlaps), 4),
            "jumlah_overlap_penuh": sum(1 for o in overlaps if o >= 0.999)},
        "target_per_genotipe": TARGET_PER_GENOTYPE,
        "seed": SEED,
    }
    (ROOT / "data" / "ringkasan_kurasi.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    print("\n== Statistik kunci ==")
    print(f"  panjang gen S : {summary['gen_S_panjang']}")
    print(f"  panjang gen P : {summary['gen_P_panjang']}")
    print(f"  overlap S di dalam P : {summary['overlap_S_dalam_P']}")


if __name__ == "__main__":
    main()

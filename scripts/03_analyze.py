"""Tahap 3 — kendali mutu (QC) dan analisis dataset.

Fokus pemeriksaan:
  1. Integritas berkas FASTA (jumlah, pasangan S-P pada dataset 3)
  2. Distribusi kelas & implikasinya untuk pelatihan JST
  3. Analisis overlap gen S di dalam gen P secara posisional (bp, bukan hanya
     fraksi) -> menguji apakah Dataset 3 benar-benar menambah informasi
     nukleotida dibanding Dataset 2
  4. Keutuhan biologis CDS (kelipatan 3, start/stop codon, stop internal)
  5. Kemiripan antar sekuens (risiko kebocoran data train/test)
"""
from __future__ import annotations

import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).parent))
from hbv_ncbi import find_cds, span_positions

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw_gb"
PROC = ROOT / "data" / "processed"
META = ROOT / "metadata"

out: dict = {}


def head(title: str) -> None:
    print(f"\n{'=' * 68}\n{title}\n{'=' * 68}")


# ---------------------------------------------------------------------------
head("1. Integritas berkas FASTA")
sets = {}
for name in ("dataset1_genS", "dataset2_genP", "dataset3_genS_genP",
             "dataset4_genP_nonS"):
    recs = list(SeqIO.parse(PROC / f"{name}.fasta", "fasta"))
    sets[name] = recs
    print(f"  {name}.fasta: {len(recs)} sekuens, "
          f"{sum(len(r.seq) for r in recs):,} bp total")

d1_ids = {r.id.split("|")[0] for r in sets["dataset1_genS"]}
d2_ids = {r.id.split("|")[0] for r in sets["dataset2_genP"]}
d3 = defaultdict(set)
for r in sets["dataset3_genS_genP"]:
    acc, _, gene = r.id.split("|")
    d3[acc].add(gene.split("=")[1])

unpaired = [a for a, genes in d3.items() if genes != {"S", "P"}]
d4_ids = {r.id.split("|")[0] for r in sets["dataset4_genP_nonS"]}
print(f"  accession identik antara dataset1 & dataset2 : {d1_ids == d2_ids}")
print(f"  accession dataset3 == dataset1               : {set(d3) == d1_ids}")
print(f"  accession dataset4 == dataset1               : {d4_ids == d1_ids}")
print(f"  accession di dataset3 tanpa pasangan S+P     : {len(unpaired)}")
dup_hdr = [k for k, v in Counter(r.id for r in sets["dataset3_genS_genP"]).items() if v > 1]
print(f"  header duplikat di dataset3                  : {len(dup_hdr)}")
out["integritas"] = {"dataset1": len(sets["dataset1_genS"]),
                     "dataset2": len(sets["dataset2_genP"]),
                     "dataset3": len(sets["dataset3_genS_genP"]),
                     "pasangan_tidak_lengkap": len(unpaired),
                     "header_duplikat": len(dup_hdr)}

# ---------------------------------------------------------------------------
head("2. Distribusi kelas (genotipe)")
dist = Counter(r.id.split("|")[1].split("=")[1] for r in sets["dataset1_genS"])
total = sum(dist.values())
for g in sorted(dist):
    n = dist[g]
    print(f"  genotipe {g}: {n:3d}  ({n / total * 100:5.2f}%)  {'#' * n}")
imb = max(dist.values()) / min(dist.values())
print(f"\n  total {total} | kelas {len(dist)} | rasio ketidakseimbangan "
      f"maks:min = {imb:.2f}")
print(f"  akurasi baseline 'tebak kelas mayoritas' = "
      f"{max(dist.values()) / total * 100:.2f}%")
out["distribusi_kelas"] = dict(sorted(dist.items()))
out["rasio_ketidakseimbangan"] = round(imb, 3)
out["baseline_mayoritas_persen"] = round(max(dist.values()) / total * 100, 2)

# ---------------------------------------------------------------------------
head("3. Overlap posisional gen S di dalam gen P")
rows = []
for acc in sorted(d1_ids):
    path = RAW / f"{acc.split('.')[0]}.gb"
    rec = next(SeqIO.parse(path, "genbank"))
    s, p = find_cds(rec, "S"), find_cds(rec, "P")
    sp, pp = span_positions(s), span_positions(p)
    rows.append({"acc": acc, "S_bp": len(sp), "P_bp": len(pp),
                 "irisan_bp": len(sp & pp), "union_bp": len(sp | pp),
                 "S_diluar_P_bp": len(sp - pp)})
ov = pd.DataFrame(rows)
print(ov[["S_bp", "P_bp", "irisan_bp", "union_bp", "S_diluar_P_bp"]].describe()
      .loc[["min", "50%", "max"]].to_string())
n_full = int((ov.S_diluar_P_bp == 0).sum())
print(f"\n  rekaman dengan gen S SELURUHNYA di dalam gen P : {n_full}/{len(ov)}")
print(f"  rekaman dengan union(S,P) == P                 : "
      f"{int((ov.union_bp == ov.P_bp).sum())}/{len(ov)}")
print(f"  rerata basa gen S yang BERADA DI LUAR gen P    : "
      f"{ov.S_diluar_P_bp.mean():.2f} bp")
print(f"\n  => Dataset 3 total bp = {ov.S_bp.sum() + ov.P_bp.sum():,}")
print(f"     informasi nukleotida unik (union)          = {ov.union_bp.sum():,}")
print(f"     redundansi Dataset 3 vs Dataset 2          = "
      f"{ov.S_bp.sum() / (ov.S_bp.sum() + ov.P_bp.sum()) * 100:.2f}% bp terduplikasi")
# Uji Dataset 4: panjang gen S + panjang region P non-S harus == panjang gen P
len_s = {r.id.split("|")[0]: len(r.seq) for r in sets["dataset1_genS"]}
len_p = {r.id.split("|")[0]: len(r.seq) for r in sets["dataset2_genP"]}
len_pn = {r.id.split("|")[0]: len(r.seq) for r in sets["dataset4_genP_nonS"]}
konsisten = sum(1 for a in len_s if len_s[a] + len_pn[a] == len_p[a])
print(f"\n  Dataset 4 -- uji |gen S| + |gen P non-S| == |gen P| : "
      f"{konsisten}/{len(len_s)}")
print(f"  panjang region gen P non-S: min {min(len_pn.values())} | "
      f"maks {max(len_pn.values())} bp")
out["dataset4"] = {"uji_panjang_konsisten": f"{konsisten}/{len(len_s)}",
                   "panjang_min": min(len_pn.values()),
                   "panjang_maks": max(len_pn.values())}

out["overlap"] = {
    "rekaman_S_penuh_di_dalam_P": n_full, "total_rekaman": len(ov),
    "union_sama_dengan_P": int((ov.union_bp == ov.P_bp).sum()),
    "rerata_S_diluar_P_bp": round(float(ov.S_diluar_P_bp.mean()), 3),
    "total_bp_dataset3": int(ov.S_bp.sum() + ov.P_bp.sum()),
    "total_bp_unik_union": int(ov.union_bp.sum()),
    "persen_bp_redundan": round(float(ov.S_bp.sum() /
                                      (ov.S_bp.sum() + ov.P_bp.sum()) * 100), 2),
}

# ---------------------------------------------------------------------------
head("4. Keutuhan biologis CDS terekstraksi")
bio = Counter()
for name, gene in (("dataset1_genS", "S"), ("dataset2_genP", "P")):
    for r in sets[name]:
        seq = str(r.seq).upper()
        bio[f"{gene}_total"] += 1
        if len(seq) % 3 == 0:
            bio[f"{gene}_kelipatan3"] += 1
        if seq.startswith("ATG"):
            bio[f"{gene}_start_ATG"] += 1
        prot = str(Seq(seq[:len(seq) // 3 * 3]).translate())
        if prot.endswith("*"):
            bio[f"{gene}_stop_diakhir"] += 1
        if "*" in prot[:-1]:
            bio[f"{gene}_ada_stop_internal"] += 1
for gene in ("S", "P"):
    t = bio[f"{gene}_total"]
    print(f"  gen {gene} (n={t}): kelipatan 3 = {bio[f'{gene}_kelipatan3']}, "
          f"start ATG = {bio[f'{gene}_start_ATG']}, "
          f"stop di akhir = {bio[f'{gene}_stop_diakhir']}, "
          f"ada stop internal = {bio[f'{gene}_ada_stop_internal']}")
out["keutuhan_cds"] = dict(bio)

# ---------------------------------------------------------------------------
head("5. Kemiripan sekuens (risiko kebocoran train/test)")
for name, gene in (("dataset1_genS", "S"), ("dataset2_genP", "P")):
    seqs = [(r.id.split("|")[0], str(r.seq)) for r in sets[name]]
    exact = len(seqs) - len({s for _, s in seqs})
    # k-mer Jaccard sebagai proksi kemiripan yang murah (k=12)
    K = 12
    profiles = [(a, {s[i:i + K] for i in range(len(s) - K + 1)}) for a, s in seqs]
    high = 0
    worst = (0.0, "", "")
    for i in range(len(profiles)):
        ai, pi = profiles[i]
        for j in range(i + 1, len(profiles)):
            aj, pj = profiles[j]
            jac = len(pi & pj) / len(pi | pj)
            if jac > 0.95:
                high += 1
            if jac > worst[0]:
                worst = (jac, ai, aj)
    print(f"  gen {gene}: sekuens identik persis = {exact}, "
          f"pasangan Jaccard-12mer > 0.95 = {high}")
    print(f"           pasangan paling mirip: {worst[1]} vs {worst[2]} "
          f"(Jaccard={worst[0]:.4f})")
    out.setdefault("kemiripan", {})[gene] = {
        "identik_persis": exact, "pasangan_jaccard_di_atas_0_95": high,
        "jaccard_maksimum": round(worst[0], 4),
        "pasangan_termirip": [worst[1], worst[2]]}

# ---------------------------------------------------------------------------
head("6. Rekapitulasi sebab penolakan (dinormalisasi)")
df = pd.read_csv(META / "metadata_rekaman.csv")
norm = Counter()
for val in df.alasan_tidak_layak.dropna():
    for part in str(val).split("; "):
        part = re.sub(r"\d+", "N", part)
        part = re.sub(r"dengan \S+", "dengan rekaman lain", part)
        norm[part] += 1
for k, v in norm.most_common():
    print(f"  {v:4d}  {k}")
print(f"\n  rekaman diperiksa : {len(df)}")
print(f"  layak             : {int((df.status_kelayakan == 'Layak').sum())}")
print(f"  masuk dataset     : {int((df.masuk_dataset_final == 'Ya').sum())}")
out["sebab_penolakan"] = dict(norm.most_common())
out["jumlah_diperiksa"] = len(df)
out["jumlah_layak"] = int((df.status_kelayakan == "Layak").sum())
out["jumlah_dataset_final"] = int((df.masuk_dataset_final == "Ya").sum())

(ROOT / "data" / "laporan_qc.json").write_text(
    json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\nhasil QC lengkap -> data/laporan_qc.json")

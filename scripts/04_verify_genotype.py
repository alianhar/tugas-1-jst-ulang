"""Tahap 4 — verifikasi genotipe secara independen (tidak memakai label NCBI).

Label genotipe pada tahap 2 diambil dari qualifier `source`->`/note`, yaitu
anotasi yang diisi penyetor rekaman. Tahap ini mengujinya dengan dua metode
yang TIDAK melihat catatan itu sama sekali:

  A. Penetapan rujukan terdekat (nearest reference) berbasis penyejajaran
     global sekuens gen S terhadap sekuens rujukan tiap genotipe.
  B. Pohon filogenetik Neighbor-Joining dari matriks jarak, lalu uji apakah
     sampel tiap genotipe mengelompok bersama rujukan genotipenya
     (uji monofili / kemurnian klad).

Set rujukan diambil dari DUA publikasi berbeda agar hasil tidak bergantung
pada satu sumber:
  set1 - Tabel rujukan terpilih berbasis panel WHO (PMC6475771)
  set2 - Tabel rujukan genotipe/subgenotipe (PMC11496132)

Catatan kejujuran metode: rujukan ini sendiri berlabel dari literatur, bukan
dari penentuan de novo. Yang diuji di sini adalah apakah label tiap sampel
konsisten dengan posisi sekuensnya relatif terhadap rujukan yang diakui
komunitas -- itu jauh lebih kuat daripada mempercayai catatan penyetor, tetapi
bukan penentuan genotipe dari nol.
"""
from __future__ import annotations

import io
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align, SeqIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

sys.path.insert(0, str(Path(__file__).parent))
from hbv_ncbi import S_LEN_RANGE, VALID_GENOTYPES, _post, curate

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
REF_CACHE = ROOT / "data" / "rujukan_genotipe.gb"

REFS = {
    "set1": dict(A="EU594392.1", B="EF473976.1", C="EU410079.1", D="EU787447.1",
                 E="AM494708.1", F="AY179734.1", G="AP007264.1", H="AY090460.1"),
    "set2": dict(A="HE974362.1", B="AB981583.1", C="LC456132.1", D="HE815465.1",
                 E="HE974384.1", F="DQ823095.1", G="AB625342.1", H="AB275308.1"),
}
KMER = 7
out: dict = {}


def head(t: str) -> None:
    print(f"\n{'=' * 68}\n{t}\n{'=' * 68}")


# ---------------------------------------------------------------------------
head("1. Penyiapan sekuens rujukan")
accs = sorted({a for s in REFS.values() for a in s.values()})
if REF_CACHE.exists():
    raw = REF_CACHE.read_text(encoding="utf-8")
else:
    raw = _post("efetch.fcgi", {"db": "nuccore", "id": ",".join(accs),
                                "rettype": "gbwithparts", "retmode": "text"})
    REF_CACHE.write_text(raw, encoding="utf-8")
refrecs = {r.id: r for r in SeqIO.parse(io.StringIO(raw), "genbank")}

ref_S: dict[str, tuple[str, str]] = {}   # nama -> (genotipe, sekuens gen S)
label_check = []
for setname, mapping in REFS.items():
    for g, acc in mapping.items():
        rec = refrecs.get(acc)
        if rec is None:
            print(f"  {setname} {g} {acc}: tidak terunduh, dilewati")
            continue
        c = curate(rec)
        label_check.append((acc, g, c.genotype))
        lo, hi = S_LEN_RANGE
        if c.S_len and lo <= c.S_len <= hi:
            ref_S[f"REF_{setname}_{g}_{acc}"] = (g, c.S_seq)
        else:
            print(f"  {setname} {g} {acc}: gen S tidak terpakai "
                  f"(panjang {c.S_len or 0} bp di luar rentang)")

cocok = sum(1 for _, g, n in label_check if g == n)
print(f"\n  rujukan dipakai untuk gen S : {len(ref_S)} dari {len(accs)}")
print(f"  genotipe literatur == catatan NCBI pada rujukan : "
      f"{cocok}/{len(label_check)}")
cov = Counter(g for g, _ in ref_S.values())
print(f"  cakupan rujukan per genotipe: "
      f"{ {g: cov[g] for g in VALID_GENOTYPES} }")
missing = [g for g in VALID_GENOTYPES if cov[g] == 0]
if missing:
    print(f"  !! genotipe tanpa rujukan gen S: {missing}")
out["rujukan"] = {"dipakai_gen_S": len(ref_S),
                  "label_literatur_vs_ncbi": f"{cocok}/{len(label_check)}",
                  "cakupan_per_genotipe": {g: cov[g] for g in VALID_GENOTYPES}}

# ---------------------------------------------------------------------------
head("2. Pemuatan sampel & uji tumpang tindih dengan rujukan")
samples = []
for r in SeqIO.parse(PROC / "dataset1_genS.fasta", "fasta"):
    acc, gt, _ = r.id.split("|")
    samples.append((acc, gt.split("=")[1], str(r.seq)))
print(f"  sampel gen S: {len(samples)}")
sample_accs = {a for a, _, _ in samples}
bentrok = sample_accs & set(accs)
print(f"  accession rujukan yang juga ada di dataset: {len(bentrok)} "
      f"{sorted(bentrok) if bentrok else ''}")
if bentrok:
    print("  !! sampel tersebut akan dikecualikan agar tidak cocok dgn dirinya")
out["bentrok_rujukan_sampel"] = sorted(bentrok)

# ---------------------------------------------------------------------------
head("3. Metode A -- rujukan terdekat via penyejajaran global gen S")
aligner = Align.PairwiseAligner(mode="global", match_score=1,
                                mismatch_score=-1, open_gap_score=-5,
                                extend_gap_score=-1)
ref_items = sorted(ref_S.items())
print(f"  {len(samples)} sampel x {len(ref_items)} rujukan = "
      f"{len(samples) * len(ref_items)} penyejajaran, mohon tunggu...")

assign, conf = [], []
for n, (acc, gt, seq) in enumerate(samples, 1):
    if acc in bentrok:
        continue
    scores = [(aligner.score(seq, rseq), rg, rname)
              for rname, (rg, rseq) in ref_items]
    scores.sort(reverse=True)
    best_g = scores[0][1]
    # margin = selisih skor terhadap rujukan terbaik dari genotipe LAIN
    other = next((s for s, g, _ in scores if g != best_g), None)
    assign.append((acc, gt, best_g, scores[0][0], scores[0][2]))
    conf.append((scores[0][0] - other) if other is not None else float("nan"))
    if n % 50 == 0:
        print(f"    {n}/{len(samples)}")

setuju = sum(1 for _, gt, pred, _, _ in assign if gt == pred)
print(f"\n  KONKORDANSI: {setuju}/{len(assign)} "
      f"({setuju / len(assign) * 100:.2f}%) label NCBI == rujukan terdekat")
margins = np.array([m for m in conf if not np.isnan(m)])
print(f"  margin skor thd genotipe lain: min={margins.min():.0f} "
      f"median={np.median(margins):.0f} maks={margins.max():.0f}")

per_g = defaultdict(lambda: [0, 0])
for _, gt, pred, _, _ in assign:
    per_g[gt][1] += 1
    if gt == pred:
        per_g[gt][0] += 1
print("\n  konkordansi per genotipe:")
for g in VALID_GENOTYPES:
    ok, tot = per_g[g]
    if tot:
        print(f"    {g}: {ok}/{tot} ({ok / tot * 100:.1f}%)")

beda = [(a, gt, pred, round(m)) for (a, gt, pred, _, _), m
        in zip(assign, conf) if gt != pred]
if beda:
    print(f"\n  {len(beda)} sampel TIDAK konkordan (accession, label, prediksi, margin):")
    for row in beda[:25]:
        print(f"    {row}")
out["metode_A"] = {
    "konkordansi": f"{setuju}/{len(assign)}",
    "persen": round(setuju / len(assign) * 100, 2),
    "per_genotipe": {g: f"{per_g[g][0]}/{per_g[g][1]}" for g in VALID_GENOTYPES
                     if per_g[g][1]},
    "margin_min": float(margins.min()), "margin_median": float(np.median(margins)),
    "tidak_konkordan": [list(map(str, r)) for r in beda],
}

# ---------------------------------------------------------------------------
head("4. Metode B -- pohon Neighbor-Joining & uji monofili")


def kmer_vec(seq: str, k: int = KMER) -> dict[str, int]:
    c = Counter(seq[i:i + k] for i in range(len(seq) - k + 1))
    return c


names, labels, vecs = [], [], []
for acc, gt, seq in samples:
    if acc in bentrok:
        continue
    names.append(acc)
    labels.append(gt)
    vecs.append(kmer_vec(seq))
for rname, (rg, rseq) in ref_items:
    names.append(rname)
    labels.append(rg)
    vecs.append(kmer_vec(rseq))

vocab = sorted({km for v in vecs for km in v})
idx = {km: i for i, km in enumerate(vocab)}
M = np.zeros((len(vecs), len(vocab)), dtype=np.float32)
for r, v in enumerate(vecs):
    for km, n in v.items():
        M[r, idx[km]] = n
M /= np.linalg.norm(M, axis=1, keepdims=True)
D = 1.0 - M @ M.T          # jarak kosinus atas frekuensi k-mer
np.fill_diagonal(D, 0.0)
D = np.clip(D, 0.0, None)
print(f"  matriks jarak {D.shape}, k-mer k={KMER}, dimensi {len(vocab)}")

dm = DistanceMatrix(names=list(names),
                    matrix=[[float(D[i][j]) for j in range(i + 1)]
                            for i in range(len(names))])
tree = DistanceTreeConstructor().nj(dm)
tree.root_at_midpoint()
label_of = dict(zip(names, labels))

print("\n  uji kemurnian klad (MRCA seluruh sampel berlabel genotipe X):")
mono = {}
for g in VALID_GENOTYPES:
    anggota = [n for n in names if label_of[n] == g]
    if len(anggota) < 2:
        continue
    terms = [t for t in tree.get_terminals() if t.name in anggota]
    clade = tree.common_ancestor(terms)
    inside = [t.name for t in clade.get_terminals()]
    asing = [n for n in inside if label_of[n] != g]
    purity = len(anggota) / len(inside)
    ref_in = [n for n in inside if n.startswith("REF_") and label_of[n] == g]
    mono[g] = {"anggota": len(anggota), "isi_klad": len(inside),
               "asing": len(asing), "kemurnian": round(purity, 4),
               "rujukan_di_dalam": len(ref_in)}
    tanda = "MONOFILETIK" if not asing else f"{len(asing)} sampel asing"
    print(f"    {g}: {len(anggota):3d} anggota, klad memuat {len(inside):3d}, "
          f"kemurnian {purity:.3f}, rujukan di dalam {len(ref_in)} -> {tanda}")
out["metode_B"] = {"k": KMER, "monofili": mono}

n_mono = sum(1 for v in mono.values() if v["asing"] == 0)
print(f"\n  genotipe yang monofiletik sempurna: {n_mono}/{len(mono)}")
print("  (catatan: metrik MRCA sangat rapuh -- SATU sampel menyimpang membuat")
print("   klad MRCA membesar hingga nyaris seluruh pohon, sehingga kemurnian")
print("   jatuh drastis tanpa berarti genotipe itu tersebar luas)")

# --- metrik tahan-outlier 1: kesepakatan label tetangga terdekat (1-NN) ---
head("5. Metrik tahan-outlier")
Dz = D.copy()
np.fill_diagonal(Dz, np.inf)
nn = Dz.argmin(axis=1)
lab = np.array(labels)
agree = lab[nn] == lab
print(f"  1-NN: tetangga terdekat berlabel sama = "
      f"{agree.sum()}/{len(lab)} ({agree.mean() * 100:.2f}%)")
per_nn = {}
for g in VALID_GENOTYPES:
    m = lab == g
    if m.sum():
        per_nn[g] = f"{int(agree[m].sum())}/{int(m.sum())}"
        print(f"    {g}: {per_nn[g]} ({agree[m].mean() * 100:.1f}%)")

# --- metrik tahan-outlier 2: klad murni terbesar yang memuat rujukan -------
print("\n  klad MURNI terbesar yang memuat rujukan genotipe tsb:")
pure = {}
for g in VALID_GENOTYPES:
    anggota = {n for n in names if label_of[n] == g}
    if len(anggota) < 2:
        continue
    best = 0
    for clade in tree.find_clades():
        isi = [t.name for t in clade.get_terminals()]
        if len(isi) < 2 or not all(label_of[n] == g for n in isi):
            continue
        if any(n.startswith("REF_") for n in isi) and len(isi) > best:
            best = len(isi)
    frac = best / len(anggota)
    pure[g] = {"klad_murni_terbesar": best, "anggota": len(anggota),
               "cakupan": round(frac, 4)}
    print(f"    {g}: {best:3d}/{len(anggota):3d} anggota ({frac * 100:5.1f}%) "
          f"berada dalam satu klad murni bersama rujukannya")

out["metode_B"]["satu_nn"] = {"total": f"{int(agree.sum())}/{len(lab)}",
                              "persen": round(float(agree.mean() * 100), 2),
                              "per_genotipe": per_nn}
out["metode_B"]["klad_murni"] = pure

head("6. Penulisan hasil per sampel")

nn_label = dict(zip(names, lab[nn]))
rows = []
for (acc, gt, pred, skor, refname), m in zip(assign, conf):
    rows.append({
        "accession_version": acc,
        "genotipe_label_ncbi": gt,
        "genotipe_rujukan_terdekat": pred,
        "konkordan": "Ya" if gt == pred else "TIDAK",
        "skor_penyejajaran": round(skor),
        "margin_thd_genotipe_lain": None if np.isnan(m) else round(m),
        "rujukan_terdekat": refname,
        "genotipe_tetangga_terdekat_1nn": nn_label.get(acc, ""),
    })
vdf = pd.DataFrame(rows)
vdf.to_csv(ROOT / "metadata" / "verifikasi_genotipe.csv", index=False,
           encoding="utf-8-sig")
xlsx = ROOT / "metadata" / "metadata_rekaman.xlsx"
with pd.ExcelWriter(xlsx, engine="openpyxl", mode="a",
                    if_sheet_exists="replace") as xw:
    vdf.to_excel(xw, sheet_name="verifikasi_genotipe", index=False)
print(f"  {len(vdf)} baris -> metadata/verifikasi_genotipe.csv")
print(f"  sheet 'verifikasi_genotipe' ditambahkan ke {xlsx.name}")

# --- daftar pengecualian berbasis DUA sinyal independen -------------------
# Sebuah rekaman hanya dikecualikan bila:
#   (1) publikasi asalnya menyebut rekombinasi -- penanda tingkat rekaman, DAN
#   (2) verifikasi rujukan-terdekat menolak label genotipenya.
# Satu sinyal saja tidak cukup: (1) sendiri terbukti terlalu kasar (ikut
# membuang rekaman yang rekombinasinya di region core), sedangkan (2) sendiri
# bersifat sirkular -- membuang apa pun yang tidak sepakat dengan rujukan akan
# menjamin konkordansi tinggi tanpa makna.
meta = pd.read_csv(ROOT / "metadata" / "metadata_rekaman.csv")
flag = dict(zip(meta.accession_version, meta.indikasi_rekombinan_publikasi))
tolak = [acc for acc, gt, pred, _, _ in assign
         if gt != pred and flag.get(acc) == "Ya"]
hanya_pub = int((meta.indikasi_rekombinan_publikasi == "Ya").sum())
print(f"\n  ditandai rekombinan oleh publikasi saja : {hanya_pub} rekaman")
print(f"  tidak konkordan saja                    : {len(beda)} rekaman")
print(f"  KEDUANYA (dasar pengecualian)           : {len(tolak)} rekaman")

# Daftar pengecualian harus DIAKUMULASI antar lintasan, bukan ditimpa.
# Membuang satu rekaman membebaskan slot deduplikasi, sehingga rekaman lain
# dari seri yang sama bisa naik menjadi layak pada lintasan berikutnya --
# persis yang terjadi pada HQ231885.1 setelah HQ231877-884 dibuang. Karena
# itu lintasan 2 dan 4 diulang sampai tidak ada penambahan baru (konvergen).
exc_path = ROOT / "data" / "pengecualian_verifikasi.json"
lama = set()
if exc_path.exists():
    lama = set(json.loads(exc_path.read_text())["accession_version"])
gabung = sorted(lama | set(tolak))
baru = sorted(set(tolak) - lama)
print(f"  sudah ada dari lintasan sebelumnya      : {len(lama)}")
print(f"  BARU pada lintasan ini                  : {len(baru)} {baru}")
print(f"  total daftar pengecualian               : {len(gabung)}")
print(f"\n  >>> {'KONVERGEN -- tidak ada penambahan' if not baru else 'BELUM konvergen -- jalankan ulang tahap 2 lalu tahap 4'}")
exc_path.write_text(
    json.dumps({"kriteria": "indikasi rekombinan di publikasi DAN tidak "
                            "konkordan dengan rujukan terdekat",
                "catatan": "akumulatif antar lintasan",
                "accession_version": gabung},
               indent=2, ensure_ascii=False), encoding="utf-8")
out["pengecualian"] = {"ditandai_publikasi": hanya_pub,
                       "tidak_konkordan": len(beda),
                       "baru_lintasan_ini": baru,
                       "dikecualikan_total": gabung,
                       "konvergen": not baru}

Path(ROOT / "data" / "laporan_verifikasi_genotipe.json").write_text(
    json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
Phylo.write(tree, ROOT / "data" / "pohon_nj_genS.newick", "newick")
print("\nhasil -> data/laporan_verifikasi_genotipe.json")
print("pohon  -> data/pohon_nj_genS.newick")

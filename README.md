# Tugas 1 JST — Dataset Sekuens DNA Virus Hepatitis B (HBV)

Pengumpulan dan kurasi dataset sekuens DNA HBV dari NCBI Nucleotide untuk
tiga skenario pemodelan klasifikasi genotipe dengan Jaringan Syaraf Tiruan.

**Laporan lengkap: [`docs/LAPORAN_PENGUMPULAN_DATASET.md`](docs/LAPORAN_PENGUMPULAN_DATASET.md)**

## Hasil ringkas

| | |
|---|---|
| Rekaman diperiksa | 2.309 |
| Rekaman layak | 902 |
| **Rekaman dataset final** | **310** (minimal 300 ✓) |
| Kelas genotipe | 8 (A–H) |
| Duplikat persis | 0 |
| Keutuhan CDS | 310/310 lolos |

Distribusi: A–F masing-masing 45, G 24, H 16 (G dan H dibatasi ketersediaan
di NCBI, bukan oleh penyaringan).

## Berkas dataset

```
data/raw_gb/*.gb                        2.309 rekaman GenBank asli (berkas utama)
data/processed/dataset1_genS.fasta        310 sekuens gen S
data/processed/dataset2_genP.fasta        310 sekuens gen P
data/processed/dataset3_genS_genP.fasta   620 sekuens (310 pasang S+P, accession sama)
metadata/metadata_rekaman.xlsx          metadata + status kelayakan tiap rekaman
```

## Menjalankan ulang

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt
.venv/Scripts/python scripts/01_collect.py    # unduh dari NCBI
.venv/Scripts/python scripts/02_curate.py     # kurasi + bentuk dataset
.venv/Scripts/python scripts/03_analyze.py    # kendali mutu
```

## Catatan penting sebelum tahap pemodelan

1. **Gen S terletak seluruhnya di dalam gen P** (310/310 rekaman, pada
   *reading frame* berbeda). Dataset 3 karena itu tidak menambah cakupan
   nukleotida dibanding Dataset 2 — 32,13% isinya duplikasi. Perbandingan
   Dataset 2 vs 3 mengukur pengaruh representasi ganda, **bukan** penambahan
   informasi genomik. Rincian dan rumusan kesimpulan yang dapat
   dipertahankan ada di §6 laporan.
2. Genotipe G (24) dan H (16) sangat tipis — gunakan validasi silang
   terstratifikasi + pembobotan kelas, dan laporkan metrik per kelas.
3. Akurasi *baseline* tebak-mayoritas = **14,52%**; itu pembanding yang
   benar, bukan 50%.
4. Masih ada 54 pasangan gen S dengan kemiripan > 0,95 — bagi data per
   klaster kemiripan agar akurasi uji tidak menggelembung.

# Tugas 1 JST — Dataset Sekuens DNA Virus Hepatitis B (HBV)

Pengumpulan dan kurasi dataset sekuens DNA HBV dari NCBI Nucleotide untuk
tiga skenario pemodelan klasifikasi genotipe dengan Jaringan Syaraf Tiruan.

**Laporan lengkap: [`docs/LAPORAN_PENGUMPULAN_DATASET.md`](docs/LAPORAN_PENGUMPULAN_DATASET.md)**

## Hasil ringkas

| | |
|---|---|
| Rekaman diperiksa | 2.309 |
| Rekaman layak | 894 |
| **Rekaman dataset final** | **302** (minimal 300 ✓) |
| Kelas genotipe | 8 (A–H) |
| Duplikat persis | 0 |
| Konkordansi verifikasi genotipe | 99,33% (§5.5 laporan) |
| Keutuhan CDS | 302/302 lolos |

Distribusi: A–F masing-masing 45, G 16, H 16. G dan H dibatasi
ketersediaan di NCBI (bukan oleh penyaringan); G juga kehilangan 9 rekaman
yang terbukti rekombinan pada verifikasi independen.

## Berkas dataset

```
data/raw_gb/*.gb                          2.309 rekaman GenBank asli (berkas utama)
data/processed/dataset1_genS.fasta          gen S                        (~1.203 bp)
data/processed/dataset2_genP.fasta          gen P                        (~2.532 bp)
data/processed/dataset3_genS_genP.fasta     pasangan S+P, accession sama (~3.735 bp)
data/processed/dataset4_genP_nonS.fasta     PEMBANDING: gen P di luar gen S (~1.329 bp)
metadata/metadata_rekaman.xlsx            metadata + status kelayakan tiap rekaman
metadata/verifikasi_genotipe.csv          hasil verifikasi genotipe independen
```

Dataset 4 bukan bagian dari arahan tugas; ia ditambahkan agar perbandingan
ketiga skenario asli dapat ditafsirkan (alasannya di bawah).

## Menjalankan ulang

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt
.venv/Scripts/python scripts/01_collect.py          # unduh dari NCBI
.venv/Scripts/python scripts/02_curate.py           # kurasi + bentuk dataset
.venv/Scripts/python scripts/03_analyze.py          # kendali mutu
.venv/Scripts/python scripts/04_verify_genotype.py  # verifikasi genotipe
.venv/Scripts/python scripts/02_curate.py           # lintasan 2: terapkan pengecualian
.venv/Scripts/python scripts/03_analyze.py
.venv/Scripts/python scripts/04_verify_genotype.py
```

Tahap 2 dan 4 dijalankan dua lintasan: tahap 4 menghasilkan daftar
pengecualian berbasis bukti yang dibaca kembali oleh tahap 2.

## Catatan penting sebelum tahap pemodelan

1. **Gen S terletak seluruhnya di dalam gen P** (302/302 rekaman, pada
   *reading frame* berbeda). Dataset 3 karena itu tidak menambah cakupan
   nukleotida dibanding Dataset 2 — 32,12% isinya duplikasi. Perbandingan
   Dataset 2 vs 3 mengukur pengaruh representasi ganda, **bukan** penambahan
   informasi genomik. Rincian dan rumusan kesimpulan yang dapat
   dipertahankan ada di §6 laporan.
2. Genotipe G (16) dan H (16) sangat tipis — gunakan validasi silang
   terstratifikasi + pembobotan kelas, dan laporkan metrik per kelas.
3. Akurasi *baseline* tebak-mayoritas = **14,90%**; itu pembanding
   yang benar, bukan 50%.
4. Masih ada 43 pasangan gen S dengan kemiripan > 0,95 — bagi data per
   klaster kemiripan agar akurasi uji tidak menggelembung.
5. Tiga rekaman tetap tidak konkordan dengan verifikasi genotipe dan sengaja
   **ditandai, bukan dibuang** (`metadata/verifikasi_genotipe.csv`).

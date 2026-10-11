# Laporan Pengumpulan Dataset Sekuens DNA Virus Hepatitis B (HBV)

Tugas 1 — Mata Kuliah Jaringan Syaraf Tiruan
Tanggal pelaksanaan: 11 Oktober 2026
Sumber data: NCBI Nucleotide (<https://www.ncbi.nlm.nih.gov/nuccore/>)

---

## 1. Ringkasan Eksekutif

| Butir | Hasil |
|---|---|
| Rekaman kandidat diunduh & diperiksa | **2.309** |
| Rekaman memenuhi seluruh kriteria kelayakan | **902** |
| Rekaman masuk dataset final | **310** (ketentuan tugas: minimal 300) |
| Jumlah kelas genotipe | **8** (A–H) |
| Rasio ketidakseimbangan kelas (maks:min) | **2,81** |
| Duplikat sekuens persis pada dataset final | **0** |
| Keutuhan biologis CDS terekstraksi | **310/310 lolos** |

Tiga berkas dataset terbentuk sesuai tiga skenario pemodelan yang diminta,
semuanya memakai **himpunan 310 rekaman yang identik** agar perbandingan
antar skenario tidak terkonfounding oleh perbedaan sampel.

---

## 2. Metode

### 2.1 Strategi pencarian

Kata kunci awal yang diarahkan tugas (`Hepatitis B virus complete genome`)
dipakai sebagai titik mulai, namun **tidak memadai bila dipakai sendirian**.
Alasannya diuji secara empiris, bukan diasumsikan:

Genotipe HBV **tidak tersimpan pada field terstruktur NCBI**. Akibatnya satu
query global menghasilkan 13.097 rekaman yang sangat didominasi genotipe C
dan D, sehingga kelas minoritas (G, H) praktis tidak akan terambil pada
pengambilan acak. Karena itu pencarian dilakukan **per genotipe**, dengan dua
term yang hasilnya digabung (union) dan diterapkan **seragam ke delapan
genotipe**:

| Label | Query |
|---|---|
| `judul` | `"Hepatitis B virus"[Organism] AND "complete genome"[Title] AND 3000:3400[SLEN]` |
| `luas` | `"Hepatitis B virus"[Organism] AND 3000:3400[SLEN] AND biomol_genomic[PROP] NOT patent[Properties]` |

Keduanya dipakai karena pengujian menunjukkan trade-off nyata di antara
keduanya (lihat §3.4): term `judul` menghasilkan mutu anotasi jauh lebih
baik, term `luas` menghasilkan jangkauan lebih besar yang menolong kelas
minoritas.

> **Catatan penting soal keseragaman kriteria.** Pelonggaran query sengaja
> diterapkan ke semua genotipe, bukan hanya ke kelas yang kekurangan data.
> Jika kriteria inklusi berbeda antar kelas, perbedaan teknis itu sendiri
> akan menjadi sinyal yang dapat dipelajari model, sehingga akurasi
> klasifikasi menjadi menggelembung secara artifisial.

Agar reprodusibel sekaligus tidak bias ke era submission tertentu (accession
lama vs baru), skrip mengambil daftar UID yang besar (hingga 900 per
genotipe per term) lalu menarik subsampel deterministik dengan seed tetap
(`SEED = 20251011`).

### 2.2 Kriteria kelayakan rekaman

Sebuah rekaman dinyatakan **Layak** hanya bila seluruh butir berikut
terpenuhi:

1. Genotipe dapat dipastikan **tunggal** dan termasuk A–H.
2. **Bukan rekombinan** antar-genotipe (label kelasnya tidak tunggal).
3. Anotasi CDS **gen S dan gen P keduanya ada**.
4. Kedua CDS **lengkap** (bukan parsial; tidak ada penanda `<` atau `>`).
5. Panjang gen S 1.100–1.310 bp dan gen P 2.400–2.700 bp.
6. Sekuens kedua gen **hanya memuat basa A/C/G/T** (tanpa `N` atau kode
   ambigu IUPAC lain).
7. Sekuens **unik** — diuji pada gen S sendiri, gen P sendiri, *dan*
   pasangannya (alasan di §5.2).

### 2.3 Definisi "gen S" dan "gen P" yang dipakai

- **Gen S** = CDS dengan `product` "large S protein". Gen S HBV adalah satu
  ORF dengan tiga kodon start in-frame (preS1/preS2/S) yang menghasilkan
  protein L, M, dan S. Region gen S yang utuh karena itu adalah CDS protein
  L (large S), bukan CDS "S protein" yang lebih pendek. Panjang modus hasil
  ekstraksi: **1.203 bp**.
- **Gen P** = CDS dengan `product` "polymerase". Panjang modus: **2.532 bp**.

Bila terdapat beberapa CDS kandidat, dipilih yang **terpanjang**, untuk
membedakan "large S" dari "middle S"/"S" yang merupakan bagian di dalamnya.

---

## 3. Temuan Empiris yang Mengubah Desain Pipeline

Bagian ini dicatat karena empat asumsi yang wajar diambil di awal ternyata
**salah**, dan keempatnya diuji lewat probe ke NCBI sebelum pipeline utama
dijalankan.

### 3.1 Qualifier `gene` tidak dapat diandalkan

Pada sampel 75 rekaman:

| Qualifier `gene` | Kemunculan | Qualifier `product` | Kemunculan |
|---|---|---|---|
| `S` | 58 / 75 (77%) | `large S protein` | **75 / 75 (100%)** |
| `P` | **29 / 75 (39%)** | `polymerase` | **74 / 75 (99%)** |

Mencari gen P lewat `gene="P"` akan **kehilangan 61% rekaman**. Ekstraksi
karena itu dibuat berbasis `product`, dengan `gene` hanya sebagai cadangan.

### 3.2 Genotipe tersimpan di `source` → `/note`, bukan di judul

Pada 100% sampel, genotipe tertulis di qualifier `note` pada feature
`source`, misalnya `genotype: C`, `subtype A1; genotype: A`,
`genotype C; subgenotype C2`.

Konsekuensi praktis: query `genotype[Title]` hanya menemukan 499 rekaman,
sedangkan pencarian `[All Fields]` menemukan 13.097. Memakai filter judul
akan membuang 96% data yang sebenarnya berlabel.

### 3.3 Genom HBV sirkular → CDS memakai *compound location*

Seluruh 310 rekaman final bertopologi **circular**. Gen P umumnya
dianotasi sebagai `join(2307..3215, 1..1623)` karena melewati titik awal
penomoran. Ekstraksi **wajib** memakai `feature.extract()`; pemotongan
sekuens biasa (`seq[start:end]`) akan menghasilkan sekuens yang salah
untuk ~40% rekaman.

### 3.4 Melonggarkan query justru menurunkan mutu anotasi

Diuji langsung dengan menjalankan pipeline dua kali:

| Genotipe | Layak dari pool `judul` | Layak dari pool `luas` |
|---|---|---|
| A | 56 dari 120 kandidat (47%) | 32 dari 150 kandidat (21%) |
| G | 30 dari 64 kandidat (47%) | 24 dari 130 kandidat (18%) |

Term `luas` memasukkan banyak rekaman tanpa anotasi gen sama sekali. Karena
itu dipilih **union kedua pool**, bukan menggantikan yang satu dengan yang
lain.

---

## 4. Hasil

### 4.1 Distribusi kelas dataset final

| Genotipe | Jumlah | Persentase | Catatan |
|---|---|---|---|
| A | 45 | 14,52% | dibatasi kuota |
| B | 45 | 14,52% | dibatasi kuota |
| C | 45 | 14,52% | dibatasi kuota |
| D | 45 | 14,52% | dibatasi kuota |
| E | 45 | 14,52% | dibatasi kuota |
| F | 45 | 14,52% | dibatasi kuota |
| G | 24 | 7,74% | **dibatasi ketersediaan NCBI** |
| H | 16 | 5,16% | **dibatasi ketersediaan NCBI** |
| **Total** | **310** | **100%** | |

Kuota 45 per genotipe dipilih agar jumlah total melampaui ketentuan 300
dengan ketidakseimbangan kelas seminimal mungkin. Enam kelas pertama dibatasi
kuota (tersedia 86–179 rekaman layak masing-masing, lihat §4.2), sedangkan
**genotipe G dan H sudah diambil seluruhnya** — keduanya memang langka di
NCBI, bukan akibat penyaringan yang terlalu ketat. Untuk genotipe H, seluruh
database NCBI hanya memuat 33–35 rekaman yang memenuhi kriteria panjang genom.

Akurasi *baseline* "selalu menebak kelas mayoritas" = **14,52%**. Model JST
harus dibandingkan terhadap angka ini, bukan terhadap 50%.

### 4.2 Rekaman layak yang tersedia (sebelum pembatasan kuota)

| Genotipe | A | B | C | D | E | F | G | H |
|---|---|---|---|---|---|---|---|---|
| Layak | 125 | 154 | 145 | 173 | 179 | 86 | 24 | 16 |

Total 902 rekaman layak. Artinya ukuran dataset masih dapat diperbesar
(lihat §7.1).

### 4.3 Sebab penolakan rekaman

Dari 2.309 rekaman yang diperiksa (satu rekaman dapat memiliki lebih dari
satu sebab):

| Sebab | Jumlah |
|---|---|
| Anotasi gen S tidak ada | 573 |
| Gen P memuat basa non-ACGT | 413 |
| Anotasi gen P tidak ada | 342 |
| Gen S memuat basa non-ACGT | 259 |
| Panjang gen S di luar rentang | 237 |
| Genotipe tidak tercantum | 146 |
| Gen S duplikat dengan rekaman lain | 143 |
| Panjang gen P di luar rentang | 38 |
| Rekombinan / antar-genotipe | 38 |
| CDS gen P parsial | 7 |
| CDS gen S parsial | 3 |
| Genotipe ambigu | 1 (A/E) |

Penyebab terbesar bukan mutu sekuens, melainkan **ketiadaan anotasi gen**.
Ini menegaskan arahan tugas untuk mengunduh format GenBank: format FASTA
saja tidak memuat anotasi sehingga penyaringan ini tidak mungkin dilakukan.

### 4.4 Statistik panjang sekuens

| | Min | Modus | Maks |
|---|---|---|---|
| Gen S | 1.149 bp | **1.203 bp** | 1.203 bp |
| Gen P | 2.478 bp | **2.532 bp** | 2.538 bp |
| Genom utuh | — | ~3.215 bp | — |

Panjang tidak seragam. Untuk JST, sekuens harus disamakan panjangnya lebih
dahulu (*padding*/pemangkasan setelah penyejajaran, atau representasi
k-mer / frekuensi yang panjangnya tetap).

---

## 5. Kendali Mutu dan Validasi

### 5.1 Keutuhan biologis sekuens terekstraksi

Uji ini memverifikasi bahwa logika ekstraksi (termasuk penanganan
*compound location*) benar:

| Uji | Gen S | Gen P |
|---|---|---|
| Panjang kelipatan 3 | 310/310 | 310/310 |
| Dimulai kodon start `ATG` | 310/310 | 310/310 |
| Diakhiri kodon stop | 310/310 | 310/310 |
| Tanpa kodon stop internal | 310/310 | 310/310 |

Seluruh 310 sekuens pada kedua gen adalah ORF yang utuh dan dapat
ditranslasikan. **Tidak ada satu pun yang gagal.**

### 5.2 Deduplikasi dan risiko kebocoran data

Deduplikasi awal yang hanya menguji pasangan (S, P) ternyata **tidak cukup**:
ia menyisakan 23 rekaman dengan gen S identik persis tetapi gen P berbeda.
Rekaman seperti itu menjadi duplikat sempurna pada Dataset 1 dan akan
membocorkan informasi antara himpunan latih dan uji. Karena ketiga dataset
harus memakai himpunan rekaman yang sama, keunikan akhirnya ditegakkan pada
**gen S sendiri, gen P sendiri, dan pasangannya sekaligus**.

Hasil setelah perbaikan:

| Metrik | Sebelum | Sesudah |
|---|---|---|
| Gen S identik persis | 23 | **0** |
| Gen P identik persis | 0 | **0** |
| Pasangan gen S sangat mirip (Jaccard 12-mer > 0,95) | 320 | **54** |
| Pasangan gen P sangat mirip | 136 | **45** |

Pasangan *sangat mirip* yang tersisa **tidak dapat dihilangkan** tanpa
merusak dataset: keragaman HBV di dalam satu genotipe memang rendah, dan
kemiripan tinggi antar isolat satu genotipe adalah sifat biologis asli, bukan
kesalahan data. Mitigasinya ada di tahap pemodelan — lihat §7.3.

### 5.3 Validasi silang label genotipe

Label genotipe diparse dari `source`→`/note`. Untuk mengujinya, dipakai dua
field lain yang **independen**:

| Uji silang | Hasil |
|---|---|
| Genotipe pada nama organisme vs hasil parse | **27/27 cocok (100%)** |
| Huruf awal subgenotipe vs hasil parse | **250/250 konsisten (100%)** |

Tidak ditemukan satu pun ketidakcocokan. Label genotipe dinilai dapat
dipercaya pada tingkat anotasi penyetor NCBI.

### 5.4 Integritas berkas dataset

| Uji | Hasil |
|---|---|
| Accession Dataset 1 == Dataset 2 | Ya |
| Accession Dataset 3 == Dataset 1 | Ya |
| Rekaman di Dataset 3 tanpa pasangan S+P lengkap | 0 |
| Header duplikat | 0 |

Syarat tugas bahwa gen S dan gen P pada Dataset 3 harus berasal dari
accession yang sama **terpenuhi untuk seluruh 310 pasang**.

---

## 6. Temuan Kritis: Gen S Seluruhnya Berada di Dalam Gen P

Ini temuan terpenting laporan ini dan **bertentangan dengan premis implisit
Dataset 3**, sehingga perlu disampaikan terbuka.

Hasil pengukuran posisional pada 310 rekaman:

| Metrik | Hasil |
|---|---|
| Rekaman dengan gen S **seluruhnya** di dalam rentang gen P | **310 / 310** |
| Rekaman dengan `union(S, P) == P` | **310 / 310** |
| Rerata basa gen S yang berada **di luar** gen P | **0,00 bp** |
| Total bp Dataset 3 | 1.154.028 |
| Total bp nukleotida **unik** (union) | 783.273 |
| Porsi bp terduplikasi pada Dataset 3 | **32,13%** |

### Artinya

Genom HBV sangat padat: keempat ORF-nya (P, S, C, X) saling bertumpang
tindih. Gen S terletak **sepenuhnya di dalam** rentang gen P, hanya pada
*reading frame* yang berbeda.

Konsekuensinya, pada tingkat **nukleotida**:

- Dataset 1 (gen S) adalah **himpunan bagian** dari Dataset 2 (gen P).
- Dataset 3 (gen S + gen P) **tidak menambah satu basa pun** dibanding
  Dataset 2. Ia hanya menuliskan ulang region gen S untuk kedua kalinya,
  sehingga 32,13% isinya adalah duplikasi.

Jadi membandingkan Dataset 2 dengan Dataset 3 **bukan** membandingkan
"informasi P" dengan "informasi S+P". Yang sebenarnya diuji adalah: apakah
memberi bobot ganda pada region gen S, serta memberi model **dua pembacaan
frame berbeda atas region yang sama**, membantu klasifikasi.

### Mengapa Dataset 3 tetap bermakna

Pertanyaan itu tetap sah secara ilmiah dan layak diteliti, karena:

1. **Reading frame berbeda menghasilkan fitur berbeda.** Bila pengkodean
   dilakukan pada tingkat kodon, asam amino, atau k-mer yang sadar frame,
   gen S dan gen P menghasilkan fitur yang benar-benar berbeda walau
   nukleotidanya sama.
2. **Pembobotan ganda region S dapat menjadi bentuk penekanan fitur.** Gen S
   (HBsAg) adalah region yang secara klasik dipakai untuk penentuan genotipe
   dan serotipe HBV, sehingga penekanan itu punya dasar biologis.
3. **Batas gen memberi informasi struktural.** Memisahkan S dan P memberi
   model informasi batas region yang hilang bila hanya diberi satu sekuens
   panjang.

### Yang harus dicantumkan di laporan penelitian

Kesimpulan penelitian **tidak boleh** dirumuskan sebagai "penggabungan gen S
dan gen P menambah informasi". Rumusan yang jujur dan dapat dipertahankan:

> "Karena gen S HBV terletak seluruhnya di dalam gen P pada *reading frame*
> yang berbeda, Dataset 3 tidak menambah cakupan nukleotida terhadap
> Dataset 2. Perbandingan Dataset 2 dengan Dataset 3 karena itu mengukur
> pengaruh representasi ganda (dua frame, bobot ganda pada region S)
> terhadap kinerja klasifikasi, bukan pengaruh penambahan informasi genomik."

Disarankan menambahkan satu skenario pembanding: **region gen P di luar gen
S saja** (≈1.329 bp). Itu akan memisahkan dengan bersih kontribusi
masing-masing region, dan menjadikan perbandingan tiga skenario tugas ini
dapat ditafsirkan secara kausal.

---

## 7. Keterbatasan dan Rekomendasi untuk Tahap Pemodelan

### 7.1 Ukuran dataset dapat diperbesar

Tersedia **902 rekaman layak**, sementara yang dipakai baru 310. Mengubah
satu konstanta (`TARGET_PER_GENOTYPE` pada `scripts/02_curate.py`) akan
memperbesar dataset:

| Kuota/genotipe | Total rekaman | Rasio ketidakseimbangan |
|---|---|---|
| 45 (sekarang) | 310 | 2,81 |
| 86 | 556 | 5,38 |
| semua | 902 | 11,19 |

Ini **trade-off**, bukan perbaikan gratis: menambah data memperbesar
ketidakseimbangan karena genotipe G dan H sudah habis terpakai. Pilihan 45
diambil karena memenuhi ketentuan 300 dengan ketidakseimbangan terkecil.

### 7.2 Genotipe G (24) dan H (16) sangat tipis

Dengan 16 sampel, genotipe H praktis tidak dapat dievaluasi secara andal:
pembagian 80/20 hanya menyisakan ~3 sampel uji. Pilihan yang tersedia:

1. Pakai **validasi silang terstratifikasi** (*stratified k-fold*) dan
   **pembobotan kelas**, lalu laporkan per-kelas *precision/recall* —
   bukan hanya akurasi total. (Direkomendasikan.)
2. Gabungkan G dan H menjadi satu kelas "genotipe langka".
3. Batasi penelitian ke genotipe A–F (6 kelas, 270 rekaman seimbang) dan
   nyatakan pembatasan itu secara eksplisit. Catatan: 270 < ketentuan 300,
   sehingga opsi ini perlu persetujuan dosen.

### 7.3 Kemiripan antar sekuens → risiko akurasi menggelembung

Masih ada 54 pasangan gen S dan 45 pasangan gen P dengan kemiripan
Jaccard > 0,95. Jika dua isolat nyaris identik terpisah ke himpunan latih dan
uji, akurasi uji akan **terlalu optimistis**. Mitigasi:

- Lakukan pengelompokan (*clustering*) sekuens pada ambang identitas
  (mis. 97%) dan bagi data **per klaster**, bukan per sekuens.
- Alternatif: laporkan juga hasil dengan pembagian berbasis negara atau
  tahun koleksi, yang keduanya tersedia di berkas metadata.

### 7.4 Label berasal dari anotasi penyetor, bukan verifikasi filogenetik

Genotipe diambil dari `/note` yang diisi penyetor rekaman. Validasi silang
(§5.3) menunjukkan konsistensi internal 100%, tetapi itu **tidak sama dengan
verifikasi independen**. Standar emas adalah penentuan genotipe secara
filogenetik (penyejajaran terhadap sekuens rujukan tiap genotipe, lalu
pembangunan pohon). Jika dosen menuntut "genotipe terverifikasi" dalam arti
kuat, langkah tambahan itu perlu dilakukan dan dapat saya tambahkan sebagai
tahap 4.

### 7.5 Bias geografis

Distribusi negara sangat tidak merata (Tiongkok 61 dari 310). Karena
genotipe HBV memang berkorelasi kuat dengan wilayah geografis, model dapat
mempelajari penanda khas galur regional alih-alih ciri genotipe yang
sesungguhnya. Kolom `negara` tersedia di metadata untuk analisis sensitivitas.

### 7.6 Panjang sekuens tidak seragam

Gen S berkisar 1.149–1.203 bp, gen P 2.478–2.538 bp. JST memerlukan dimensi
masukan tetap, sehingga perlu salah satu dari: penyejajaran berganda lalu
*padding*, pemangkasan ke panjang minimum, atau representasi berdimensi
tetap (frekuensi k-mer, one-hot dengan *padding*).

---

## 8. Reprodusibilitas

Seluruh hasil dapat dibangun ulang dari nol:

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt
.venv/Scripts/python scripts/01_collect.py    # unduh dari NCBI
.venv/Scripts/python scripts/02_curate.py     # kurasi + bentuk dataset
.venv/Scripts/python scripts/03_analyze.py    # kendali mutu
```

Pengambilan subsampel memakai seed tetap (`SEED = 20251011`), sehingga
pemilihan rekaman bersifat deterministik. Query lengkap, jumlah kandidat per
genotipe, dan seed tercatat di `data/manifest_pencarian.json`.

Catatan: jumlah rekaman di NCBI bertambah seiring waktu, sehingga
menjalankan ulang `01_collect.py` di masa depan dapat menghasilkan kumpulan
kandidat yang berbeda. Berkas `.gb` mentah disertakan di repositori agar
tahap 2 dan 3 tetap dapat direproduksi secara persis.

Akses NCBI dilakukan tanpa API key, dengan pembatasan laju 3 permintaan per
detik sesuai ketentuan NCBI E-utilities.

---

## 9. Daftar Berkas yang Dikumpulkan

### Berkas utama (deliverable)

| Berkas | Isi |
|---|---|
| `data/raw_gb/*.gb` | 2.309 rekaman GenBank asli dari NCBI (berkas utama — memuat sekuens **dan** anotasi) |
| `data/processed/dataset1_genS.fasta` | 310 sekuens gen S |
| `data/processed/dataset2_genP.fasta` | 310 sekuens gen P |
| `data/processed/dataset3_genS_genP.fasta` | 620 sekuens (310 pasang S+P dari accession yang sama) |
| `metadata/metadata_rekaman.xlsx` | Metadata 2.309 rekaman; sheet `dataset_final` memuat 310 rekaman terpilih |
| `metadata/metadata_rekaman.csv` | Versi CSV dari metadata yang sama |
| `docs/LAPORAN_PENGUMPULAN_DATASET.md` | Laporan ini |

Format header FASTA mengikuti contoh pada arahan tugas, dengan accession
berversi sebagai identitas sampel:

```
>AB241114.1|Genotype=A|Gene=S
ATGGGAGGTTGGTCATCAAAACCTCGCAAAGGCATGGGGACGAACCTTTCTGTTCCCAAC
...
```

### Kolom metadata (memuat seluruh 8 butir yang diminta tugas)

| Butir yang diminta tugas | Kolom |
|---|---|
| 1. Accession number / identitas unik | `accession`, `accession_version` |
| 2. Nama organisme | `nama_organisme` |
| 3. Genotipe HBV terverifikasi | `genotipe`, `subgenotipe` |
| 4. Panjang genom | `panjang_genom_bp` |
| 5. Ketersediaan anotasi gen S | `anotasi_gen_S`, `panjang_gen_S_bp`, `lokasi_gen_S` |
| 6. Ketersediaan anotasi gen P | `anotasi_gen_P`, `panjang_gen_P_bp`, `lokasi_gen_P` |
| 7. Sumber / tautan rekaman NCBI | `tautan_ncbi` |
| 8. Status kelayakan rekaman | `status_kelayakan`, `alasan_tidak_layak` |

Kolom tambahan: `topologi`, `fraksi_S_overlap_P`, `negara`,
`tanggal_koleksi`, `masuk_dataset_final`.

Metadata sengaja memuat **seluruh 2.309 rekaman yang diperiksa**, bukan hanya
310 yang lolos, agar proses penyaringan dapat diaudit dan setiap penolakan
dapat ditelusuri alasannya.

### Berkas pendukung

| Berkas | Isi |
|---|---|
| `scripts/hbv_ncbi.py` | Modul akses NCBI + logika kurasi |
| `scripts/01_collect.py` | Tahap 1 — pencarian & pengunduhan |
| `scripts/02_curate.py` | Tahap 2 — penyaringan & pembentukan dataset |
| `scripts/03_analyze.py` | Tahap 3 — kendali mutu |
| `data/manifest_pencarian.json` | Query, seed, jumlah kandidat per genotipe |
| `data/ringkasan_kurasi.json` | Statistik hasil kurasi |
| `data/laporan_qc.json` | Hasil kendali mutu lengkap |
| `requirements.txt` | Versi paket yang dipakai |

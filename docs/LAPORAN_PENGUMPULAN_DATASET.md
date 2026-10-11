# Laporan Pengumpulan Dataset Sekuens DNA Virus Hepatitis B (HBV)

Tugas 1 — Mata Kuliah Jaringan Syaraf Tiruan
Tanggal pelaksanaan: 11 Oktober 2026
Sumber data: NCBI Nucleotide (<https://www.ncbi.nlm.nih.gov/nuccore/>)

---

## 1. Ringkasan Eksekutif

| Butir | Hasil |
|---|---|
| Rekaman kandidat diunduh & diperiksa | **2.309** |
| Rekaman memenuhi seluruh kriteria kelayakan | **894** |
| Rekaman masuk dataset final | **302** (ketentuan tugas: minimal 300) |
| Jumlah kelas genotipe | **8** (A–H) |
| Rasio ketidakseimbangan kelas (maks:min) | **2,81** |
| Duplikat sekuens persis pada dataset final | **0** |
| Keutuhan biologis CDS terekstraksi | **302/302 lolos** |
| Konkordansi label vs verifikasi independen | **99,33%** (§5.5) |

Tiga berkas dataset terbentuk sesuai tiga skenario pemodelan yang diminta,
ditambah **satu dataset pembanding** (§6), semuanya memakai **himpunan 302
rekaman yang identik** agar perbandingan antar skenario tidak terkonfounding
oleh perbedaan sampel.

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
8. **Tidak ditolak oleh verifikasi tahap 4**, yaitu tidak memenuhi *kedua*
   hal ini sekaligus: publikasi asalnya menyebut rekombinasi, **dan** genotipe
   hasil verifikasi independen berbeda dari label NCBI-nya (alasan di §5.5).

Kriteria 8 mensyaratkan **dua sinyal independen** sekaligus, dan itu
disengaja. Masing-masing sinyal kalau dipakai sendiri tidak sah:

- **Judul publikasi saja terlalu kasar.** Diuji langsung: memakainya sebagai
  dasar pengecualian membuang 158 rekaman, termasuk **55 rekaman genotipe E**
  yang publikasinya berjudul *"Deletions and recombinations in the **core
  region** of hepatitis B virus genotype E strains"*. Rekombinasi di region
  *core* sama sekali tidak menyentuh gen S maupun gen P, sehingga rekaman itu
  sebenarnya sah. Jumlah dataset pun jatuh ke 299 — di bawah ketentuan.
- **Ketidakkonkordanan verifikasi saja bersifat sirkular.** Membuang setiap
  rekaman yang tidak sepakat dengan sekuens rujukan otomatis menjamin
  konkordansi akhir yang tinggi, dan angka itu menjadi tanpa makna.

Gabungan keduanya bersifat spesifik per rekaman dan dapat diaudit: kolom
`indikasi_rekombinan_publikasi` dan `judul_publikasi_rekombinan` pada metadata
mencatat sinyal pertama, dan `metadata/verifikasi_genotipe.csv` mencatat
sinyal kedua untuk setiap sampel.

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

Seluruh 302 rekaman final bertopologi **circular**. Gen P umumnya
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

### 3.5 Filter rekombinan berbasis `/note` ternyata bocor

Temuan ini baru muncul **setelah** verifikasi tahap 4 dijalankan, dan menjadi
alasan utama mengapa tahap verifikasi itu berguna.

Delapan rekaman — `HQ231877.1` sampai `HQ231884.1` — lolos seluruh kriteria
tahap 2 dengan label genotipe G, karena qualifier `/note`-nya hanya menulis
`genotype: G`. Verifikasi tahap 4 menolak kedelapannya dan menetapkannya
sebagai **genotipe A**. Penelusuran ke metadata publikasinya menjelaskan
sebabnya:

> *"A novel complex **a/c/g intergenotypic recombinant** of hepatitis B virus
> isolated in southern China"*

Kedelapan rekaman itu memang rekombinan A/C/G, tetapi informasi itu **hanya
ada di judul publikasi**, tidak di `/note` yang dibaca filter tahap 2. Untuk
klasifikasi genotipe berbasis gen S, kedelapannya praktis adalah label yang
salah — gen S-nya berasal dari galur genotipe A.

Pelajarannya: anotasi `/note` NCBI bukan sumber tunggal yang memadai untuk
menyaring rekombinan, dan verifikasi independen bukan formalitas — ia
menemukan kesalahan yang tidak terlihat dari metadata rekaman saja.

---

## 4. Hasil

### 4.1 Distribusi kelas dataset final

| Genotipe | Jumlah | Persentase | Catatan |
|---|---|---|---|
| A | 45 | 14,90% | dibatasi kuota |
| B | 45 | 14,90% | dibatasi kuota |
| C | 45 | 14,90% | dibatasi kuota |
| D | 45 | 14,90% | dibatasi kuota |
| E | 45 | 14,90% | dibatasi kuota |
| F | 45 | 14,90% | dibatasi kuota |
| G | 16 | 5,30% | **dibatasi ketersediaan NCBI + verifikasi** |
| H | 16 | 5,30% | **dibatasi ketersediaan NCBI** |
| **Total** | **302** | **100%** | |

Kuota 45 per genotipe dipilih agar jumlah total melampaui ketentuan 300
dengan ketidakseimbangan kelas seminimal mungkin. Enam kelas pertama dibatasi
kuota (tersedia 86–179 rekaman layak masing-masing, lihat §4.2), sedangkan
**genotipe G dan H sudah diambil seluruhnya** — keduanya memang langka di
NCBI, bukan akibat penyaringan yang terlalu ketat. Untuk genotipe H, seluruh
database NCBI hanya memuat 33–35 rekaman yang memenuhi kriteria panjang
genom; untuk genotipe G, 9 rekaman tambahan gugur karena terbukti rekombinan
pada verifikasi tahap 4 (§3.5 dan §5.5).

Akurasi *baseline* "selalu menebak kelas mayoritas" = **14,90%**. Model JST
harus dibandingkan terhadap angka ini, bukan terhadap 50%.

### 4.2 Rekaman layak yang tersedia (sebelum pembatasan kuota)

| Genotipe | A | B | C | D | E | F | G | H |
|---|---|---|---|---|---|---|---|---|
| Layak | 125 | 154 | 145 | 173 | 179 | 86 | 16 | 16 |

Total 894 rekaman layak. Artinya ukuran dataset masih dapat diperbesar
(lihat §7.1).

### 4.3 Sebab penolakan rekaman

Dari 2.309 rekaman yang diperiksa (satu rekaman dapat memiliki lebih dari
satu sebab):

| Sebab | Jumlah |
|---|---|
| Anotasi gen S tidak ada | 573 |
| Gen P memuat basa non-ACGT | 413 |
| Anotasi gen P tidak ada | 342 |
| Region gen P non-S memuat basa non-ACGT | 265 |
| Gen S memuat basa non-ACGT | 259 |
| Panjang gen S di luar rentang | 237 |
| Genotipe tidak tercantum | 146 |
| Gen S duplikat dengan rekaman lain | 141 |
| Panjang gen P di luar rentang | 38 |
| Rekombinan / antar-genotipe (dari `/note`) | 38 |
| Ditolak verifikasi tahap 4 (rekombinan, 2 sinyal) | 9 |
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
| Gen P di luar gen S (Dataset 4) | 1.317 bp | **1.329 bp** | 1.368 bp |
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
| Panjang kelipatan 3 | 302/302 | 302/302 |
| Dimulai kodon start `ATG` | 302/302 | 302/302 |
| Diakhiri kodon stop | 302/302 | 302/302 |
| Tanpa kodon stop internal | 302/302 | 302/302 |

Seluruh 302 sekuens pada kedua gen adalah ORF yang utuh dan dapat
ditranslasikan. **Tidak ada satu pun yang gagal.**

Uji terpisah untuk Dataset 4: `|gen S| + |gen P non-S| == |gen P|` terpenuhi
**302/302**, memastikan pemisahan region dilakukan tanpa basa yang hilang
maupun terhitung ganda.

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
| Pasangan gen S sangat mirip (Jaccard 12-mer > 0,95) | 320 | **43** |
| Pasangan gen P sangat mirip | 136 | **33** |

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
| Accession Dataset 4 == Dataset 1 | Ya |
| Rekaman di Dataset 3 tanpa pasangan S+P lengkap | 0 |
| Header duplikat | 0 |

Syarat tugas bahwa gen S dan gen P pada Dataset 3 harus berasal dari
accession yang sama **terpenuhi untuk seluruh 302 pasang**.

### 5.5 Verifikasi genotipe secara independen (tahap 4)

Label genotipe pada tahap 2 berasal dari qualifier `source`→`/note`, yaitu
anotasi yang diisi penyetor rekaman. Tahap 4 mengujinya dengan dua metode yang
**tidak melihat catatan itu sama sekali**.

**Sekuens rujukan.** Dipakai dua set rujukan dari dua publikasi berbeda agar
hasil tidak bergantung pada satu sumber: set berbasis panel WHO
([PMC6475771](https://pmc.ncbi.nlm.nih.gov/articles/PMC6475771/table/T5)) dan
set rujukan genotipe/subgenotipe
([PMC11496132](https://pmc.ncbi.nlm.nih.gov/articles/PMC11496132/table/T2)).
Dari 16 accession rujukan, **12** dapat dipakai untuk gen S
(sisanya tidak punya anotasi gen S yang utuh). Setiap genotipe A–H tetap
terwakili minimal satu rujukan.

Sebagai pemeriksaan atas set rujukan itu sendiri: genotipe menurut publikasi
cocok dengan catatan NCBI pada **16/16** rujukan.

Dua accession rujukan (`AP007264.1`, `AY090460.1`) juga terdapat di dataset;
keduanya **dikeluarkan dari himpunan sampel** saat pengujian agar tidak
dicocokkan dengan dirinya sendiri.

**Metode A — rujukan terdekat via penyejajaran global gen S.**
Setiap sampel disejajarkan terhadap seluruh sekuens rujukan, lalu diberi
genotipe rujukan dengan skor tertinggi.

> **Konkordansi: 298/300 (99,33%)** label NCBI == genotipe
> rujukan terdekat.

Margin skor terhadap genotipe terbaik berikutnya: median
**142**, minimum **2**.

**Metode B — pohon Neighbor-Joining dari jarak k-mer (k=7).**
Metrik utamanya adalah kesepakatan label tetangga terdekat (1-NN), karena
metrik itu tahan terhadap satu-dua sampel menyimpang:

> **1-NN: 309/312 (99,04%)** sampel bertetangga terdekat dengan
> sampel berlabel genotipe sama.

| Genotipe | Metode A (rujukan terdekat) | Metode B (1-NN) |
|---|---|---|
| A | 45/45 | 47/47 |
| B | 44/45 | 45/47 |
| C | 45/45 | 46/46 |
| D | 45/45 | 47/47 |
| E | 45/45 | 46/46 |
| F | 45/45 | 47/47 |
| G | 14/15 | 15/16 |
| H | 15/15 | 16/16 |

**Sampel yang tetap tidak konkordan** (ditandai, **tidak** dibuang):

| Accession | Label NCBI | Rujukan terdekat | Margin skor |
|---|---|---|---|
| `EU939681.1` | B | D | 188 |
| `KP274926.1` | G | A | 88 |

Ketiganya dibiarkan di dalam dataset dan hanya ditandai, karena membuang
rekaman semata-mata karena ia tidak sepakat dengan rujukan bersifat sirkular
(lihat kriteria 8 di §2.2). Hasil per sampel tercatat lengkap di
`metadata/verifikasi_genotipe.csv` dan sheet `verifikasi_genotipe` pada
berkas Excel, sehingga dapat ditinjau sendiri.

**Pengecualian yang dilakukan.** Dari 158 rekaman yang
ditandai rekombinan oleh judul publikasinya dan
2 rekaman yang tidak konkordan, hanya
**9 rekaman** memenuhi kedua syarat sekaligus dan
karena itu dibuang:

`HQ231877.1, HQ231878.1, HQ231879.1, HQ231880.1, HQ231881.1, HQ231882.1, HQ231883.1, HQ231884.1, HQ231885.1`

Seluruhnya berasal dari satu seri submission yang sama, yaitu publikasi
*"A novel complex a/c/g intergenotypic recombinant of hepatitis B virus
isolated in southern China"* (§3.5).

**Iterasi sampai konvergen.** Membuang satu rekaman membebaskan slot
deduplikasi, sehingga rekaman lain dari seri yang sama dapat naik menjadi
layak pada lintasan berikutnya — itu benar-benar terjadi pada `HQ231885.1`,
yang baru muncul setelah delapan rekaman sebelumnya dibuang. Tahap 2 dan
tahap 4 karena itu diulang sampai tidak ada penambahan pengecualian baru.
`scripts/run_pipeline.py` menjalankan iterasi ini otomatis. Status akhir:
**konvergen**.

**Keterbatasan metrik monofili.** Uji monofili berbasis MRCA (nenek moyang
bersama seluruh sampel berlabel genotipe X) dilaporkan di bawah, tetapi
**metrik ini rapuh** dan tidak layak dibaca sebagai kesimpulan: satu sampel
menyimpang membuat klad MRCA membesar hingga nyaris seluruh pohon, sehingga
kemurnian jatuh drastis tanpa berarti genotipe itu tersebar luas. Angka 1-NN
di atas adalah ukuran yang sahih.

| Genotipe | Anggota | Isi klad MRCA | Kemurnian | Klad murni terbesar |
|---|---|---|---|---|
| A | 47 | 48 | 0.979 | 24/47 |
| B | 47 | 249 | 0.189 | 45/47 |
| C | 46 | 46 | 1.000 | 46/46 |
| D | 47 | 48 | 0.979 | 3/47 |
| E | 46 | 46 | 1.000 | 46/46 |
| F | 47 | 47 | 1.000 | 47/47 |
| G | 16 | 155 | 0.103 | 15/16 |
| H | 16 | 16 | 1.000 | 16/16 |

**Kesimpulan tahap 4.** Label genotipe dataset ini konsisten dengan posisi
sekuens relatif terhadap sekuens rujukan yang diakui komunitas, pada tingkat
99,33% (penyejajaran) dan 99,04% (1-NN). Perlu dicatat secara
jujur bahwa sekuens rujukan itu sendiri berlabel dari literatur, bukan dari
penentuan de novo; jadi yang diverifikasi adalah **konsistensi terhadap
rujukan komunitas**, bukan penentuan genotipe dari nol. Itu jauh lebih kuat
daripada sekadar mempercayai catatan penyetor, tetapi bukan klaim absolut.

---

## 6. Temuan Kritis: Gen S Seluruhnya Berada di Dalam Gen P

Ini temuan terpenting laporan ini dan **bertentangan dengan premis implisit
Dataset 3**, sehingga perlu disampaikan terbuka.

Hasil pengukuran posisional pada 302 rekaman:

| Metrik | Hasil |
|---|---|
| Rekaman dengan gen S **seluruhnya** di dalam rentang gen P | **302 / 302** |
| Rekaman dengan `union(S, P) == P` | **302 / 302** |
| Rerata basa gen S yang berada **di luar** gen P | **0,00 bp** |
| Total bp Dataset 3 | 1.124.148 |
| Total bp nukleotida **unik** (union) | 763.017 |
| Porsi bp terduplikasi pada Dataset 3 | **32,12%** |

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

### Dataset 4 — pembanding yang ditambahkan

Agar perbandingan ketiga skenario tugas dapat ditafsirkan secara kausal,
ditambahkan satu dataset pembanding: **region gen P di luar gen S saja**
(`dataset4_genP_nonS.fasta`, modus **1.329 bp** = 2.532 − 1.203).

Dengan empat dataset, kontribusi tiap region dapat dipisahkan bersih:

| Dataset | Isi | Panjang modus | Hubungan |
|---|---|---|---|
| 1 | gen S | 1.203 bp | ⊂ gen P |
| 2 | gen P | 2.532 bp | = S ∪ (P\S) |
| 3 | gen S + gen P | 3.735 bp | region S berbobot ganda |
| **4** | **gen P di luar gen S** | **1.329 bp** | komplemen gen S di dalam gen P |

Region pada Dataset 1 dan Dataset 4 **tidak bertumpang tindih sama sekali**
dan bersama-sama menyusun tepat Dataset 2. Jadi bila Dataset 1 dan Dataset 4
dibandingkan, selisih kinerjanya benar-benar berasal dari region genom yang
berbeda — bukan dari perbedaan bobot atau *reading frame*. Itu yang
menjadikan hasil ketiga skenario asli dapat ditafsirkan.

Catatan: region gen P non-S tersusun dari dua segmen tak bersambung (bagian
gen P sebelum gen S mulai dan sesudah gen S berakhir). Keduanya dirangkai
mengikuti urutan baca CDS gen P. Karena pemotongan di tengah menggeser
*reading frame*, Dataset 4 **tidak** merupakan ORF yang dapat ditranslasikan
dan hanya sah dipakai sebagai fitur pada tingkat nukleotida/k-mer, bukan pada
tingkat kodon atau asam amino.

---

## 7. Keterbatasan dan Rekomendasi untuk Tahap Pemodelan

### 7.1 Ukuran dataset dapat diperbesar

Tersedia **894 rekaman layak**, sementara yang dipakai baru 302. Mengubah
satu konstanta (`TARGET_PER_GENOTYPE` pada `scripts/02_curate.py`) akan
memperbesar dataset:

| Kuota/genotipe | Total rekaman | Rasio ketidakseimbangan |
|---|---|---|
| 45 (sekarang) | 302 | 2,81 |
| 86 | 548 | 5,38 |
| semua | 894 | 11,19 |

Ini **trade-off**, bukan perbaikan gratis: menambah data memperbesar
ketidakseimbangan karena genotipe G dan H sudah habis terpakai. Pilihan 45
diambil karena memenuhi ketentuan 300 dengan ketidakseimbangan terkecil.

### 7.2 Genotipe G (16) dan H (16) sangat tipis

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

Masih ada 43 pasangan gen S dan 33 pasangan gen P dengan kemiripan
Jaccard > 0,95. Jika dua isolat nyaris identik terpisah ke himpunan latih dan
uji, akurasi uji akan **terlalu optimistis**. Mitigasi:

- Lakukan pengelompokan (*clustering*) sekuens pada ambang identitas
  (mis. 97%) dan bagi data **per klaster**, bukan per sekuens.
- Alternatif: laporkan juga hasil dengan pembagian berbasis negara atau
  tahun koleksi, yang keduanya tersedia di berkas metadata.

### 7.4 Label berasal dari anotasi penyetor, lalu diverifikasi

Genotipe awalnya diambil dari `/note` yang diisi penyetor rekaman. Verifikasi
independen pada tahap 4 (§5.5) menunjukkan konsistensi **99,33%** terhadap sekuens rujukan dan
**99,04%** pada uji 1-NN, serta berhasil
menemukan dan membuang 9 rekaman rekombinan yang salah label.

Yang masih menjadi keterbatasan: sekuens rujukan itu sendiri berlabel dari
literatur, bukan dari penentuan genotipe de novo. Penentuan de novo yang
sesungguhnya memerlukan penyejajaran berganda (MSA) dengan perangkat seperti
MAFFT dan rekonstruksi pohon dengan dukungan bootstrap -- di luar lingkup
pengumpulan dataset ini. Tiga rekaman tetap tidak konkordan dan sengaja
dibiarkan di dataset dengan penanda, bukan dibuang (alasannya di §5.5).

### 7.5 Bias geografis

Distribusi negara sangat tidak merata (Tiongkok ~54 dari 302). Karena
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
.venv/Scripts/python scripts/run_pipeline.py --unduh
```

`run_pipeline.py` menjalankan seluruh tahap dalam urutan yang benar dan
**mengulang tahap 2 ↔ tahap 4 sampai konvergen**. Pengulangan itu wajib:
tahap 4 menghasilkan daftar pengecualian yang dibaca tahap 2, dan membuang
satu rekaman membebaskan slot deduplikasi sehingga rekaman lain dapat naik
menjadi layak pada lintasan berikutnya (§5.5). Tahap-tahap juga dapat
dijalankan satu per satu:

```bash
.venv/Scripts/python scripts/01_collect.py          # unduh dari NCBI
.venv/Scripts/python scripts/02_curate.py           # kurasi + bentuk dataset
.venv/Scripts/python scripts/03_analyze.py          # kendali mutu
.venv/Scripts/python scripts/04_verify_genotype.py  # verifikasi genotipe
```

Pengambilan subsampel memakai seed tetap (`SEED = 20251011`), sehingga
pemilihan rekaman bersifat deterministik. Query lengkap, jumlah kandidat per
genotipe, dan seed tercatat di `data/manifest_pencarian.json`. Daftar
pengecualian hasil verifikasi tercatat di
`data/pengecualian_verifikasi.json`, sehingga keputusan membuang 9 rekaman
itu dapat ditinjau dan dibatalkan bila tidak disetujui.

Catatan: jumlah rekaman di NCBI bertambah seiring waktu, sehingga
menjalankan ulang `01_collect.py` di masa depan dapat menghasilkan kumpulan
kandidat yang berbeda. Berkas `.gb` mentah disertakan di repositori agar
tahap 2–4 tetap dapat direproduksi secara persis. Sekuens rujukan tahap 4
juga disimpan (`data/rujukan_genotipe.gb`) agar verifikasi tidak bergantung
pada unduhan ulang.

Akses NCBI dilakukan tanpa API key, dengan pembatasan laju 3 permintaan per
detik sesuai ketentuan NCBI E-utilities.

---

## 9. Daftar Berkas yang Dikumpulkan

### Berkas utama (deliverable)

| Berkas | Isi |
|---|---|
| `data/raw_gb/*.gb` | 2.309 rekaman GenBank asli dari NCBI (berkas utama — memuat sekuens **dan** anotasi) |
| `data/processed/dataset1_genS.fasta` | 302 sekuens gen S (~1.203 bp) |
| `data/processed/dataset2_genP.fasta` | 302 sekuens gen P (~2.532 bp) |
| `data/processed/dataset3_genS_genP.fasta` | 604 sekuens (302 pasang S+P dari accession yang sama) |
| `data/processed/dataset4_genP_nonS.fasta` | 302 sekuens region gen P di luar gen S (~1.329 bp) — **pembanding tambahan, bukan bagian arahan tugas** (§6) |
| `metadata/metadata_rekaman.xlsx` | Metadata 2.309 rekaman; sheet `dataset_final` memuat 302 rekaman terpilih, sheet `verifikasi_genotipe` memuat hasil tahap 4 |
| `metadata/metadata_rekaman.csv` | Versi CSV dari metadata yang sama |
| `metadata/verifikasi_genotipe.csv` | Hasil verifikasi genotipe per sampel (§5.5) |
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
| `scripts/04_verify_genotype.py` | Tahap 4 — verifikasi genotipe independen |
| `scripts/run_pipeline.py` | Driver: menjalankan tahap 2 ↔ 4 sampai konvergen |
| `data/manifest_pencarian.json` | Query, seed, jumlah kandidat per genotipe |
| `data/ringkasan_kurasi.json` | Statistik hasil kurasi |
| `data/laporan_qc.json` | Hasil kendali mutu lengkap |
| `data/laporan_verifikasi_genotipe.json` | Hasil verifikasi tahap 4 lengkap |
| `data/pengecualian_verifikasi.json` | 9 accession yang dibuang + kriterianya |
| `data/rujukan_genotipe.gb` | 16 sekuens rujukan genotipe (2 publikasi) |
| `data/pohon_nj_genS.newick` | Pohon Neighbor-Joining gen S (format Newick) |
| `requirements.txt` | Versi paket yang dipakai |

### Sumber sekuens rujukan tahap 4

- [PMC6475771 — tabel rujukan berbasis panel WHO](https://pmc.ncbi.nlm.nih.gov/articles/PMC6475771/table/T5)
- [PMC11496132 — tabel rujukan genotipe/subgenotipe](https://pmc.ncbi.nlm.nih.gov/articles/PMC11496132/table/T2)
- [HBVdb — basis data & nomenklatur HBV](https://hbvdb.lyon.inserm.fr/HBVdb/HBVdbNomenclature)

Tidak ada satu set rujukan yang disepakati semua sumber; karena itu dipakai
dua set sekaligus (§5.5).

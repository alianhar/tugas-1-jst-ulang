## **PENGUMPULAN DATASET SEKUENS DNA VIRUS HEPATITIS B \(HBV\)**
### 1.1.1 **Tujuan Penugasan**
### Mengumpulkan, memeriksa, dan mendokumentasikan data sekuens DNA virus Hepatitis B \(HBV\) yang akan digunakan dalam penelitian pemodelan Jaringan Syaraf Tiruan \(JST\)\. Kumpulkan dataset untuk tiga skenario pemodelan:
### 1.1.2 Dataset yang menggunakan sekuens gen S saja untuk klasifikasi genotipe HBV\.
### 1.1.3 Dataset yang menggunakan sekuens gen P saja untuk klasifikasi genotipe HBV\.
### 1.1.4 Dataset yang menggunakan gabungan sekuens gen S dan gen P untuk klasifikasi genotipe HBV\.
### 1.1.5 **Sumber Data**
Pencarian data dapat dilakukan pada  sumber public NCBI Nucleotide: [<u>https://www\.ncbi\.nlm\.nih\.gov/nuccore/</u>](https://www.ncbi.nlm.nih.gov/nuccore/)\. Gunakan kata kunci pencarian awal berikut :
- Hepatitis B virus complete genome
Pilih rekaman genom HBV lengkap yang identitas rekaman, informasi organisme, dan anotasi gen\-nya dapat diperiksa\.
### 1.1.6 **Data yang Harus Dikumpulkan**
Untuk setiap rekaman genom HBV, pastikan memuat data :
1. Accession number atau identitas unik rekaman\.
2. Nama organisme\.
3. Genotipe HBV yang terverifikasi\.
4. Panjang genom\.
5. Ketersediaan anotasi gen S\.
6. Ketersediaan anotasi gen P\.
7. Sumber atau tautan rekaman NCBI\.
8. Status kelayakan rekaman untuk penelitian\.
Unduh rekaman dalam GenBank \(\.gb\) atau FASTA \(\.fasta, \.fa\) sebagai berkas utama karena format ini menyimpan sekuens dan anotasi\.
1. **Prinsip Penyusunan Dataset**
Data sekuens DNA virus Hepatitis B \(HBV\) dikumpulkan dari rekaman genom yang memiliki accession number dan sumber yang dapat ditelusuri\. Setiap rekaman perlu diperiksa untuk memastikan identitas organisme, genotipe, serta anotasi gen S dan gen P\.
Data akan dipersiapkan untuk tiga skenario pemodelan:
1. **Dataset 1:** menggunakan sekuens gen S saja\.
2. **Dataset 2:** menggunakan sekuens gen P saja\.
3. **Dataset 3:** menggunakan gabungan sekuens gen S dan gen P\.
Pada Dataset 3, sekuens gen S dan gen P yang dipasangkan sebagai satu sampel harus berasal dari rekaman genom yang sama\. Sebagai contoh, jika gen S berasal dari accession REC\_A, gen P yang dipasangkan dengannya juga harus berasal dari REC\_A, bukan dari REC\_B\. Accession REC\_A, REC\_B, dan seterusnya pada contoh berikut hanya merupakan kode ilustrasi, bukan accession NCBI yang telah diverifikasi\. Sekuens DNA yang digunakan sebagai contoh juga harus diganti dengan sekuens hasil pengunduhan dan ekstraksi yang sebenarnya\.
## 1.2 **Contoh Bentuk Dataset Sekuens**
Data sekuens disimpan dalam format **GenBank \(\.gb\)** atau FASTA \(\.fasta atau \.fa\)\. Setiap rekaman GenBank atau FASTA terdiri atas baris header yang diawali tanda > dan diikuti oleh sekuens DNA\. Header asli dari NCBI sebaiknya dipertahankan pada berkas sumber\. Pada berkas hasil pengolahan, informasi identitas sampel, accession, gen, dan genotipe dapat dicantumkan secara konsisten\. Contoh di bawah menunjukkan struktur berkas, bukan sekuens HBV yang dapat langsung digunakan untuk penelitian\.
### **Dataset 1 — Sekuens Gen S Saja**
Berkas: dataset1\_genS\.fasta
Dataset ini hanya memuat sekuens gen S dari setiap rekaman yang memenuhi kriteria\.


<table><tr><td><p>>REC_A|Genotype=A|Gene=S</p><p>[sekuens DNA gen S hasil ekstraksi dari REC_A]</p><p>>REC_B|Genotype=B|Gene=S</p><p>[sekuens DNA gen S hasil ekstraksi dari REC_B]</p><p>>REC_C|Genotype=C|Gene=S</p><p>[sekuens DNA gen S hasil ekstraksi dari REC_C]</p><p>>REC_D|Genotype=A|Gene=S</p><p>[sekuens DNA gen S hasil ekstraksi dari REC_D]</p></td></tr></table>


### **Dataset 2 — Sekuens Gen P Saja**
Berkas: dataset2\_genP\.fasta
Dataset ini hanya memuat sekuens gen P dari rekaman yang memenuhi kriteria\.


<table><tr><td><p>>REC_A|Genotype=A|Gene=P</p><p>[sekuens DNA gen P hasil ekstraksi dari REC_A]</p><p>>REC_B|Genotype=B|Gene=P</p><p>[sekuens DNA gen P hasil ekstraksi dari REC_B]</p><p>>REC_C|Genotype=C|Gene=P</p><p>[sekuens DNA gen P hasil ekstraksi dari REC_C]</p><p>>REC_D|Genotype=A|Gene=P</p><p>[sekuens DNA gen P hasil ekstraksi dari REC_D]</p></td></tr></table>


### **Dataset 3 — Gabungan Sekuens Gen S dan Gen P**
Berkas: dataset3\_genS\_genP\.fasta
Pada tahap penyimpanan awal, kedua gen dapat dituliskan sebagai rekaman FASTA terpisah dalam satu berkas\. Identitas sampel dan accession digunakan untuk menghubungkan gen S dan gen P yang berasal dari genom yang sama\.


<table><tr><td><p>>REC_A|Genotype=A|Gene=S</p><p>[sekuens DNA gen S hasil ekstraksi dari REC_A]</p><p>>REC_A|Genotype=A|Gene=P</p><p>[sekuens DNA gen P hasil ekstraksi dari REC_A]</p><p>>REC_B|Genotype=B|Gene=S</p><p>[sekuens DNA gen S hasil ekstraksi dari REC_B]</p><p>>REC_B|Genotype=B|Gene=P</p><p>[sekuens DNA gen P hasil ekstraksi dari REC_B]</p></td></tr></table>


Struktur tersebut menunjukkan bahwa kedua gen berasal dari accession yang sama\. Namun, berkas ini masih berupa penyimpanan dua sekuens terpisah; belum merupakan matriks fitur gabungan yang siap dimasukkan ke JST\. Pada tahap pemodelan, sekuens gen S dan gen P harus dikodekan secara konsisten, kemudian fitur keduanya digabungkan berdasarkan ID sampel\.
## 1.3 **Berkas yang Dikumpulkan Mahasiswa**
## Unduh rekaman dalam format **GenBank** sebagai berkas utama karena format ini menyimpan sekuens dan anotasi\. Format FASTA dapat digunakan sebagai format tambahan untuk menyimpan sekuens DNA\.
## 1.4 **Ketentuan Jumlah Data**
## Kumpulkan minimal 300 rekaman genom lengkap HBV yang valid dan unik sebagai dataset sumber\. 
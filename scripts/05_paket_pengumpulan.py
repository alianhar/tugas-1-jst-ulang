"""Tahap 5 — susun folder pengumpulan yang siap diserahkan + arsip ZIP.

Masalah yang diselesaikan: `data/raw_gb/` memuat 2.309 berkas GenBank yang
bercampur — 302 di antaranya adalah dataset final, 2.007 sisanya kandidat
yang ditolak penyaringan. Diserahkan apa adanya, pemeriksa tidak dapat
membedakan keduanya dan berisiko salah menghitung jumlah rekaman.

Struktur keluaran:
    pengumpulan/
      01_genbank_final/     302 berkas .gb dataset final (berkas utama)
      02_dataset_fasta/     4 berkas .fasta
      03_metadata/          metadata + hasil verifikasi genotipe
      04_laporan/           laporan dalam Markdown dan PDF
      05_lampiran_audit/    2.007 berkas .gb kandidat yang ditolak
      06_skrip/             skrip pipeline (laporan merujuk ke berkas ini)
      README.txt
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "pengumpulan"
PENANDA = "# berkas ini dibuat otomatis oleh scripts/05_paket_pengumpulan.py"
CHROME = Path(r"C:\Program Files\Google\Chrome\Application\chrome.exe")

CSS = """
@page { size: A4; margin: 18mm 16mm; }
body { font-family: "Segoe UI", Arial, sans-serif; font-size: 10.5pt;
       line-height: 1.5; color: #111; }
h1 { font-size: 19pt; border-bottom: 2px solid #333; padding-bottom: 6px; }
h2 { font-size: 14pt; margin-top: 22px; border-bottom: 1px solid #bbb;
     padding-bottom: 3px; page-break-after: avoid; }
h3 { font-size: 12pt; margin-top: 16px; page-break-after: avoid; }
h4 { font-size: 11pt; margin-top: 12px; page-break-after: avoid; }
table { border-collapse: collapse; width: 100%; margin: 10px 0;
        font-size: 9.5pt; page-break-inside: avoid; }
th, td { border: 1px solid #999; padding: 4px 7px; text-align: left;
         vertical-align: top; }
th { background: #eee; font-weight: 600; }
code { background: #f2f2f2; padding: 1px 4px; border-radius: 3px;
       font-family: Consolas, monospace; font-size: 9pt; }
pre { background: #f6f6f6; border: 1px solid #ddd; padding: 8px;
      border-radius: 4px; page-break-inside: avoid; }
pre code { background: none; padding: 0; }
blockquote { border-left: 4px solid #888; margin-left: 0; padding-left: 12px;
             color: #333; background: #fafafa; }
a { color: #1a4f8a; }
"""


def bersihkan_keluaran() -> None:
    """Hapus folder keluaran lama, tetapi HANYA bila jelas milik skrip ini."""
    if not OUT.exists():
        return
    readme = OUT / "README.txt"
    if not (readme.exists() and PENANDA in readme.read_text(encoding="utf-8")):
        raise SystemExit(
            f"'{OUT}' sudah ada tetapi bukan hasil skrip ini (penanda tidak\n"
            f"ditemukan di README.txt). Dihentikan agar tidak menimpa berkas\n"
            f"Anda. Pindahkan atau hapus folder itu lebih dahulu.")
    print(f"  folder lama dihapus: {OUT.name}/")
    shutil.rmtree(OUT)


def buat_pdf(md_path: Path, pdf_path: Path) -> bool:
    import markdown
    html = markdown.markdown(
        md_path.read_text(encoding="utf-8"),
        extensions=["tables", "fenced_code", "sane_lists", "attr_list"])
    tmp = OUT / "_laporan_sementara.html"
    tmp.write_text(
        f"<!doctype html><html lang='id'><head><meta charset='utf-8'>"
        f"<style>{CSS}</style></head><body>{html}</body></html>",
        encoding="utf-8")
    if not CHROME.exists():
        print(f"  !! Chrome tidak ditemukan di {CHROME}; PDF dilewati")
        tmp.unlink()
        return False
    hasil = subprocess.run(
        [str(CHROME), "--headless=new", "--disable-gpu", "--no-sandbox",
         "--no-pdf-header-footer", f"--print-to-pdf={pdf_path}",
         tmp.resolve().as_uri()],
        capture_output=True, text=True, timeout=180)
    tmp.unlink(missing_ok=True)
    if pdf_path.exists() and pdf_path.stat().st_size > 20_000:
        print(f"  PDF dibuat: {pdf_path.name} "
              f"({pdf_path.stat().st_size / 1024:.0f} KB)")
        return True
    print(f"  !! PDF gagal (exit {hasil.returncode}): "
          f"{(hasil.stderr or '')[-300:]}")
    return False


def main() -> None:
    meta_csv = ROOT / "metadata" / "metadata_rekaman.csv"
    df = pd.read_csv(meta_csv)
    final = df[df.masuk_dataset_final == "Ya"].accession.astype(str).tolist()
    final_set = set(final)
    raw = ROOT / "data" / "raw_gb"
    semua = sorted(raw.glob("*.gb"))
    ditolak = [p for p in semua if p.stem not in final_set]

    qc = json.loads((ROOT / "data" / "laporan_qc.json").read_text(encoding="utf-8"))
    ver = json.loads((ROOT / "data" / "laporan_verifikasi_genotipe.json")
                     .read_text(encoding="utf-8"))

    print("== Tahap 5: menyusun folder pengumpulan ==")
    print(f"  rekaman dataset final : {len(final)}")
    print(f"  kandidat ditolak      : {len(ditolak)}")
    bersihkan_keluaran()

    d1 = OUT / "01_genbank_final"
    d2 = OUT / "02_dataset_fasta"
    d3 = OUT / "03_metadata"
    d4 = OUT / "04_laporan"
    d5 = OUT / "05_lampiran_audit"
    d6 = OUT / "06_skrip"
    for d in (d1, d2, d3, d4, d5, d6):
        d.mkdir(parents=True)

    hilang = []
    for acc in final:
        src = raw / f"{acc}.gb"
        if src.exists():
            shutil.copy2(src, d1 / src.name)
        else:
            hilang.append(acc)
    if hilang:
        raise SystemExit(f"berkas .gb tidak ditemukan untuk: {hilang[:10]}")
    print(f"  01_genbank_final  : {len(list(d1.glob('*.gb')))} berkas")

    for p in sorted((ROOT / "data" / "processed").glob("*.fasta")):
        shutil.copy2(p, d2 / p.name)
    print(f"  02_dataset_fasta  : {len(list(d2.glob('*.fasta')))} berkas")

    for p in sorted((ROOT / "metadata").iterdir()):
        if p.is_file():
            shutil.copy2(p, d3 / p.name)
    print(f"  03_metadata       : {len(list(d3.iterdir()))} berkas")

    laporan = ROOT / "docs" / "LAPORAN_PENGUMPULAN_DATASET.md"
    shutil.copy2(laporan, d4 / laporan.name)
    buat_pdf(laporan, d4 / "LAPORAN_PENGUMPULAN_DATASET.pdf")

    for p in ditolak:
        shutil.copy2(p, d5 / p.name)
    print(f"  05_lampiran_audit : {len(list(d5.glob('*.gb')))} berkas")

    for p in sorted((ROOT / "scripts").glob("*.py")):
        shutil.copy2(p, d6 / p.name)
    shutil.copy2(ROOT / "requirements.txt", d6 / "requirements.txt")
    print(f"  06_skrip          : {len(list(d6.iterdir()))} berkas")

    dist = qc["distribusi_kelas"]
    (OUT / "README.txt").write_text(f"""{PENANDA}

PENGUMPULAN TUGAS 1 - JARINGAN SYARAF TIRUAN
Dataset Sekuens DNA Virus Hepatitis B (HBV) untuk Klasifikasi Genotipe
Sumber data: NCBI Nucleotide (https://www.ncbi.nlm.nih.gov/nuccore/)

====================================================================
RINGKASAN
====================================================================
Rekaman kandidat diperiksa    : {qc['jumlah_diperiksa']}
Rekaman memenuhi kriteria     : {qc['jumlah_layak']}
REKAMAN DATASET FINAL         : {qc['jumlah_dataset_final']}   (ketentuan tugas: minimal 300)
Jumlah kelas genotipe         : 8 (A-H)
Distribusi genotipe           : {', '.join(f'{g}={n}' for g, n in dist.items())}
Duplikat sekuens persis       : 0
Keutuhan ORF gen S dan gen P  : {qc['jumlah_dataset_final']}/{qc['jumlah_dataset_final']} lolos
Konkordansi verifikasi        : {ver['metode_A']['persen']}% terhadap sekuens rujukan

====================================================================
ISI FOLDER
====================================================================
01_genbank_final/     {len(final)} berkas .gb -- INI DATASET YANG DINILAI.
                      Format GenBank sebagai berkas utama sesuai butir 1.3
                      arahan tugas, karena memuat sekuens DAN anotasi gen.
                      Satu berkas = satu rekaman genom, dinamai menurut
                      accession number-nya.

02_dataset_fasta/     Empat berkas FASTA sebagai format tambahan:
                        dataset1_genS.fasta       {len(final)} sekuens gen S
                        dataset2_genP.fasta       {len(final)} sekuens gen P
                        dataset3_genS_genP.fasta  {len(final) * 2} sekuens ({len(final)} pasang S+P
                                                  dari accession yang sama)
                        dataset4_genP_nonS.fasta  {len(final)} sekuens region gen P
                                                  di luar gen S -- TAMBAHAN,
                                                  bukan bagian arahan tugas.
                                                  Alasannya di bagian 6 laporan.
                      Format header: >ACCESSION|Genotype=X|Gene=Y

03_metadata/          metadata_rekaman.xlsx -- memuat seluruh 8 butir data
                      yang diminta arahan tugas (butir 1.1.6), untuk SEMUA
                      {qc['jumlah_diperiksa']} rekaman yang diperiksa, bukan hanya yang lolos,
                      agar proses penyaringan dapat diaudit.
                        sheet "semua_rekaman"       semua kandidat + alasan tolak
                        sheet "dataset_final"       {qc['jumlah_dataset_final']} rekaman terpilih
                        sheet "verifikasi_genotipe" hasil verifikasi genotipe
                      Kolom kunci: status_kelayakan, alasan_tidak_layak,
                      masuk_dataset_final.

04_laporan/           Laporan lengkap (PDF untuk dibaca, Markdown untuk
                      versi sumbernya). Memuat metode, temuan, kendali mutu,
                      verifikasi genotipe, dan keterbatasan.

05_lampiran_audit/    {len(ditolak)} berkas .gb kandidat yang DITOLAK penyaringan.
                      BUKAN bagian dataset. Disertakan sebagai jejak audit
                      agar setiap penolakan dapat ditelusuri lewat kolom
                      alasan_tidak_layak di metadata.

06_skrip/             Skrip pipeline yang dirujuk di dalam laporan, agar
                      rujukan seperti "scripts/02_curate.py" dapat dibuka.

====================================================================
CATATAN PENTING UNTUK PEMERIKSA
====================================================================
1. Jumlah rekaman dataset adalah {qc['jumlah_dataset_final']}, yaitu isi 01_genbank_final/.
   Folder 05_lampiran_audit/ berisi kandidat yang DITOLAK dan tidak boleh
   ikut dihitung.

2. Gen S HBV terletak SELURUHNYA di dalam gen P pada {qc['jumlah_dataset_final']}/{qc['jumlah_dataset_final']} rekaman,
   hanya berbeda reading frame. Akibatnya Dataset 3 (gen S + gen P) tidak
   menambah cakupan nukleotida dibanding Dataset 2 ({qc['overlap']['persen_bp_redundan']}% isinya
   duplikasi). Implikasinya terhadap perumusan kesimpulan dibahas di
   bagian 6 laporan. Dataset 4 ditambahkan untuk mengatasi hal ini.

3. Genotipe G ({dist.get('G', 0)}) dan H ({dist.get('H', 0)}) jumlahnya sedikit karena memang
   langka di NCBI, bukan karena penyaringan yang terlalu ketat. Rinciannya
   di bagian 4.1 dan 7.2 laporan.

4. Sembilan rekaman dibuang karena terbukti rekombinan antar-genotipe
   lewat verifikasi independen; tiga rekaman lain yang tidak konkordan
   sengaja DITANDAI dan tidak dibuang. Alasannya di bagian 5.5 laporan.
""", encoding="utf-8")
    print("  README.txt dibuat")

    arsip = ROOT / "pengumpulan_tugas1_jst_dataset_hbv"
    print("\n  membuat arsip ZIP...")
    zip_path = Path(shutil.make_archive(str(arsip), "zip", root_dir=OUT))
    print(f"  {zip_path.name}: {zip_path.stat().st_size / 1024 / 1024:.1f} MB")

    total = sum(p.stat().st_size for p in OUT.rglob("*") if p.is_file())
    print(f"\n  folder pengumpulan/: {total / 1024 / 1024:.1f} MB "
          f"({sum(1 for p in OUT.rglob('*') if p.is_file())} berkas)")
    print(f"  arsip ZIP          : {zip_path.stat().st_size / 1024 / 1024:.1f} MB")


if __name__ == "__main__":
    main()

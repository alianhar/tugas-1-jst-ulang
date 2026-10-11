"""Driver pipeline: jalankan seluruh tahap sampai daftar pengecualian konvergen.

Tahap 2 (kurasi) dan tahap 4 (verifikasi) saling bergantung: tahap 4
menghasilkan daftar pengecualian berbasis bukti yang dibaca tahap 2, dan
membuang satu rekaman dapat membuat rekaman lain naik menjadi layak (slot
deduplikasi terbebaskan). Karena itu keduanya diulang sampai tahap 4 tidak
lagi menambahkan pengecualian baru.

Pemakaian:
    .venv/Scripts/python scripts/run_pipeline.py            # tahap 2-4
    .venv/Scripts/python scripts/run_pipeline.py --unduh    # termasuk tahap 1
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY = sys.executable
MAX_LINTASAN = 6


def jalankan(skrip: str) -> str:
    print(f"\n{'#' * 70}\n# {skrip}\n{'#' * 70}")
    hasil = subprocess.run([PY, "-I", str(ROOT / "scripts" / skrip)],
                           capture_output=True, text=True, cwd=ROOT)
    if hasil.returncode != 0:
        print(hasil.stdout[-3000:])
        print(hasil.stderr[-3000:])
        raise SystemExit(f"{skrip} gagal (exit {hasil.returncode})")
    return hasil.stdout


def ringkas(keluaran: str, kunci: tuple[str, ...]) -> None:
    for line in keluaran.splitlines():
        if any(k in line for k in kunci):
            print(line)


def main() -> None:
    if "--unduh" in sys.argv:
        ringkas(jalankan("01_collect.py"), ("selesai:", "total UID"))

    laporan = ROOT / "data" / "laporan_verifikasi_genotipe.json"
    for lintasan in range(1, MAX_LINTASAN + 1):
        print(f"\n\n{'=' * 70}\nLINTASAN {lintasan}\n{'=' * 70}")
        ringkas(jalankan("02_curate.py"),
                ("dikecualikan oleh", "rekaman layak sebelum", "terpilih",
                 "dataset", "PERINGATAN"))
        ringkas(jalankan("04_verify_genotype.py"),
                ("KONKORDANSI", "1-NN:", "BARU pada lintasan",
                 "total daftar", ">>>"))
        data = json.loads(laporan.read_text())
        if data.get("pengecualian", {}).get("konvergen"):
            print(f"\n>>> Konvergen pada lintasan {lintasan}.")
            break
    else:
        print(f"\n!! belum konvergen setelah {MAX_LINTASAN} lintasan")

    print(f"\n{'#' * 70}\n# kendali mutu akhir\n{'#' * 70}")
    keluaran = jalankan("03_analyze.py")
    ringkas(keluaran, ("dataset", "total ", "rekaman dengan", "Dataset 4 --",
                       "gen S (n=", "gen P (n=", "identik persis",
                       "masuk dataset", "layak  "))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""fetch_phenology_data.py -- put the Harvard Forest and Hyytiala phenology data in data/.

The data are not in the repository. Sources (all open; cite them if you use the data):

Harvard Forest, Massachusetts (42.54 N), Harvard Forest Data Archive:
  HF001  Boose & Gould, Fisher Meteorological Station: daily mean air temperature.
  HF003  O'Keefe, Phenology of woody species: percent leaf fall of tagged trees.
  HF069  Munger & Wofsy, EMS tower: litter baskets, by species and date.
  MODIS  MCD15A3H leaf area index (Myneni et al.), the ORNL DAAC fixed subset for the site; the
         median of the 3 x 3 pixels around the tower is kept.
Hyytiala, Finland (61.85 N), ICOS ETC Level 2 archive for FI-Hyy (CC BY 4.0):
  the daily FLUXNET file (air temperature TA_F) and the ancillary file (needle litter).

Usage:
  python fetch_phenology_data.py              # downloads ~180 MB (the ICOS archive), keeps ~20 MB
"""
import argparse
import http.cookiejar
import io
import json
import os
import statistics
import time
import urllib.error
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
HF = "https://harvardforest.fas.harvard.edu/data"
HF_FILES = {"hf001-06-daily-m.csv": f"{HF}/p00/hf001/hf001-06-daily-m.csv",
            "hf003-04-fall.csv":    f"{HF}/p00/hf003/hf003-04-fall.csv",
            "hf069-05-litter.csv":  f"{HF}/p06/hf069/hf069-05-litter.csv"}
MODIS = ("https://modis.ornl.gov/rst/api/v1/MCD15A3H/us_massachusetts_harvard_forest/"
         "subsetFiltered?band=Lai_500m&startDate=A{y}001&endDate=A{y}366")
MODIS_YEARS = range(2003, 2024)
ICOS_ARCHIVE = "XkewKEuf9Bv592orXDy4QpFi"          # ICOSETC_FI-Hyy_ARCHIVE_INTERIM_L2.zip
ICOS = "https://data.icos-cp.eu/licence_accept?ids=%5B%22{id}%22%5D"
ICOS_MEMBERS = {"ICOSETC_FI-Hyy_FLUXNET_DD_INTERIM_L2.csv": "hyytiala_fluxnet_dd.csv",
                "ICOSETC_FI-Hyy_ANCILLARY_INTERIM_L2.csv":  "hyytiala_ancillary.csv"}
#----- The Harvard Forest archive refuses urllib's default User-Agent, so name the client. --------#
HEADERS = {"User-Agent": "MEDS-example02 (urllib)"}


def get(url, opener=urllib.request.urlopen, **headers):
    """Open a URL, retrying a few times on server errors (the MODIS service has transient 502s)."""
    request = urllib.request.Request(url, headers={**HEADERS, **headers})
    for attempt in range(1, 5):
        try:
            return opener(request, timeout=1800)
        except urllib.error.HTTPError as err:
            if err.code < 500 or attempt == 4:
                raise
            print(f"  server error {err.code}; retrying in {15 * attempt} s")
            time.sleep(15 * attempt)


def modis_lai(dest):
    """One row per MODIS date: the median LAI of the 3 x 3 pixels at the centre of the subset."""
    rows = []
    for y in MODIS_YEARS:
        with get(MODIS.format(y=y), Accept="application/json") as r:
            subset = json.load(r)["subset"]
        for rec in subset:
            n = int(round(len(rec["data"]) ** 0.5))
            grid = [rec["data"][i * n:(i + 1) * n] for i in range(n)]
            c = n // 2
            vals = [float(v) for row in grid[c - 1:c + 2] for v in row[c - 1:c + 2]
                    if isinstance(v, (int, float))]
            rows.append((rec["calendar_date"], statistics.median(vals) if vals else ""))
        print(f"  MODIS {y}: {len(subset)} dates")
    with open(dest, "w") as fh:
        fh.write("date,lai\n")
        fh.writelines(f"{d},{v}\n" for d, v in sorted(rows))


def icos_files(data_dir):
    """The ICOS archive needs its licence accepted in a cookie session; keep two files of it."""
    opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar()))
    with get(ICOS.format(id=ICOS_ARCHIVE), opener.open) as r:
        archive = zipfile.ZipFile(io.BytesIO(r.read()))
    for member, name in ICOS_MEMBERS.items():
        with archive.open(member) as src, open(os.path.join(data_dir, name), "wb") as dst:
            dst.write(src.read())


def main(argv=None):
    ap = argparse.ArgumentParser(description="Fetch the example02 phenology data into data/.")
    ap.add_argument("--data-dir", default=os.path.join(HERE, "data"))
    args = ap.parse_args(argv)
    os.makedirs(args.data_dir, exist_ok=True)
    path = lambda name: os.path.join(args.data_dir, name)          # noqa: E731
    for name, url in HF_FILES.items():
        if not os.path.exists(path(name)):
            print(f"{name}: downloading from the Harvard Forest Data Archive")
            with get(url) as r, open(path(name), "wb") as fh:
                fh.write(r.read())
    if not os.path.exists(path("harvard_modis_lai.csv")):
        print("harvard_modis_lai.csv: MODIS MCD15A3H from the ORNL DAAC")
        modis_lai(path("harvard_modis_lai.csv"))
    if not all(os.path.exists(path(n)) for n in ICOS_MEMBERS.values()):
        print("hyytiala_*.csv: the ICOS FI-Hyy archive (~160 MB download)")
        icos_files(args.data_dir)
    print(f"data in {args.data_dir}")


if __name__ == "__main__":
    main()

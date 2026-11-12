#!/usr/bin/env python3
"""January weather for NaturalConvectionPCM from measured data: PVGIS TMY.

Source: PVGIS 5.3 typical meteorological year (European Commission, Joint
Research Centre; ISO 15927-4 month selection; satellite irradiance + ERA5
reanalysis), downloaded from https://re.jrc.ec.europa.eu/api/v5_3/tmy.
Every month of a TMY is a real month of a real year (the year is written in
the header). Default site: Santiago, Quinta Normal (-33.445, -70.683).

Writes the -pcm_weather file of paine_weather.py: hour (solar time), I [W/m^2]
on the plane of the roof, T_a [C] (2 m), w [m/s] (10 m), every 30 min.
Roof irradiance: DNI cos(theta) + isotropic diffuse + ground reflection
(albedo 0.2), sun position at each record in solar time.

  python3 pvgis_weather.py -o santiago_enero.csv              # mean January day
  python3 pvgis_weather.py --dia caluroso -o santiago_cal.csv # real day with max T
  python3 pvgis_weather.py --dia 15 -o santiago_15ene.csv     # 15 January of the TMY
  python3 pvgis_weather.py --input tmy.json ...               # already downloaded JSON

  curl -o tmy.json "https://re.jrc.ec.europa.eu/api/v5_3/tmy?lat=-33.445&lon=-70.683&outputformat=json"
"""
import argparse
import json
import math
import urllib.request

import numpy as np

URL = "https://re.jrc.ec.europa.eu/api/{v}/tmy?lat={lat}&lon={lon}&outputformat=json"


def load(a):
    if a.input:
        with open(a.input) as f:
            return json.load(f), a.input
    err = None
    for v in ("v5_3", "v5_2"):
        url = URL.format(v=v, lat=a.lat, lon=a.lon)
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                data = json.load(r)
            if a.save:
                with open(a.save, "w") as f:
                    json.dump(data, f)
            return data, url
        except Exception as e:  # noqa: BLE001
            err = e
    raise SystemExit("No se pudo descargar PVGIS (%s). Descarga el JSON con curl y usa --input." % err)


def field(rec, *names):
    for n in names:
        if n in rec:
            return float(rec[n])
    raise KeyError(names)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lat", type=float, default=-33.445)
    ap.add_argument("--lon", type=float, default=-70.683)
    ap.add_argument("--name", default="Santiago (Quinta Normal)")
    ap.add_argument("--month", type=int, default=1)
    ap.add_argument("--dia", default="promedio", help="promedio | caluroso | despejado | day of the month")
    ap.add_argument("--tilt", type=float, default=25.0)
    ap.add_argument("--facing", choices=["north", "south", "east", "west", "flat"], default="north")
    ap.add_argument("--albedo", type=float, default=0.2)
    ap.add_argument("--input", help="PVGIS TMY JSON already downloaded")
    ap.add_argument("--save", help="also save the downloaded JSON here")
    ap.add_argument("-o", "--output", default="santiago_enero.csv")
    a = ap.parse_args()

    data, src = load(a)
    recs = data["outputs"]["tmy_hourly"]
    lat = float(data.get("inputs", {}).get("location", {}).get("latitude", a.lat))
    lon = float(data.get("inputs", {}).get("location", {}).get("longitude", a.lon))
    year = "?"
    for m in data["outputs"].get("months_selected", []):
        if int(m.get("month", 0)) == a.month:
            year = m.get("year", "?")
    db = data.get("inputs", {}).get("meteo_data", {}).get("radiation_db", "PVGIS")

    d2r = math.pi / 180
    phi = lat * d2r
    beta = 0.0 if a.facing == "flat" else a.tilt * d2r
    gamma = {"south": 0.0, "north": 180.0, "east": -90.0, "west": 90.0, "flat": 0.0}[a.facing] * d2r

    rows = []
    for r in recs:
        ts = r.get("time(UTC)") or r.get("time")
        mo, dd, hh, mi = int(ts[4:6]), int(ts[6:8]), int(ts[9:11]), int(ts[11:13])
        if mo != a.month:
            continue
        rows.append((dd, hh + mi / 60.0, field(r, "G(h)"), field(r, "Gb(n)"), field(r, "Gd(h)"),
                     field(r, "T2m"), field(r, "WS10m")))
    if not rows:
        raise SystemExit("month %d not in the TMY" % a.month)
    R = np.array(rows)
    day, utc, ghi, dni, dhi, T, W = R.T
    nday = np.array([sum([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31][: a.month - 1]) + int(d) for d in day])
    B = 2 * math.pi * (nday - 81) / 364.0
    eot = 9.87 * np.sin(2 * B) - 7.53 * np.cos(B) - 1.5 * np.sin(B)  # min
    delta = 23.45 * d2r * np.sin(2 * math.pi * (284 + nday) / 365.0)

    def geometry(shift):
        solar = utc + shift + lon / 15.0 + eot / 60.0  # h, may leave [0, 24)
        w = (solar - 12.0) * 15.0 * d2r
        cz = np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.cos(w)
        ct = (np.sin(delta) * np.sin(phi) * np.cos(beta)
              - np.sin(delta) * np.cos(phi) * np.sin(beta) * np.cos(gamma)
              + np.cos(delta) * np.cos(phi) * np.cos(beta) * np.cos(w)
              + np.cos(delta) * np.sin(phi) * np.sin(beta) * np.cos(gamma) * np.cos(w)
              + np.cos(delta) * np.sin(beta) * np.sin(gamma) * np.sin(w))
        return solar, cz, ct

    # Time stamp convention (instant vs. start/end of the hour): the one where
    # GHI = DNI cos z + DHI holds best.
    best = None
    for shift in (-0.5, -1.0 / 6, 0.0, 1.0 / 6, 0.5):
        _, cz, _ = geometry(shift)
        res = float(np.sum(np.abs(ghi - (dni * np.maximum(cz, 0) + dhi))))
        if best is None or res < best[0]:
            best = (res, shift)
    shift = best[1]
    solar, cz, ct = geometry(shift)
    closure = best[0] / max(float(np.sum(ghi)), 1.0)
    poa = dni * np.maximum(ct, 0) + dhi * (1 + math.cos(beta)) / 2 + ghi * a.albedo * (1 - math.cos(beta)) / 2
    poa[cz <= 0] = np.maximum(dhi[cz <= 0], 0) * (1 + math.cos(beta)) / 2

    # Continuous solar-time axis over the month.
    t = (day - 1) * 24.0 + solar
    order = np.argsort(t)
    t, poa, T, W, ghi = t[order], poa[order], T[order], W[order], ghi[order]
    grid = np.arange(48) * 0.5
    days = sorted(set(int(d) for d in day))
    ok = [d for d in days if (d - 1) * 24.0 >= t[0] and (d - 1) * 24.0 + 23.5 <= t[-1]]

    def at(d, y):
        return np.interp((d - 1) * 24.0 + grid, t, y)

    if a.dia == "promedio":
        I = np.mean([at(d, poa) for d in ok], axis=0)
        Ta = np.mean([at(d, T) for d in ok], axis=0)
        w = np.mean([at(d, W) for d in ok], axis=0)
        label = "mean of %d days" % len(ok)
    else:
        if a.dia == "caluroso":
            d = max(ok, key=lambda d: at(d, T).max())
        elif a.dia == "despejado":
            d = max(ok, key=lambda d: at(d, ghi).sum())
        else:
            d = int(a.dia)
            if d not in ok:
                raise SystemExit("day %d not complete in solar time (days %d-%d)" % (d, ok[0], ok[-1]))
        I, Ta, w = at(d, poa), at(d, T), at(d, W)
        label = "day %d" % d

    daily_poa = float(np.sum(I) * 0.5 / 1000)
    with open(a.output, "w") as f:
        f.write("# %s, %s month %d of PVGIS TMY (%s, JRC; real year %s), %s; lat %.3f lon %.3f;"
                " tilt %.1f deg facing %s; plane of roof %.2f kWh/m2/day; T_a 2 m, w 10 m\n"
                % (a.name, label, a.month, db, year, src.split("?")[0], lat, lon, a.tilt, a.facing, daily_poa))
        f.write("# hour, I [W/m^2] (plane of roof), T_a [C], w [m/s]\n")
        for row in zip(grid, I, Ta, w):
            f.write("%.2f,%.1f,%.2f,%.2f\n" % row)
    print("%s: %s, year %s; roof %.2f kWh/m2/day, peak %.0f W/m2; T_a %.1f-%.1f C; w %.1f-%.1f m/s;"
          " stamp shift %+.2f h (GHI closure %.1f%%)"
          % (a.output, label, year, daily_poa, I.max(), Ta.min(), Ta.max(), w.min(), w.max(),
             shift, 100 * closure))


if __name__ == "__main__":
    main()

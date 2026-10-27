#!/usr/bin/env python3
"""Typical clear January day in Paine (RM, Chile) for NaturalConvectionPCM.

Writes the -pcm_weather file: hour (solar time), I [W/m^2] on the plane of the
roof, T_a [C], w [m/s], every 30 min. Periodic in 24 h.

Site: Paine, -33.81, -70.74, 370 m (Santiago basin, southern end).
  - Irradiance: Hottel (1976) clear-sky beam transmittance for 0.37 km,
    mid-latitude summer, Liu-Jordan diffuse; scaled so that the daily global
    horizontal irradiation is --ghi kWh/m^2/day (default 8.2: January, clear
    sky in the Santiago basin; annual 1750-1900 kWh/m^2 per Explorador Solar).
    Transposed to the roof (tilt --tilt, azimuth --facing) with the isotropic
    sky model, ground albedo 0.2.
  - T_a: asymmetric cosine between --tmin at sunrise and --tmax at 15 h
    (January means near Paine: ~31 / ~13 C).
  - w: valley breeze, --wnight at night to --wday at 16 h.

    python3 paine_weather.py --tilt 25 --facing north -o paine_enero.csv
"""
import argparse
import math


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lat", type=float, default=-33.81)
    ap.add_argument("--alt", type=float, default=0.37, help="km")
    ap.add_argument("--day", type=int, default=15, help="day of the year")
    ap.add_argument("--tilt", type=float, default=25.0, help="roof slope, deg")
    ap.add_argument("--facing", choices=["north", "south", "east", "west", "flat"], default="north")
    ap.add_argument("--ghi", type=float, default=8.2, help="daily GHI, kWh/m^2/day")
    ap.add_argument("--albedo", type=float, default=0.2)
    ap.add_argument("--tmin", type=float, default=13.0)
    ap.add_argument("--tmax", type=float, default=31.0)
    ap.add_argument("--wnight", type=float, default=1.5)
    ap.add_argument("--wday", type=float, default=4.0)
    ap.add_argument("-o", "--output", default="paine_enero.csv")
    a = ap.parse_args()

    d2r = math.pi / 180.0
    phi = a.lat * d2r
    n = a.day
    delta = 23.45 * d2r * math.sin(2 * math.pi * (284 + n) / 365.0)
    gon = 1367.0 * (1 + 0.033 * math.cos(2 * math.pi * n / 365.0))
    A = a.alt
    a0 = 0.97 * (0.4237 - 0.00821 * (6 - A) ** 2)
    a1 = 0.99 * (0.5055 + 0.00595 * (6.5 - A) ** 2)
    k = 1.02 * (0.2711 + 0.01858 * (2.5 - A) ** 2)

    beta = 0.0 if a.facing == "flat" else a.tilt * d2r
    # Duffie-Beckman surface azimuth: 0 facing the equator... of the north
    # hemisphere (south), so a north-facing roof is 180 deg.
    gamma = {"south": 0.0, "north": 180.0, "east": -90.0, "west": 90.0, "flat": 0.0}[a.facing] * d2r

    ws = math.acos(-math.tan(phi) * math.tan(delta))  # sunset hour angle
    sunrise = 12.0 - ws / d2r / 15.0

    def raw(hour):
        w = (hour - 12.0) * 15.0 * d2r
        cz = math.sin(phi) * math.sin(delta) + math.cos(phi) * math.cos(delta) * math.cos(w)
        if cz <= 0.0:
            return 0.0, 0.0, 0.0
        tb = a0 + a1 * math.exp(-k / cz)
        td = 0.271 - 0.294 * tb
        dni = gon * tb
        dhi = gon * cz * td
        ghi = dni * cz + dhi
        ct = (math.sin(delta) * math.sin(phi) * math.cos(beta)
              - math.sin(delta) * math.cos(phi) * math.sin(beta) * math.cos(gamma)
              + math.cos(delta) * math.cos(phi) * math.cos(beta) * math.cos(w)
              + math.cos(delta) * math.sin(phi) * math.sin(beta) * math.cos(gamma) * math.cos(w)
              + math.cos(delta) * math.sin(beta) * math.sin(gamma) * math.sin(w))
        poa = dni * max(ct, 0.0) + dhi * (1 + math.cos(beta)) / 2 + ghi * a.albedo * (1 - math.cos(beta)) / 2
        return ghi, poa, cz

    # Scale to the target daily GHI (fine quadrature, 1 min).
    m = 1440
    daily = sum(raw(24.0 * i / m)[0] for i in range(m)) * 24.0 / m / 1000.0
    scale = a.ghi / daily

    def ta(hour):
        tpk = 15.0
        h = hour if hour >= sunrise else hour + 24.0
        if h <= tpk:  # rising branch
            s = (h - sunrise) / (tpk - sunrise)
            return a.tmin + (a.tmax - a.tmin) * 0.5 * (1 - math.cos(math.pi * s))
        s = (h - tpk) / (24.0 + sunrise - tpk)
        return a.tmax - (a.tmax - a.tmin) * 0.5 * (1 - math.cos(math.pi * s))

    def wind(hour):
        return a.wnight + (a.wday - a.wnight) * max(0.0, math.cos(math.pi * (hour - 16.0) / 12.0)) ** 2

    rows = []
    poa_day = 0.0
    for i in range(48):
        h = 0.5 * i
        _, poa, _ = raw(h)
        rows.append((h, scale * poa, ta(h), wind(h)))
    for i in range(m):
        poa_day += scale * raw(24.0 * i / m)[1]
    poa_day *= 24.0 / m / 1000.0

    with open(a.output, "w") as f:
        f.write("# Paine (RM) clear January day: lat %.2f, day %d, tilt %.1f deg facing %s;"
                " GHI %.2f kWh/m2/day, plane of roof %.2f kWh/m2/day\n"
                % (a.lat, n, a.tilt, a.facing, a.ghi, poa_day))
        f.write("# hour, I [W/m^2] (plane of roof), T_a [C], w [m/s]\n")
        for r in rows:
            f.write("%.2f,%.1f,%.2f,%.2f\n" % r)
    print("%s: sunrise %.2f h, peak I %.0f W/m^2, roof %.2f kWh/m^2/day (GHI %.2f)"
          % (a.output, sunrise, max(r[1] for r in rows), poa_day, a.ghi))


if __name__ == "__main__":
    main()

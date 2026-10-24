#!/usr/bin/env python3
"""The Bolund roughness map for bc.surface_z0_file, from the Bolund terrain
as a terrain.type = xyz file (bolund.xyz of windninja2019.py prepare): one
"x y z0" line per terrain point; lines that do not hold three numbers are skipped.

The values and rules are those of the Bolund blind comparison (Bechmann et
al. 2011; Bolund_Validation.pdf in Bolund_Measurements.zip, Section 2 and
Figure 2): water z0 = 0.0003 m where the terrain is at the water level
(z = 0.75 m in the topography files) or below; land (the grass of Bolund)
z0 = 0.015 m elsewhere; and land roughness for x > 327 m in Bolund's frame,
where the coast begins east of Bolund (the extended terrain is flat water
there).

Usage: bolund_z0_map.py IN.xyz OUT.xyz [--water-level 0.75] [--z0-water 0.0003]
                        [--z0-land 0.015] [--land-x-min 327 | --land-x-min none]
"""
import argparse
import sys


def classify(z, x, water_level, land_x_min):
    """True for water: at or below the water level, and not in the eastern land."""
    if land_x_min is not None and x > land_x_min:
        return False
    return z <= water_level


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("inp", help="terrain points, 'x y z' per line ('#' comments)")
    ap.add_argument("out", help="roughness points, 'x y z0' per line")
    ap.add_argument("--water-level", type=float, default=0.75, help="water surface height (m), default 0.75")
    ap.add_argument("--z0-water", type=float, default=0.0003, help="water roughness (m), default 0.0003")
    ap.add_argument("--z0-land", type=float, default=0.015, help="land roughness (m), default 0.015")
    ap.add_argument("--land-x-min", default="327",
                    help="land roughness for x above this (m; the coast east of Bolund), or 'none'; default 327")
    a = ap.parse_args(argv)
    land_x_min = None if a.land_x_min == "none" else float(a.land_x_min)
    for name, v in (("--z0-water", a.z0_water), ("--z0-land", a.z0_land)):
        if not v > 0:
            sys.exit(f"bolund_z0_map.py: {name} = {v}: must be > 0")
    n = [0, 0]
    with open(a.inp) as fi, open(a.out, "w") as fout:
        fout.write(f"# Bolund roughness from {a.inp}: z0 = {a.z0_water} m (water, z <= {a.water_level} m"
                   + (f", x <= {land_x_min} m" if land_x_min is not None else "") + f"), {a.z0_land} m elsewhere"
                   " (Bechmann et al. 2011)\n")
        for line in fi:
            tok = line.split("#")[0].replace(",", " ").split()
            if len(tok) != 3:
                continue
            try:
                x, y, z = (float(t) for t in tok)
            except ValueError:
                continue
            water = classify(z, x, a.water_level, land_x_min)
            n[0 if water else 1] += 1
            fout.write(f"{tok[0]} {tok[1]} {a.z0_water if water else a.z0_land}\n")
    if n[0] + n[1] == 0:
        sys.exit(f"bolund_z0_map.py: {a.inp} holds no 'x y z' point")
    print(f"{a.out}: {n[0]} water points (z0 {a.z0_water} m), {n[1]} land points (z0 {a.z0_land} m)")


if __name__ == "__main__":
    main(sys.argv[1:])

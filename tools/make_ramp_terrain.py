#!/usr/bin/env python3
"""A smooth periodic plateau for the 2-D terrain comparison:

    h(x) = (H/4) [1 + tanh((x + x_r)/L)] [1 - tanh((x - x_r)/L)],

a tanh ramp up from flat ground, a flat top and a ramp down (H = 100 m,
x_r = 750 m, L = 250 m: the maximum slope 0.2 of ERF's Witch of Agnesi), written
in ERF's structured terrain-file format (nx, ny, the x values, the y values,
then nx*ny heights contiguous in y) for `terrain.type = erf_grid` and ERF's
`erf.terrain_file_name`, on 129 x 129 points over [-2560, 2560] m.

    python3 tools/make_ramp_terrain.py ramp2d.txt
"""
import math
import sys

N, H, XR, L, X0, X1 = 129, 100.0, 750.0, 250.0, -2560.0, 2560.0


def height(x):
    return H / 4 * (1 + math.tanh((x + XR) / L)) * (1 - math.tanh((x - XR) / L))


def main(path):
    xs = [X0 + (X1 - X0) * i / (N - 1) for i in range(N)]
    with open(path, "w") as fh:
        fh.write("%d\n%d\n" % (N, N))
        for v in xs:
            fh.write("%.6f\n" % v)
        for v in xs:
            fh.write("%.6f\n" % v)
        for x in xs:
            h = height(x)
            for _ in xs:
                fh.write("%.6f\n" % h)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "ramp2d.txt")

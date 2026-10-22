#!/usr/bin/env python3
"""A single-band GeoTIFF (float32/float64, uncompressed or LZW, strips, no
predictor) as a GUST terrain.type = xyz file, without GDAL: one
"x y z" line per pixel centre (every STRIDE-th pixel), in the TIFF's own
frame (ModelTiepoint + PixelScale), z shifted by -ZSHIFT. Pixels at or below
-1e30 (no data) are skipped.

Usage: tiff_to_xyz.py IN.tif OUT.xyz [STRIDE] [ZSHIFT]
"""
import struct
import sys

import numpy as np


def lzw_decode(data):
    out = bytearray()
    table = [bytes([i]) for i in range(256)] + [b"", b""]
    width, bitpos, prev = 9, 0, None
    nbits = len(data) * 8
    while bitpos + width <= nbits:
        byte, bit = divmod(bitpos, 8)
        chunk = int.from_bytes(data[byte:byte + 3].ljust(3, b"\0"), "big")
        code = (chunk >> (24 - bit - width)) & ((1 << width) - 1)
        bitpos += width
        if code == 256:  # clear
            table = table[:258]
            width, prev = 9, None
            continue
        if code == 257:  # end of information
            break
        if prev is None:
            entry = table[code]
        elif code < len(table):
            entry = table[code]
            table.append(prev + entry[:1])
        else:
            entry = prev + prev[:1]
            table.append(entry)
        out += entry
        prev = entry
        if len(table) + 1 >= (1 << width) and width < 12:  # early change
            width += 1
    return bytes(out)


def read_tiff(path):
    f = open(path, "rb").read()
    assert f[:2] == b"II", "little-endian TIFF only"
    off = struct.unpack("<I", f[4:8])[0]
    n = struct.unpack("<H", f[off:off + 2])[0]
    size = {3: 2, 4: 4, 12: 8}
    fmt = {3: "H", 4: "I", 12: "d"}
    tags = {}
    for i in range(n):
        tag, typ, cnt, val = struct.unpack("<HHII", f[off + 2 + 12 * i:off + 14 + 12 * i])
        if typ not in fmt:
            continue
        nb = size[typ] * cnt
        raw = f[val:val + nb] if nb > 4 else f[off + 10 + 12 * i:off + 10 + 12 * i + nb]
        tags[tag] = struct.unpack("<" + fmt[typ] * cnt, raw)
    nx, ny, bps = tags[256][0], tags[257][0], tags[258][0]
    comp, pred = tags[259][0], tags.get(317, (1,))[0]
    assert comp in (1, 5) and pred == 1 and bps in (32, 64), (comp, pred, bps)
    rows = []
    for o, c in zip(tags[273], tags[279]):
        chunk = f[o:o + c]
        rows.append(lzw_decode(chunk) if comp == 5 else chunk)
    a = np.frombuffer(b"".join(rows), dtype="<f8" if bps == 64 else "<f4")[:nx * ny].reshape(ny, nx)
    sx, sy = tags[33550][0], tags[33550][1]
    tie = tags[33922]
    x0, y0 = tie[3] - tie[0] * sx, tie[4] + tie[1] * sy  # the corner of pixel (0, 0)
    return a, x0, y0, sx, sy


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    stride = int(argv[2]) if len(argv) > 2 else 1
    zshift = float(argv[3]) if len(argv) > 3 else 0.0
    a, x0, y0, sx, sy = read_tiff(argv[0])
    ny, nx = a.shape
    n = 0
    with open(argv[1], "w") as out:
        out.write(f"# {argv[0]}: {nx} x {ny} pixels of {sx} x {sy} m, stride {stride}, z - {zshift}\n")
        for j in range(0, ny, stride):
            y = y0 - (j + 0.5) * sy
            for i in range(0, nx, stride):
                z = a[j, i]
                if z > -1e30:
                    out.write(f"{x0 + (i + 0.5) * sx:.3f} {y:.3f} {z - zshift:.3f}\n")
                    n += 1
    good = a[a > -1e30]
    print(f"{argv[1]}: {n} points; x {x0 + 0.5 * sx:.2f} .. {x0 + (nx - 0.5) * sx:.2f}, "
          f"y {y0 - (ny - 0.5) * sy:.2f} .. {y0 - 0.5 * sy:.2f}, z {good.min():.2f} .. {good.max():.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

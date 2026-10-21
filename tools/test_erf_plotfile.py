"""Self-tests for extern_erf/erf_plotfile.py: the size and byte order of a FAB's data.

    cd tools && python3 -m unittest -v test_erf_plotfile

Each test writes a one-box AMReX plotfile by hand (Header, Level_0/Cell_H,
Level_0/Cell_D_00000) with known values, in double or float precision and
little- or big-endian byte order, and reads it back with read_fields and
planar_averages. The FAB header carries the real descriptor as AMReX writes
it: the floating-point format first (8 fields for float and double alike),
then the byte order, one entry per byte. Before the local change the reader
took the format's field count as the size, so float data was read as doubles.
"""
import os
import pathlib
import struct
import sys
import tempfile
import unittest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "extern_erf"))
import erf_plotfile  # noqa: E402

FORMAT = {8: "(64 11 52 0 1 12 0 1023)", 4: "(32 8 23 0 1 9 0 127)"}
NX, NY, NZ = 3, 2, 4
NAMES = ["u", "v"]


def value(c, i, j, k):
    """A field value that differs at every point and component, exact in float."""
    return (c + 1) * 100.0 + i + 0.25 * j + 0.5 * k


def write_plotfile(root, nbytes, little):
    root = pathlib.Path(root)
    (root / "Level_0").mkdir(parents=True)
    hi = f"({NX - 1},{NY - 1},{NZ - 1})"
    (root / "Header").write_text(
        "HyperCLaw-V1.1\n" + f"{len(NAMES)}\n" + "".join(n + "\n" for n in NAMES)
        + "3\n0\n0\n0 0 0\n" + f"{NX * 10.0} {NY * 10.0} {NZ * 10.0}\n\n"
        + f"((0,0,0) {hi} (0,0,0))\n0\n10 10 10\n0\n0\n")
    order = " ".join(str(b) for b in (range(nbytes, 0, -1) if little else range(1, nbytes + 1)))
    head = f"FAB (({8}, {FORMAT[nbytes]}),({nbytes}, ({order})))((0,0,0) {hi} (0,0,0)) {len(NAMES)}\n"
    fmt = ("<" if little else ">") + ("d" if nbytes == 8 else "f")
    data = b"".join(struct.pack(fmt, value(c, i, j, k))
                    for c in range(len(NAMES)) for k in range(NZ) for j in range(NY) for i in range(NX))
    (root / "Level_0" / "Cell_D_00000").write_bytes(head.encode() + data)
    (root / "Level_0" / "Cell_H").write_text(
        f"1\n1\n{len(NAMES)}\n0\n(1 0\n((0,0,0) {hi} (0,0,0))\n)\n1\nFabOnDisk: Cell_D_00000 0\n")
    return str(root)


class FabSizeAndOrder(unittest.TestCase):
    def check(self, nbytes, little):
        with tempfile.TemporaryDirectory() as d:
            pf = write_plotfile(os.path.join(d, "plt00000"), nbytes, little)
            _, out = erf_plotfile.read_fields(pf, NAMES)
            for c, name in enumerate(NAMES):
                for i in range(NX):
                    for j in range(NY):
                        for k in range(NZ):
                            self.assertEqual(out[name][i][j][k], value(c, i, j, k), (name, i, j, k))
            _, avg, _ = erf_plotfile.planar_averages(pf, ["v"])
            for k in range(NZ):
                want = sum(value(1, i, j, k) for i in range(NX) for j in range(NY)) / (NX * NY)
                self.assertAlmostEqual(avg["v"][k], want, places=12)

    def test_double_little_endian(self):
        self.check(8, True)

    def test_double_big_endian(self):
        self.check(8, False)

    def test_float_little_endian(self):
        self.check(4, True)

    def test_float_big_endian(self):
        self.check(4, False)


if __name__ == "__main__":
    unittest.main()

"""Self-tests for extern_erf/erf_plotfile.py: the size and byte order of a FAB's data.

    cd tools && python3 -m unittest -v test_erf_plotfile

Each test writes an AMReX plotfile by hand (Header, Level_0/Cell_H and two
FABs in Level_0/Cell_D_00000, the second at an offset, optionally with ghost
cells) with known values, in double or float precision and little- or
big-endian byte order, and reads it back with read_fields and
planar_averages; other byte orders must be refused. The FAB header carries the real descriptor as AMReX writes
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


def write_plotfile(root, nbytes, little, order=None, ng=0):
    """Two boxes split in x (i < 2 and i >= 2), both FABs in one data file, the
    second at a non-zero offset; with ng > 0 each FAB is grown by ng cells (its
    header box differs from the Cell_H box, as AMReX writes ghost cells)."""
    root = pathlib.Path(root)
    (root / "Level_0").mkdir(parents=True)
    hi = f"({NX - 1},{NY - 1},{NZ - 1})"
    (root / "Header").write_text(
        "HyperCLaw-V1.1\n" + f"{len(NAMES)}\n" + "".join(n + "\n" for n in NAMES)
        + "3\n0\n0\n0 0 0\n" + f"{NX * 10.0} {NY * 10.0} {NZ * 10.0}\n\n"
        + f"((0,0,0) {hi} (0,0,0))\n0\n10 10 10\n0\n0\n")
    if order is None:
        order = list(range(nbytes, 0, -1)) if little else list(range(1, nbytes + 1))
    fmt = ("<" if little else ">") + ("d" if nbytes == 8 else "f")
    boxes = [((0, 0, 0), (1, NY - 1, NZ - 1)), ((2, 0, 0), (NX - 1, NY - 1, NZ - 1))]
    blob, offsets = b"", []
    for lo, bhi in boxes:
        flo = [lo[d] - ng for d in range(3)]
        fhi = [bhi[d] + ng for d in range(3)]
        head = (f"FAB (({8}, {FORMAT[nbytes]}),({nbytes}, ({' '.join(str(b) for b in order)})))"
                f"(({flo[0]},{flo[1]},{flo[2]}) ({fhi[0]},{fhi[1]},{fhi[2]}) (0,0,0)) {len(NAMES)}\n")
        offsets.append(len(blob))
        blob += head.encode() + b"".join(
            struct.pack(fmt, value(c, i, j, k)) for c in range(len(NAMES)) for k in range(flo[2], fhi[2] + 1)
            for j in range(flo[1], fhi[1] + 1) for i in range(flo[0], fhi[0] + 1))
    (root / "Level_0" / "Cell_D_00000").write_bytes(blob)
    box_lines = "".join(f"(({lo[0]},{lo[1]},{lo[2]}) ({bhi[0]},{bhi[1]},{bhi[2]}) (0,0,0))\n" for lo, bhi in boxes)
    (root / "Level_0" / "Cell_H").write_text(
        f"1\n1\n{len(NAMES)}\n{ng}\n(2 0\n{box_lines})\n2\n"
        + "".join(f"FabOnDisk: Cell_D_00000 {o}\n" for o in offsets))
    return str(root)


class FabSizeAndOrder(unittest.TestCase):
    def check(self, nbytes, little, ng=0):
        with tempfile.TemporaryDirectory() as d:
            pf = write_plotfile(os.path.join(d, "plt00000"), nbytes, little, ng=ng)
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

    def test_float_with_ghost_cells(self):
        self.check(4, True, ng=1)

    def test_double_with_ghost_cells(self):
        self.check(8, False, ng=1)

    def test_other_byte_orders_are_refused(self):
        for nbytes, order in ((8, [2, 1, 4, 3, 6, 5, 8, 7]), (4, [3, 4, 1, 2])):
            with tempfile.TemporaryDirectory() as d:
                pf = write_plotfile(os.path.join(d, "plt00000"), nbytes, True, order=order)
                with self.assertRaises(ValueError) as e:
                    erf_plotfile.read_fields(pf, NAMES)
                self.assertIn(" ".join(str(b) for b in order), str(e.exception))


if __name__ == "__main__":
    unittest.main()

"""Self-tests for input_reference.py, on a throw-away docs/ and Source/.

    cd tools && python3 -m unittest -v test_input_reference
"""
import pathlib
import tempfile
import unittest

import input_reference as ir

DOC = """\
Input reference
===============

Quick reference
---------------

.. list-table::
   :header-rows: 1

   * - Input
     - Default
     - Description
   * - **Mesh**
     -
     -
   * - ``mesh.n_cell``
     - required
     - Cells
   * - ``bc.zlo``
     - ``slip_wall``
     - Ground

``mesh``
--------

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - Input
     - Default
     - Range and meaning
   * - ``mesh.n_cell``
     - required
     - Three integers ≥ 1.

``bc``
------

See :doc:`boundary_conditions` and :ref:`the wall <wall-law>`.

.. list-table::
   :header-rows: 1

   * - Input
     - Default
     - Range and meaning
   * - ``bc.zlo``
     - ``slip_wall``
     - The ground, :math:`f = 4\\pi\\sin\\varphi/T`, :math:`(1 + a_m\\zeta)`.
   * - ``bc.<face>_wall_velocity``
     - ``0 0 0``
     - Three values.
   * - ``bc.inflow_speed``, ``bc.inflow_zref``
     - 10 m/s, 10 m
     - power_law only.
   * - ``bc.rayleigh_zdamp``
     - 0 (off)
     - Height.
"""

INPUTS_H = """\
inline const std::vector<std::string>& owned_namespaces ()
{
    static const std::vector<std::string> ns{"mesh", "bc"};
    return ns;
}
"""

READER = """\
const char* const kFaceKeys[6] = {"xlo", "ylo", "zlo", "xhi", "yhi", "zhi"};
void read_mesh ()
{
    amrex::ParmParse pp("mesh");
    amrex::Vector<int> n;
    pp.queryarr("n_cell", n);
}
void read_bc ()
{
    amrex::ParmParse pp("bc");
    const char* keys[] = {"zlo", "xlo_wall_velocity", "xhi_wall_velocity", "ylo_wall_velocity",
                          "yhi_wall_velocity", "zlo_wall_velocity", "zhi_wall_velocity", "inflow_speed", "inflow_zref",
                          "rayleigh_zdamp"};
    (void)keys;
}
"""

LABELS = ".. _wall-law:\n\nWall law\n========\n"


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(text.encode())


class InputReference(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = pathlib.Path(self.tmp.name)
        write(self.root / "docs" / "inputs.rst", DOC)
        write(self.root / "docs" / "boundary_conditions.rst", LABELS)
        write(self.root / "Source" / "Inputs.H", INPUTS_H)
        write(self.root / "Source" / "Reader.cpp", READER)
        self.refresh()

    def tearDown(self):
        self.tmp.cleanup()

    def refresh(self):
        doc = (self.root / "docs" / "inputs.rst").read_text(encoding="utf-8")
        write(self.root / "Tests" / "inputs_reference", ir.render(doc, ir.labels(self.root)))

    def errors(self):
        return ir.check(self.root)[0]

    def edit(self, path, old, new):
        p = self.root / path
        s = p.read_text(encoding="utf-8")
        self.assertIn(old, s)
        write(p, s.replace(old, new, 1))

    def test_clean_tree_passes(self):
        errors, n = ir.check(self.root)
        self.assertEqual(errors, [])
        self.assertEqual(n, 11)   # n_cell, zlo, six wall velocities, speed, zref, zdamp

    def test_rendering(self):
        text = (self.root / "Tests" / "inputs_reference").read_text(encoding="utf-8")
        self.assertIn("# bc.zhi_wall_velocity = 0 0 0\n", text)
        self.assertIn("# bc.inflow_speed = 10    # m/s\n", text)
        self.assertIn("# bc.inflow_zref = 10    # m\n", text)
        self.assertIn("# bc.rayleigh_zdamp = 0    # off\n", text)
        self.assertIn("# mesh.n_cell =    # required\n", text)
        self.assertIn("f = 4π sin φ/T", text)
        self.assertIn("(1 + a_m ζ)", text)
        self.assertIn("the wall, docs/boundary_conditions.rst", text)
        self.assertIn("docs/boundary_conditions.rst and", text)

    def test_stale_reference(self):
        self.edit("Tests/inputs_reference", "# bc.zlo = slip_wall", "# bc.zlo = no_slip_wall")
        self.assertTrue(any("differs" in e for e in self.errors()))

    def test_crlf_checkout(self):
        for p in (self.root / "docs" / "inputs.rst", self.root / "Source" / "Reader.cpp"):
            p.write_bytes(p.read_bytes().replace(b"\n", b"\r\n"))
        ref = self.root / "Tests" / "inputs_reference"
        ref.write_bytes(ref.read_bytes().replace(b"\n", b"\r\n"))
        self.assertEqual(self.errors(), [])

    def test_key_read_but_not_documented(self):
        self.edit("Source/Reader.cpp", 'pp.queryarr("n_cell", n);', 'pp.queryarr("n_cell", n);\n    pp.query("dz0", n);')
        self.assertIn("mesh.dz0: read in Source/ but not in docs/inputs.rst", self.errors())

    def test_key_in_a_key_array(self):
        self.edit("Source/Reader.cpp", '"rayleigh_zdamp"};', '"rayleigh_zdamp", "sponge_width"};')
        self.assertIn("bc.sponge_width: read in Source/ but not in docs/inputs.rst", self.errors())

    def test_other_variable_and_temporary(self):
        self.edit("Source/Reader.cpp", "void read_bc ()",
                  'void other ()\n{\n    amrex::ParmParse q{"bc"};\n    q.query ("top", x);\n'
                  '    amrex::ParmParse("mesh").contains("lo");\n}\nvoid read_bc ()')
        errors = self.errors()
        self.assertIn("bc.top: read in Source/ but not in docs/inputs.rst", errors)
        self.assertIn("mesh.lo: read in Source/ but not in docs/inputs.rst", errors)

    def test_loop_list_takes_the_namespace_its_body_names(self):
        self.edit("Source/Reader.cpp", "void read_bc ()",
                  'void check ()\n{\n    for (const char* k : {"idw_power"}) {\n'
                  '        err(std::string("mesh.") + k + " has no effect");\n    }\n}\nvoid read_bc ()')
        self.assertIn("mesh.idw_power: read in Source/ but not in docs/inputs.rst", self.errors())

    def test_documented_key_not_read(self):
        self.edit("docs/inputs.rst", "   * - ``bc.rayleigh_zdamp``", "   * - ``bc.lusgs``\n     - ``1``\n     - x.\n"
                  "   * - ``bc.rayleigh_zdamp``")
        self.refresh()
        self.assertIn("bc.lusgs: in docs/inputs.rst but not read in Source/", self.errors())

    def test_placeholder_keys_must_be_read(self):
        self.edit("Source/Reader.cpp", '"zhi_wall_velocity", ', "")
        self.assertIn("bc.zhi_wall_velocity: in docs/inputs.rst but not read in Source/", self.errors())

    def test_malformed_tables(self):
        self.edit("docs/inputs.rst", "   * - ``bc.rayleigh_zdamp``", "    * - ``bc.rayleigh_zdamp``")
        self.assertTrue(any("not a row, cell or continuation" in e for e in self.errors()))

    def test_row_without_a_key_or_twice(self):
        self.edit("docs/inputs.rst", "   * - ``bc.rayleigh_zdamp``", "   * - rayleigh_zdamp")
        self.assertTrue(any("names no key" in e for e in self.errors()))
        self.setUp()
        self.edit("docs/inputs.rst", "   * - ``bc.rayleigh_zdamp``", "   * - ``bc.zlo``")
        self.assertTrue(any("in two rows" in e for e in self.errors()))

    def test_quick_reference(self):
        self.edit("docs/inputs.rst", "   * - ``bc.zlo``\n     - ``slip_wall``\n     - Ground",
                  "   * - ``bc.zlo``\n     - ``no_slip_wall``\n     - Ground")
        self.assertIn("bc.zlo: default 'no_slip_wall' in the quick-reference table, 'slip_wall' in its section",
                      self.errors())

    def test_namespace_without_section(self):
        self.edit("Source/Inputs.H", '{"mesh", "bc"}', '{"mesh", "bc", "turb"}')
        self.assertIn("namespace turb: no section in docs/inputs.rst", self.errors())

    def test_split_defaults(self):
        self.assertEqual(ir.split_defaults(2, "``1``, ``2``", {}), [("1", ""), ("2", "")])
        self.assertEqual(ir.split_defaults(1, "0 (off)", {}), [("0", "(off)")])
        self.assertEqual(ir.split_defaults(1, "1E-3 m", {}), [("1E-3", "m")])
        self.assertEqual(ir.split_defaults(1, "required", {}), [(None, "required")])


if __name__ == "__main__":
    unittest.main()

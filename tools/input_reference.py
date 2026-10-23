#!/usr/bin/env python3
"""The master input reference, Tests/inputs_reference, built from docs/inputs.rst.

Every key of the detailed tables of docs/inputs.rst (one section per
namespace) is written as a commented line ``# namespace.key = default`` with
its range and meaning below it, and each section's notes above its keys. The
file is a reference to copy key lines from, not a deck: every line is a
comment.

  input_reference.py --write   rewrite Tests/inputs_reference
  input_reference.py --check   (CTest) fail when the file differs from what
                               --write would give, when docs/inputs.rst does
                               not parse, or when its keys and the keys the
                               sources read disagree

What the check compares:

* The detailed tables. Every data row must name at least one key (a
  ``<face>`` placeholder stands for xlo, xhi, ylo, yhi, zlo and zhi); a key
  may appear once; a line inside a table that is not a row, a cell or a
  continuation is an error.
* The quick-reference table: every key in it is a key of the detailed tables,
  with the same default where both give a single literal.
* The sources (Source/, recursively). A key is read when a reader names it:
  ``ParmParse v("ns")`` (or ``{"ns"}``) and then ``v.query("key")``, ``v.get*``,
  ``v.contains``, ``v.count*`` on that variable; a temporary
  ``ParmParse("ns").query("key")``; a quoted word in a key array
  (``const char* kFooKeys[] = {...}`` or ``keys[]``) or a loop list (``for
  (const char* k : {...})``), in the namespace the loop's body names
  (``"ns." + k``), else that of the nearest ``ParmParse`` before it (or the
  file's only one); and a ``"ns.key"`` literal. Every key read must be in the
  docs and every key in the docs must be read. A key built by string
  concatenation (``kFaceKeys[o] + "_wall_velocity"``) is found only if its
  parts are listed as above.
* The owned namespaces (Source/Inputs.H, owned_namespaces) and every other
  namespace a ParmParse in Source/ reads (refavg) each have a detailed section.
"""

import argparse
import os
import pathlib
import re
import sys
import textwrap

WIDTH = 116
PLACEHOLDERS = {"face": ["xlo", "xhi", "ylo", "yhi", "zlo", "zhi"]}
KEY = r"[a-z_]+\.[A-Za-z0-9_<>]+"

HEADER = """\
# GUST master input reference (all inputs, with defaults, ranges and meaning).
#
# Generated from docs/inputs.rst by tools/input_reference.py --write; do not
# edit by hand. The CTest input_reference_current checks that this file, the
# docs and the keys read in Source/ agree.
#
# Every line is a comment. Copy the key lines you need into a deck, remove the
# leading "# " and set the values. A key line reads "# namespace.key = default",
# with a note after it where the default has a unit or a condition;
# "# required" means the key has no default, "# none" that it is off or empty
# unless given. The indented text below each key is its range and meaning. A key
# that has no effect with the other inputs is refused, as are unknown and
# repeated keys (docs/inputs.rst). Units are SI (m, s, K) unless the text says
# otherwise.
"""

GREEK = {"alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "Delta": "Δ", "epsilon": "ε", "varepsilon": "ε",
         "theta": "θ", "kappa": "κ", "lambda": "λ", "mu": "μ", "nu": "ν", "pi": "π", "rho": "ρ", "sigma": "σ",
         "tau": "τ", "phi": "φ", "varphi": "φ", "psi": "ψ", "omega": "ω", "Omega": "Ω", "zeta": "ζ", "ell": "ℓ",
         "ne": "≠", "le": "≤", "ge": "≥", "ldots": "…", "times": "×", "cdot": "·", "approx": "≈", "infty": "∞",
         "pm": "±", "to": "→", "lt": "<", "gt": ">"}


class DocError(Exception):
    """docs/inputs.rst does not parse as the reference expects."""


def labels(root):
    """{label: doc name} from the ``.. _label:`` lines of docs/*.rst."""
    out = {}
    for f in sorted((root / "docs").glob("*.rst")):
        for m in re.finditer(r"^\.\. _([^:]+):\s*$", read(f), re.M):
            out[m.group(1)] = f.stem
    return out


def math(s):
    s = re.sub(r"\\(?:mathrm|text|mathit|mathbf)\{([^}]*)\}", r"\1", s)
    s = re.sub(r"\\sqrt\{([^}]*)\}", r"√(\1)", s)
    s = re.sub(r"\^\{([^}]*)\}", lambda g: "^" + (g.group(1) if len(g.group(1)) == 1 else "(" + g.group(1) + ")"), s)

    def macro(g):
        name = g.group(1)
        sym = GREEK.get(name, name)
        before = s[g.start() - 1] if g.start() > 0 else ""
        after = s[g.end()] if g.end() < len(s) else ""
        # A space between a symbol and a letter or another macro (4π sin φ, a_m ζ).
        if before.isalpha() or (before.isdigit() and name not in GREEK):
            sym = " " + sym
        if after.isalpha() or after == "\\":
            sym = sym + " "
        return sym

    s = re.sub(r"\\([A-Za-z]+)", macro, s)
    s = re.sub(r"_\{([^}]*)\}", lambda g: "_" + g.group(1).replace(" ", ""), s)
    return s.replace("{", "").replace("}", "").replace("\\,", " ").replace("\\ ", " ")


def plain(text, refs):
    """docs RST inline markup to plain text."""
    text = re.sub(r":doc:`([^`<]*?)\s*<([^>]+)>`", r"\1 (docs/\2.rst)", text)
    text = re.sub(r":doc:`([^`]+)`", r"docs/\1.rst", text)

    def ref(m):
        title, label = (m.group(1), m.group(2)) if m.group(2) else (m.group(3), m.group(3))
        doc = refs.get(label)
        return f"{title.strip()}, docs/{doc}.rst" if doc else title.strip()

    text = re.sub(r":ref:`(?:([^`<]*?)\s*<([^>]+)>|([^`<]+))`", ref, text)
    text = re.sub(r":math:`([^`]+)`", lambda m: math(m.group(1)), text)
    text = re.sub(r":[a-z]+:`([^`]+)`", r"\1", text)
    text = text.replace("``", "").replace("\\|", "|").replace("**", "")
    text = re.sub(r"(?<![\w*])\*([^*\s][^*]*?)\*(?![\w*])", r"\1", text)
    text = re.sub(r"`([^`]+)`_", r"\1", text)
    return re.sub(r"\s+", " ", text).strip()


def expand(key):
    m = re.search(r"<([a-z]+)>", key)
    if not m:
        return [key]
    if m.group(1) not in PLACEHOLDERS:
        raise DocError(f"{key}: unknown placeholder <{m.group(1)}>")
    return [key[:m.start()] + v + key[m.end():] for v in PLACEHOLDERS[m.group(1)]]


def list_table(body, i, where):
    """The rows (lists of cell texts) of the list-table whose directive is body[i - 1]; and the next index."""
    header_rows = 1
    rows, row = [], None
    while i < len(body) and (body[i].startswith("   ") or body[i] == ""):
        line = body[i]
        m = re.match(r"   (\*| ) -(?: (.*))?$", line)
        opt = re.match(r"   :([a-z-]+):\s*(.*)$", line)
        if opt:
            if opt.group(1) == "header-rows":
                header_rows = int(opt.group(2))
        elif m:
            if m.group(1) == "*":
                if row is not None:
                    rows.append(row)
                row = [m.group(2) or ""]
            elif row is None:
                raise DocError(f"{where}: a cell before the first row: {line.strip()!r}")
            else:
                row.append(m.group(2) or "")
        elif line.startswith("       ") and row is not None:
            row[-1] += " " + line.strip()
        elif line != "":
            raise DocError(f"{where}: not a row, cell or continuation in a list-table: {line!r}")
        i += 1
    if row is not None:
        rows.append(row)
    return rows[header_rows:], i


def sections(doc, refs):
    """[(title, [notes], [(keys, default, meaning)])] for the detailed tables."""
    lines = doc.split("\n")
    heads = [i for i in range(1, len(lines)) if re.fullmatch(r"-{3,}", lines[i]) and lines[i - 1].startswith("``")]
    out, seen = [], set()
    for n, h in enumerate(heads):
        end = heads[n + 1] - 1 if n + 1 < len(heads) else len(lines)
        body = lines[h + 1:end]
        title = plain(lines[h - 1], refs)
        notes, rows, para, i = [], [], [], 0
        while i < len(body):
            line = body[i]
            if line.startswith(".. list-table::"):
                if para:
                    notes.append(plain(" ".join(para), refs))
                    para = []
                cells, i = list_table(body, i + 1, f"section {title}")
                for c in cells:
                    if len(c) != 3:
                        raise DocError(f"section {title}: a row with {len(c)} cells: {c[0]!r}")
                    names = re.findall(rf"``({KEY})``", c[0])
                    if not names:
                        raise DocError(f"section {title}: a row that names no key: {c[0]!r}")
                    keys = [k for name in names for k in expand(name)]
                    for k in keys:
                        if k in seen:
                            raise DocError(f"{k}: in two rows of the detailed tables")
                        seen.add(k)
                    rows.append((names, keys, c[1], c[2]))
                continue
            if line.strip() == "":
                if para:
                    notes.append(plain(" ".join(para), refs))
                    para = []
            elif not line.startswith(".. "):
                para.append(line)
            i += 1
        if para:
            notes.append(plain(" ".join(para), refs))
        out.append((title, notes, rows))
    return out


def quick_reference(doc):
    """[(keys, default)] of the quick-reference table (the first list-table)."""
    lines = doc.split("\n")
    start = next(i for i, line in enumerate(lines) if line.startswith(".. list-table::"))
    cells, _ = list_table(lines, start + 1, "quick reference")
    out = []
    for c in cells:
        names = re.findall(rf"``({KEY})``", c[0])
        if names:
            out.append(([k for name in names for k in expand(name)], c[1] if len(c) > 1 else ""))
    return out


NUMBER = re.compile(r"(-?(?:\d+\.?\d*|\.\d+)(?:e[-+]?\d+)?)(?:\s+(.*))?$", re.I)


def split_defaults(n, default, refs):
    """One (value, note) per key of a row of n keys. ``a``, ``b`` gives a and b;
    a number with a unit or a remark (0.4 m/s, 0 (off)) gives the number and the
    rest as a note; anything else (required, none, ...) gives no value."""
    vals = re.findall(r"``([^`]*)``", default)
    if n > 1 and len(vals) == n:
        return [(v, "") for v in vals]
    if len(vals) == 1 and re.fullmatch(r"\s*``[^`]*``\s*", default):
        return [(vals[0], "")] * n
    text = plain(default, refs)
    parts = [p.strip() for p in text.split(",")]
    if len(parts) != n:
        parts = [text] * n
    out = []
    for p in parts:
        m = NUMBER.match(p)
        out.append((m.group(1), (m.group(2) or "").strip()) if m else (None, p))
    return out


def render(doc, refs):
    w = textwrap.TextWrapper(width=WIDTH, initial_indent="#     ", subsequent_indent="#     ",
                             break_long_words=False, break_on_hyphens=False)
    wn = textwrap.TextWrapper(width=WIDTH, initial_indent="# ", subsequent_indent="# ",
                              break_long_words=False, break_on_hyphens=False)
    out = [HEADER]
    for title, notes, rows in sections(doc, refs):
        out += ["", "# " + "=" * (WIDTH - 2), "# " + title, "# " + "=" * (WIDTH - 2)]
        for note in notes:
            out += wn.wrap(note) + ["#"]
        if out[-1] == "#":
            out.pop()
        for names, keys, default, meaning in rows:
            out.append("#")
            # A placeholder row gives every key the row's one default.
            vals = split_defaults(len(names), default, refs)
            per_key = [v for name, v in zip(names, vals) for _ in expand(name)]
            for k, (v, note) in zip(keys, per_key):
                if note.startswith("(") and note.endswith(")") and note.count("(") == 1:
                    note = note[1:-1]
                if v is None:
                    out.append(f"# {k} =    # {note}")
                elif note:
                    out.append(f"# {k} = {v}    # {note}")
                else:
                    out.append(f"# {k} = {v}")
            out += w.wrap(plain(meaning, refs))
    return "\n".join(out) + "\n"


def read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def namespaces(root):
    """The owned namespaces (Source/Inputs.H) and every other namespace a ParmParse in Source/ reads
    (refavg, the reference-averaging tool's)."""
    m = re.search(r"owned_namespaces \(\)\s*\{.*?\{(.*?)\};", read(root / "Source" / "Inputs.H"), re.S)
    if not m:
        raise DocError("Source/Inputs.H: owned_namespaces() not found")
    out = re.findall(r'"([a-z_]+)"', m.group(1))
    for f in sorted(p for p in (root / "Source").rglob("*") if p.suffix in (".cpp", ".H")):
        for ns in re.findall(r'ParmParse\s*(?:\w+\s*)?[({]\s*"([a-z_]+)"\s*[)}]', read(f)):
            if ns not in out:
                out.append(ns)
    return out


READER = r"(?:query\w*|get\w*|contains|count\w*)"


def source_keys(root, spaces):
    """The keys the readers in Source/ name, as namespace.key."""
    keys = set()
    lit = re.compile(r'"((?:%s)\.[A-Za-z0-9_]*[A-Za-z0-9])"' % "|".join(spaces))
    files = sorted(p for p in (root / "Source").rglob("*") if p.suffix in (".cpp", ".H"))
    for f in files:
        s = read(f)
        keys.update(lit.findall(s))
        defs = [(m.start(), m.group(1), m.group(2))
                for m in re.finditer(r'ParmParse\s+(\w+)\s*[({]\s*"([a-z_]+)"\s*[)}]', s)]
        for m in re.finditer(r'ParmParse\s*\(\s*"([a-z_]+)"\s*\)\s*\.\s*%s\s*\(\s*"(\w+)"' % READER, s):
            keys.add(m.group(1) + "." + m.group(2))
        for m in re.finditer(r'\b(\w+)\s*\.\s*%s\s*\(\s*"(\w+)"' % READER, s):
            ns = [d[2] for d in defs if d[0] < m.start() and d[1] == m.group(1)]
            if ns:
                keys.add(ns[-1] + "." + m.group(2))
        file_ns = {d[2] for d in defs}
        # Key arrays (named k...Keys or keys, the readers' convention) and loop
        # lists. A loop list takes the namespace its body names ("ns." + k),
        # else the nearest ParmParse before it; an array the nearest ParmParse
        # before it, else the file's only namespace.
        lists = [(m, False) for m in re.finditer(
            r'const\s+char\s*\*\s*(?:const\s+)?(?:k\w*Keys|keys)\s*\[\s*\w*\s*\]\s*=\s*\{([^}]*)\}', s)]
        lists += [(m, True) for m in re.finditer(r'for\s*\(\s*const\s+char\s*\*\s*\w+\s*:\s*\{([^}]*)\}', s)]
        for m, loop in lists:
            named = re.search(r'"(%s)\."' % "|".join(spaces), s[m.end():m.end() + 400]) if loop else None
            before = [d[2] for d in defs if d[0] < m.start()]
            if named:
                ns = named.group(1)
            elif before:
                ns = before[-1]
            else:
                ns = next(iter(file_ns)) if len(file_ns) == 1 else None
            if ns:
                keys.update(ns + "." + word for word in re.findall(r'"(\w+)"', m.group(1)))
    return keys


def check(root):
    """The errors of --check, and the number of keys."""
    refs = labels(root)
    doc = read(root / "docs" / "inputs.rst")
    try:
        secs = sections(doc, refs)
        quick = quick_reference(doc)
        text = render(doc, refs)
    except DocError as e:
        return [f"docs/inputs.rst: {e}"], 0
    errors = []
    out = root / "Tests" / "inputs_reference"
    if not out.exists() or read(out) != text:
        errors.append("Tests/inputs_reference differs from docs/inputs.rst: run tools/input_reference.py --write")
    detail = {}
    for _, _, rows in secs:
        for names, keys, default, _ in rows:
            vals = split_defaults(len(names), default, refs)
            for k, v in zip(keys, [v for name, v in zip(names, vals) for _ in expand(name)]):
                detail[k] = v
    for keys, default in quick:
        vals = re.findall(r"``([^`]*)``", default)
        for n, k in enumerate(keys):
            if k not in detail:
                errors.append(f"{k}: in the quick-reference table but in no detailed table")
            elif len(vals) == len(keys) and detail[k][0] is not None and detail[k][1] == "" and vals[n] != detail[k][0]:
                errors.append(f"{k}: default {vals[n]!r} in the quick-reference table, {detail[k][0]!r} in its section")
    spaces = namespaces(root)
    titles = {t.split()[0] for t, _, _ in secs}
    for ns in spaces:
        if ns not in titles:
            errors.append(f"namespace {ns}: no section in docs/inputs.rst")
    read_keys = source_keys(root, spaces)
    for k in sorted(read_keys - set(detail)):
        errors.append(f"{k}: read in Source/ but not in docs/inputs.rst")
    for k in sorted(set(detail) - read_keys):
        errors.append(f"{k}: in docs/inputs.rst but not read in Source/")
    return errors, len(detail)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--write", action="store_true")
    g.add_argument("--check", action="store_true")
    ap.add_argument("--root", type=pathlib.Path, default=pathlib.Path(__file__).resolve().parent.parent,
                    help="the repository (default: this script's)")
    a = ap.parse_args()
    root = a.root
    if a.write:
        text = render(read(root / "docs" / "inputs.rst"), labels(root))
        with open(root / "Tests" / "inputs_reference", "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
        print("wrote Tests/inputs_reference")
        return 0
    errors, n = check(root)
    for e in errors:
        print("input_reference: " + e)
    print(f"input_reference: {n} keys, {'FAIL' if errors else 'OK'}")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())

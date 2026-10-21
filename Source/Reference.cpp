#include "Reference.H"

#include "FileMagic.H"
#include "Reductions.H"
#include "Sha256.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <cctype>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <math.h>   // fma: host and CUDA device overloads for float and double
#include <set>
#include <sstream>
#include <system_error>

namespace fs = std::filesystem;

namespace gust {

namespace {

constexpr const char* kMagic = "gust-reference";
constexpr int kVersion = 1;

// Every file-system call below runs on the I/O rank only, uses the
// std::error_code overloads (no exceptions), and records failures in a list
// that is broadcast before anyone acts on it, so all ranks take the same
// branch and no rank is left waiting in a collective.

std::string hexfloat (double x)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%a", x);
    return buf;
}

bool parse_double (const std::string& tok, double& x)
{
    if (tok.empty()) { return false; }
    char* end = nullptr;
    x = std::strtod(tok.c_str(), &end);
    return *end == '\0';
}

bool parse_int (const std::string& tok, int& x)
{
    if (tok.empty()) { return false; }
    char* end = nullptr;
    errno = 0;
    const long v = std::strtol(tok.c_str(), &end, 10);
    if (*end != '\0' || errno != 0 || v < INT_MIN || v > INT_MAX) { return false; }
    x = static_cast<int>(v);
    return true;
}

bool is_sha (const std::string& s)
{
    if (s.size() != 40) { return false; }
    for (const char c : s) {
        if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) { return false; }
    }
    return true;
}

bool is_word (const std::string& s)
{
    if (s.empty()) { return false; }
    for (const char c : s) {
        if (std::isspace(static_cast<unsigned char>(c))) { return false; }
    }
    return true;
}

[[noreturn]] void fail (const std::string& msg)
{
    amrex::Abort("gust: " + msg);
    std::abort();   // not reached; amrex::Abort does not return
}

[[noreturn]] void fail_all (const std::string& head, const std::vector<std::string>& err)
{
    std::string all = head;
    for (const auto& e : err) { all += "\n  " + e; }
    fail(all);
}

// The last path component, ignoring trailing separators ("a/b/" -> "b").
std::string base_name (const std::string& path)
{
    fs::path p(path);
    if (p.filename().empty()) { p = p.parent_path(); }
    return p.filename().string();
}

// Every regular file under dir, relative to dir with '/' separators. I/O rank only.
std::set<std::string> files_under (const std::string& dir, std::vector<std::string>& err)
{
    std::set<std::string> out;
    std::error_code ec;
    fs::recursive_directory_iterator it(dir, ec), end;
    for (; !ec && it != end; it.increment(ec)) {
        if (it->is_regular_file(ec)) { out.insert(fs::relative(it->path(), dir, ec).generic_string()); }
    }
    if (ec) { err.emplace_back(dir + ": cannot be listed (" + ec.message() + ")"); }
    return out;
}

// The I/O rank's lines, on every rank.
std::vector<std::string> broadcast_lines (const std::vector<std::string>& lines)
{
    std::string joined;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        for (const auto& l : lines) { joined += l + "\n"; }
    }
    long n = static_cast<long>(joined.size());
    const int root = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&n, 1, root);
    joined.resize(static_cast<std::size_t>(n));
    if (n > 0) { amrex::ParallelDescriptor::Bcast(joined.data(), static_cast<std::size_t>(n), root); }
    std::vector<std::string> out;
    std::istringstream is(joined);
    for (std::string l; std::getline(is, l);) { out.push_back(l); }
    return out;
}

struct Manifest
{
    ReferenceMeta meta;
    std::map<std::string, std::string> files;   // relative path -> sha256
};

// MANIFEST lines as read on the calling rank (I/O rank use only).
std::vector<std::string> read_lines (const std::string& path)
{
    std::vector<std::string> lines;
    std::ifstream is(path);
    for (std::string l; std::getline(is, l);) { lines.push_back(l); }
    return lines;
}

// Parse MANIFEST text (one entry per line); problems are appended to err.
Manifest parse_manifest (const std::vector<std::string>& lines, std::vector<std::string>& err)
{
    Manifest m;
    if (lines.empty()) {
        err.emplace_back("MANIFEST: missing or empty");
        return m;
    }
    {
        std::istringstream hs(lines[0]);
        std::string magic, ver;
        hs >> magic >> ver;
        int version = 0;
        if (!is_file_magic(magic, kMagic) || !parse_int(ver, version) || version != kVersion) {
            err.emplace_back("MANIFEST: not a GUST reference manifest of version " + std::to_string(kVersion));
            return m;
        }
    }
    int n_snap = -1;
    std::set<std::string> seen;
    for (std::size_t li = 1; li < lines.size(); ++li) {
        const std::string& line = lines[li];
        if (line.empty()) { continue; }
        std::istringstream ls(line);
        std::string key;
        ls >> key;
        std::vector<std::string> v;
        for (std::string t; ls >> t;) { v.push_back(t); }
        const std::string where = "MANIFEST line " + std::to_string(li + 1) + ": ";
        auto want = [&] (std::size_t n) {
            if (v.size() != n) { err.emplace_back(where + "'" + key + "' needs " + std::to_string(n) + " value(s)"); }
            return v.size() == n;
        };
        auto num = [&] (const std::string& tok, double& x) {
            if (!parse_double(tok, x)) { err.emplace_back(where + "bad number '" + tok + "'"); }
        };
        auto inum = [&] (const std::string& tok, int& x) {
            if (!parse_int(tok, x)) { err.emplace_back(where + "bad integer '" + tok + "'"); }
        };
        const bool repeatable = key == "file" || key == "snapshot" || key == "field";
        if (!repeatable && !seen.insert(key).second) {
            err.emplace_back(where + "'" + key + "' is given twice");
            continue;
        }
        if (key == "case") {
            if (want(1)) { m.meta.case_name = v[0]; }
        } else if (key == "erf_sha") {
            if (want(1)) { m.meta.erf_sha = v[0]; }
        } else if (key == "machine") {
            if (want(1) && !parse_machine(v[0], m.meta.machine)) { err.emplace_back(where + "unknown machine '" + v[0] + "'"); }
        } else if (key == "nodes") {
            if (want(1)) { inum(v[0], m.meta.nodes); }
        } else if (key == "ranks") {
            if (want(1)) { inum(v[0], m.meta.ranks); }
        } else if (key == "wall_seconds") {
            if (want(1)) { num(v[0], m.meta.wall_seconds); }
        } else if (key == "au") {
            if (want(1)) { num(v[0], m.meta.au); }
        } else if (key == "snapshots") {
            if (want(1)) { inum(v[0], n_snap); }
        } else if (key == "snapshot") {
            if (want(3)) {
                Snapshot s;
                num(v[0], s.time);
                inum(v[1], s.step);
                s.plotfile = v[2];
                m.meta.snapshots.push_back(s);
            }
        } else if (key == "field") {
            if (want(6)) {
                FieldStat f;
                f.name = v[0];
                num(v[1], f.sum); num(v[2], f.min); num(v[3], f.max);
                num(v[4], f.half_change); num(v[5], f.max_abs);
                m.meta.fields.push_back(f);
            }
        } else if (key == "file") {
            if (want(2) && !m.files.emplace(v[0], v[1]).second) {
                err.emplace_back(where + "file '" + v[0] + "' is listed twice");
            }
        } else {
            err.emplace_back(where + "unknown key '" + key + "'");
        }
    }
    for (const char* k : {"case", "erf_sha", "machine", "nodes", "ranks", "wall_seconds", "au", "snapshots"}) {
        if (seen.count(k) == 0) { err.emplace_back(std::string("MANIFEST: no '") + k + "' line"); }
    }
    if (n_snap != static_cast<int>(m.meta.snapshots.size())) {
        err.emplace_back("MANIFEST: 'snapshots' says " + std::to_string(n_snap) + " but "
                         + std::to_string(m.meta.snapshots.size()) + " snapshot lines follow");
    }
    if (m.meta.fields.empty()) { err.emplace_back("MANIFEST: no 'field' lines"); }
    return m;
}

// Statistics of component n of mf (collective, deterministic for a fixed BoxArray).
FieldStat field_stat (const std::string& name, const amrex::MultiFab& mf, int n)
{
    FieldStat f;
    f.name = name;
    f.sum = static_cast<double>(deterministic_sum(mf, n));
    f.min = static_cast<double>(mf.min(n, 0));
    f.max = static_cast<double>(mf.max(n, 0));
    return f;
}

} // namespace

double au_per_node_hour (Machine m)
{
    switch (m) {
    case Machine::KestrelCpu: return 10.0;
    case Machine::KestrelGpu: return 100.0;
    case Machine::Local: return 0.0;
    }
    return 0.0;
}

std::string machine_name (Machine m)
{
    switch (m) {
    case Machine::KestrelCpu: return "kestrel-cpu";
    case Machine::KestrelGpu: return "kestrel-gpu";
    case Machine::Local: return "local";
    }
    return "local";
}

bool parse_machine (const std::string& s, Machine& m)
{
    for (const auto c : {Machine::KestrelCpu, Machine::KestrelGpu, Machine::Local}) {
        if (s == machine_name(c)) { m = c; return true; }
    }
    return false;
}

double allocation_units (Machine m, int nodes, double wall_seconds)
{
    return static_cast<double>(nodes) * (wall_seconds / 3600.0) * au_per_node_hour(m);
}

std::vector<std::string> validate (const ReferenceMeta& m)
{
    std::vector<std::string> err;
    if (!is_word(m.case_name)) { err.emplace_back("refavg.case must be one word, not empty"); }
    if (!is_sha(m.erf_sha)) {
        err.emplace_back("refavg.erf_sha must be a full 40-character lowercase hex commit");
    }
    if (m.nodes < 1) { err.emplace_back("refavg.nodes must be >= 1"); }
    if (m.ranks < 1) { err.emplace_back("refavg.ranks must be >= 1"); }
    if (!std::isfinite(m.wall_seconds) || !(m.wall_seconds >= 0)) {
        err.emplace_back("refavg.wall_seconds must be finite and >= 0");
    }
    return err;
}

std::vector<double> trapezoid_weights (const std::vector<double>& t)
{
    const std::size_t n = t.size();
    if (n < 2) { fail("a time average needs at least 2 snapshots"); }
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(t[i])) { fail("snapshot " + std::to_string(i) + " has a non-finite time"); }
        if (i > 0 && !(t[i] > t[i - 1])) {
            fail("snapshot times must increase strictly (snapshot " + std::to_string(i) + ")");
        }
    }
    const double span = t[n - 1] - t[0];
    std::vector<double> w(n);
    w[0] = 0.5 * (t[1] - t[0]) / span;
    w[n - 1] = 0.5 * (t[n - 1] - t[n - 2]) / span;
    for (std::size_t i = 1; i + 1 < n; ++i) { w[i] = 0.5 * (t[i + 1] - t[i - 1]) / span; }
    return w;
}

Averaged average_plotfiles (const std::vector<std::string>& plotfiles)
{
    const std::size_t ns = plotfiles.size();
    if (ns < 3) {
        fail("refavg.plotfiles: a time average needs at least 3 plotfiles (two half windows "
             "for the stationarity check)");
    }

    Averaged avg;
    std::vector<double> times;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;

    for (std::size_t s = 0; s < ns; ++s) {
        const std::string& name = plotfiles[s];
        std::vector<std::string> missing;
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::error_code ec;
            if (!fs::exists(fs::path(name) / "Header", ec)) {
                missing.emplace_back("refavg.plotfiles: " + name + " is not a plotfile");
            }
        }
        missing = broadcast_lines(missing);
        if (!missing.empty()) { fail(missing.front()); }

        const std::string base = base_name(name);
        if (!is_word(base)) { fail("refavg.plotfiles: '" + name + "' has no usable directory name"); }

        amrex::PlotFileData pf(name);   // Header read once and broadcast by AMReX
        if (pf.spaceDim() != 3) { fail("refavg.plotfiles: " + name + " is not 3-D"); }
        if (pf.finestLevel() != 0) {
            fail("refavg.plotfiles: " + name + " has " + std::to_string(pf.finestLevel() + 1)
                 + " levels; references use single-level runs");
        }
        if (s == 0) {
            avg.names = pf.varNames();
            for (const auto& v : avg.names) {
                if (!is_word(v)) { fail("refavg.plotfiles: variable '" + v + "' contains white space"); }
            }
            avg.domain = pf.probDomain(0);
            avg.prob_lo = pf.probLo();
            avg.prob_hi = pf.probHi();
            ba = pf.boxArray(0);
            dm = pf.DistributionMap(0);
        } else {
            if (pf.varNames() != avg.names) { fail("refavg.plotfiles: " + name + " has other variables than " + plotfiles[0]); }
            if (pf.probDomain(0) != avg.domain || pf.probLo() != avg.prob_lo || pf.probHi() != avg.prob_hi
                || pf.boxArray(0) != ba) {
                fail("refavg.plotfiles: " + name + " has another domain or box layout than " + plotfiles[0]);
            }
        }
        times.push_back(static_cast<double>(pf.time()));
        avg.snapshots.push_back({base, static_cast<double>(pf.time()), pf.levelStep(0)});
    }

    // Whole window, first half [t_0, t_m] and second half [t_m, t_n].
    const std::size_t m = (ns - 1) / 2;
    const std::vector<double> w = trapezoid_weights(times);
    const std::vector<double> w1 = trapezoid_weights(
        std::vector<double>(times.begin(), times.begin() + static_cast<std::ptrdiff_t>(m) + 1));
    const std::vector<double> w2 = trapezoid_weights(
        std::vector<double>(times.begin() + static_cast<std::ptrdiff_t>(m), times.end()));

    const int nv = static_cast<int>(avg.names.size());
    avg.fields.define(ba, dm, nv, 0);
    amrex::MultiFab h1(ba, dm, nv, 0);
    amrex::MultiFab h2(ba, dm, nv, 0);
    avg.fields.setVal(0.0);
    h1.setVal(0.0);
    h2.setVal(0.0);

    // Accumulate in snapshot order with an explicit fma, so the sums are the
    // same on every compiler and rank count.
    for (std::size_t s = 0; s < ns; ++s) {
        amrex::PlotFileData pf(plotfiles[s]);
        for (int n = 0; n < nv; ++n) {
            amrex::MultiFab f = pf.get(0, avg.names[n]);
            if (f.contains_nan(0, 1, 0)) {
                fail("refavg.plotfiles: " + plotfiles[s] + " has NaN in " + avg.names[n]);
            }
            if (f.DistributionMap() != dm) {
                amrex::MultiFab g(ba, dm, 1, 0);
                g.ParallelCopy(f, 0, 0, 1);
                f = std::move(g);
            }
            const auto wa = static_cast<amrex::Real>(w[s]);
            const auto wf = static_cast<amrex::Real>(s <= m ? w1[s] : 0.0);
            const auto ws = static_cast<amrex::Real>(s >= m ? w2[s - m] : 0.0);
            for (amrex::MFIter mfi(avg.fields, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const amrex::Box bx = mfi.tilebox();
                const auto a = avg.fields.array(mfi, n);
                const auto a1 = h1.array(mfi, n);
                const auto a2 = h2.array(mfi, n);
                const auto t = f.const_array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    a(i, j, k) = fma(wa, t(i, j, k), a(i, j, k));
                    a1(i, j, k) = fma(wf, t(i, j, k), a1(i, j, k));
                    a2(i, j, k) = fma(ws, t(i, j, k), a2(i, j, k));
                });
            }
        }
    }

    amrex::MultiFab::Subtract(h2, h1, 0, 0, nv, 0);
    for (int n = 0; n < nv; ++n) {
        FieldStat f = field_stat(avg.names[n], avg.fields, n);
        f.half_change = static_cast<double>(h2.norm0(n));
        f.max_abs = static_cast<double>(avg.fields.norm0(n));
        avg.stats.push_back(f);
    }
    return avg;
}

void write_reference (const std::string& dir, const Averaged& avg, ReferenceMeta meta,
                      const std::vector<std::string>& inputs_files)
{
    meta.snapshots = avg.snapshots;
    meta.fields = avg.stats;
    meta.au = allocation_units(meta.machine, meta.nodes, meta.wall_seconds);

    // Metadata checks (same on every rank), then file-system checks and the
    // directory creation on the I/O rank, broadcast.
    std::vector<std::string> err = validate(meta);
    std::vector<std::string> fs_err;
    if (err.empty() && amrex::ParallelDescriptor::IOProcessor()) {
        std::error_code ec;
        if (fs::exists(dir, ec) || ec) { fs_err.emplace_back("refavg.output = " + dir + " already exists"); }
        std::set<std::string> bases;
        for (const auto& f : inputs_files) {
            if (!fs::is_regular_file(f, ec)) {
                fs_err.emplace_back("refavg.inputs_files: " + f + " is not a file");
            } else if (!bases.insert(base_name(f)).second) {
                fs_err.emplace_back("refavg.inputs_files: two files are named " + base_name(f));
            }
        }
        if (fs_err.empty() && !fs::create_directories(fs::path(dir) / "inputs", ec)) {
            fs_err.emplace_back("refavg.output = " + dir + ": cannot be created (" + ec.message() + ")");
        }
    }
    for (const auto& e : broadcast_lines(fs_err)) { err.push_back(e); }
    if (!err.empty()) { fail_all("invalid reference:", err); }

    const amrex::RealBox rb(avg.prob_lo.data(), avg.prob_hi.data());
    const amrex::Geometry geom(avg.domain, rb, 0, {0, 0, 0});
    const Snapshot& last = avg.snapshots.back();
    amrex::WriteSingleLevelPlotfile(dir + "/average", avg.fields, avg.names, geom,
                                    static_cast<amrex::Real>(last.time), last.step);
    amrex::ParallelDescriptor::Barrier();

    std::vector<std::string> werr;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::error_code ec;
        for (const auto& f : inputs_files) {
            fs::copy_file(f, fs::path(dir) / "inputs" / base_name(f), ec);
            if (ec) { werr.emplace_back("refavg.inputs_files: " + f + " cannot be copied (" + ec.message() + ")"); }
        }
        std::ofstream os(dir + "/MANIFEST");
        os << kMagic << " " << kVersion << "\n";
        os << "case " << meta.case_name << "\n";
        os << "erf_sha " << meta.erf_sha << "\n";
        os << "machine " << machine_name(meta.machine) << "\n";
        os << "nodes " << meta.nodes << "\n";
        os << "ranks " << meta.ranks << "\n";
        os << "wall_seconds " << hexfloat(meta.wall_seconds) << "\n";
        os << "au " << hexfloat(meta.au) << "\n";
        os << "snapshots " << meta.snapshots.size() << "\n";
        for (const auto& s : meta.snapshots) {
            os << "snapshot " << hexfloat(s.time) << " " << s.step << " " << s.plotfile << "\n";
        }
        for (const auto& f : meta.fields) {
            os << "field " << f.name << " " << hexfloat(f.sum) << " " << hexfloat(f.min) << " "
               << hexfloat(f.max) << " " << hexfloat(f.half_change) << " " << hexfloat(f.max_abs) << "\n";
        }
        for (const auto& f : files_under(dir, werr)) {
            if (f == "MANIFEST") { continue; }
            const std::string sha = sha256_file(dir + "/" + f);
            if (sha.empty()) { werr.emplace_back(f + ": cannot be read to compute its SHA-256"); }
            os << "file " << f << " " << sha << "\n";
        }
        os.close();
        if (!os) { werr.emplace_back(dir + "/MANIFEST: write failed"); }
    }
    werr = broadcast_lines(werr);
    if (!werr.empty()) { fail_all("reference " + dir + " was not written completely:", werr); }
}

std::vector<std::string> verify_reference (const std::string& dir)
{
    std::vector<std::string> err;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::error_code ec;
        if (!fs::is_directory(dir, ec)) {
            err.emplace_back(dir + ": not a directory");
        } else {
            Manifest m = parse_manifest(read_lines(dir + "/MANIFEST"), err);
            if (err.empty()) {
                for (const auto& e : validate(m.meta)) { err.emplace_back("MANIFEST: " + e); }
                if (m.meta.au != allocation_units(m.meta.machine, m.meta.nodes, m.meta.wall_seconds)) {
                    err.emplace_back("MANIFEST: au does not equal nodes x wall hours x the "
                                     + machine_name(m.meta.machine) + " rate");
                }
                if (m.files.empty()) { err.emplace_back("MANIFEST: lists no files"); }
                const auto present = files_under(dir, err);
                for (const auto& [f, sha] : m.files) {
                    if (present.count(f) == 0) {
                        err.emplace_back(f + ": listed in MANIFEST but missing");
                    } else if (sha256_file(dir + "/" + f) != sha) {
                        err.emplace_back(f + ": SHA-256 differs from MANIFEST");
                    }
                }
                for (const auto& f : present) {
                    if (f != "MANIFEST" && m.files.count(f) == 0) {
                        err.emplace_back(f + ": not listed in MANIFEST");
                    }
                }
                if (present.count("average/Header") == 0) { err.emplace_back("average: no averaged plotfile"); }
            }
        }
    }
    return broadcast_lines(err);
}

Reference load_reference (const std::string& dir)
{
    const auto err = verify_reference(dir);
    if (!err.empty()) { fail_all("reference " + dir + " does not verify:", err); }

    // One reading of MANIFEST, on the I/O rank, broadcast: every rank parses
    // the same text.
    std::vector<std::string> lines;
    if (amrex::ParallelDescriptor::IOProcessor()) { lines = read_lines(dir + "/MANIFEST"); }
    lines = broadcast_lines(lines);
    std::vector<std::string> perr;
    Manifest m = parse_manifest(lines, perr);
    if (!perr.empty()) { fail_all("reference " + dir + ": MANIFEST changed after verification:", perr); }

    Reference r;
    r.meta = m.meta;
    amrex::PlotFileData pf(dir + "/average");
    r.names = pf.varNames();
    r.domain = pf.probDomain(0);
    r.prob_lo = pf.probLo();
    r.prob_hi = pf.probHi();
    const amrex::BoxArray& ba = pf.boxArray(0);
    const amrex::DistributionMapping& dm = pf.DistributionMap(0);
    const int nv = static_cast<int>(r.names.size());
    r.fields.define(ba, dm, nv, 0);
    for (int n = 0; n < nv; ++n) {
        r.fields.ParallelCopy(pf.get(0, r.names[n]), 0, n, 1);
    }

    // The fields as read must have the recorded sum, minimum and maximum,
    // exactly. This ties the data to the manifest whatever the number of
    // files the averaged plotfile was written in, and catches a change made
    // after the file checksums were checked.
    std::vector<std::string> ferr;
    if (static_cast<int>(r.meta.fields.size()) != nv) {
        ferr.emplace_back("MANIFEST lists " + std::to_string(r.meta.fields.size())
                          + " fields but average has " + std::to_string(nv));
    } else {
        for (int n = 0; n < nv; ++n) {
            const FieldStat& want = r.meta.fields[static_cast<std::size_t>(n)];
            const FieldStat got = field_stat(r.names[n], r.fields, n);
            if (want.name != got.name) {
                ferr.emplace_back("field " + std::to_string(n) + ": MANIFEST says " + want.name
                                  + ", average has " + got.name);
            } else if (want.sum != got.sum || want.min != got.min || want.max != got.max) {
                ferr.emplace_back("field " + got.name + ": sum, minimum or maximum differs from MANIFEST");
            }
        }
    }
    if (!ferr.empty()) { fail_all("reference " + dir + " does not verify:", ferr); }
    return r;
}

} // namespace gust

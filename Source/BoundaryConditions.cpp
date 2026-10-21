#include "BoundaryConditions.H"

#include "MassFlux_K.H"

#include "CVFluxes.H"
#include "Heights_K.H"
#include "ImplicitOperator.H"
#include "Stencil_K.H"


#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>

namespace gust {

namespace {

constexpr int kNumTypes = 7;
const char* const kNames[kNumTypes] = {"periodic", "inflow", "outflow", "slip_wall", "no_slip_wall",
                                       "surface_layer", "geostrophic"};
const char* const kErfNames[kNumTypes] = {"Periodic", "Inflow", "Outflow", "SlipWall", "NoSlipWall",
                                          "surface_layer", "SlipWall"};
const char* const kFaceKeys[6] = {"xlo", "ylo", "zlo", "xhi", "yhi", "zhi"};

//! How a ghost is filled along a face normal (see the header).
enum class Rule { even, odd, value, profile, gradient };

//! The fields the fills act on, in the order of State::fields.
enum Field : int { f_u = 0, f_v, f_w, f_p, f_theta, f_tke, n_fields };

//! The rule for field f at a face of type t normal to direction d.
Rule rule_for (Field f, FaceType t, int d, bool lo)
{
    const bool normal_vel = (f == f_u && d == 0) || (f == f_v && d == 1) || (f == f_w && d == 2);
    const bool vel = (f == f_u || f == f_v || f == f_w);
    switch (t) {
    case FaceType::inflow:
        return (f == f_p) ? Rule::even : Rule::profile;
    case FaceType::outflow:
        return (f == f_p) ? Rule::odd : Rule::even;
    case FaceType::slip_wall:
        if (normal_vel) { return Rule::value; }
        // theta's gradient bc.zhi_theta_grad is the top's; a slip-wall ground is even.
        return (f == f_theta && d == 2 && !lo) ? Rule::gradient : Rule::even;
    case FaceType::no_slip_wall:
        // Every velocity component takes the wall's value: 0, or the moving
        // wall's tangential velocity (bc.<face>_wall_velocity).
        return vel ? Rule::value : Rule::even;
    case FaceType::surface_layer:
        if (normal_vel) { return Rule::value; }
        return (vel || f == f_theta) ? Rule::gradient : Rule::even;
    case FaceType::geostrophic:
        if (normal_vel) { return Rule::value; }
        if (f == f_u || f == f_v) { return Rule::value; }
        return (f == f_theta) ? Rule::gradient : Rule::even;
    default:
        return Rule::even;
    }
}

//! The inflow profile's value for field f at height h (0 when homogeneous).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real profile_value (const InflowProfile::View& p, int f, amrex::Real h) noexcept
{
    amrex::Real u, v, t, k, nut;
    p.eval(h, u, v, t, k, nut);
    switch (f) {
    case f_u: return u;
    case f_v: return v;
    case f_theta: return t;
    case f_tke: return k;
    default: return 0;
    }
}

//! The ghost slab of layer L (1 or 2) outside face (d, lo) of a box, grown
//! over the other directions as fill_bc_ghosts: fully over the directions
//! already filled (e < d), and over later ones only where values exist
//! (inside the domain, or across periodic faces).
amrex::Box ghost_slab (const amrex::Box& vb, const amrex::Box& dom, const amrex::Geometry& geom, int d, bool lo,
                       int layer, int ng)
{
    const bool nodal_d = vb.type(d) == amrex::IndexType::NODE;
    const int edge = lo ? dom.smallEnd(d) : (nodal_d ? dom.bigEnd(d) + 1 : dom.bigEnd(d));
    amrex::Box gb = vb;
    for (int e = 0; e < AMREX_SPACEDIM; ++e) {
        if (e == d) { continue; }
        gb.grow(e, ng);
        if (e > d && !geom.isPeriodic(e)) {
            const bool nodal_e = vb.type(e) == amrex::IndexType::NODE;
            gb.setSmall(e, std::max(gb.smallEnd(e), dom.smallEnd(e)));
            gb.setBig(e, std::min(gb.bigEnd(e), nodal_e ? dom.bigEnd(e) + 1 : dom.bigEnd(e)));
        }
    }
    const int g = lo ? edge - layer : edge + layer;
    gb.setSmall(d, g);
    gb.setBig(d, g);
    return gb;
}

//! Whether box vb touches the domain face (d, lo).
bool touches (const amrex::Box& vb, const amrex::Box& dom, int d, bool lo)
{
    const bool nodal_d = vb.type(d) == amrex::IndexType::NODE;
    return lo ? vb.smallEnd(d) == dom.smallEnd(d) : vb.bigEnd(d) == (nodal_d ? dom.bigEnd(d) + 1 : dom.bigEnd(d));
}

} // namespace

const char* face_type_name (FaceType t) { return kNames[static_cast<int>(t)]; }
const char* erf_bc_name (FaceType t) { return kErfNames[static_cast<int>(t)]; }

FaceType face_type_from_name (const std::string& name, const std::string& key)
{
    for (int n = 0; n < kNumTypes; ++n) {
        if (name == kNames[n]) { return static_cast<FaceType>(n); }
    }
    std::string all;
    for (int n = 0; n < kNumTypes; ++n) { all += std::string(n ? ", " : "") + kNames[n]; }
    amrex::Abort("gust: " + key + " = " + name + ": unknown face type (one of " + all + ")");
    return FaceType::periodic;
}

const char* const kSurfaceKeys[] = {"surface_z0", "surface_zref", "surface_stability", "surface_heat_flux", "surface_wall_k",
                                    "surface_denominator_floor"};

void use_closure_constants (BoundaryInputs& bin, SurfaceLayerOptions& surf, amrex::Real Cmu0, amrex::Real kappa,
                            amrex::Real max_geom_lscale)
{
    surf.Cmu0 = Cmu0;
    surf.kappa = kappa;
    bin.inflow.Cmu0 = Cmu0;
    bin.inflow.kappa = kappa;
    bin.inflow.l_max = max_geom_lscale;
}

BoundaryInputs read_boundary_inputs ()
{
    BoundaryInputs in;
    amrex::ParmParse pp("bc");
    const char* keys[] = {"xlo", "xhi", "ylo", "yhi", "zlo", "zhi", "inflow_profile", "inflow_ustar", "inflow_z0",
                          "inflow_speed", "inflow_zref", "inflow_alpha", "inflow_direction", "inflow_file", "inflow_tke_length",
                          "theta_profile_z", "theta_profile", "zhi_theta_grad", "rayleigh_zdamp",
                          "rayleigh_dampcoef", "rayleigh_fields", "rayleigh_height", "sponge_width", "sponge_strength",
                          "xlo_wall_velocity", "xhi_wall_velocity", "ylo_wall_velocity", "yhi_wall_velocity",
                          "zlo_wall_velocity", "zhi_wall_velocity"};
    std::vector<std::string> errs;
    for (const char* key : keys) {
        if (!pp.contains(key)) { continue; }
        in.given.emplace_back(key);
        if (pp.countname(key) > 1) { errs.emplace_back(std::string("bc.") + key + " is given more than once"); }
        const std::string k(key);
        const bool array = k.rfind("theta_profile", 0) == 0;
        const bool three = k.size() > 14 && k.compare(k.size() - 14, 14, "_wall_velocity") == 0;
        if (k == "rayleigh_fields") {
            if (pp.countval(key) < 1) { errs.emplace_back("bc.rayleigh_fields takes at least one field name"); }
            continue;
        }
        if (three && pp.countval(key) != 3) { errs.emplace_back(std::string("bc.") + key + " takes exactly three values (u v w)"); }
        if (!array && !three && pp.countval(key) != 1) { errs.emplace_back(std::string("bc.") + key + " takes exactly one value"); }
    }
    for (const char* key : kSurfaceKeys) {
        if (pp.contains(key)) { in.given.emplace_back(key); }   // read by read_surface_layer_options
    }
    if (!errs.empty()) { amrex::Abort("gust: invalid inputs:\n  " + errs.front()); }
    for (int o = 0; o < 6; ++o) {
        std::string name;
        if (pp.query(kFaceKeys[o], name)) { in.face[o] = face_type_from_name(name, std::string("bc.") + kFaceKeys[o]); }
    }
    std::string kind = "log_law";
    pp.query("inflow_profile", kind);
    using K = InflowProfileOptions::Kind;
    if (kind == "log_law") { in.inflow.kind = K::log_law; }
    else if (kind == "power_law") { in.inflow.kind = K::power_law; }
    else if (kind == "file") { in.inflow.kind = K::file; }
    else { amrex::Abort("gust: bc.inflow_profile = " + kind + ": must be log_law, power_law or file"); }
    pp.query("inflow_ustar", in.inflow.ustar);
    pp.query("inflow_z0", in.inflow.z0);
    pp.query("inflow_speed", in.inflow.speed_ref);
    pp.query("inflow_zref", in.inflow.z_ref);
    pp.query("inflow_alpha", in.inflow.alpha);
    pp.query("inflow_direction", in.inflow.direction);
    pp.query("inflow_file", in.inflow.file);
    std::string length;
    if (pp.query("inflow_tke_length", length)) {
        if (length == "closure") { in.inflow.length = InflowProfileOptions::Length::closure; }
        else if (length == "unbounded") { in.inflow.length = InflowProfileOptions::Length::unbounded; }
        else { amrex::Abort("gust: bc.inflow_tke_length = " + length + ": must be closure or unbounded"); }
    }
    amrex::Vector<amrex::Real> tz, tt;
    pp.queryarr("theta_profile_z", tz);
    pp.queryarr("theta_profile", tt);
    if (!(tz.empty() && tt.empty())) {
        in.theta_ref.z = tz;
        in.theta_ref.theta = tt;
    }
    pp.query("zhi_theta_grad", in.zhi_theta_grad);
    for (int o = 0; o < 6; ++o) {
        amrex::Vector<amrex::Real> wv;
        if (pp.queryarr((std::string(kFaceKeys[o]) + "_wall_velocity").c_str(), wv) && wv.size() == 3) {
            for (int c = 0; c < 3; ++c) { in.wall_velocity[o][c] = wv[c]; }
        }
    }
    pp.query("rayleigh_zdamp", in.rayleigh_zdamp);
    pp.query("rayleigh_dampcoef", in.rayleigh_dampcoef);
    {
        std::string rh = "absolute";
        if (pp.query("rayleigh_height", rh)) {
            if (rh == "absolute") { in.rayleigh_above_ground = false; }
            else if (rh == "above_ground") { in.rayleigh_above_ground = true; }
            else { amrex::Abort("gust: bc.rayleigh_height = " + rh + ": must be absolute or above_ground"); }
        }
        std::vector<std::string> names;
        if (pp.queryarr("rayleigh_fields", names)) {
            in.rayleigh_fields = {false, false, false, false, false};
            for (const auto& n : names) {
                const int f = (n == "u") ? 0 : (n == "v") ? 1 : (n == "w") ? 2 : (n == "theta") ? 3 : (n == "tke") ? 4 : -1;
                if (f < 0) { amrex::Abort("gust: bc.rayleigh_fields = " + n + ": must be among u, v, w, theta, tke"); }
                in.rayleigh_fields[static_cast<std::size_t>(f)] = true;
            }
        }
    }
    pp.query("sponge_width", in.sponge_width);
    pp.query("sponge_strength", in.sponge_strength);
    return in;
}

std::vector<std::string> check_boundary_inputs (const BoundaryInputs& in, const amrex::Geometry& geom,
                                                const CoriolisOptions& coriolis)
{
    return check_boundary_inputs(in, geom, coriolis, geom.ProbHi(2));
}

std::vector<std::string> check_boundary_inputs (const BoundaryInputs& in, const amrex::Geometry& geom,
                                                const CoriolisOptions& coriolis, amrex::Real ztop)
{
    std::vector<std::string> err;
    bool any_inflow = false, any_outflow = false;
    for (int o = 0; o < 6; ++o) {
        const amrex::Orientation ori(o % 3, o < 3 ? amrex::Orientation::low : amrex::Orientation::high);
        const int d = ori.coordDir();
        const FaceType t = in.face[o];
        const std::string key = std::string("bc.") + kFaceKeys[o];
        if (geom.isPeriodic(d) != (t == FaceType::periodic)) {
            err.emplace_back(key + " = " + face_type_name(t) + ": periodic exactly when mesh.is_periodic is 1 in that direction");
            continue;
        }
        if (d < 2) {
            if (t == FaceType::surface_layer || t == FaceType::geostrophic) {
                err.emplace_back(key + " = " + face_type_name(t) + ": only the ground (zlo) or the top (zhi) take this type");
            }
        } else if (ori.isLow()) {
            if (t != FaceType::surface_layer && t != FaceType::slip_wall && t != FaceType::no_slip_wall) {
                err.emplace_back(key + " = " + face_type_name(t) + ": the ground must be surface_layer, slip_wall or no_slip_wall");
            }
        } else if (t != FaceType::slip_wall && t != FaceType::no_slip_wall && t != FaceType::geostrophic) {
            err.emplace_back(key + " = " + face_type_name(t) + ": the top must be slip_wall, no_slip_wall or geostrophic");
        }
        if (t == FaceType::geostrophic && !(coriolis.enabled && coriolis.geostrophic)) {
            err.emplace_back(key + " = geostrophic needs Coriolis with geostrophic forcing (the wind it holds)");
        }
        any_inflow = any_inflow || t == FaceType::inflow;
        any_outflow = any_outflow || t == FaceType::outflow;
        const std::string wkey = key + "_wall_velocity";
        const bool wgiven = std::find(in.given.begin(), in.given.end(), std::string(kFaceKeys[o]) + "_wall_velocity") != in.given.end();
        if (wgiven && t != FaceType::no_slip_wall) {
            err.emplace_back(wkey + " has no effect without " + key + " = no_slip_wall");
        } else if (wgiven) {
            const auto& wv = in.wall_velocity[o];
            if (!(std::isfinite(wv[0]) && std::isfinite(wv[1]) && std::isfinite(wv[2]))) {
                err.emplace_back(wkey + ": every component must be finite (m/s)");
            } else if (wv[d] != 0) {
                err.emplace_back(wkey + " = " + std::to_string(wv[0]) + " " + std::to_string(wv[1]) + " "
                                 + std::to_string(wv[2]) + ": the normal component must be 0 (the wall moves in its plane)");
            }
        }
    }
    auto given = [&] (const char* k) { return std::find(in.given.begin(), in.given.end(), k) != in.given.end(); };
    const bool needs_profile = any_inflow || in.sponge_width > 0 || in.rayleigh_zdamp > 0 || in.profile_needed;
    if (needs_profile) {
        for (const auto& e : check_inflow_profile_options(in.inflow)) { err.push_back(e); }
        // The keys of the other kinds would be read and ignored.
        using K = InflowProfileOptions::Kind;
        const char* kind_name = (in.inflow.kind == K::log_law) ? "log_law" : (in.inflow.kind == K::power_law) ? "power_law" : "file";
        auto only_for = [&] (const char* k, K kind) {
            if (given(k) && in.inflow.kind != kind) {
                err.emplace_back(std::string("bc.") + k + " has no effect with bc.inflow_profile = " + kind_name);
            }
        };
        only_for("inflow_ustar", K::log_law);
        only_for("inflow_speed", K::power_law);
        only_for("inflow_zref", K::power_law);
        only_for("inflow_alpha", K::power_law);
        only_for("inflow_file", K::file);
    } else {
        for (const char* k : {"inflow_profile", "inflow_ustar", "inflow_z0", "inflow_speed", "inflow_zref",
                              "inflow_alpha", "inflow_direction", "inflow_file", "inflow_tke_length"}) {
            if (given(k)) {
                err.emplace_back(std::string("bc.") + k
                                 + " has no effect without an inflow face, a damping layer or init.first_guess = profile or column");
            }
        }
    }
    for (const auto& e : check_reference_profile(in.theta_ref)) { err.push_back("bc.theta_profile: " + e); }
    if (in.face[2] != FaceType::surface_layer) {
        for (const char* k : kSurfaceKeys) {
            if (given(k)) { err.emplace_back(std::string("bc.") + k + " has no effect without bc.zlo = surface_layer"); }
        }
    }
    if (in.theta_ref.z.size() != in.theta_ref.theta.size()) {
        err.emplace_back("bc.theta_profile_z and bc.theta_profile need the same number of values");
    }
    if (!std::isfinite(in.zhi_theta_grad)) {
        err.emplace_back("bc.zhi_theta_grad = " + std::to_string(in.zhi_theta_grad) + ": must be finite (K/m)");
    } else if (given("zhi_theta_grad") && in.face[5] == FaceType::no_slip_wall) {
        err.emplace_back("bc.zhi_theta_grad has no effect with bc.zhi = no_slip_wall");
    }
    if (std::isnan(in.rayleigh_zdamp) || !(in.rayleigh_zdamp >= 0 && in.rayleigh_zdamp < ztop)) {
        char top[64];
        std::snprintf(top, sizeof(top), "%.1f", static_cast<double>(ztop));
        err.emplace_back("bc.rayleigh_zdamp = " + std::to_string(in.rayleigh_zdamp) + ": must be in [0, the top level "
                         + top + " m) (0: off)");
    }
    if (std::isnan(in.rayleigh_dampcoef) || !(in.rayleigh_dampcoef >= 0 && std::isfinite(in.rayleigh_dampcoef))) {
        err.emplace_back("bc.rayleigh_dampcoef = " + std::to_string(in.rayleigh_dampcoef) + ": must be finite and >= 0 (1/s)");
    } else if (given("rayleigh_dampcoef") && !(in.rayleigh_zdamp > 0)) {
        err.emplace_back("bc.rayleigh_dampcoef has no effect without bc.rayleigh_zdamp > 0");
    } else if (given("rayleigh_fields") && !(in.rayleigh_zdamp > 0)) {
        err.emplace_back("bc.rayleigh_fields has no effect without bc.rayleigh_zdamp > 0");
    } else if (given("rayleigh_height") && !(in.rayleigh_zdamp > 0)) {
        err.emplace_back("bc.rayleigh_height has no effect without bc.rayleigh_zdamp > 0");
    } else if (in.rayleigh_zdamp > 0 && !(in.rayleigh_dampcoef > 0)) {
        err.emplace_back("bc.rayleigh_zdamp > 0 needs bc.rayleigh_dampcoef > 0");
    }
    if (std::isnan(in.sponge_width) || !(in.sponge_width >= 0 && std::isfinite(in.sponge_width))) {
        err.emplace_back("bc.sponge_width = " + std::to_string(in.sponge_width) + ": must be finite and >= 0 (0: off)");
    } else if (in.sponge_width > 0 && !(any_inflow || any_outflow)) {
        err.emplace_back("bc.sponge_width > 0 needs an inflow or outflow face to damp towards");
    }
    if (std::isnan(in.sponge_strength) || !(in.sponge_strength >= 0 && std::isfinite(in.sponge_strength))) {
        err.emplace_back("bc.sponge_strength = " + std::to_string(in.sponge_strength) + ": must be finite and >= 0 (1/s)");
    } else if (given("sponge_strength") && !(in.sponge_width > 0)) {
        err.emplace_back("bc.sponge_strength has no effect without bc.sponge_width > 0");
    } else if (in.sponge_width > 0 && !(in.sponge_strength > 0)) {
        err.emplace_back("bc.sponge_width > 0 needs bc.sponge_strength > 0");
    }
    return err;
}

PoissonBC poisson_bc_from_faces (const BoundaryInputs& in, const amrex::Geometry& geom)
{
    PoissonBC bc;
    for (int o = 0; o < 6; ++o) {
        const amrex::Orientation ori(o % 3, o < 3 ? amrex::Orientation::low : amrex::Orientation::high);
        if (geom.isPeriodic(ori.coordDir())) {
            bc.face[o] = PoissonFace::periodic;
        } else if (in.face[o] == FaceType::outflow) {
            bc.face[o] = PoissonFace::dirichlet;
        } else {
            bc.face[o] = PoissonFace::neumann;
        }
    }
    return bc;
}

namespace {

//! The inputs, after the checks (every message, each naming its key); the
//! inflow profile is built only from checked inputs.
const BoundaryInputs& checked (const BoundaryInputs& in, const amrex::Geometry& geom, const Mesh& mesh,
                               const CoriolisOptions& coriolis)
{
    if (mesh.zlevels().empty()) { amrex::Abort("gust: BoundaryConditions: the mesh has no levels (not built)"); }
    const auto err = check_boundary_inputs(in, geom, coriolis, mesh.zlevels().back());
    if (!err.empty()) {
        std::string all = "gust: invalid boundary conditions:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    return in;
}

} // namespace

BoundaryConditions::BoundaryConditions (const BoundaryInputs& in, const amrex::Geometry& geom, const Mesh& mesh,
                                        const CoriolisOptions& coriolis)
    : m_in(checked(in, geom, mesh, coriolis)), m_geom(geom), m_mesh(mesh), m_coriolis(coriolis), m_profile(m_in.inflow, m_in.theta_ref)
{
    // Known only once the table is read: with its own k column the length
    // changes nothing (the fill, the damping and the first guess take k, not
    // the profile's nu_t), so the key would be read and ignored.
    if (m_profile.has_k() && std::find(m_in.given.begin(), m_in.given.end(), "inflow_tke_length") != m_in.given.end()) {
        amrex::Abort("gust: invalid boundary conditions:\n  bc.inflow_tke_length has no effect with a k column in "
                     "bc.inflow_file = " + m_in.inflow.file);
    }
}

bool BoundaryConditions::has_inflow () const
{
    for (const FaceType t : m_in.face) { if (t == FaceType::inflow) { return true; } }
    return false;
}

namespace {

//! Fill the ghosts of one field at one face for one rule (not gradient).
template <int T>
void fill_face (amrex::MultiFab& q, const amrex::MultiFab& z_nd, const amrex::Geometry& geom, int d, bool lo,
                int f, Rule rule, amrex::Real value, const InflowProfile::View& prof, bool homogeneous)
{
    const amrex::Box& dom = geom.Domain();
    const int klo = dom.smallEnd(2);
    const bool nodal_d = stencil::nodal(T, d);
    const int edge = lo ? dom.smallEnd(d) : (nodal_d ? dom.bigEnd(d) + 1 : dom.bigEnd(d));
    const amrex::Real V = homogeneous ? amrex::Real(0) : value;
    const int sgn = (rule == Rule::odd) ? -1 : 1;
    for (amrex::MFIter mfi(q, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        if (!touches(vb, dom, d, lo)) { continue; }
        const auto a = q.array(mfi);
        const amrex::Box zb = z_nd[mfi].box();
        const NodeZ nz{z_nd.const_array(mfi), amrex::lbound(zb), amrex::ubound(zb)};
        // The face point of a normal velocity takes the value itself.
        if (nodal_d && (rule == Rule::value || rule == Rule::profile)) {
            amrex::Box fb = ghost_slab(vb, dom, geom, d, lo, 0, q.nGrow());
            fb &= mfi.fabbox();
            const bool is_profile = (rule == Rule::profile);
            amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                a(i, j, k) = (is_profile && !homogeneous) ? profile_value(prof, f, height_above_ground<T>(nz, i, j, k, klo)) : V;
            });
        }
        for (int layer = 1; layer <= q.nGrow(); ++layer) {
            amrex::Box gb = ghost_slab(vb, dom, geom, d, lo, layer, q.nGrow());
            gb &= mfi.fabbox();
            if (gb.isEmpty()) { continue; }
            // The mirror point about the face (nodal) or the face plane (cell).
            const int mirror = nodal_d ? 2 * edge : (lo ? 2 * edge - 1 : 2 * edge + 1);
            const int ex = (d == 0) ? 1 : 0, ey = (d == 1) ? 1 : 0, ez = (d == 2) ? 1 : 0;
            if (rule == Rule::profile) {
                amrex::ParallelFor(gb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    a(i, j, k) = homogeneous ? amrex::Real(0) : profile_value(prof, f, height_above_ground<T>(nz, i, j, k, klo));
                });
            } else if (rule == Rule::value) {
                amrex::ParallelFor(gb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const int m = mirror - (d == 0 ? i : d == 1 ? j : k);
                    a(i, j, k) = 2 * V - a(i + ex * (m - i), j + ey * (m - j), k + ez * (m - k));
                });
            } else {
                amrex::ParallelFor(gb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const int m = mirror - (d == 0 ? i : d == 1 ? j : k);
                    a(i, j, k) = sgn * a(i + ex * (m - i), j + ey * (m - j), k + ez * (m - k));
                });
            }
        }
    }
}

//! The gradient fill along z (the ground or the top) for a field at T with
//! the gradient G(i, j) of a per-column slab (one component), or a constant.
template <int T>
void fill_gradient_z (amrex::MultiFab& q, const amrex::MultiFab& z_nd, const amrex::Geometry& geom, bool lo,
                      const amrex::MultiFab* grad, int comp, amrex::Real constant, bool homogeneous)
{
    const amrex::Box& dom = geom.Domain();
    const int klo = dom.smallEnd(2);   // the ground slab's k
    const int edge = lo ? klo : dom.bigEnd(2);
    const int c1 = lo ? edge + 1 : edge - 1;
    for (amrex::MFIter mfi(q, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        if (!touches(vb, dom, 2, lo)) { continue; }
        const auto a = q.array(mfi);
        const amrex::Box zb = z_nd[mfi].box();
        const NodeZ nz{z_nd.const_array(mfi), amrex::lbound(zb), amrex::ubound(zb)};
        amrex::Array4<amrex::Real const> g;
        if (grad != nullptr) { g = grad->const_array(mfi); }
        const bool per_column = (grad != nullptr);
        amrex::Box gb = ghost_slab(vb, dom, geom, 2, lo, 1, q.nGrow());
        gb &= mfi.fabbox();
        const int g1 = lo ? edge - 1 : edge + 1, g2 = lo ? edge - 2 : edge + 2;
        // The slab g holds the two ghost columns the slab gb spans; a face
        // point of a nodal direction averages the two cells beside it, and
        // at the outermost ghost faces (one beyond the columns) the one it has.
        const amrex::Dim3 glo = per_column ? amrex::lbound(g) : amrex::Dim3{0, 0, 0};
        const amrex::Dim3 ghi = per_column ? amrex::ubound(g) : amrex::Dim3{0, 0, 0};
        amrex::ParallelFor(gb, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            amrex::Real G = constant;
            if (per_column) {
                const int ia = amrex::max(glo.x, amrex::min(ghi.x, i)), ja = amrex::max(glo.y, amrex::min(ghi.y, j));
                const int ib = amrex::max(glo.x, amrex::min(ghi.x, stencil::nodal(T, 0) ? i - 1 : i));
                const int jb = amrex::max(glo.y, amrex::min(ghi.y, stencil::nodal(T, 1) ? j - 1 : j));
                G = amrex::Real(0.5) * (g(ib, jb, klo, comp) + g(ia, ja, klo, comp));
            }
            if (homogeneous) { G = amrex::Real(0); }
            const amrex::Real zc1 = point_z<T>(nz, i, j, c1), zg1 = point_z<T>(nz, i, j, g1), zg2 = point_z<T>(nz, i, j, g2);
            a(i, j, g1) = a(i, j, c1) + (zg1 - zc1) * G;
            a(i, j, g2) = a(i, j, g1) + (zg2 - zg1) * G;
        });
    }
}

template <int T>
void fill_field (amrex::MultiFab& q, int f, int d, const BoundaryInputs& in, const amrex::MultiFab& z_nd,
                 const amrex::Geometry& geom, const InflowProfile::View& prof, const CoriolisOptions& cor,
                 const amrex::MultiFab* ground_grad, bool homogeneous)
{
    for (const bool lo : {true, false}) {
            const FaceType t = in.face[amrex::Orientation(d, lo ? amrex::Orientation::low : amrex::Orientation::high)];
            const Rule rule = rule_for(static_cast<Field>(f), t, d, lo);
            amrex::Real value = 0;
            if (t == FaceType::geostrophic) { value = (f == f_u) ? cor.ug : (f == f_v) ? cor.vg : amrex::Real(0); }
            if (t == FaceType::no_slip_wall && f <= f_w) {
                value = in.wall_velocity[amrex::Orientation(d, lo ? amrex::Orientation::low : amrex::Orientation::high)][f];
            }
            if (rule == Rule::gradient) {
                if (t == FaceType::surface_layer) {
                    fill_gradient_z<T>(q, z_nd, geom, lo, ground_grad, f == f_u ? 0 : f == f_v ? 1 : 2, 0, homogeneous);
                } else {
                    fill_gradient_z<T>(q, z_nd, geom, lo, nullptr, 0, in.zhi_theta_grad, homogeneous);
                }
            } else {
                fill_face<T>(q, z_nd, geom, d, lo, f, rule, value, prof, homogeneous);
            }
    }
}

} // namespace

namespace {
void fill_no_flux_w (State& q, const BoundaryInputs& in, const amrex::Geometry& geom);
} // namespace

void BoundaryConditions::fill_ghosts (State& q, const SurfaceLayerOptions* surface_opt, const BuoyancyOptions& buoyancy,
                                      bool homogeneous, SurfaceState* surface_out) const
{
    const bool ground_wall = (m_in.face[2] == FaceType::surface_layer);
    if (ground_wall && !homogeneous && surface_opt == nullptr) {
        amrex::Abort("gust: fill_ghosts: bc.zlo = surface_layer needs the surface-layer options");
    }
    amrex::MultiFab* fields[n_fields] = {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    for (auto* mf : fields) {
        if (mf->nGrow() < 2) { amrex::Abort("gust: fill_ghosts: the state needs two ghost cells"); }
        mf->FillBoundary(m_geom.periodicity());
    }
    const InflowProfile::View prof = m_profile.view();
    SurfaceState local;
    SurfaceState* ss = (surface_out != nullptr) ? surface_out : &local;
    amrex::MultiFab ground_grad;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (d == 2 && ground_wall && !homogeneous) {
            // The surface layer reads the lateral ghost columns filled above.
            compute_surface_state(q, *surface_opt, buoyancy, m_mesh, m_geom, *ss);
            // Two lateral ghost columns, as the z fills span (the state's
            // ngrow). The first is computed from the surface state's own
            // ghost column; FillBoundary gives every column that images a
            // domain column (the neighbouring boxes and the periodic shifts);
            // what is left (outside a non-periodic face beyond the first
            // ghost column, and the corners there) takes the gradient of the
            // nearest column inside the domain along the non-periodic
            // directions (the state's ghost columns there are the profile's
            // or mirrors; no valid row reads those corner ghosts).
            const amrex::Real unset = std::numeric_limits<amrex::Real>::max();
            ground_grad.define(ss->ustar.boxArray(), ss->ustar.DistributionMap(), 3, amrex::IntVect(2, 2, 0));
            ground_grad.setVal(unset);
            // The wall law's gradient per column: u*/(kappa (d0 + z0)) along the
            // tangential wind for u and v (components 0, 1), -q/(kappa u* (d0 +
            // z0)) for theta (component 2, 0 when neutral).
            const amrex::Real kappa = surface_opt->kappa, z0 = surface_opt->z0;
            const amrex::Real qs = (surface_opt->stability == SurfaceLayerOptions::Stability::heat_flux)
                                       ? surface_opt->heat_flux : amrex::Real(0);
            for (amrex::MFIter mfi(ground_grad, false); mfi.isValid(); ++mfi) {
                const auto g = ground_grad.array(mfi);
                const auto us = ss->ustar.const_array(mfi);
                const auto d0 = ss->d0.const_array(mfi);
                const auto ta = ss->tau.const_array(mfi);
                const auto il = ss->inv_L.const_array(mfi);
                amrex::ParallelFor(mfi.growntilebox(amrex::IntVect(1, 1, 0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // MOST's flux-gradient relations, phi_m and phi_h of
                    // zeta = d0/L (1 when neutral), with the wall model's clamp.
                    const amrex::Real zeta = amrex::max(amrex::min(d0(i, j, k) * il(i, j, k), amrex::Real(100)), -amrex::Real(100));
                    amrex::Real phi_m, phi_h;
                    most_phi(zeta, phi_m, phi_h);
                    const amrex::Real s = kappa * (d0(i, j, k) + z0);
                    const amrex::Real tmag = std::sqrt(ta(i, j, k, 0) * ta(i, j, k, 0) + ta(i, j, k, 1) * ta(i, j, k, 1)
                                                       + ta(i, j, k, 2) * ta(i, j, k, 2));
                    // The unit stress; with no stress (tmag = 0) the gradient is 0.
                    // No select or max: the compiler evaluates both sides of one
                    // and moves it past the division, and 0/0 would raise the
                    // invalid trap although discarded. The added term is below
                    // the resolution of any stress (bit-identical quotient).
                    const amrex::Real tsafe = tmag + std::numeric_limits<amrex::Real>::min();
                    const amrex::Real dir_x = ta(i, j, k, 0) / tsafe;
                    const amrex::Real dir_y = ta(i, j, k, 1) / tsafe;
                    g(i, j, k, 0) = us(i, j, k) / s * dir_x * phi_m;
                    g(i, j, k, 1) = us(i, j, k) / s * dir_y * phi_m;
                    g(i, j, k, 2) = -qs / (us(i, j, k) * s) * phi_h;
                });
            }
            ground_grad.FillBoundary(m_geom.periodicity());
            {
                const amrex::Box& dom = m_geom.Domain();
                const bool perx = m_geom.isPeriodic(0), pery = m_geom.isPeriodic(1);
                for (amrex::MFIter mfi(ground_grad, false); mfi.isValid(); ++mfi) {
                    const auto g = ground_grad.array(mfi);
                    amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        if (g(i, j, k, 0) != unset) { return; }
                        const int ic = perx ? i : amrex::max(dom.smallEnd(0), amrex::min(dom.bigEnd(0), i));
                        const int jc = pery ? j : amrex::max(dom.smallEnd(1), amrex::min(dom.bigEnd(1), j));
                        for (int n = 0; n < 3; ++n) { g(i, j, k, n) = g(ic, jc, k, n); }
                    });
                }
            }
        }
        if (m_geom.isPeriodic(d)) { continue; }
        const amrex::MultiFab* gg = ground_grad.empty() ? nullptr : &ground_grad;
        fill_field<stencil::XF>(q.u, f_u, d, m_in, q.z_phys_nd, m_geom, prof, m_coriolis, gg, homogeneous);
        fill_field<stencil::YF>(q.v, f_v, d, m_in, q.z_phys_nd, m_geom, prof, m_coriolis, gg, homogeneous);
        fill_field<stencil::ZF>(q.w, f_w, d, m_in, q.z_phys_nd, m_geom, prof, m_coriolis, gg, homogeneous);
        if (d == 2) { fill_no_flux_w(q, m_in, m_geom); }
        fill_field<stencil::CC>(q.p, f_p, d, m_in, q.z_phys_nd, m_geom, prof, m_coriolis, gg, homogeneous);
        fill_field<stencil::CC>(q.theta, f_theta, d, m_in, q.z_phys_nd, m_geom, prof, m_coriolis, gg, homogeneous);
        fill_field<stencil::CC>(q.tke, f_tke, d, m_in, q.z_phys_nd, m_geom, prof, m_coriolis, gg, homogeneous);
    }
}

namespace {

//! Whether the normal velocity of direction d has a Dirichlet value on face
//! (d, lo), and that value (profile: NaN marks "the inflow profile").
bool dirichlet_normal (const BoundaryInputs& in, int d, bool lo, FaceType& t)
{
    t = in.face[amrex::Orientation(d, lo ? amrex::Orientation::low : amrex::Orientation::high)];
    return t == FaceType::inflow || t == FaceType::slip_wall || t == FaceType::no_slip_wall
        || t == FaceType::surface_layer || t == FaceType::geostrophic;
}

//! The ground and top faces of w take the no-flux value Omega = 0,
//! w = z_x u_f + z_y v_f (MassFlux_K.H; 0 over flat ground), with the two
//! ghost layers mirrored about it (2 V - the interior), after the generic
//! value fill of w in z. u and v are filled (their lateral ghosts) before w.
//! The same formula serves the homogeneous fill (it is linear).
// With the homogeneous fill the value is linear in the perturbation's u, v
// given their homogeneous ghosts (the wall law's gradient of the nonlinear
// fill is not differentiated here; no production caller fills homogeneously).
void fill_no_flux_w (State& q, const BoundaryInputs& in, const amrex::Geometry& geom)
{
    const amrex::Box& dom = geom.Domain();
    const auto dxinv = geom.InvCellSizeArray();
    for (const bool lo : {true, false}) {
        FaceType t;
        if (!dirichlet_normal(in, 2, lo, t)) { continue; }
        const int face = lo ? dom.smallEnd(2) : dom.bigEnd(2) + 1;
        for (amrex::MFIter mfi(q.w, false); mfi.isValid(); ++mfi) {
            amrex::Box fb = mfi.validbox();
            if (fb.smallEnd(2) > face || fb.bigEnd(2) < face) { continue; }
            fb.setSmall(2, face);
            fb.setBig(2, face);
            fb.grow(0, q.w.nGrow());
            fb.grow(1, q.w.nGrow());
            fb &= mfi.fabbox();
            // Only the faces whose four u, v points and four nodal heights
            // exist in this box's fabs (the heights have one lateral ghost
            // layer, w two: the outermost lateral ghost column is set below
            // from its periodic image, and on a non-periodic side keeps the
            // lateral rule's value, which only the Dirichlet rows read).
            const amrex::Box ub = q.u[mfi].box(), vb = q.v[mfi].box(), zb = q.z_phys_nd[mfi].box();
            fb.setSmall(0, amrex::max(fb.smallEnd(0), ub.smallEnd(0), vb.smallEnd(0), zb.smallEnd(0)));
            fb.setBig(0, amrex::min(fb.bigEnd(0), ub.bigEnd(0) - 1, vb.bigEnd(0), zb.bigEnd(0) - 1));
            fb.setSmall(1, amrex::max(fb.smallEnd(1), ub.smallEnd(1), vb.smallEnd(1), zb.smallEnd(1)));
            fb.setBig(1, amrex::min(fb.bigEnd(1), ub.bigEnd(1), vb.bigEnd(1) - 1, zb.bigEnd(1) - 1));
            if (fb.isEmpty()) { continue; }
            const auto w = q.w.array(mfi);
            const auto uu = q.u.const_array(mfi), vv = q.v.const_array(mfi);
            const auto zn = q.z_phys_nd.const_array(mfi);
            const int ng = q.w.nGrow();
            const int sgn = lo ? -1 : 1;
            amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const amrex::Real V = massflux::no_flux_w(zn, uu, vv, i, j, k, dxinv[0], dxinv[1]);
                w(i, j, k) = V;
                for (int m = 1; m <= ng; ++m) { w(i, j, k + sgn * m) = 2 * V - w(i, j, k - sgn * m); }
            });
        }
    }
    // The periodic images of the face rows (the lateral ghost columns the
    // heights do not reach).
    q.w.FillBoundary(geom.periodicity());
}

} // namespace

void BoundaryConditions::apply_rows (const State& q, ResidualFields& R, bool homogeneous) const
{
    const amrex::Box& dom = m_geom.Domain();
    const int klo = dom.smallEnd(2);
    const InflowProfile::View prof = m_profile.view();
    const amrex::MultiFab* qf[3] = {&q.u, &q.v, &q.w};
    amrex::MultiFab* rf[3] = {&R.u, &R.v, &R.w};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (m_geom.isPeriodic(d)) { continue; }
        for (const bool lo : {true, false}) {
            FaceType t;
            if (!dirichlet_normal(m_in, d, lo, t)) { continue; }
            const bool is_profile = (t == FaceType::inflow) && !homogeneous;
            const int face = lo ? dom.smallEnd(d) : dom.bigEnd(d) + 1;
            for (amrex::MFIter mfi(*rf[d], false); mfi.isValid(); ++mfi) {
                amrex::Box fb = mfi.validbox();
                if (fb.smallEnd(d) > face || fb.bigEnd(d) < face) { continue; }
                fb.setSmall(d, face);
                fb.setBig(d, face);
                const auto r = rf[d]->array(mfi);
                const auto a = qf[d]->const_array(mfi);
                const amrex::Box zb = q.z_phys_nd[mfi].box();
                const NodeZ nz{q.z_phys_nd.const_array(mfi), amrex::lbound(zb), amrex::ubound(zb)};
                const int f = d;
                if (d == 0) {
                    amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        const amrex::Real V = is_profile ? profile_value(prof, f, height_above_ground<stencil::XF>(nz, i, j, k, klo)) : amrex::Real(0);
                        r(i, j, k) = a(i, j, k) - V;
                    });
                } else if (d == 1) {
                    amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        const amrex::Real V = is_profile ? profile_value(prof, f, height_above_ground<stencil::YF>(nz, i, j, k, klo)) : amrex::Real(0);
                        r(i, j, k) = a(i, j, k) - V;
                    });
                } else {
                    // The ground and the top: w minus the no-flux value
                    // (Omega = 0; 0 over flat ground). Linear in (u, v, w):
                    // the homogeneous fill uses the same formula.
                    const auto uu = q.u.const_array(mfi), vv = q.v.const_array(mfi);
                    const auto zn = q.z_phys_nd.const_array(mfi);
                    const auto dxinv = m_geom.InvCellSizeArray();
                    amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        r(i, j, k) = a(i, j, k) - massflux::no_flux_w(zn, uu, vv, i, j, k, dxinv[0], dxinv[1]);
                    });
                }
            }
        }
    }
}

void BoundaryConditions::dirichlet_rows (Equation eq, Loc loc, const Mesh& mesh, amrex::MultiFab& mask) const
{
    if (mask.empty()) {
        const amrex::MultiFab* like = (loc == Loc::xface) ? &mesh.ax : (loc == Loc::yface) ? &mesh.ay : (loc == Loc::zface) ? &mesh.az : &mesh.detJ;
        mask.define(like->boxArray(), mesh.detJ.DistributionMap(), 1, 0);
    }
    check_operator_field("dirichlet_rows", "mask", mask, loc, mesh, 0);
    mask.setVal(0.0);
    const amrex::Box& dom = m_geom.Domain();
    const int f = (eq == Equation::u) ? f_u : (eq == Equation::v) ? f_v : (eq == Equation::w) ? f_w
                : (eq == Equation::theta) ? f_theta : f_tke;
    const int t = static_cast<int>(loc);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (m_geom.isPeriodic(d)) { continue; }
        for (const bool lo : {true, false}) {
            const FaceType ft = m_in.face[amrex::Orientation(d, lo ? amrex::Orientation::low : amrex::Orientation::high)];
            const Rule rule = rule_for(static_cast<Field>(f), ft, d, lo);
            const bool nodal_d = stencil::nodal(t, d);
            if (!(nodal_d && (rule == Rule::value || rule == Rule::profile))) { continue; }
            const int edge = lo ? dom.smallEnd(d) : dom.bigEnd(d) + 1;
            for (amrex::MFIter mfi(mask, false); mfi.isValid(); ++mfi) {
                amrex::Box fb = mfi.validbox();
                if (fb.smallEnd(d) > edge || fb.bigEnd(d) < edge) { continue; }
                fb.setSmall(d, edge);
                fb.setBig(d, edge);
                mask[mfi].template setVal<amrex::RunOn::Device>(amrex::Real(1), fb, 0, 1);
            }
        }
    }
}

void BoundaryConditions::apply_to_implicit (Equation eq, Loc loc, const amrex::MultiFab& z_nd, amrex::MultiFab& A) const
{
    // The damping first: the Dirichlet face points below become identity rows.
    add_damping_to_implicit(eq, loc, z_nd, A);
    const amrex::Box& dom = m_geom.Domain();
    const int f = (eq == Equation::u) ? f_u : (eq == Equation::v) ? f_v : (eq == Equation::w) ? f_w
                : (eq == Equation::theta) ? f_theta : f_tke;
    const int t = static_cast<int>(loc);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (m_geom.isPeriodic(d)) { continue; }
        for (const bool lo : {true, false}) {
            const FaceType ft = m_in.face[amrex::Orientation(d, lo ? amrex::Orientation::low : amrex::Orientation::high)];
            const Rule rule = rule_for(static_cast<Field>(f), ft, d, lo);
            // A field nodal in d is the normal velocity, with a valid point on
            // the face: a Dirichlet value makes that row the identity; an even
            // mirror about the point (outflow) couples it to its other neighbour.
            const bool nodal_d = stencil::nodal(t, d);
            const bool dirichlet_pt = nodal_d && (rule == Rule::value || rule == Rule::profile);
            const int edge = lo ? dom.smallEnd(d) : (nodal_d ? dom.bigEnd(d) + 1 : dom.bigEnd(d));
            const int cm = (d == 0) ? coef_XM : (d == 1) ? coef_YM : coef_ZM;
            const int cp = (d == 0) ? coef_XP : (d == 1) ? coef_YP : coef_ZP;
            const int ghost_coef = lo ? cm : cp, other_coef = lo ? cp : cm;
            // Fold sign: +1 value/odd, -1 even, 0 profile; gradient moves the
            // coefficient to the opposite neighbour.
            amrex::Real fold = 0;
            bool to_other = false;
            if (nodal_d) { to_other = true; }
            else if (rule == Rule::even) { fold = -1; }
            else if (rule == Rule::odd || rule == Rule::value) { fold = 1; }
            else if (rule == Rule::gradient) { to_other = true; }
            for (amrex::MFIter mfi(A, false); mfi.isValid(); ++mfi) {
                amrex::Box fb = mfi.validbox();
                if (fb.smallEnd(d) > edge || fb.bigEnd(d) < edge) { continue; }
                fb.setSmall(d, edge);
                fb.setBig(d, edge);
                const auto a = A.array(mfi);
                if (dirichlet_pt) {
                    // Dirichlet face point: the identity row.
                    amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        for (int n = 0; n < n_implicit_coefs; ++n) { a(i, j, k, n) = 0; }
                        a(i, j, k, coef_P) = 1;
                    });
                } else {
                    amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        const amrex::Real ag = a(i, j, k, ghost_coef);
                        a(i, j, k, ghost_coef) = 0;
                        if (to_other) { a(i, j, k, other_coef) += ag; }
                        else { a(i, j, k, coef_P) += fold * ag; }
                    });
                }
            }
        }
    }
}

namespace {

//! The damping coefficient c(x, y, z) of both layers (1/s), device-usable.
struct DampingCoef
{
    bool rayleigh, sponge;
    amrex::Real zd, rc, ztop, sw, inv_w, ss;
    amrex::GpuArray<amrex::Real, 3> plo, phi;
    amrex::GpuArray<int, 4> lateral;   //!< xlo, xhi, ylo, yhi hold an inflow or outflow
    //! The Rayleigh part alone (1/s).
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real rayleigh_part (amrex::Real z) const noexcept
    {
        const amrex::Real pi = amrex::Real(3.14159265358979323846);
        if (rayleigh && z > zd) {
            const amrex::Real s = std::sin(amrex::Real(0.5) * pi * (z - zd) / (ztop - zd));
            return rc * s * s;
        }
        return 0;
    }
    //! The sponge part alone (1/s).
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real sponge_part (amrex::Real x, amrex::Real y) const noexcept
    {
        amrex::Real c = 0;
        if (sponge) {
            const amrex::Real dist[4] = {x - plo[0], phi[0] - x, y - plo[1], phi[1] - y};
            for (int n = 0; n < 4; ++n) {
                if (lateral[n] && dist[n] < sw) {
                    const amrex::Real r = 1 - dist[n] * inv_w;   // inv_w: 1/sw, finite when off
                    c += ss * r * r;
                }
            }
        }
        return c;
    }
    //! Both parts; with_rayleigh false leaves the Rayleigh part out (a
    //! field bc.rayleigh_fields does not name).
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (amrex::Real x, amrex::Real y, amrex::Real z, bool with_rayleigh = true) const noexcept
    {
        return (with_rayleigh ? rayleigh_part(z) : amrex::Real(0)) + sponge_part(x, y);
    }
};

//! ztop: the mesh's top level, where the Rayleigh ramp ends.
DampingCoef damping_coef (const BoundaryInputs& in, const amrex::Geometry& geom, amrex::Real ztop)
{
    DampingCoef d{};
    d.zd = in.rayleigh_zdamp;
    d.rc = in.rayleigh_dampcoef;
    d.ztop = ztop;
    d.sw = in.sponge_width;
    d.ss = in.sponge_strength;
    d.rayleigh = d.zd > 0;
    d.sponge = d.sw > 0;
    d.inv_w = d.sponge ? 1 / d.sw : amrex::Real(0);   // no x/0 for the compiler to evaluate early
    d.plo = geom.ProbLoArray();
    d.phi = geom.ProbHiArray();
    for (int dd = 0; dd < 2; ++dd) {
        for (const bool lo : {true, false}) {
            const FaceType t = in.face[amrex::Orientation(dd, lo ? amrex::Orientation::low : amrex::Orientation::high)];
            d.lateral[2 * dd + (lo ? 0 : 1)] = (t == FaceType::inflow || t == FaceType::outflow) ? 1 : 0;
        }
    }
    return d;
}

//! Height above the ground and absolute height of the point (i, j, k) of type t.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void point_heights (int t, const NodeZ& nz, int i, int j, int k, int klo, amrex::Real& h, amrex::Real& zc) noexcept
{
    if (t == stencil::XF) { h = height_above_ground<stencil::XF>(nz, i, j, k, klo); zc = point_z<stencil::XF>(nz, i, j, k); }
    else if (t == stencil::YF) { h = height_above_ground<stencil::YF>(nz, i, j, k, klo); zc = point_z<stencil::YF>(nz, i, j, k); }
    else if (t == stencil::ZF) { h = height_above_ground<stencil::ZF>(nz, i, j, k, klo); zc = point_z<stencil::ZF>(nz, i, j, k); }
    else { h = height_above_ground<stencil::CC>(nz, i, j, k, klo); zc = point_z<stencil::CC>(nz, i, j, k); }
}

} // namespace

void BoundaryConditions::add_damping (const State& q, ResidualFields& R) const
{
    if (!damping_on()) { return; }
    const int klo = m_geom.Domain().smallEnd(2);
    const InflowProfile::View prof = m_profile.view();
    const DampingCoef coef = damping_coef(m_in, m_geom, m_mesh.zlevels().back());
    const auto plo = m_geom.ProbLoArray();
    const auto dx = m_geom.CellSizeArray();
    const amrex::MultiFab* qf[5] = {&q.u, &q.v, &q.w, &q.theta, &q.tke};
    amrex::MultiFab* rf[5] = {&R.u, &R.v, &R.w, &R.theta, &R.tke};
    const int ff[5] = {f_u, f_v, f_w, f_theta, f_tke};
    for (int n = 0; n < 5; ++n) {
        const int f = ff[n];
        const bool with_rayleigh = m_in.rayleigh_fields[static_cast<std::size_t>(n)];
        const bool above_ground = m_in.rayleigh_above_ground;
        for (amrex::MFIter mfi(*rf[n], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.tilebox();
            const auto r = rf[n]->array(mfi);
            const auto a = qf[n]->const_array(mfi);
            const amrex::Box zb = q.z_phys_nd[mfi].box();
            const NodeZ nz{q.z_phys_nd.const_array(mfi), amrex::lbound(zb), amrex::ubound(zb)};
            const int tt = (n < 3) ? (1 << n) : 0;
            const bool nx = (tt & 1) != 0, ny = (tt & 2) != 0;
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                amrex::Real h, zc;
                point_heights(tt, nz, i, j, k, klo, h, zc);
                const amrex::Real x = plo[0] + (i + (nx ? 0 : amrex::Real(0.5))) * dx[0];
                const amrex::Real y = plo[1] + (j + (ny ? 0 : amrex::Real(0.5))) * dx[1];
                // The Rayleigh layer relaxes towards the free atmosphere:
                // theta towards the sounding at the point's absolute height
                // (the terrain-following profile height would impose a warm
                // anomaly over a crest), u and v towards the profile, w to 0.
                // The sponge relaxes towards the inflow profile (the height
                // above the ground, as the inflow itself).
                const amrex::Real zr = above_ground ? h : zc;   // the layer's height coordinate
                const amrex::Real cr = with_rayleigh ? coef.rayleigh_part(zr) : amrex::Real(0);
                const amrex::Real cs = coef.sponge_part(x, y);
                if (cr > 0) {
                    const amrex::Real ref = (f == f_w) ? amrex::Real(0) : profile_value(prof, f, (f == f_theta) ? zr : h);
                    r(i, j, k) += cr * (a(i, j, k) - ref);
                }
                if (cs > 0) {
                    const amrex::Real ref = (f == f_w) ? amrex::Real(0) : profile_value(prof, f, h);
                    r(i, j, k) += cs * (a(i, j, k) - ref);
                }
            });
        }
    }
}

void BoundaryConditions::add_damping_to_implicit (Equation eq, Loc loc, const amrex::MultiFab& z_nd, amrex::MultiFab& A) const
{
    if (!damping_on()) { return; }
    const int klo = m_geom.Domain().smallEnd(2);
    const DampingCoef coef = damping_coef(m_in, m_geom, m_mesh.zlevels().back());
    const int n = (eq == Equation::u) ? 0 : (eq == Equation::v) ? 1 : (eq == Equation::w) ? 2 : (eq == Equation::theta) ? 3 : 4;
    const bool with_rayleigh = m_in.rayleigh_fields[static_cast<std::size_t>(n)];
    const bool above_ground = m_in.rayleigh_above_ground;
    const auto plo = m_geom.ProbLoArray();
    const auto dx = m_geom.CellSizeArray();
    const int tt = static_cast<int>(loc);
    const bool nx = (tt & 1) != 0, ny = (tt & 2) != 0;
    for (amrex::MFIter mfi(A, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto a = A.array(mfi);
        const amrex::Box zb = z_nd[mfi].box();
        const NodeZ nz{z_nd.const_array(mfi), amrex::lbound(zb), amrex::ubound(zb)};
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            amrex::Real h, zc;
            point_heights(tt, nz, i, j, k, klo, h, zc);
            const amrex::Real x = plo[0] + (i + (nx ? 0 : amrex::Real(0.5))) * dx[0];
            const amrex::Real y = plo[1] + (j + (ny ? 0 : amrex::Real(0.5))) * dx[1];
            a(i, j, k, coef_P) += (with_rayleigh ? coef.rayleigh_part(above_ground ? h : zc) : amrex::Real(0)) + coef.sponge_part(x, y);
        });
    }
}

} // namespace gust

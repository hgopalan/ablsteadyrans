#ifndef ABLSTEADY_INFLOWPROFILE_H_
#define ABLSTEADY_INFLOWPROFILE_H_

// Inflow profiles by height above the ground (sub-phase 1.8): the wind,
// potential temperature and turbulence that enter through an inflow face,
// and that the damping layers relax towards.
//
// Wind speed S(z) above the ground, along a fixed direction (u, v) = S (cos a,
// sin a), a the wind direction in degrees from the x axis:
//   log_law    S = (u*/kappa) ln((z + z0)/z0)                 (the wall law of SurfaceLayer.H)
//   power_law  S = S_ref (z / z_ref)^alpha  (its slope below z0 is the slope at z0)
//   file       S, direction, theta and (optionally) k tabulated in a file:
//              "z u v theta [k]" per line, z strictly increasing from 0,
//              linear in between, constant beyond the ends.
// theta(z) is the hydrostatic reference profile (ReferenceProfile) for the
// analytic kinds.
//
// Turbulence in equilibrium with the wind: with nu_t = kappa u*_l (z + z0)
// and the local friction velocity u*_l = kappa (z + z0) dS/dz (shear
// production equals dissipation in the wall law),
//   k = u*_l^2 / Cmu0^2,   nu_t = kappa u*_l (z + z0).
// For the log law u*_l = u* at every height, so k = u*^2/Cmu0^2, the wall
// value. A file with a k column uses it as given.

#include "Forcing.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <string>
#include <vector>

namespace ablsteady {

struct InflowProfileOptions
{
    enum class Kind { log_law, power_law, file };
    Kind kind{Kind::log_law};
    //! log_law: friction velocity (m/s, > 0) and roughness (m, > 0).
    amrex::Real ustar{amrex::Real(0.4)};
    amrex::Real z0{amrex::Real(0.1)};
    //! power_law: speed S_ref (>= 0) at height z_ref (> 0), exponent alpha (> 0).
    amrex::Real speed_ref{10};
    amrex::Real z_ref{10};
    amrex::Real alpha{amrex::Real(0.14)};
    //! Wind direction (degrees from the x axis), log_law and power_law.
    amrex::Real direction{0};
    //! file: the table's path.
    std::string file;
    //! For k and nu_t: Cmu0 and kappa (> 0; the closure's values).
    amrex::Real Cmu0{amrex::Real(0.5562)};
    amrex::Real kappa{amrex::Real(0.41)};
};

//! Problems with the options, naming each (the file is not read here).
std::vector<std::string> check_inflow_profile_options (const InflowProfileOptions& opt);

//! The profile, ready to evaluate on the device: the analytic kinds from
//! their parameters, a file as a table copied to the device.
class InflowProfile
{
public:
    //! Build from the options and the reference theta profile (used by the
    //! analytic kinds). Aborts, naming bc.inflow_file, on a file that cannot
    //! be read or is not "z u v theta [k]" with z increasing from 0.
    InflowProfile (const InflowProfileOptions& opt, const ReferenceProfile& theta_ref);

    //! A device-usable view.
    struct View
    {
        int kind;                  //!< 0 log_law, 1 power_law, 2 file
        amrex::Real ustar, z0, speed_ref, z_ref, alpha, cosd, sind, Cmu0, kappa;
        const amrex::Real* tz;     //!< file table: heights (n)
        const amrex::Real* tu;     //!< file table: u, v, theta, k, each n long (k: NaN if absent)
        const amrex::Real* tv;
        const amrex::Real* tt;
        const amrex::Real* tk;
        int n;
        bool has_k;
        const amrex::Real* rz;     //!< reference theta profile (m)
        const amrex::Real* rt;
        int rn;

        //! The values at height z (>= 0) above the ground.
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void eval (amrex::Real z, amrex::Real& u, amrex::Real& v, amrex::Real& theta, amrex::Real& k,
                   amrex::Real& nut) const noexcept;
    };
    [[nodiscard]] View view () const;
    //! The same view over host copies of the tables (for host-side tests and
    //! tools; the device view's pointers are device memory).
    [[nodiscard]] View host_view () const;

    //! Table sizes (for tests).
    [[nodiscard]] int table_size () const { return static_cast<int>(m_z.size()); }
    [[nodiscard]] bool has_k () const { return m_has_k; }

private:
    InflowProfileOptions m_opt;
    bool m_has_k{false};
    amrex::Gpu::DeviceVector<amrex::Real> m_z, m_u, m_v, m_t, m_k, m_rz, m_rt;
    amrex::Vector<amrex::Real> m_hz, m_hu, m_hv, m_ht, m_hk, m_hrz, m_hrt;   //!< host copies
    View view_with (const amrex::Real* tz, const amrex::Real* tu, const amrex::Real* tv, const amrex::Real* tt,
                    const amrex::Real* tk, const amrex::Real* rz, const amrex::Real* rt) const;
};

//! Linear interpolation in a table of n increasing abscissae (constant
//! beyond the ends).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real table_interp (const amrex::Real* x, const amrex::Real* y, int n, amrex::Real xq) noexcept
{
    if (n == 1 || xq <= x[0]) { return y[0]; }
    if (xq >= x[n - 1]) { return y[n - 1]; }
    int m = 1;
    while (m < n - 1 && !(xq < x[m])) { ++m; }
    const amrex::Real w = (xq - x[m - 1]) / (x[m] - x[m - 1]);
    return y[m - 1] + w * (y[m] - y[m - 1]);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void InflowProfile::View::eval (amrex::Real z, amrex::Real& u, amrex::Real& v, amrex::Real& theta,
                                amrex::Real& k, amrex::Real& nut) const noexcept
{
    z = amrex::max(z, amrex::Real(0));
    amrex::Real ustar_l;
    if (kind == 2) {
        u = table_interp(tz, tu, n, z);
        v = table_interp(tz, tv, n, z);
        theta = table_interp(tz, tt, n, z);
        // Local friction velocity from the table's shear, kappa (z + z0) dS/dz
        // with z0 the options' (the speed's own gradient between the table
        // points, one-sided at the ends).
        int m = 1;
        while (m < n - 1 && !(z < tz[m])) { ++m; }
        const amrex::Real s0 = std::sqrt(tu[m - 1] * tu[m - 1] + tv[m - 1] * tv[m - 1]);
        const amrex::Real s1 = std::sqrt(tu[m] * tu[m] + tv[m] * tv[m]);
        const amrex::Real dsdz = (n > 1) ? (s1 - s0) / (tz[m] - tz[m - 1]) : amrex::Real(0);
        ustar_l = kappa * (z + z0) * amrex::max(dsdz, amrex::Real(0));
        if (has_k) {
            k = table_interp(tz, tk, n, z);
            ustar_l = Cmu0 * std::sqrt(amrex::max(k, amrex::Real(0)));
        } else {
            k = ustar_l * ustar_l / (Cmu0 * Cmu0);
        }
        nut = kappa * ustar_l * (z + z0);
        return;
    }
    amrex::Real S;
    if (kind == 0) {
        S = ustar / kappa * std::log((z + z0) / z0);
        ustar_l = ustar;
    } else {
        S = speed_ref * std::pow(z / z_ref, alpha);
        // dS/dz = alpha S / z diverges at the ground: below z0 the slope at
        // z0 is used, so u*_l -> kappa alpha S(z0) there instead of 0 or infinity.
        // max(z, z0) in arithmetic: a select would let the compiler evaluate
        // S/z at z = 0 (0/0, the invalid trap) although unused.
        const amrex::Real dsdz = alpha * S / (amrex::Real(0.5) * (z + z0 + std::abs(z - z0)));
        ustar_l = kappa * (z + z0) * dsdz;
    }
    u = S * cosd;
    v = S * sind;
    theta = table_interp(rz, rt, rn, z);
    k = ustar_l * ustar_l / (Cmu0 * Cmu0);
    nut = kappa * ustar_l * (z + z0);
}

} // namespace ablsteady

#endif

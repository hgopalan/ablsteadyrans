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
// Turbulence in equilibrium with the wind (shear production equals
// dissipation of the k-L closure, Closure_K.H, neutral): with the length l and
// the local friction velocity u*_l = l dS/dz,
//   k = u*_l^2 / Cmu0^2,   nu_t = l u*_l.
// bc.inflow_tke_length = closure (the default) takes the closure's own
// geometric length, kappa (z + z0) capped by turb.max_geom_lscale:
//   1/l = 1/(kappa (z + z0)) + 1/l_max   (Closure_K.H geom_length);
// unbounded takes l = kappa (z + z0), the wall-law length (the form before
// it), for which the log law gives u*_l = u* and k = u*^2/Cmu0^2 at every
// height. A file with a k column uses its k as given (bc.inflow_tke_length is
// then refused: the fill, the damping and the first guess use k, not the
// profile's nu_t, which eval reports for tests and tools).

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
    //! bc.inflow_tke_length: the length of the equilibrium k and nu_t; closure
    //! caps kappa (z + z0) at l_max (turb.max_geom_lscale, > 0), unbounded does not.
    enum class Length { closure, unbounded };
    Length length{Length::closure};
    amrex::Real l_max{30};
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
    //! A tabulated profile given in memory (sub-phase 2.1b, the converged
    //! column of ColumnGuess.H): the file kind with the table z u v theta k
    //! (z strictly increasing from 0, at least two rows, theta > 0, k >= 0;
    //! aborts naming "column guess" otherwise) and the sounding.
    InflowProfile (const amrex::Vector<amrex::Real>& z, const amrex::Vector<amrex::Real>& u,
                   const amrex::Vector<amrex::Real>& v, const amrex::Vector<amrex::Real>& theta,
                   const amrex::Vector<amrex::Real>& k, const ReferenceProfile& theta_ref);

    //! A device-usable view.
    struct View
    {
        int kind;                  //!< 0 log_law, 1 power_law, 2 file
        amrex::Real ustar, z0, speed_ref, z_ref, alpha, cosd, sind, Cmu0, kappa;
        amrex::Real inv_lmax;      //!< 1/l_max, or 0 for bc.inflow_tke_length = unbounded
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

        //! The cap's factor at height z, l / l_raw = 1 / (1 + l_raw / l_max)
        //! with l_raw = kappa (z + z0) (geom_length of Closure_K.H), written
        //! without a select: exactly 1 for unbounded (inv_lmax = 0), so the
        //! old expressions times this factor are the old values bit for bit.
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::Real cap_factor (amrex::Real z) const noexcept
        {
            return amrex::Real(1) / (amrex::Real(1) + kappa * (z + z0) * inv_lmax);
        }
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
        const amrex::Real c = cap_factor(z);
        ustar_l = kappa * (z + z0) * amrex::max(dsdz, amrex::Real(0)) * c;
        if (has_k) {
            k = table_interp(tz, tk, n, z);
            ustar_l = Cmu0 * std::sqrt(amrex::max(k, amrex::Real(0)));
        } else {
            k = ustar_l * ustar_l / (Cmu0 * Cmu0);
        }
        nut = kappa * ustar_l * (z + z0) * c;   // l u*_l (diagnostic: the fill and the first guess take k)
        return;
    }
    amrex::Real S;
    const amrex::Real c = cap_factor(z);
    if (kind == 0) {
        S = ustar / kappa * std::log((z + z0) / z0);
        // dS/dz = u* / (kappa (z + z0)): u*_l = l dS/dz = u* l / l_raw.
        ustar_l = ustar * c;
    } else {
        S = speed_ref * std::pow(z / z_ref, alpha);
        // dS/dz = alpha S / z diverges at the ground: below z0 the slope at
        // z0 is used, so u*_l -> kappa alpha S(z0) there instead of 0 or infinity.
        // max(z, z0) in arithmetic: a select would let the compiler evaluate
        // S/z at z = 0 (0/0, the invalid trap) although unused.
        const amrex::Real dsdz = alpha * S / (amrex::Real(0.5) * (z + z0 + std::abs(z - z0)));
        ustar_l = kappa * (z + z0) * dsdz * c;
    }
    u = S * cosd;
    v = S * sind;
    theta = table_interp(rz, rt, rn, z);
    k = ustar_l * ustar_l / (Cmu0 * Cmu0);
    nut = kappa * ustar_l * (z + z0) * c;   // l u*_l (diagnostic: the fill and the first guess take k)
}

} // namespace ablsteady

#endif

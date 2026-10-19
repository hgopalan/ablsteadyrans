#ifndef ABLSTEADY_LINESWEEPS_H_
#define ABLSTEADY_LINESWEEPS_H_

// The line LU-SGS sweeps (sub-phases 2.1 and 2.2): symmetric line
// Gauss-Seidel passes over the columns of each box, an exact z-line (Thomas)
// solve per column with the x/y neighbours already updated in the pass's
// direction, in place. One sweep is a forward and a backward pass.
//
// What a box reads across its faces (steady.line_ghosts):
//
//  * lagged (the default; the path of sub-phase 2.1, bit-identical to it):
//    before every pass delta.FillBoundary copies the neighbouring boxes'
//    values of the previous pass into the ghosts, a block Jacobi across
//    boxes with one exchange per pass. Each box sweeps its whole valid box,
//    so the two copies of a face shared by two boxes (and of a periodic seam
//    face) are both solved; delta.OverrideSync keeps the copy of AMReX's
//    owner (the lowest box index) at the end. No overlap.
//
//  * zero (sub-phase 2.2, box-local): no communication inside the sweeps.
//    Every box sweeps the region it owns, grown by n_o = steady.line_overlap
//    cells in x and y (clipped to the domain on a non-periodic axis; not
//    grown along an axis the domain is one cell wide, as y in a 2-D deck), with
//    delta = 0 outside that region: the ghost delta is zeroed, never lagged.
//    The system (A, a_P, rhs) on the grown region comes from one
//    FillBoundary of copies with n_o ghost cells (periodic images through
//    the geometry); the swept values outside the owned region are dropped.
//    This is a restricted additive Schwarz step per pseudo-time iteration
//    (the outer defect correction recomputes the residual, so the converged
//    state does not depend on n_o or on the box layout; the iterates do).
//    With n_o = 0 the sweeps are confined to the owned region. Along a
//    periodic axis one cell wide (y in a 2-D deck) the two neighbours are
//    the column itself, not another box: their coefficients are moved to
//    the diagonal (line_self_axes), which solves that coupling exactly.
//
// C-grid face ownership (zero): a box owns the faces of its valid box in its
// nodal direction except its high face, which belongs to the next box (as
// its low face); a high face on the domain's high side belongs to the box,
// on a periodic axis too, and the domain's low face on a periodic axis is
// then that face's image, owned by the last box, not the first. Every face
// is swept by exactly one box per pass (line_owned_box) and the non-owned
// copies take the owner's value at the end (OverrideSync with this owner
// mask; AMReX's default mask, the lowest box index, would keep a copy that
// was not swept). Cells and z-faces are owned by their box.
//
// The ordering of the columns inside a box (steady.line_order):
//
//  * sequential: j outer, i inner, ascending in the forward pass and
//    descending in the backward pass (one thread per box on a GPU);
//  * diagonal: the wavefronts i + j = s, one kernel launch per plane s
//    (ascending s forward, descending backward), one column per thread
//    inside a plane. A column's x/y neighbours lie on the planes s - 1 and
//    s + 1, so a forward plane reads the updated lower neighbours and the
//    old upper ones exactly as the sequential order does: the two orderings
//    perform the same operations on the same values and give bit-identical
//    results (UT-2.2.4).

#include "CVFluxes.H"
#include "Mesh.H"

#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>
#include <vector>

namespace ablsteady {

//! The line solver of the sweeps (steady.sweep).
enum class SweepKind : int { zjacobi = 0, lusgs = 1 };
//! What a box's LU-SGS sweep reads across its faces (steady.line_ghosts).
enum class LineGhosts : int { lagged = 0, zero = 1 };
//! The order of the columns inside a box (steady.line_order).
enum class LineOrder : int { sequential = 0, diagonal = 1 };

//! The widest overlap steady.line_overlap accepts (cells); the layout also
//! bounds it by the narrowest box (check_line_overlap_layout).
inline constexpr int max_line_overlap = 8;

//! The sweep options (steady.sweep, line_ghosts, line_overlap, line_order).
//! Implicit from a SweepKind: the other options keep their defaults (the
//! path of sub-phase 2.1).
struct LineSweepOptions
{
    LineSweepOptions (SweepKind k = SweepKind::zjacobi) noexcept : kind(k) {}   // NOLINT(google-explicit-constructor)
    SweepKind kind;
    LineGhosts ghosts{LineGhosts::lagged};
    int overlap{0};
    LineOrder order{LineOrder::sequential};
};

//! The names of the options ("zjacobi", "lusgs"; "lagged", "zero";
//! "sequential", "diagonal").
const char* sweep_name (SweepKind kind) noexcept;
const char* line_ghosts_name (LineGhosts g) noexcept;
const char* line_order_name (LineOrder o) noexcept;

//! The range checks of the LU-SGS options, each message naming its
//! steady.<key>: line_overlap an integer >= 0 and <= max_line_overlap, and 0
//! unless line_ghosts = zero. The kind is not checked here.
std::vector<std::string> check_line_options (const LineSweepOptions& opt);

//! 1 along x or y when that axis is periodic and the domain one cell wide.
amrex::IntVect line_self_axes (const amrex::Geometry& geom);

//! The overlap per axis: overlap in x and y, except along an axis the
//! domain (any index type) is one cell wide; 0 in z.
amrex::IntVect line_overlap_width (int overlap, const amrex::Box& domain);

//! Empty, or the message naming steady.line_overlap when the overlap is wider
//! than the narrowest box of ba along an axis it applies to
//! (line_overlap_width of ba's minimal box: the grown region would reach
//! past the neighbouring box; on a periodic axis past one period).
std::string check_line_overlap_layout (const amrex::BoxArray& ba, int overlap);

//! The points of vbox (a valid box of a field at loc) its box owns (the
//! header's face ownership); the cell box itself for cells and z-faces.
amrex::Box line_owned_box (Loc loc, const amrex::Box& vbox, const amrex::Geometry& geom);

//! The region a box sweeps with steady.line_ghosts = zero: owned grown by
//! line_overlap_width in x and y, clipped to the domain (at loc) on a non-periodic axis.
amrex::Box line_sweep_box (Loc loc, const amrex::Box& owned, int overlap, const amrex::Geometry& geom);

//! The checks shared by the line solvers (zline_jacobi, zline_lusgs): sweeps
//! >= 1, the layouts, every box spanning the column; defines delta if empty
//! (at loc, one ghost cell) and zeroes it.
void check_line_system (const char* op, Loc loc, const amrex::MultiFab& A, const amrex::MultiFab& ap,
                        const amrex::MultiFab& rhs, int sweeps, const Mesh& mesh, const amrex::Geometry& geom,
                        amrex::MultiFab& delta);

//! The line LU-SGS solve: sweeps symmetric passes of the system of
//! zline_jacobi (Predictor.H; the same arguments) with opt's ghosts,
//! overlap and ordering (opt.kind is not read). Returns max |change over
//! the last backward pass| (over the owned points with ghosts = zero).
//! Aborts, naming steady.line_overlap, on an overlap check_line_options or
//! check_line_overlap_layout refuses.
amrex::Real zline_lusgs (Loc loc, const amrex::MultiFab& A, const amrex::MultiFab& ap, const amrex::MultiFab& rhs,
                         int sweeps, const Mesh& mesh, const amrex::Geometry& geom, amrex::MultiFab& delta,
                         const LineSweepOptions& opt = LineSweepOptions(SweepKind::lusgs));

} // namespace ablsteady

#endif

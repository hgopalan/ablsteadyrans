#ifndef ABLSTEADY_ERF_SHIM_INDEX_DEFINES_H_
#define ABLSTEADY_ERF_SHIM_INDEX_DEFINES_H_

// Shim for ERF's ERF_IndexDefines.H, which the files in extern_erf/ include.
// They use only the boundary-condition enum below; it is ERF's ERF_BC
// (Source/ERF_IndexDefines.H at the commit in external/VERSIONS), copied with
// its enumerators in ERF's order. ERF's licence: extern_erf/LICENSE.md.

enum struct ERF_BC {
    symmetry,
    inflow,
    outflow,
    inflow_outflow,
    ho_outflow,
    open,
    no_slip_wall,
    slip_wall,
    periodic,
    surface_layer,
    undefined
};

#endif

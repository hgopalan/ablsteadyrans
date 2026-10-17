#ifndef ABLSTEADY_ERF_SHIM_UTILS_H_
#define ABLSTEADY_ERF_SHIM_UTILS_H_

// Shim for ERF's ERF_Utils.H. The real header pulls in most of ERF (data
// structures, surface layer, fill patchers); the files in extern_erf/ need
// only these declarations, whose definitions are in extern_erf/.

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

// extern_erf/ERF_ChopGrids.cpp
amrex::BoxArray ERFPostProcessBaseGrids (const amrex::Box& domain, bool decompose_in_z);
void ChopGrids2D (amrex::BoxArray& ba, const amrex::Box& domain, int target_size);

// extern_erf/ERF_TerrainMetrics.cpp
void make_J (const amrex::Geometry& geom, amrex::MultiFab& z_phys_nd, amrex::MultiFab& detJ_cc);
void make_areas (const amrex::Geometry& geom, amrex::MultiFab& z_phys_nd,
                 amrex::MultiFab& ax, amrex::MultiFab& ay, amrex::MultiFab& az);
void make_zcc (const amrex::Geometry& geom, amrex::MultiFab& z_phys_nd, amrex::MultiFab& z_phys_cc);

#endif

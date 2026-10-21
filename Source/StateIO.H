#ifndef GUST_STATEIO_H_
#define GUST_STATEIO_H_

// Named StateIO.H, not IO.H: Windows file systems ignore case, so an IO.H on
// the include path would stand in for the system <io.h>.

#include "State.H"

#include <string>

namespace gust {

//! Plotfile variable names, in order. Velocities are face averages to cell
//! centres; z is the mean of the cell's eight node heights.
const amrex::Vector<std::string>& plot_var_names ();

//! Write <prefix>NNNNN (step from state). Besides the cell-centred variables,
//! the nodal heights go to <plotfile>/Level_0/z_phys_nd as a VisMF MultiFab,
//! ghost nodes included. With detJ (the mesh's cell volume factor), it is
//! added as an eighth variable "detJ". Returns the directory name.
std::string write_plotfile (const std::string& prefix, const State& state,
                            const amrex::Geometry& geom, const amrex::MultiFab* detJ = nullptr);

//! Write <prefix>NNNNN: a text Header (format version, bytes per real, step,
//! mesh in hex floats, BoxArray) and every State field as a VisMF MultiFab
//! with ghosts.
//! Returns the directory name.
std::string write_checkpoint (const std::string& prefix, State& state,
                              const MeshInputs& mesh);

//! Read a checkpoint into a State already defined on the current layout.
//! Aborts, naming io.restart, if the Header is corrupt, was written with
//! another precision, or its mesh differs from mesh. Valid values are
//! restored exactly; of the ghost cells, only z_phys_nd's are restored. The
//! checkpoint may come from a different box layout or rank count.
void read_checkpoint (const std::string& dir, State& state, const MeshInputs& mesh);

} // namespace gust

#endif

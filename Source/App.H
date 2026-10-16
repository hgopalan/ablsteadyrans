#ifndef ABLSTEADY_APP_H_
#define ABLSTEADY_APP_H_

namespace ablsteady {

//! Read and check the inputs, build the mesh and state (fresh or from
//! io.restart), and write the plotfile and checkpoint.
void run ();

} // namespace ablsteady

#endif

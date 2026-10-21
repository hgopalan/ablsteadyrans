#ifndef GUST_APP_H_
#define GUST_APP_H_

namespace gust {

//! Read and check the inputs, build the mesh and state (fresh or from
//! io.restart), and write the plotfile and checkpoint.
void run ();

} // namespace gust

#endif

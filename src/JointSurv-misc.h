#ifndef JOINTSURV_MISC_H
#define JOINTSURV_MISC_H

#include "VA-joint-config.h"

namespace survival {

/**
 * struct to hold nodes, weights, and the number of nodes. In this namespace,
 * the quadrature that are used are some quadrature rule on the interval (0, 1)
 */
struct node_weight {
  /// the nodes should be between zero and one
  double const * ns, * ws;
  vajoint_uint n_nodes;
};

}

#endif

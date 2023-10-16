#pragma once

#include "defines.h"

/**
 * @brief Create a Random Binomial Graph with n vertices and m edges
 *        using the Erdös-Rény model
 * 
 * @param n         Number of vertices
 * @param m         Number of edges
 * @return Graph    Random graph
 */
Graph createRandomBinomialGraph(int n, int m);
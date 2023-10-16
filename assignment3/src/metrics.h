#pragma once

#include "defines.h"

/**
 * @brief Computes exactly the closeness centrality of the graph
 * 
 * @param g         Input graph
 * @return double   Closeness Centrality $ = \frac{1}{N}\sum_{i}\frac{1}{N-1}\sum_{j\neq i} d_{i,j}$
 */
double exactClosenessCentrality(const Graph& g);

/**
 * @brief Estimates the closeness centrality of the graph
 * 
 * @param g         Input graph
 * @param T         Number of repetitions
 * @param delta     Percentage of vertices to take into account. 0 <= delta <= 1
 * @return double   Closeness Centrality
 */
double montecarloClosenessCentrality(const Graph& g, int T, double delta);

/**
 * @brief Returns a lower and upper bound of the closeness centrality
 * 
 * @param g                         Input graph
 * @return std::pair<double,double> The first element contain a 
 *                                  lowerbound and the second one 
 *                                  is the upperbound
 */
std::pair<double,double> boundedClosenessCentrality(const Graph& g);
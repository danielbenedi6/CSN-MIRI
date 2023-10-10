#include "input.h"
#include <fstream>
#include <iostream>
#include <unordered_map>

bool read_file(std::string filename, Graph& g) {
    std::ifstream file(filename);
    if(!file.is_open()) {
        std::cerr << "Could not open file graph: " << filename << std::endl;
        return false;
    }

    int num_vertices = 0, max_vertices = 0, max_edges = 0;
    std::unordered_map<std::string,int> equivalence;
    file >> max_vertices >> max_edges;
    
    while(!file.eof()) {
        std::string vertex1, vertex2;
        int id1 = -1, id2 = -1;
        
        // Read label of both vertices
        file >> vertex1 >> vertex2;

        // Check if label already readed If so, use id. 
        // Otherwise, create new and assign
        if(equivalence.count(vertex1) == 0) {
            equivalence.insert({vertex1, num_vertices++});
            boost::add_vertex(g);
        }
        id1 = equivalence[vertex1];
        if(equivalence.count(vertex2) == 0) {
            equivalence.insert({vertex2, num_vertices++});
            boost::add_vertex(g);
        }
        id2 = equivalence[vertex2];

        // If no self-loop, add edge to graph g
        if(id1 != id2) {
            boost::add_edge(id1, id2, g);
        }
    }

    if(boost::num_vertices(g) != max_vertices) {
        std::cerr << "Wrong number of vertices read. File may be corrupted!" << std::endl;
        std::cerr << "Expected: " << max_vertices << ". Read: " << boost::num_vertices(g) << std::endl;
        return false;
    }

    if(boost::num_edges(g) != max_edges) {
        std::cerr << "Wrong number of edges read. File may be corrupted!" << std::endl;
        std::cerr << "Expected: " << max_edges << ". Read: " << boost::num_edges(g) << std::endl;
        return false;
    }

    return true;
}
#include "input.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

bool read_file(std::string filename, Graph& g) {
    std::ifstream file(filename);
    if(!file.is_open()) {
        std::cerr << "Could not open file graph: " << filename << std::endl;
        return false;
    }

    int num_vertices = 0, num_edges = 0, max_vertices = 0, max_edges = 0;
    std::unordered_map<std::string,int> equivalence;
    file >> max_vertices >> max_edges;
    file.get(); // Discard newline 
    
    while(!file.eof()) {
        std::string vertex1, vertex2, line;
        int id1 = -1, id2 = -1;
        
        // Read label of both vertices
        std::getline(file, line);
        if(line == "")
            break;
        std::istringstream ss(line);
        std::getline(ss, vertex1, ' ');
        std::getline(ss, vertex2, ' ');
        num_edges++;

        // If self-loop, discard even before creating vertex
        if(vertex1 == vertex2)
            continue;
        

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

        // Add edge to graph g
        boost::add_edge(id1, id2, g);
    }

    if(boost::num_vertices(g) != max_vertices) {
        std::cerr << "Wrong number of vertices read. File may be corrupted!" << std::endl;
        std::cerr << "Expected: " << max_vertices << ". Read: " << boost::num_vertices(g) << std::endl;
        return false;
    }

    if(num_edges != max_edges) {
        std::cerr << "Wrong number of edges read. File may be corrupted!" << std::endl;
        std::cerr << "Expected: " << max_edges << ". Read: " << boost::num_edges(g) << std::endl;
        return false;
    }

    return true;
}
//
// Created by 罗程阳 on 2024/3/16.
//

#ifndef MLCS_PROJECTED_GRAPH_H
#define MLCS_PROJECTED_GRAPH_H

#include <vector>
#include <string>
#include <map>
#include "multilayer_graph.h"
class ProjectedGraph {
public:
    int number_of_vertex;
    int maximum_vertex;
    int maximum_degree;
    std::vector<int> vertex_iterator;
    std::vector<std::vector<int>> adjacency_list;
    std::vector<std::map<int, int>> adjacency_list_loc;

public:
    ProjectedGraph();
    ProjectedGraph(MultilayerGraph multilayer_graph);
    void delete_vertex(int vertex);
};
#endif //MLCS_PROJECTED_GRAPH_H

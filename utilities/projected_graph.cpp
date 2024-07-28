//
// Created by 罗程阳 on 2024/3/16.
//

#include "projected_graph.h"
#include <set>
#include <algorithm>
ProjectedGraph::ProjectedGraph() {
}

ProjectedGraph::ProjectedGraph(MultilayerGraph multilayer_graph) {
    this->maximum_vertex = multilayer_graph.maximum_vertex;
    this->number_of_vertex = multilayer_graph.number_of_vertex;
    this->vertex_iterator = multilayer_graph.vertex_iterator;
    this->adjacency_list.resize(maximum_vertex + 1);
    this->adjacency_list_loc.resize(maximum_vertex + 1);
    int number_of_layers = multilayer_graph.number_of_layers;
    maximum_degree = 0;
    for (auto vertex : vertex_iterator) {
        std::set<int> neighbors;
        for (int layer = 0; layer < number_of_layers; layer++) {
            for (auto &neighbor : multilayer_graph.adjacency_list[vertex][layer])
            neighbors.insert(neighbor);
        }
        this->adjacency_list[vertex].reserve(neighbors.size());
        maximum_degree = std::max(maximum_degree, (int)neighbors.size());
        int loc = 0;
        for (auto neighbor : neighbors) {
            adjacency_list[vertex].push_back(neighbor);
            adjacency_list_loc[vertex][neighbor] = loc;
            loc++;
        }
    }
}

void ProjectedGraph::delete_vertex(int vertex) {
    for (auto &neighbor : adjacency_list[vertex]) {
        if (adjacency_list_loc[neighbor].find(vertex) != adjacency_list_loc[neighbor].end()) {
            int loc = adjacency_list_loc[neighbor][vertex];
            int last_vertex = adjacency_list[neighbor].back();
            adjacency_list[neighbor][loc] = last_vertex;
            adjacency_list_loc[neighbor][last_vertex] = loc;
            adjacency_list[neighbor].pop_back();
            adjacency_list_loc[neighbor].erase(vertex);
        }
    }
    adjacency_list[vertex].clear();
    adjacency_list_loc[vertex].clear();

//    for (auto &neighbor : adjacency_list[vertex]) {
//        if (adjacency_list_loc[neighbor].find(vertex) != adjacency_list_loc[neighbor].end()) {
//            auto iter = std::find(adjacency_list[neighbor].begin(), adjacency_list[neighbor].end(), vertex);
//            if (iter != adjacency_list[neighbor].end()) {
//                adjacency_list[neighbor].erase(iter);
//            }
//        }
//    }
//    adjacency_list[vertex].clear();
//    adjacency_list_loc[vertex].clear();
}
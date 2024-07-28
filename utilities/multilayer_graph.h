//
// Created by 罗程阳 on 2024/3/15.
//

#ifndef MLCS_MULTILAYER_GRAPH_H
#define MLCS_MULTILAYER_GRAPH_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>

using namespace std;
class MultilayerGraph {

public:
    MultilayerGraph();
    MultilayerGraph(std::string dataset_path);
    void load_dataset(std::string dataset_path);
    void add_edge(int from_vertex, int to_vertex, int layer);
    std::unordered_map<int, int> order_layers(std::ifstream& dataset_file);
    std::vector<int> get_vertex();
    double get_number_of_edges(int layer = -1);
    std::unordered_map<int, double> get_number_of_edges_layer_by_layer();
    int get_layer_mapping(int layer);
    void delete_vertex(int vertex);

public:
    int number_of_layers;
    std::vector<int> layers_iterator;
    std::unordered_map<int, int> layers_map;

    int number_of_vertex;
    int maximum_vertex;
    std::vector<int> vertex_iterator;
    std::vector<std::vector<std::vector<int>>> adjacency_list;
    std::vector<std::vector<std::map<int, int>>> adjacency_list_loc;
    std::string dataset_path;
};

#endif //MLCS_MULTILAYER_GRAPH_H

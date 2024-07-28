//
// Created by 罗程阳 on 2024/3/20.
//

#ifndef MLCS_SEARCH_SPACE_H
#define MLCS_SEARCH_SPACE_H

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include "multilayer_graph.h"
#include "projected_graph.h"
using namespace std;

class BitMap {
public:
    vector<unsigned int> value;
    int num_one;
    int num_total;
    int size_total;
public:
    BitMap();
    BitMap(int layer_num);
    void set_one(int pos);
    int count_one();
    void print_value();
    bool compare_other_bitmap(vector<unsigned int> compare_value);
};

class SearchNode {
public:
    vector<int> candidate_vertex;
    unordered_set<int> candidate_vertex_set;
    map<int, int> vertex_pos;
    int d;
    int vertex_num;
    int layer_num;
    map<int, vector<vector<int>>> adjacency_list;
    map<int, vector<int>> adjacency_list_projected;
    map<int, vector<map<int, int>>> adjacency_list_loc;
    map<int, map<int, int>> adjacency_list_projected_loc;
    vector<int> layer_composition;
    unordered_map<int, int> top_d_degree;
    unordered_map<int, vector<int>> vertex_layer_degree;
    unordered_map<int, int> vertex_projected_degree;
    bool is_prune1;
    bool is_construct_graph;
    unordered_map<int, BitMap> vertex_bitmap;
    bool is_search;
public:
    SearchNode();
    SearchNode(vector<int> &vertex_vec, int num_of_layer, int d);
    void initialize_vec(vector<int> &vertex_vec, int num_of_layer, int d);
    void get_position();
    void get_top_d_degree_delete();
    void delete_vertex(int v);
    int get_com_degree(int v);

    // do not delete vertex from origin graph
    void get_degree(MultilayerGraph &multilayer_graph, ProjectedGraph &projected_graph);
    void get_top_d_degree();
    void get_degree_composition(MultilayerGraph &multilayer_graph, ProjectedGraph &projected_graph);
    void get_degree_composition_vec(MultilayerGraph &multilayer_graph, ProjectedGraph &projected_graph, vector<bool> &vertex_exist);
    void get_top_d_degree_composition();
    int update_vertex_degree(int updated_vertex, int deleted_vertex, MultilayerGraph &multilayer_graph, ProjectedGraph &projected_graph);
    void cal_bitmap(int k);
};


#endif //MLCS_SEARCH_SPACE_H

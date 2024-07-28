//
// Created by 罗程阳 on 2024/3/20.
//

#include "search_space.h"
#include <algorithm>

SearchNode::SearchNode() {
    is_prune1 = false;
    is_search = false;
    is_construct_graph = false;
}

SearchNode::SearchNode(vector<int> &vertex_vec, int num_of_layer, int d) {
    is_prune1 = false;
    is_construct_graph = false;
    is_search = false;
    candidate_vertex = vertex_vec;
    for (auto vertex : candidate_vertex) {
        candidate_vertex_set.insert(vertex);
    }
    vertex_num = candidate_vertex.size();
    this->d = d;
    this->layer_num = num_of_layer;
}

void SearchNode::initialize_vec(vector<int> &vertex_vec, int num_of_layer, int d) {
    is_prune1 = false;
    is_construct_graph = false;
    is_search = false;
    candidate_vertex = vertex_vec;
    vertex_num = candidate_vertex.size();
    this->d = d;
    this->layer_num = num_of_layer;
}
void SearchNode::get_position() {
    for (int i = 0; i < candidate_vertex.size(); i++) {
        vertex_pos[candidate_vertex[i]] = i;
    }
}

void SearchNode::get_top_d_degree_delete() {
    for (auto vertex : candidate_vertex) {
        vector<int> degrees;
        for (int layer = 0; layer < layer_num; layer++) {
            degrees.push_back(adjacency_list[vertex][layer].size());
        }
        std::sort(degrees.begin(), degrees.end(), [](int a, int b) {
            return a > b;
        });
        top_d_degree[vertex] = degrees[d-1];
    }
}

int SearchNode::get_com_degree(int v) {
    vector<int> degrees;
    for (int layer = 0; layer < layer_num; layer++) {
        degrees.push_back(adjacency_list[v][layer].size());
    }
    std::sort(degrees.begin(), degrees.end(), [](int a, int b) {
        return a > b;
    });
    top_d_degree[v] = degrees[d-1];
    int min_deg = min(top_d_degree[v], (int)adjacency_list_projected[v].size() - 1);
    if (min_deg < 0)
        min_deg = 0;
    return min_deg;
//    return top_d_degree[v];
}
void SearchNode::delete_vertex(int vertex) {
    // delete vertex from projected graph
    for (auto &neighbor : adjacency_list_projected[vertex]) {
        auto iter = std::find(adjacency_list_projected[neighbor].begin(), adjacency_list_projected[neighbor].end(), vertex);
        if (iter != adjacency_list_projected[neighbor].end()) {
            adjacency_list_projected[neighbor].erase(iter);
        }
    }
    adjacency_list_projected.erase(vertex);
    adjacency_list_projected_loc.erase(vertex);

    // delete vertex from multilayer graph
    for (int i = 0; i < layer_num; i++) {
        for (auto &neighbor: adjacency_list[vertex][i]) {
            auto iter = std::find(adjacency_list[neighbor][i].begin(), adjacency_list[neighbor][i].end(), vertex);
            if (iter != adjacency_list[neighbor][i].end()) {
                adjacency_list[neighbor][i].erase(iter);
            }
        }
        adjacency_list[vertex][i].clear();
        adjacency_list_loc[vertex][i].clear();
    }
    adjacency_list.erase(vertex);
    adjacency_list_loc.erase(vertex);

//    // delete vertex from projected graph
//    for (auto &neighbor : adjacency_list_projected[vertex]) {
//        if (adjacency_list_projected_loc[neighbor].find(vertex) != adjacency_list_projected_loc[neighbor].end()) {
//            int loc = adjacency_list_projected_loc[neighbor][vertex];
//            int last_vertex = adjacency_list_projected[neighbor].back();
//            adjacency_list_projected[neighbor][loc] = last_vertex;
//            adjacency_list_projected_loc[neighbor][last_vertex] = loc;
//            adjacency_list_projected[neighbor].pop_back();
//            adjacency_list_projected_loc[neighbor].erase(vertex);
//        }
//    }
//    adjacency_list_projected.erase(vertex);
//    adjacency_list_projected_loc.erase(vertex);
//
//    // delete vertex from multilayer graph
//    for (int i = 0; i < layer_num; i++) {
//        for (auto &neighbor : adjacency_list[vertex][i]) {
//            if (adjacency_list_loc[neighbor][i].find(vertex) != adjacency_list_loc[neighbor][i].end()) {
//                int loc = adjacency_list_loc[neighbor][i][vertex];
//                int last_vertex = adjacency_list[neighbor][i].back();
//                adjacency_list[neighbor][i][loc] = last_vertex;
//                adjacency_list_loc[neighbor][i][last_vertex] = loc;
//                adjacency_list[neighbor][i].pop_back();
//                adjacency_list_loc[neighbor][i].erase(vertex);
//            }
//        }
//    }
//    adjacency_list.erase(vertex);
//    adjacency_list_loc.erase(vertex);
}

void SearchNode::get_degree(MultilayerGraph &multilayer_graph, ProjectedGraph &projected_graph) {
    // get projected degree of each vertex
    for (auto &vertex : candidate_vertex_set) {
        int temp_deg = 0;
        for (auto &neighbor : projected_graph.adjacency_list[vertex]) {
            if (candidate_vertex_set.find(neighbor) != candidate_vertex_set.end()) {
                temp_deg++;
            }
        }
        vertex_projected_degree[vertex] = temp_deg;
    }

    for (auto &vertex : candidate_vertex_set) {
        vector<int> temp_layer_degree;
        for (int layer = 0; layer < layer_num; layer++) {
            int temp_deg = 0;
            for (auto &neighbor: multilayer_graph.adjacency_list[vertex][layer]) {
                if (candidate_vertex_set.find(neighbor) != candidate_vertex_set.end()) {
                    temp_deg++;
                }
            }
            temp_layer_degree.push_back(temp_deg);
        }
        vertex_layer_degree[vertex] = temp_layer_degree;
    }
}

void SearchNode::get_top_d_degree() {
    for (auto vertex : candidate_vertex) {
        vector<int> degrees = vertex_layer_degree[vertex];
        std::sort(degrees.begin(), degrees.end(), [](int a, int b) {
            return a > b;
        });
        top_d_degree[vertex] = degrees[d-1];
    }
}

void SearchNode::get_degree_composition(MultilayerGraph &multilayer_graph, ProjectedGraph &projected_graph) {
    // get projected degree of each vertex
    for (auto &vertex : candidate_vertex) {
        int temp_deg = 0;
        for (auto &neighbor : projected_graph.adjacency_list[vertex]) {
            if (candidate_vertex_set.find(neighbor) != candidate_vertex_set.end()) {
                temp_deg++;
            }
        }
        vertex_projected_degree[vertex] = temp_deg;
    }

    for (auto &vertex : candidate_vertex_set) {
        vector<int> temp_layer_degree(layer_num, -1);
        for (int layer : layer_composition) {
            int temp_deg = 0;
            for (auto &neighbor: multilayer_graph.adjacency_list[vertex][layer]) {
                if (candidate_vertex_set.find(neighbor) != candidate_vertex_set.end()) {
                    temp_deg++;
                }
            }
            temp_layer_degree[layer] = temp_deg;
        }
        vertex_layer_degree[vertex] = temp_layer_degree;
    }
}

void SearchNode::get_degree_composition_vec(MultilayerGraph &multilayer_graph, ProjectedGraph &projected_graph, vector<bool> &vertex_exist) {
    // get projected degree of each vertex
    for (auto &vertex : candidate_vertex) {
        int temp_deg = 0;
        for (auto &neighbor : projected_graph.adjacency_list[vertex]) {
            if (vertex_exist[neighbor]) {
                temp_deg++;
            }
        }
        vertex_projected_degree[vertex] = temp_deg;
    }

    for (auto &vertex : candidate_vertex) {
        vector<int> temp_layer_degree(layer_num, -1);
        for (int layer : layer_composition) {
            int temp_deg = 0;
            for (auto &neighbor: multilayer_graph.adjacency_list[vertex][layer]) {
                if (vertex_exist[neighbor]) {
                    temp_deg++;
                }
            }
            temp_layer_degree[layer] = temp_deg;
        }
        vertex_layer_degree[vertex] = temp_layer_degree;
    }
}

void SearchNode::get_top_d_degree_composition() {
    for (auto vertex : candidate_vertex_set) {
        vector<int> degrees = vertex_layer_degree[vertex];
        std::sort(degrees.begin(), degrees.end(), [](int a, int b) {
            return a > b;
        });
        top_d_degree[vertex] = degrees[d-1];
    }
}

int SearchNode::update_vertex_degree(int updated_vertex, int deleted_vertex, MultilayerGraph &multilayer_graph, ProjectedGraph &projected_graph) {
//    for (int layer = 0; layer < layer_num; layer++) {
//        if (std::find(multilayer_graph.adjacency_list[updated_vertex][layer].begin(),
//                      multilayer_graph.adjacency_list[updated_vertex][layer].end(), deleted_vertex) != multilayer_graph.adjacency_list[updated_vertex][layer].end()) {
//            vertex_layer_degree[updated_vertex][layer]--;
//        }
//    }
    vector<int> degrees = vertex_layer_degree[updated_vertex];
    std::sort(degrees.begin(), degrees.end(), [](int a, int b) {
        return a > b;
    });
    top_d_degree[updated_vertex] = degrees[d-1];

    // update degree of projected graph
//    vertex_projected_degree[updated_vertex]--;
//    int temp_deg = 0;
//    for (auto &neighbor : projected_graph.adjacency_list[updated_vertex]) {
//        if (candidate_vertex_set.find(neighbor) != candidate_vertex_set.end()) {
//            temp_deg++;
//        }
//    }
//    if (vertex_projected_degree[updated_vertex] - temp_deg > 2)
//        cout << "?";
//    vertex_projected_degree[updated_vertex] = temp_deg;

    int new_deg = min(top_d_degree[updated_vertex], vertex_projected_degree[updated_vertex] - 1);
    if (new_deg < 0)
        new_deg = 0;
    return new_deg;
}

void SearchNode::cal_bitmap(int k) {
    for (auto vertex : candidate_vertex) {
        BitMap temp_bitmap(layer_num);
        for (int i = 0; i < layer_num; i++) {
            if (vertex_layer_degree[vertex][i] >= k) {
                temp_bitmap.set_one(i);
            }
        }
        vertex_bitmap[vertex] = temp_bitmap;
    }
}


BitMap::BitMap() {

}

BitMap::BitMap(int layer_num) {
    num_total = layer_num;
    for (int i = 0; i < layer_num / 32 + 1; i++) {
        value.push_back(0);
    }
    size_total = value.size();
}

void BitMap::set_one(int pos) {
    if (pos >= num_total) {
        cout << "fail to set one" << endl;
    } else {
        int pos1 = pos / 32;
        int pos2 = pos % 32;
        unsigned int temp_value = 0x01;
        temp_value = temp_value << pos2;
        value[pos1] |= temp_value;
    }
}

int BitMap::count_one() {
    int count = 0;
    for (int i = 0; i < size_total; i++) {
        unsigned int temp_value = value[i];
        for (int j = 0; j < 32; j++) {
            count += temp_value & 0x01;
            temp_value = temp_value >> 1;
        }
    }
    return count;
}

void BitMap::print_value() {
    for (int i = 0; i < num_total; i++) {
        int pos1 = i / 32;
        int pos2 = i % 32;
        unsigned int temp_value = 0x01;
        temp_value = temp_value << pos2;
        unsigned int bit = value[pos1] & temp_value;
        bit = bit >> pos2;
        cout << bit;
    }
}

bool BitMap::compare_other_bitmap(vector<unsigned int> compare_value) {
    for (int i = 0; i < value.size(); i++) {
        if ((value[i] & compare_value[i]) != compare_value[i]) {
            return false;
        }
    }
    return true;
}
//
// Created by 罗程阳 on 2024/4/7.
//


#include <iostream>
#include <algorithm>
#include <queue>
#include <map>
#include <set>
#include <time.h>
#include <unordered_set>
#include "utilities/multilayer_graph.h"
#include "utilities/projected_graph.h"
#include "utilities/search_space.h"
#include <fstream>
#include <iomanip>


using namespace std;

//string dataset_name = "higgs";
//string dataset_name = "dblp";
//string dataset_name = "slashdot";
//string dataset_name = "obamainisrael";
//string dataset_name = "yeast";
//string dataset_name = "homo";
//string dataset_name = "sacchcere";
//string dataset_name = "stack";
string dataset_name = "wiki";


int number_of_layers;
int number_of_edge;
MultilayerGraph multilayer_graph;
ProjectedGraph projected_graph;
vector<map<vector<int>, int>> vertex_core_number;

vector<vector<vector<vector<vector<int>>>>> PMC_index;

double avg_dominated_layers;

void generate_layer_composition(int layer, int d, int select, vector<vector<int>> &layer_composition,
                                vector<int> &layer_array, vector<bool> &layer_visited) {
    if (select == d) {
        vector<int> temp_layer_composition;
        for (int i = 0; i < number_of_layers; ++i) {
            if (layer_visited[i]) {
                temp_layer_composition.push_back(i);
            }
        }
        layer_composition.push_back(temp_layer_composition);
        return;
    }
    if (layer == number_of_layers)
        return;
    layer_visited[layer] = true;
    generate_layer_composition(layer + 1, d, select + 1, layer_composition, layer_array, layer_visited);
    layer_visited[layer] = false;
    generate_layer_composition(layer + 1, d, select, layer_composition, layer_array, layer_visited);
}


void decomposition_top_down_prune_vec2(vector<int> &layers, vector<int> &vertex_upper_bound,
                                       vector<vector<int>> &vertex_support, vector<int> &vertex_support_projected,
                                       vector<bool> &vertex_exist, vector<bool> &updated_vertex_flag) {
    // calculate support of each vertex
    for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
        if (!vertex_exist[v])
            continue;
        for (auto layer: layers) {
            int support_num = 0;
            for (auto neighbor: multilayer_graph.adjacency_list[v][layer]) {
                if (vertex_exist[neighbor] && vertex_upper_bound[neighbor] >= vertex_upper_bound[v]) {
                    support_num++;
                }
            }
            vertex_support[v][layer] = support_num;
        }
        int support_projected_num = 0;
        for (auto neighbor: projected_graph.adjacency_list[v]) {
            if (vertex_exist[neighbor] && vertex_upper_bound[neighbor] >= vertex_upper_bound[v]) {
                support_projected_num++;
            }
        }
        vertex_support_projected[v] = support_projected_num - 1;
        if (vertex_support_projected[v] < 0) {
            vertex_support_projected[v] = 0;
        }
    }
    queue<int> updated_vertex;
    std::fill(updated_vertex_flag.begin(), updated_vertex_flag.end(), false);
//    std::fill_n(updated_vertex_flag.begin(), multilayer_graph.maximum_vertex + 1, false);
    for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
        bool change = false;
        for (auto l: layers) {
            if (vertex_support[v][l] < vertex_upper_bound[v]) {
                change = true;
                break;
            }
        }
        if (vertex_support_projected[v] < vertex_upper_bound[v]) {
            change = true;
        }
        if (change) {
            updated_vertex.push(v);
//            updated_vertex_set.insert(v);
            updated_vertex_flag[v] = true;
        }
    }
    while (!updated_vertex.empty()) {
        vector<int> zero_vertex;
        int v = updated_vertex.front();
        updated_vertex.pop();
        bool change = false;
        for (auto l: layers) {
            if (vertex_support[v][l] < vertex_upper_bound[v]) {
                change = true;
                break;
            }
        }
        if (vertex_support_projected[v] < vertex_upper_bound[v])
            change = true;
        if (change) {
            int origin_upper_bound = vertex_upper_bound[v];
            // multilayer graph
            for (auto l: layers) {
                int upper_bound = multilayer_graph.adjacency_list[v][l].size();
                int lower_bound = 0;
                int mid = (upper_bound + lower_bound) / 2 + 1;
                while (lower_bound < upper_bound) {
                    int neighbor_num = 0;
                    for (auto neighbor: multilayer_graph.adjacency_list[v][l]) {
                        if (vertex_exist[neighbor] && vertex_upper_bound[neighbor] >= mid) {
                            neighbor_num++;
                        }
                    }
                    if (neighbor_num >= mid) {
                        lower_bound = mid;
                    } else {
                        upper_bound = mid - 1;
                    }
                    mid = (upper_bound + lower_bound) / 2 + 1;
                }
                int max_k_l = lower_bound;
                if (vertex_upper_bound[v] > max_k_l) {
                    vertex_upper_bound[v] = max_k_l;
                }
//                int max_k_l = multilayer_graph.adjacency_list[v].size();
//                for (; max_k_l > 0; max_k_l--) {
//                    int neighbor_num = 0;
//                    for (auto neighbor: multilayer_graph.adjacency_list[v][l]) {
//                        if (vertex_exist[neighbor] && vertex_upper_bound[neighbor] >= max_k_l) {
//                            neighbor_num++;
//                        }
//                    }
//                    if (neighbor_num >= max_k_l) {
//                        break;
//                    }
//                }
//                if (vertex_upper_bound[v] > max_k_l) {
//                    vertex_upper_bound[v] = max_k_l;
//                }
            }

            // projected graph
            int upper_bound = projected_graph.adjacency_list[v].size() - 1;
            int lower_bound = 0;
            int mid = (upper_bound + lower_bound) / 2 + 1;
            while (lower_bound < upper_bound) {
                int neighbor_num = 0;
                for (auto neighbor: projected_graph.adjacency_list[v]) {
                    if (vertex_exist[neighbor] && vertex_upper_bound[neighbor] >= mid) {
                        neighbor_num++;
                    }
                }
                if (neighbor_num >= mid + 1) {
                    lower_bound = mid;
                } else {
                    upper_bound = mid - 1;
                }
                mid = (upper_bound + lower_bound) / 2 + 1;
            }
            int max_k = lower_bound;
            if (vertex_upper_bound[v] > max_k) {
                vertex_upper_bound[v] = max_k;
            }
//            int max_k = projected_graph.adjacency_list[v].size() - 1;
//            for (; max_k > 0; max_k--) {
//                int neighbor_num = 0;
//                for (auto neighbor: projected_graph.adjacency_list[v]) {
//                    if (vertex_exist[neighbor] && vertex_upper_bound[neighbor] >= max_k) {
//                        neighbor_num++;
//                    }
//                }
//                if (neighbor_num >= max_k + 1) {
//                    break;
//                }
//            }
//            if (vertex_upper_bound[v] > max_k) {
//                vertex_upper_bound[v] = max_k;
//            }

            // multilayer graph
            for (auto l: layers) {
                int support_neighbor_num = 0;
                for (auto neighbor: multilayer_graph.adjacency_list[v][l]) {
                    if (vertex_exist[neighbor] && vertex_upper_bound[neighbor] >= vertex_upper_bound[v]) {
                        support_neighbor_num++;
                    }
                }
                vertex_support[v][l] = support_neighbor_num;
                for (auto neighbor: multilayer_graph.adjacency_list[v][l]) {
                    if (vertex_exist[neighbor] && vertex_upper_bound[neighbor] <= origin_upper_bound &&
                        vertex_upper_bound[neighbor] > vertex_upper_bound[v]) {
                        vertex_support[neighbor][l]--;
                        if (!updated_vertex_flag[neighbor]) {
                            updated_vertex.push(neighbor);
                            updated_vertex_flag[neighbor] = true;
                        }
                    }
                }
            }

            // projected graph
            int support_neighbor_num = 0;
            for (auto neighbor: projected_graph.adjacency_list[v]) {
                if (vertex_exist[neighbor] && vertex_upper_bound[neighbor] >= vertex_upper_bound[v]) {
                    support_neighbor_num++;
                }
            }
            vertex_support_projected[v] = support_neighbor_num - 1;
            if (vertex_support_projected[v] < 0) {
                vertex_support_projected[v] = 0;
            }
            for (auto neighbor: projected_graph.adjacency_list[v]) {
                if (vertex_exist[neighbor] && vertex_upper_bound[neighbor] <= origin_upper_bound &&
                    vertex_upper_bound[neighbor] > vertex_upper_bound[v]) {
                    vertex_support_projected[neighbor]--;
//                    if (updated_vertex_set.find(neighbor) == updated_vertex_set.end()) {
//                        updated_vertex_set.insert(neighbor);
//                        updated_vertex.push(neighbor);
//                    }
                    if (!updated_vertex_flag[neighbor]) {
                        updated_vertex_flag[neighbor] = true;
                        updated_vertex.push(neighbor);
                    }
                }
            }

            // record zero v
            if (vertex_upper_bound[v] <= 0) {
                vertex_exist[v] = false;
//                updated_vertex_set.erase(v);
                updated_vertex_flag[v] = false;
            } else {
                updated_vertex.push(v);
            }
        } else {
//            updated_vertex_set.erase(v);
            updated_vertex_flag[v] = false;
        }
    }

    for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
        if (vertex_upper_bound[v] != 0) {
            vertex_core_number[v][layers] = vertex_upper_bound[v];
//            if (vertex_upper_bound[v] > number_of_layers)
//                cout << "!!!!";
        }
    }
}

void decomposition() {
    vector<int> vertex_upper_bound(multilayer_graph.maximum_vertex + 1, 0);
    vector<int> vertex_support_projected(multilayer_graph.maximum_vertex + 1);
    vector<vector<int>> vertex_support(multilayer_graph.maximum_vertex + 1);
    vector<bool> vertex_exist(multilayer_graph.maximum_vertex + 1);
    vector<bool> updated_vertex_flag(multilayer_graph.maximum_vertex + 1, false);
    for (int i = 0; i < vertex_support.size(); i++) {
        vertex_support[i].resize(number_of_layers);
    }
    vertex_core_number.resize(multilayer_graph.maximum_vertex + 1);
    for (int layer_num = 1; layer_num <= number_of_layers; layer_num++) {
//    for (int layer_num = 1; layer_num <= 3; layer_num++) {
        cout << layer_num << endl;
        // generate layer composition
        vector<vector<int>> layer_composition;
        vector<int> layer_array = multilayer_graph.layers_iterator;
        vector<bool> layer_visited(multilayer_graph.number_of_layers, false);
        generate_layer_composition(0, layer_num, 0, layer_composition, layer_array, layer_visited);
        for (auto layers: layer_composition) {
            for (auto layer_id: layers) {
                cout << layer_id << " ";
            }
            cout << endl;
//            vertex_upper_bound.clear();
//            vertex_support.clear();
//            unordered_set<int> vertex_set;
//            set<int> vertex_set;

            bool non_zero_flag = false;
            if (layers.size() == 1) {
                for (auto &vertex: multilayer_graph.vertex_iterator) {
                    vertex_upper_bound[vertex] = multilayer_graph.adjacency_list[vertex][layers[0]].size();
                    int vertex_projected_degree = projected_graph.adjacency_list[vertex].size();
                    if (vertex_projected_degree - 1 < vertex_upper_bound[vertex]) {
                        vertex_upper_bound[vertex] = vertex_projected_degree - 1;
                        if (vertex_upper_bound[vertex] < 0) {
                            vertex_upper_bound[vertex] = 0;
                        }
                    }
                    if (vertex_upper_bound[vertex] != 0) {
//                        vertex_set.insert(vertex);
                        vertex_exist[vertex] = true;
                        non_zero_flag = true;
                    } else {
                        vertex_exist[vertex] = false;
                    }
                }
            } else {
                // generate sub layers
                vector<vector<int>> sub_layer_composition;
                for (int i = 0; i < layers.size(); i++) {
                    vector<int> temp_sub_layer_composition;
                    for (int j = 0; j < layers.size(); j++) {
                        if (j != i) {
                            temp_sub_layer_composition.push_back(layers[j]);
                        }
                    }
                    sub_layer_composition.push_back(temp_sub_layer_composition);
                }
                for (int vertex = 0; vertex < multilayer_graph.maximum_vertex + 1; vertex++) {
                    int min_core_number = INT32_MAX;
                    for (auto sub_layer: sub_layer_composition) {
                        if (vertex_core_number[vertex].find(sub_layer) != vertex_core_number[vertex].end()) {
                            min_core_number = min(min_core_number, vertex_core_number[vertex][sub_layer]);
                        }
                    }
                    if (min_core_number != INT32_MAX && min_core_number > 0) {
                        vertex_upper_bound[vertex] = min_core_number;
//                        vertex_set.insert(vertex);
                        vertex_exist[vertex] = true;
                        non_zero_flag = true;
                    } else {
                        vertex_upper_bound[vertex] = 0;
                        vertex_exist[vertex] = false;
                    }
                }
            }
            if (non_zero_flag) {
//                decomposition_top_down(layers, vertex_upper_bound, vertex_support, vertex_support_projected);
//                decomposition_top_down_prune_set(layers, vertex_upper_bound, vertex_support, vertex_support_projected, vertex_set);
//                decomposition_top_down_prune_vec(layers, vertex_upper_bound, vertex_support, vertex_support_projected, vertex_exist);
                decomposition_top_down_prune_vec2(layers, vertex_upper_bound, vertex_support, vertex_support_projected, vertex_exist, updated_vertex_flag);
            } else {
                cout << "All Zero!" << endl;
            }
        }
    }
}

void decomposition_top_down_prune_no_projected(vector<int> &layers, vector<int> &vertex_upper_bound,
                                       vector<vector<int>> &vertex_support, vector<int> &vertex_support_projected,
                                       vector<bool> &vertex_exist, vector<bool> &updated_vertex_flag) {
    // calculate support of each vertex
    for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
        if (!vertex_exist[v])
            continue;
        for (auto layer: layers) {
            int support_num = 0;
            for (auto neighbor: multilayer_graph.adjacency_list[v][layer]) {
                if (vertex_exist[neighbor] && vertex_upper_bound[neighbor] >= vertex_upper_bound[v]) {
                    support_num++;
                }
            }
            vertex_support[v][layer] = support_num;
        }
    }
    queue<int> updated_vertex;
    std::fill(updated_vertex_flag.begin(), updated_vertex_flag.end(), false);
//    std::fill_n(updated_vertex_flag.begin(), multilayer_graph.maximum_vertex + 1, false);
    for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
        bool change = false;
        for (auto l: layers) {
            if (vertex_support[v][l] < vertex_upper_bound[v]) {
                change = true;
                break;
            }
        }
        if (change) {
            updated_vertex.push(v);
//            updated_vertex_set.insert(v);
            updated_vertex_flag[v] = true;
        }
    }
    while (!updated_vertex.empty()) {
        vector<int> zero_vertex;
        int v = updated_vertex.front();
        updated_vertex.pop();
        bool change = false;
        for (auto l: layers) {
            if (vertex_support[v][l] < vertex_upper_bound[v]) {
                change = true;
                break;
            }
        }
        if (change) {
            int origin_upper_bound = vertex_upper_bound[v];
            // multilayer graph
            for (auto l: layers) {
                int upper_bound = multilayer_graph.adjacency_list[v][l].size();
                int lower_bound = 0;
                int mid = (upper_bound + lower_bound) / 2 + 1;
                while (lower_bound < upper_bound) {
                    int neighbor_num = 0;
                    for (auto neighbor: multilayer_graph.adjacency_list[v][l]) {
                        if (vertex_exist[neighbor] && vertex_upper_bound[neighbor] >= mid) {
                            neighbor_num++;
                        }
                    }
                    if (neighbor_num >= mid) {
                        lower_bound = mid;
                    } else {
                        upper_bound = mid - 1;
                    }
                    mid = (upper_bound + lower_bound) / 2 + 1;
                }
                int max_k_l = lower_bound;
                if (vertex_upper_bound[v] > max_k_l) {
                    vertex_upper_bound[v] = max_k_l;
                }
//                int max_k_l = multilayer_graph.adjacency_list[v][l].size();
//                for (; max_k_l > 0; max_k_l--) {
//                    int neighbor_num = 0;
//                    for (auto neighbor: multilayer_graph.adjacency_list[v][l]) {
//                        if (vertex_exist[neighbor] && vertex_upper_bound[neighbor] >= max_k_l) {
//                            neighbor_num++;
//                        }
//                    }
//                    if (neighbor_num >= max_k_l) {
//                        break;
//                    }
//                }
//                if (vertex_upper_bound[v] > max_k_l) {
//                    vertex_upper_bound[v] = max_k_l;
//                }
            }


            // multilayer graph
            for (auto l: layers) {
                int support_neighbor_num = 0;
                for (auto neighbor: multilayer_graph.adjacency_list[v][l]) {
                    if (vertex_exist[neighbor] && vertex_upper_bound[neighbor] >= vertex_upper_bound[v]) {
                        support_neighbor_num++;
                    }
                }
                vertex_support[v][l] = support_neighbor_num;
                for (auto neighbor: multilayer_graph.adjacency_list[v][l]) {
                    if (vertex_exist[neighbor] && vertex_upper_bound[neighbor] <= origin_upper_bound &&
                        vertex_upper_bound[neighbor] > vertex_upper_bound[v]) {
                        vertex_support[neighbor][l]--;
                        if (!updated_vertex_flag[neighbor]) {
                            updated_vertex.push(neighbor);
                            updated_vertex_flag[neighbor] = true;
                        }
                    }
                }
            }

            // record zero v
            if (vertex_upper_bound[v] <= 0) {
                vertex_exist[v] = false;
                updated_vertex_flag[v] = false;
            } else {
                updated_vertex.push(v);
            }
            updated_vertex.push(v);
        } else {
//            updated_vertex_set.erase(v);
            updated_vertex_flag[v] = false;
        }
    }

    for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
        if (vertex_upper_bound[v] != 0) {
            vertex_core_number[v][layers] = vertex_upper_bound[v];
            if (vertex_upper_bound[v] > number_of_layers)
                cout << "!!!!";
        }
    }
}

void construct_index() {
    PMC_index.resize(multilayer_graph.maximum_vertex + 1);
    for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
//    for (int v = 0; v < 10; v++) {
//        cout << v << endl;
        PMC_index[v].push_back({});
        for (int k = 1; k <= projected_graph.adjacency_list[v].size(); k++) {
            vector<vector<vector<int>>> temp_layers;
            temp_layers.resize(number_of_layers + 1);
            int layers_count = 0;
            for (auto core_number: vertex_core_number[v]) {
                vector<int> layers = core_number.first;
                int temp_k = core_number.second;
                if (temp_k >= k) {
                    temp_layers[layers.size()].push_back(layers);
                    layers_count++;
                }
            }
            if (layers_count == 0) {
                break;
            }
            for (int l = number_of_layers; l > 0; l--) {
                for (auto layers: temp_layers[l]) {
                    for (int ll = l - 1; ll > 0; ll--) {
                        vector<vector<int>> new_layers;
                        for (auto compare_layers: temp_layers[ll]) {
                            bool dominate = true;
                            for (int i = 0; i < compare_layers.size(); i++) {
                                if (std::find(layers.begin(), layers.end(), compare_layers[i]) == layers.end()) {
                                    dominate = false;
                                    break;
                                }
                            }
                            if (!dominate) {
                                new_layers.push_back(compare_layers);
                            }
                        }
                        temp_layers[ll] = new_layers;
                    }
                }
            }
            PMC_index[v].push_back(temp_layers);

//            cout << "xx";
//            for (int i = 1; i < temp_layers.size(); i++) {
//                cout << i << endl;
//                for (auto layers : temp_layers[i]) {
//                    for (auto layer : layers) {
//                        cout << layer << " ";
//                    }
//                    cout << endl;
//                }
//            }
        }
    }

    int total_dominate_layers = 0;
    //// write index
    ofstream index_file("datasets/" + dataset_name + "/" + dataset_name + "_index.txt");
    for (int v = 0; v < multilayer_graph.maximum_vertex + 1; v++) {
        cout << "v " << v << endl;
        index_file << "v " << v << endl;
        for (int k = 1; k < PMC_index[v].size(); k++) {
            cout << "k " << k << endl;
            index_file << "k " << k << endl;
            for (int d = 1; d < PMC_index[v][k].size(); d++) {
                if (PMC_index[v][k][d].size() == 0)
                    continue;
//                cout << "d=" << d << endl;
//                index_file << "d=" << d << endl;
                for (auto layers : PMC_index[v][k][d]) {
                    for (auto layer : layers) {
                        cout << layer << " ";
                        index_file << layer << " ";
                    }
                    total_dominate_layers++;
                    cout << endl;
                    index_file << endl;
                }
            }
        }
    }
    index_file.close();
    avg_dominated_layers = (double)total_dominate_layers / (double)multilayer_graph.maximum_vertex;
    cout << "Average dominated layers:" << (double)total_dominate_layers / (double)multilayer_graph.maximum_vertex << endl;
}


int main(int argc, char* argv[]) {
    srand(time(0));
    clock_t start_time, end_time;
    dataset_name = argv[1];
    string dataset_path = "datasets/" + dataset_name + "/" + dataset_name + ".txt";
    cout << "Loading graph...." << endl;
    MultilayerGraph multilayer_graph_dataset(dataset_path);
    multilayer_graph = multilayer_graph_dataset;
    ProjectedGraph projected_graph_dataset(multilayer_graph);
    projected_graph = projected_graph_dataset;
    number_of_layers = multilayer_graph.number_of_layers;

    number_of_edge = 0;
    for (auto v = 0; v < multilayer_graph.adjacency_list.size(); v++) {
        for (auto l = 0; l < number_of_layers; l++) {
            number_of_edge += multilayer_graph.adjacency_list[v][l].size();
        }
    }
    number_of_edge = number_of_edge / 2.0;
    cout << "Number of Edges : " << number_of_edge << endl;

    //// decomposition
    start_time = clock();
//    decomposition_no_projected();
    ofstream out_file("datasets/" + dataset_name + "/" + dataset_name + "_decomposition.txt");
    decomposition();

    construct_index();

    end_time = clock();
    cout << "Decomposition Time: " << (double) (end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
    out_file << "Decomposition Time: " << (double) (end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
    out_file << "Average dominated layers:" << avg_dominated_layers << endl;
    out_file.close();
    return 0;
}
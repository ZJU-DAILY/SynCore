//
// Created by 罗程阳 on 2024/3/15.
//
#include <iostream>
#include <algorithm>
#include "utilities/multilayer_graph.h"
#include "utilities/projected_graph.h"
#include <time.h>
#include <sstream>
#include <chrono>
#include <unordered_set>

using namespace std;

string dataset_name = "stack";
//string dataset_name = "dblp";
//string dataset_name = "higgs";
//string dataset_name = "homo";
int k = 0;
int d = 0;
int number_of_layers;
int k_1_num = 0;
vector<int> k_1_vec;
vector<int> result;

vector<int> peeling_projected_graph(ProjectedGraph &projected_graph) {
    vector<int> bin;  // store vertices
    bin.reserve(projected_graph.maximum_vertex + 1);
    vector<int> deg(projected_graph.maximum_vertex + 1);  // the degree of each vertex
    vector<int> pos(projected_graph.maximum_vertex + 1);  // the position of each vertex in bin
    vector<int> bou(projected_graph.maximum_degree + 1);  // the boundary of each degree value in bin

    //// initialize arrays
    vector<vector<int>> degree_set(projected_graph.maximum_degree + 1);
    for (auto &vertex : projected_graph.vertex_iterator) {
        int temp_deg = projected_graph.adjacency_list[vertex].size();
        degree_set[temp_deg].push_back(vertex);
        deg[vertex] = temp_deg;
    }
    int pos_count = 0;
    for (int i = 0; i < projected_graph.maximum_degree + 1; ++i) {
        bou[i] = pos_count;
        for (auto &vertex : degree_set[i]) {
            bin.push_back(vertex);
            pos[vertex] = pos_count;
            pos_count++;
        }
    }

    //// start peeling
    for (int i = 0; i < projected_graph.maximum_vertex + 1; i++) {
        int v = bin[i];
        if (i >= bou[k + 1]) {
            break;
        }
        for (auto &u : projected_graph.adjacency_list[v]) {
            if (deg[u] > deg[v]) {
                int deg_u = deg[u];
                int pos_u = pos[u];
                int pos_w = bou[deg_u];  // w is the first vertex in the group that u is in;
                int w = bin[pos_w];
                if (u != w) {
                    // swap u and w;
                    bin[pos_u] = w;
                    bin[pos_w] = u;
                    pos[u] = pos_w;
                    pos[w] = pos_u;
                }
                bou[deg_u]++;
                deg[u]--;
            }
        }
    }

    //// get projected (k+1)-core
    vector<int> projected_result;
    for (int i = bou[k+1]; i < bin.size(); i++) {
        projected_result.push_back(bin[i]);
    }
    cout << "num of (k+1)-core:";
    cout << projected_result.size() << endl;
    k_1_num = projected_result.size();
    k_1_vec = projected_result;
    return projected_result;
}

void generate_layer_composition(int layer, int select, vector<vector<int>> &layer_composition, vector<int> &layer_array, vector<bool> &layer_visited) {
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
    generate_layer_composition(layer + 1, select + 1, layer_composition, layer_array, layer_visited);
    layer_visited[layer] = false;
    generate_layer_composition(layer + 1, select, layer_composition, layer_array, layer_visited);
}

void MLCS(MultilayerGraph &multilayer_graph) {
    //// construct the projected graph of the multilayer graph
    ProjectedGraph projected_graph(multilayer_graph);

    //// search connected (k+1)-core in the projected graph
    vector<int>projected_k1_core = peeling_projected_graph(projected_graph);

    //// get the composition of d layers
    vector<vector<int>> layer_composition;
    vector<int> layer_array = multilayer_graph.layers_iterator;
    vector<bool> layer_visited(multilayer_graph.number_of_layers, false);
    generate_layer_composition(0, 0, layer_composition, layer_array, layer_visited);

    unordered_set<int> candidate_vertex;
    vector<bool> vertex_exist(multilayer_graph.maximum_vertex + 1, false);
    for (auto &vertex : projected_k1_core) {
        candidate_vertex.insert(vertex);
        vertex_exist[vertex] = true;
    }
    for (auto &composition : layer_composition) {
        for (auto &num : composition) {
            cout << num << " ";
        }
        cout << endl;
        auto temp_candidate_vertex = candidate_vertex;
        auto temp_vertex_exist = vertex_exist;

        // compute minimum degree of vertex in L layers
        map<int, vector<int>> degree_set;
        vector<int> deg(multilayer_graph.maximum_vertex + 1);
        int max_min_deg = 0;
//        unordered_map<int, vector<int>>
        for (auto &vertex : candidate_vertex) {
            int min_degree = multilayer_graph.maximum_vertex - 1;
            for (auto layer : composition) {
                int temp_layer_deg = 0;
                for (auto neighbor : multilayer_graph.adjacency_list[vertex][layer]) {
//                    if (temp_candidate_vertex.find(neighbor) != temp_candidate_vertex.end()) {
//                        temp_layer_deg++;
//                    }
                    if (vertex_exist[neighbor]) {
                        temp_layer_deg++;
                    }
                }
                min_degree = min(min_degree, temp_layer_deg);
            }
            int temp_projected_degree = 0;
            for (auto neighbor : projected_graph.adjacency_list[vertex]) {
//                if (temp_candidate_vertex.find(neighbor) != temp_candidate_vertex.end()) {
//                    temp_projected_degree++;
//                }
                if (vertex_exist[neighbor]) {
                    temp_projected_degree++;
                }
            }
            min_degree = min(min_degree, temp_projected_degree - 1);
            if (min_degree < 0)
                min_degree = 0;
            deg[vertex] = min_degree;
            degree_set[min_degree].push_back(vertex);
            max_min_deg = max(max_min_deg, min_degree);
        }

        vector<int> bin, pos, bou;
        pos.resize(multilayer_graph.maximum_vertex + 1);  // the index of v in the bin
        bou.resize(max_min_deg + 1);  // the start position of bou structure
        int pos_count = 0;
        for (int i = 0; i < max_min_deg + 1; i++) {
            bou[i] = pos_count;
            if (degree_set.find(i) != degree_set.end()) {
                for (auto vertex : degree_set[i]) {
                    bin.push_back(vertex);
                    pos[vertex] = pos_count;
                    pos_count++;
                }
            }
        }

        // start peeling
        int i = 0;
        for (i = 0; i < bin.size(); i++) {
            if (temp_candidate_vertex.size() < result.size())
                break;
//            if (i % 1000 == 0) {
//                cout << i << "/" << bin.size() << endl;
//            }
            int vertex = bin[i];
            if (deg[vertex] < k) {
                unordered_set<int> affect_neighbor;
//                for (auto layer : composition) {
//                    for (auto neighbor : multilayer_graph.adjacency_list[vertex][layer]) {
//                        if (temp_candidate_vertex.find(neighbor) != temp_candidate_vertex.end()) {
//                            affect_neighbor.insert(neighbor);
//                        }
//                    }
//                }
                for (auto neighbor : projected_graph.adjacency_list[vertex]) {
//                    if (temp_candidate_vertex.find(neighbor) != temp_candidate_vertex.end()) {
//                        affect_neighbor.insert(neighbor);
//                    }
                    if (temp_vertex_exist[neighbor]) {
                        affect_neighbor.insert(neighbor);
                    }
                }

                bou[deg[vertex]]++;
                if (deg[vertex] > 0) {
                    bou[deg[vertex] - 1] = bou[deg[vertex]];
                }
                temp_candidate_vertex.erase(vertex);
                temp_vertex_exist[vertex] = false;
                for (auto u : affect_neighbor) {
                    int min_degree = multilayer_graph.number_of_vertex - 1;
                    for (auto layer : composition) {
                        int temp_layer_deg = 0;
                        for (auto neighbor : multilayer_graph.adjacency_list[u][layer]) {
//                            if (temp_candidate_vertex.find(neighbor) != temp_candidate_vertex.end()) {
//                                temp_layer_deg++;
//                            }
                            if (temp_vertex_exist[neighbor]) {
                                temp_layer_deg++;
                            }
                        }
                        min_degree = min(min_degree, temp_layer_deg);
                    }
                    int temp_projected_degree = 0;
                    for (auto neighbor : projected_graph.adjacency_list[u]) {
//                        if (temp_candidate_vertex.find(neighbor) != temp_candidate_vertex.end()) {
//                            temp_projected_degree++;
//                        }
                        if (temp_vertex_exist[neighbor]) {
                            temp_projected_degree++;
                        }
                    }
                    min_degree = min(min_degree, temp_projected_degree - 1);
                    if (min_degree < 0)
                        min_degree = 0;
                    if (deg[u] - min_degree > 1)
                        cout << "!!!" << endl;
                    if (min_degree != deg[u]) {
                        int w = bin[bou[deg[u]]];
                        int pos_w = pos[w];
                        int pos_u = pos[u];
                        bin[pos_u] = w;
                        bin[pos_w] = u;
                        pos[u] = pos_w;
                        pos[w] = pos_u;
                        bou[deg[u]]++;
                        deg[u] = min_degree;
                    }
                }
            } else {
                break;
            }
        }
//        cout << "(k+1)-core: " << temp_result.size() - i<< endl;
        cout << "temp result: " << temp_candidate_vertex.size() << endl;
        cout << "current result: " << result.size() << endl;
        if (temp_candidate_vertex.size() > result.size()) {
            result.clear();
            for (i; i < bin.size(); i++) {
                result.push_back(bin[i]);
            }
        }
    }
}

void exp(MultilayerGraph multilayer_graph) {
    ifstream test_case_file("datasets/" + dataset_name + "/" + dataset_name + "_test_case.txt");
    ofstream exp_file("datasets/" + dataset_name + "/" + dataset_name + "_baseline_exp.txt");

    string line;
    while (getline(test_case_file, line)) {
        if (line.size() == 0)
            continue;
        istringstream iss(line);
        vector<string> words;
        string word;
        while (iss >> word) {
            words.push_back(word);
        }
        if (word.size() == 0)
            continue;
        result.clear();
        k = stoi(words[0]);
        d = stoi(words[1]);
        int query_vertex = stoi(words[2]);
        exp_file << "k:" << k << " " << "d:" << d << endl;
        cout << "Baseline on " << dataset_name << endl;
        cout << "k:" << k << " " << "d:" << d << endl;
        cout << "Query vertex :" << query_vertex << endl;
        exp_file << "Query vertex :" << query_vertex << endl;
        auto start = std::chrono::high_resolution_clock::now();
        MLCS(multilayer_graph);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        cout << "Result Num:" << result.size() << endl;
        exp_file << "Result Num:" << result.size() << endl;
        exp_file << "k+1 core num:" << k_1_num << endl;

        cout << "Search Time (s): " << duration.count() / 1000.0 << std::endl;
        exp_file << "Search Time (s): " << duration.count() / 1000.0 <<  endl << endl;

    }
    exp_file.close();
}

int main(int argc, char* argv[]) {
    srand(time(0));
    clock_t start_time, end_time;
    dataset_name = argv[1];
    string dataset_path = "datasets/" + dataset_name + "/" + dataset_name + ".txt";
    cout << "Loading graph...." << endl;
    // load dataset
    MultilayerGraph multilayer_graph(dataset_path);
    number_of_layers = multilayer_graph.number_of_layers;

    exp(multilayer_graph);
    return 0;
}